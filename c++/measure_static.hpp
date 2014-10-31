/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <vector>
#include <utility>
#include <boost/serialization/complex.hpp>
#include <boost/mpi/collectives.hpp>

#include <triqs/operators/many_body_operator.hpp>
#include "./triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp"
#include "./sorted_spaces.hpp"

namespace cthyb {

using triqs::arrays::matrix;

// Measure a static observable
struct measure_static {
 using mc_sign_type = std::complex<double>;

 qmc_data const& data;
 double & result;
 std::vector<std::pair<long,matrix<double>>> observable_matrices;   // pairs (subspace number,matrix elements)
 mc_sign_type z;
 int64_t num;

 measure_static(triqs::utility::many_body_operator<double> const& observable, double & result,
                qmc_data const& data, fundamental_operator_set const& fops) :
 result(result), data(data), z(0), num(0)
 {
  observable_matrices.reserve(data.sosp.n_subspaces());
  imperative_operator<hilbert_space> op(observable,fops);

  // Iterate over all subspaces;
  for(long spn = 0; spn < data.sosp.n_subspaces(); ++spn){
    auto M = data.sosp.make_op_matrix(op,spn,spn);
    if(M.second>0) observable_matrices.emplace_back(spn,M.first); // At least one non-zero element in this subspace
  }
 }

 void accumulate(mc_sign_type s) {

  num += 1;
  if (num < 0) TRIQS_RUNTIME_ERROR << " Overflow of counter ";

  auto corr = real(data.imp_trace.full_trace_over_estimator());
  if (!std::isfinite(corr)) TRIQS_RUNTIME_ERROR << " measure_static :corr not finite" << corr;

  z += s * corr;

  double numerator(0), denominator(0);

  for(auto const& M : observable_matrices){
   auto const& block_map = data.imp_trace.get_block_table()[M.first];
   auto const& trace_matrix = data.imp_trace.get_trace_matrices()[block_map];

   auto prod_matrix = trace_matrix * M.second;
   for(int n = 0; n < first_dim(trace_matrix); ++n){
    numerator += prod_matrix(n,n);
    denominator += trace_matrix(n,n);
   }
  }

  result += real(s) * (numerator/denominator) * corr;
 }

 void collect_results(boost::mpi::communicator const& c) {

  int64_t total_num;
  mc_sign_type total_z;
  boost::mpi::all_reduce(c, z, total_z, std::c14::plus<>());
  boost::mpi::all_reduce(c, num, total_num, std::c14::plus<>());

  double result_out;
  boost::mpi::all_reduce(c, result, result_out, std::c14::plus<>());
  result = result_out / real(total_z);
 }

};
}
