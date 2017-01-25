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
#include <triqs/gfs.hpp>
#include "./qmc_data.hpp"

namespace cthyb {
using namespace triqs::gfs;


struct measure_g_pp_tau {

  qmc_data const& data;
  std::vector<block_gf<imtime>>& g_pp_tau;
  mc_weight_t z;
  int64_t num;
  mc_weight_t average_sign;
  int n_g_blocks;
  double beta;
  block_gf<imtime> g_tmp;
  std::vector<matrix_t> const& density_matrix;
  size_t i_pp;

  measure_g_pp_tau(std::vector<block_gf<imtime>>& g_pp_tau, qmc_data const& data,
		   block_gf<imtime> const& g_tau, std::vector<matrix_t> const& density_matrix)
    : g_pp_tau(g_pp_tau), data(data), g_tmp(g_tau), density_matrix(density_matrix) {
    g_pp_tau.resize(data.h_diag.get_full_hilbert_space_dim(), g_tau);
    for (auto& g_i: g_pp_tau) g_i() = 0.0;
    z = 0;
    num = 0;
    n_g_blocks = n_blocks(g_pp_tau[0]);
  }

  void accumulate(mc_weight_t s) {
    num += 1;
    if (num < 0) TRIQS_RUNTIME_ERROR << " Overflow of counter ";
    data.imp_trace.compute();
    z += s * data.atomic_reweighting;
    s /= data.atomic_weight;
    g_tmp() = 0.0;
    for (int i_block = 0; i_block < n_g_blocks; ++i_block) {
      foreach(data.dets[i_block], [i_block, s, this](std::pair<time_pt, int> const& x, std::pair<time_pt, int> const& y, det_scalar_t M) {
	  this->g_tmp[i_block][closest_mesh_pt(double(y.first - x.first))](y.second, x.second) += (y.first >= x.first ? s : -s) * M;});
      for (size_t i_subspace = 0; i_subspace < density_matrix.size(); ++i_subspace) {
	if (data.imp_trace.get_density_matrix()[i_subspace].is_valid) {
	  for (size_t i_state = 0; i_state < density_matrix[i_subspace].shape()[0]; ++i_state) {
	    i_pp = data.h_diag.flatten_block_index(i_subspace, i_state);
	    g_pp_tau[i_pp][i_block] += g_tmp[i_block] * data.imp_trace.get_density_matrix()[i_subspace].mat(i_state, i_state);
	  }
	}
      }
    }
  }

  void collect_results(triqs::mpi::communicator const& c) {
    z = mpi_all_reduce(z,c);
    for (auto& g_i: g_pp_tau) {
      for (auto& g_i_block: g_i) {
	g_i_block[0] = g_i_block[0] * 2;// Multiply first and last bins by 2 to account for full bins
	g_i_block[g_i_block.mesh().size() - 1] = g_i_block[g_i_block.mesh().size() - 1] * 2;
	g_i_block = mpi_all_reduce(g_i_block, c);
	g_i_block = g_i_block / (-real(z) * data.config.beta() * g_i_block.mesh().delta());
      }
    }
    for (size_t i_subspace = 0; i_subspace < density_matrix.size(); ++i_subspace) {
      for (size_t i_state = 0; i_state < density_matrix[i_subspace].shape()[0]; ++i_state) {
	int i_pp = data.h_diag.flatten_block_index(i_subspace, i_state);
	for (auto& g_block: g_pp_tau[i_pp]) {
	  g_block.singularity()(1) = density_matrix[i_subspace](i_state, i_state);
	}
      }
    }
  }
};
}
