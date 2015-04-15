/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#include <algorithm>
#include "./qmc_data.hpp"

namespace cthyb {

// Removal of C, C^dagger operator
class move_remove_c_c_cdag_cdag {

 qmc_data& data;
 configuration& config;
 mc_tools::random_generator& rng;
 int block_index1, block_index2, block_size1, block_size2;
 qmc_data::trace_t new_trace;
 time_pt tau1, tau2, tau3, tau4;

 public:
 //----------------------------------

 move_remove_c_c_cdag_cdag(int block_index1, int block_index2, int block_size1, int block_size2, qmc_data& data, mc_tools::random_generator& rng)
    : data(data),
      config(data.config),
      rng(rng),
      block_index1(block_index1),
      block_size1(block_size1),
      block_index2(block_index2),
      block_size2(block_size2) {}

 //----------------
 
 mc_weight_type attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "In config " << config.id << std::endl;
  std::cerr << "* Attempt for move_remove_c_c_cdag_cdag (block " << block_index1 << ", " << block_index2 << ")" << std::endl;
//  std::cerr << "* Configuration before:" << std::endl << config;
//  data.imp_trace.tree.graphviz(std::ofstream("tree_before"));
#endif

  auto& det1 = data.dets[block_index1];
  auto& det2 = data.dets[block_index2];
  double det_ratio; //FIXME

  // Pick two pairs of C, Cdagger to remove at random
  // Remove the operators from the traces
  int det1_size = det1.size();
  int det2_size = det2.size();
  if (block_index1 == block_index2) {
   if (det1_size < 2) return 0; // det does not contain two operators
  } else {
   if ((det1_size == 0) || (det2_size == 0)) return 0; // one of two dets is empty
  }
  int num_c_dag1 = rng(det1_size), num_c1 = rng(det1_size);
  int num_c_dag2 = rng(det2_size), num_c2 = rng(det2_size);
  if ((block_index1 == block_index2) && ((num_c_dag1 == num_c_dag2) || (num_c1 == num_c2))) return 0; // picked the same operator twice
#ifdef EXT_DEBUG
  std::cerr << "* Proposing to remove: ";
  std::cerr << num_c_dag1 << "-th Cdag(" << block_index1 << ",...), ";
  std::cerr << num_c1 << "-th C(" << block_index1 << ",...)" << std::endl;
  std::cerr << " and ";
  std::cerr << num_c_dag2 << "-th Cdag(" << block_index2 << ",...), ";
  std::cerr << num_c2 << "-th C(" << block_index2 << ",...)" << std::endl;
#endif

  // now mark 2 nodes for deletion
  tau1 = data.imp_trace.try_delete(num_c1, block_index1, false);
  tau2 = data.imp_trace.try_delete(num_c_dag1, block_index1, true);
  tau3 = data.imp_trace.try_delete(num_c2, block_index2, false);
  tau4 = data.imp_trace.try_delete(num_c_dag2, block_index2, true);

  if (block_index1 == block_index2) {
   det_ratio = det1.try_remove2(num_c_dag1, num_c_dag2, num_c1, num_c2);
  } else { // block_index1 != block_index2
   auto det_ratio1 = det1.try_remove(num_c_dag1, num_c1);
   auto det_ratio2 = det2.try_remove(num_c_dag2, num_c2);
   det_ratio = det_ratio1 * det_ratio2;
  }

  // proposition probability
  double t_ratio;
  // Note: Must use the size of the det before the try_delete!
  if (block_index1 == block_index2) {
   // Here, we use the fact that the two cdag/c proposed to be removed in the det can be at the same 
   // positions in the det, and thus remove prob is NOT (detsize+2)*(detsize+1)
   t_ratio = std::pow(block_size1 * config.beta() / double(det1_size), 4);
  } else {
   t_ratio = std::pow(block_size1 * config.beta() / double(det1_size), 2) * std::pow(block_size2 * config.beta() / double(det2_size), 2);
  }
  
  // For quick abandon 
  double random_number = rng.preview();
  if (random_number == 0.0) return 0;
  double p_yee = std::abs(det_ratio / t_ratio / data.trace);

  // recompute the trace
  new_trace = data.imp_trace.estimate(p_yee, random_number);
  if (new_trace == 0.0) {
#ifdef EXT_DEBUG
   std::cout << "trace == 0" << std::endl;
#endif
   return 0;
  }
  auto trace_ratio = new_trace / data.trace;
  if (!std::isfinite(trace_ratio)) TRIQS_RUNTIME_ERROR << "trace_ratio not finite " << new_trace << " " << data.trace << " " << new_trace/data.trace << " in config " << config.id;

  mc_weight_type p = trace_ratio * det_ratio;

#ifdef EXT_DEBUG
  std::cerr << "Trace ratio: " << trace_ratio << '\t';
  std::cerr << "Det ratio: " << det_ratio << '\t';
  std::cerr << "Prefactor: " << t_ratio << '\t';
  std::cerr << "Weight: " << p/ t_ratio << std::endl;
#endif

  if (!std::isfinite(p)) TRIQS_RUNTIME_ERROR << "(remove) p not finite :" << p << " in config " << config.id;
  if (!std::isfinite(p / t_ratio)) TRIQS_RUNTIME_ERROR << "p / t_ratio not finite p : " << p << " t_ratio :  " << t_ratio << " in config " << config.id;
  return p / t_ratio;
 }

 //----------------

 mc_weight_type accept() {

  config.id++; // increment the config id

  // remove from the tree
  data.imp_trace.confirm_delete();

  // remove from the configuration
  config.erase(tau1);
  config.erase(tau2);
  config.erase(tau3);
  config.erase(tau4);
  
  // remove from the determinants
  if (block_index1 == block_index2) {
   data.dets[block_index1].complete_operation();
  } else {
   data.dets[block_index1].complete_operation();
   data.dets[block_index2].complete_operation();
  }
  data.update_sign();
  data.trace = new_trace;

#ifdef EXT_DEBUG
//  std::cerr << "* Configuration after: " << config.id << std::endl;
//  std::cerr << config;
  check_det_sequence(data.dets[block_index1],config.id);
  check_det_sequence(data.dets[block_index2],config.id);
#endif
#ifdef PRINT_CONF_DEBUG
  config.print_to_h5();
#endif

  return data.current_sign / data.old_sign;
 }

 //----------------

 void reject() {

  config.id++; // increment the config id

  data.imp_trace.cancel_delete();

#ifdef EXT_DEBUG
//  std::cerr << "* Configuration after: " << config.id << std::endl;
//  std::cerr << config;
  check_det_sequence(data.dets[block_index1],config.id);
  check_det_sequence(data.dets[block_index2],config.id);
#endif
#ifdef PRINT_CONF_DEBUG
  config.print_to_h5();
#endif

 }
};
}