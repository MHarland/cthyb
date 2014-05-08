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
#include <vector>
#include <triqs/gfs.hpp>

namespace cthyb {
    
class importance {

 gf<imtime> delta_norm; // Delta(\tau) / \int_0^\beta d\tau \Delta(\tau)
 const double beta;
 double imp;
 mutable double new_imp;
 int order;
 mutable int new_order;
    
public:
 
 importance(gf_const_view<imtime> delta) : 
  imp(1.0),
  order(0),
  beta(delta.domain().beta),
  delta_norm(delta)
 {
  // normalize delta_norm   
   auto shape = delta_norm.data().shape().front_pop();
   triqs::arrays::array<double,2> norm(shape);
   norm() = 0;

   auto dtau = delta_norm.mesh().delta();
   
   using mesh_point_t = typename gf<imtime>::mesh_t::mesh_point_t;
   foreach(delta_norm.mesh(), [&norm,this,dtau](mesh_point_t const& p){
    norm += delta_norm[p.index()]*dtau;
   });
   
   for(int n1=0; n1<shape[0];++n1)
   for(int n2=0; n2<shape[1];++n2){
    if(norm(n1,n2)==.0) continue;
    foreach(delta_norm.mesh(),[this,norm,n1,n2](mesh_point_t const& p){delta_norm[p.index()](n1,n2) /= norm(n1,n2);});
   }
 }
 
 double try_insert(std::pair<time_pt,int> const& x, std::pair<time_pt,int> const& y) const {
  // x <-> y
  double r = beta * delta_norm[closest_mesh_pt(double(y.first - x.first))](y.second, x.second);
  new_imp = std::pow(std::pow(imp,order) * r, 1.0/(order+1));
  new_order = order+1;
  return new_imp/imp;
 }
 
 double try_remove(std::pair<time_pt,int> const& x, std::pair<time_pt,int> const& y) const {
  // x <-> y
  double r = beta * delta_norm[closest_mesh_pt(double(y.first - x.first))](y.second, x.second);
  new_imp = std::pow(std::pow(imp,order) / r, order==1 ? .0 : 1.0/(order-1));
  new_order = order-1;
  return new_imp/imp;
 }
 
 void complete_operation() { imp = new_imp; order = new_order; }
 
 double operator()() const { return imp; }
};
    
}

