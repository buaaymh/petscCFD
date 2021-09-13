/// @file euler.h
///
/// Definition of the class related to linear models.
///
//*****************************************************************************
//
//  Copyright Minghao Yang, CFD Solver project
//  Created August 15, 2021
//  Last modification: August 15, 2021
//
//=============================================================================

#ifndef INCLUDE_EULER_H_
#define INCLUDE_EULER_H_

#include "defs.h"

namespace cfd {

struct Euler
{
  static constexpr Equations type = Equations::Euler;
  static constexpr int nEqual = 4;
  using Flux = Matrix<Real, nEqual, 1>;
  using State = Flux;

  Euler() = default;

  static unordered_map<string, int> CreateFieldDiscription() {
    unordered_map<string, int>  field_desc;
    field_desc.emplace("Pressure", 1);
    field_desc.emplace("Y_Velocity", 1);
    field_desc.emplace("X_Velocity", 1);
    field_desc.emplace("Density", 1);
    return field_desc;
  }
  static void PrimToCons(const Real* primitive, Real* conserved);
  static void ConsToPrim(const Real* conserved, Real* primitive);
  static Flux GetFlux(const Real* normal, const Real* pv_l, const Real* pv_r);
  static Flux GetFlux(const Real* normal, const Real* pv);
  static Flux GetWallFlux(const Real* normal, const Real* pv);
};

void Euler::PrimToCons(const Real* primitive, Real* conserved)
{
  conserved[0] = primitive[0];
  // V^2 = u^2 + v^2
  Real V2  = primitive[1] * primitive[1] + primitive[2] * primitive[2];
  // xmom = rho * uvel
  conserved[1] = primitive[1] * primitive[0];
  // ymom = rho * vvel
  conserved[2] = primitive[2] * primitive[0];
  // energy = p/(gamma - 1) + 0.5*rho*V^2
  conserved[3] = primitive[3] * One_GamM1;
  conserved[3] += 0.5 * primitive[0] * V2;
}

void Euler::ConsToPrim(const Real* conserved, Real* primitive)
{
  primitive[0] = conserved[0];
  // uvel = xmom / rho
  primitive[1] = conserved[1] / conserved[0];
  // vvel = ymom / rho
  primitive[2] = conserved[2] / conserved[0];
  // energy = p/(gamma - 1) + 0.5*rho*V^2
  Real V2  = primitive[1]*primitive[1] + primitive[2]*primitive[2];
  primitive[3] = conserved[3] - 0.5 * conserved[0] * V2;
  primitive[3] *= GamM1;
}

Euler::Flux Euler::GetFlux(const Real* normal, const Real* pv)
{
  const Real& r  = pv[0];
  const Real& u  = pv[1];
  const Real& v  = pv[2];
  const Real& p  = pv[3];
  Real        h   = Gam_GamM1*p/r + 0.5*(u*u+v*v);

  Flux flux;
  Real qs = u*normal[0] + v*normal[1];
  flux(0) = qs*r;
  flux(1) = flux(0)*u + p*normal[0];
  flux(2) = qs*r*v    + p*normal[1];
  flux(3) = qs*r*h;
  return flux;
}

Euler::Flux Euler::GetFlux(const Real* normal, const Real* pv_l, const Real* pv_r)
{
  // left & right state
  const Real& rl   = pv_l[0];
  const Real& ul   = pv_l[1];
  const Real& vl   = pv_l[2];
  const Real& pl   = pv_l[3];
  Real        hl   = Gam_GamM1*pl/rl + 0.5*(ul*ul+vl*vl);

  const Real& rr   = pv_r[0];
  const Real& ur   = pv_r[1];
  const Real& vr   = pv_r[2];
  const Real& pr   = pv_r[3];
  Real        hr   = Gam_GamM1*pr/rr + 0.5*(ur*ur+vr*vr);

  if (rl < 0 || rr < 0) { cout << "Rho_l : " << rl << " Rho_r : " << rr << endl; }
  if (pl < 0 || pr < 0) { cout << "P_l : "   << pl << " P_r : "   << pr << endl; }
  
  // left and right V*n
  Real qsl = ul*normal[0] + vl*normal[1];
  Real qsr = ur*normal[0] + vr*normal[1];

  Flux flux;
  Real pav = 0.5 * (pl + pr);
  flux(0)  = 0.5 * (qsl*rl    + qsr*rr);
  flux(1)  = 0.5 * (qsl*rl*ul + qsr*rr*ur) + pav*normal[0];
  flux(2)  = 0.5 * (qsl*rl*vl + qsr*rr*vr) + pav*normal[1];
  flux(3)  = 0.5 * (qsl*rl*hl + qsr*rr*hr);

  // Roe's average
  Real rav  = Sqrt(rl*rr);
  Real dd   = rav/rl;
  Real dd1  = 1.0/(1.0+dd);
  Real uav  = (ul+dd*ur)*dd1;
  Real vav  = (vl+dd*vr)*dd1;
  Real hav  = (hl+dd*hr)*dd1;
  Real q2a  = 0.5*(uav*uav+vav*vav);
  Real c2a  = GamM1*(hav-q2a);
  Real cav  = Sqrt(c2a);
  Real uv   =     uav*normal[0] +     vav*normal[1];
  Real du   = (ur-ul)*normal[0] + (vr-vl)*normal[1];

  // eigenvalues
  Real h1    = Abs(uv - cav);
  Real h2    = Abs(uv);
  Real h4    = Abs(uv + cav);
  Real delta = EPS_Entr * cav;

  Real eabs1 = EntropyCorr( h1,delta );
  Real eabs2 = EntropyCorr( h2,delta );
  Real eabs4 = EntropyCorr( h4,delta );

  // upwind fluxes

  h1      = rav*cav*du;
  h2      = eabs1*(pr-pl - h1)/(2.0*c2a);
  Real h3 = eabs2*(rr-rl - (pr-pl)/c2a);
  h4      = eabs2*rav;
  Real h5 = eabs4*(pr-pl + h1)/(2.0*c2a);

  Flux dissp;
  dissp(0) = h2                     + h3                                          + h5;
  dissp(1) = h2*(uav-cav*normal[0]) + h3*uav + h4*(ur-ul-du*normal[0])            + h5*(uav+cav*normal[0]);
  dissp(2) = h2*(vav-cav*normal[1]) + h3*vav + h4*(vr-vl-du*normal[1])            + h5*(vav+cav*normal[1]);
  dissp(3) = h2*(hav-cav*uv)        + h3*q2a + h4*(uav*(ur-ul)+vav*(vr-vl)-uv*du) + h5*(hav+cav*uv);

  return flux - 0.5*dissp;
}

Euler::Flux Euler::GetWallFlux(const Real* normal, const Real* pv)
{
  Flux flux = Flux::Zero();
  flux(1) += pv[3] * normal[0];
  flux(2) += pv[3] * normal[1];
  return flux;
}

}  // cfd

#endif // INCLUDE_EULER_H_