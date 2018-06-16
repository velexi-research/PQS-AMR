/*
 * spatial_derivatives3d.h
 *
 * Header file for Fortran 77 3D ENO/WENO routines.
 *
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE.  This file is part of the PQS package.  It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution.  No part of the PQS
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

#ifndef INCLUDED_PQS_SPATIAL_DERIVATIVES_3D_H
#define INCLUDED_PQS_SPATIAL_DERIVATIVES_3D_H

#include "PQS_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \file pqs_spatial_derivatives3d.h
 *
 * \brief
 * @ref pqs_spatial_derivatives3d.h provides support for computing spatial
 * derivatives in three space dimensions using high-order ENO and WENO
 * discretizations.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                      name in
 *      C/C++ code                   Fortran code
 *      ----------                   ------------
 */
#define PQS3D_HJ_ENO1                pqs3dhjeno1_
#define PQS3D_HJ_ENO2                pqs3dhjeno2_
#define PQS3D_HJ_ENO3                pqs3dhjeno3_
#define PQS3D_HJ_WENO5               pqs3dhjweno5_
#define PQS3D_UPWIND_HJ_ENO1         pqs3dupwindhjeno1_
#define PQS3D_UPWIND_HJ_ENO2         pqs3dupwindhjeno2_
#define PQS3D_UPWIND_HJ_ENO3         pqs3dupwindhjeno3_
#define PQS3D_UPWIND_HJ_WENO5        pqs3dupwindhjweno5_
#define PQS3D_CENTRAL_GRAD_ORDER2    pqs3dcentralgradorder2_
#define PQS3D_CENTRAL_GRAD_ORDER4    pqs3dcentralgradorder4_
#define PQS3D_LAPLACIAN_ORDER2       pqs3dlaplacianorder2_
#define PQS3D_PHI_UPWIND_GRAD_F      pqs3dphiupwindgradf_
#define PQS3D_AVERAGE_GRAD_PHI       pqs3daveragegradphi_
#define PQS3D_GRADIENT_MAGNITUDE     pqs3dgradientmagnitude_
#define PQS3D_CENTRAL_HESSIAN        pqs3dcentralhessian_


/*!
 * PQS3D_HJ_ENO1() computes the forward (plus) and backward (minus)
 * first-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - dx, dy, dz (in):    grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *      
 * Return value:          none
 *
 * NOTES:
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *  - phi_x_plus and phi_x_minus are face-centered data (i.e. their
 *    indices are of the form (i+1/2) and (i-1/2)).  For phi_x_plus,
 *    the data array position corresponding to the (i+1/2) is i (i.e. 
 *    the index shifted down to the nearest integer index).  For 
 *    phi_x_minus, the data array position corresponding to the (i-1/2) 
 *    is i (i.e.  the index shifted up to the nearest integer index).  
 *    Analogous conventions hold for phi_y_plus, phi_y_minus, phi_z_plus,
 *    and phi_z_minus.
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void PQS3D_HJ_ENO1(
  PQS_REAL *phi_x_plus,
  PQS_REAL *phi_y_plus,
  PQS_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  PQS_REAL *phi_x_minus,
  PQS_REAL *phi_y_minus,
  PQS_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  PQS_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*!
 * PQS3D_HJ_ENO2() computes the forward (plus) and backward (minus)
 * second-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - D2 (in):            scratch space for holding undivided second-differences
 *  - dx, dy, dz (in):    grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *      
 * Return value:          none
 *
 * NOTES:
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *  - phi_x_plus and phi_x_minus are face-centered data (i.e. their
 *    indices are of the form (i+1/2) and (i-1/2)).  For phi_x_plus,
 *    the data array position corresponding to the (i+1/2) is i (i.e. 
 *    the index shifted down to the nearest integer index).  For 
 *    phi_x_minus, the data array position corresponding to the (i-1/2) 
 *    is i (i.e.  the index shifted up to the nearest integer index).  
 *    Analogous conventions hold for phi_y_plus, phi_y_minus, phi_z_plus,
 *    and phi_z_minus.
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void PQS3D_HJ_ENO2(
  PQS_REAL *phi_x_plus,
  PQS_REAL *phi_y_plus,
  PQS_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  PQS_REAL *phi_x_minus,
  PQS_REAL *phi_y_minus,
  PQS_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  PQS_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  PQS_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const int *klo_D2_gb,
  const int *khi_D2_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*!
 * PQS3D_HJ_ENO3() computes the forward (plus) and backward (minus)
 * third-order Hamilton-Jacobi ENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - D2 (in):            scratch space for holding undivided second-differences
 *  - D3 (in):            scratch space for holding undivided third-differences
 *  - dx, dy, dz (in):    grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *      
 * Return value:          none
 *
 * NOTES:
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *  - phi_x_plus and phi_x_minus are face-centered data (i.e. their
 *    indices are of the form (i+1/2) and (i-1/2)).  For phi_x_plus,
 *    the data array position corresponding to the (i+1/2) is i (i.e. 
 *    the index shifted down to the nearest integer index).  For 
 *    phi_x_minus, the data array position corresponding to the (i-1/2) 
 *    is i (i.e.  the index shifted up to the nearest integer index).  
 *    Analogous conventions hold for phi_y_plus, phi_y_minus, phi_z_plus,
 *    and phi_z_minus.
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void PQS3D_HJ_ENO3(
  PQS_REAL *phi_x_plus,
  PQS_REAL *phi_y_plus,
  PQS_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  PQS_REAL *phi_x_minus,
  PQS_REAL *phi_y_minus,
  PQS_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  PQS_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  PQS_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const int *klo_D2_gb,
  const int *khi_D2_gb,
  PQS_REAL *D3,
  const int *ilo_D3_gb,
  const int *ihi_D3_gb,
  const int *jlo_D3_gb,
  const int *jhi_D3_gb,
  const int *klo_D3_gb,
  const int *khi_D3_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*!
 * PQS3D_HJ_WENO5() computes the forward (plus) and backward (minus)
 * fifth-order Hamilton-Jacobi WENO approximations to the gradient of
 * \f$ \phi \f$.
 *           
 * Arguments:
 *  - phi_*_plus (out):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (out):  components of \f$ \nabla \phi \f$ in minus direction
 *  - phi (in):           \f$ \phi \f$
 *  - D1 (in):            scratch space for holding undivided first-differences
 *  - dx, dy, dz (in):    grid spacing
 *  - *_gb (in):          index range for ghostbox
 *  - *_fb (in):          index range for fillbox
 *      
 * Return value:          none
 *
 * NOTES:
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *  - phi_x_plus and phi_x_minus are face-centered data (i.e. their
 *    indices are of the form (i+1/2) and (i-1/2)).  For phi_x_plus,
 *    the data array position corresponding to the (i+1/2) is i (i.e. 
 *    the index shifted down to the nearest integer index).  For 
 *    phi_x_minus, the data array position corresponding to the (i-1/2) 
 *    is i (i.e.  the index shifted up to the nearest integer index).  
 *    Analogous conventions hold for phi_y_plus, phi_y_minus, phi_z_plus,
 *    and phi_z_minus.
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void PQS3D_HJ_WENO5(
  PQS_REAL *phi_x_plus,
  PQS_REAL *phi_y_plus,
  PQS_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  PQS_REAL *phi_x_minus,
  PQS_REAL *phi_y_minus,
  PQS_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  PQS_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*! 
 * PQS3D_UPWIND_HJ_ENO1() computes the first-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_* (out):      components of \f$ \nabla \phi \f$
 *  - phi (in):         \f$ \phi \f$
 *  - vel_* (in):       components of the velocity 
 *  - D1 (in):          scratch space for holding undivided first-differences
 *  - dx, dy, dz (in):  grid cell size
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void PQS3D_UPWIND_HJ_ENO1(
  PQS_REAL *phi_x,
  PQS_REAL *phi_y,
  PQS_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const PQS_REAL *vel_x,
  const PQS_REAL *vel_y,
  const PQS_REAL *vel_z,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  const int *klo_vel_gb,
  const int *khi_vel_gb,
  PQS_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*! 
 * PQS3D_UPWIND_HJ_ENO2() computes the second-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_* (out):      components of \f$ \nabla \phi \f$
 *  - phi (in):         \f$ \phi \f$
 *  - vel_* (in):       components of the velocity 
 *  - D1 (in):          scratch space for holding undivided first-differences
 *  - D2 (in):          scratch space for holding undivided second-differences
 *  - dx, dy, dz (in):  grid cell size
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void PQS3D_UPWIND_HJ_ENO2(
  PQS_REAL *phi_x,
  PQS_REAL *phi_y,
  PQS_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const PQS_REAL *vel_x,
  const PQS_REAL *vel_y,
  const PQS_REAL *vel_z,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  const int *klo_vel_gb,
  const int *khi_vel_gb,
  PQS_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  PQS_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const int *klo_D2_gb,
  const int *khi_D2_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*! 
 * PQS3D_UPWIND_HJ_ENO3() computes the third-order Hamilton-Jacobi ENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_* (out):      components of \f$ \nabla \phi \f$
 *  - phi (in):         \f$ \phi \f$
 *  - vel_* (in):       components of the velocity 
 *  - D1 (in):          scratch space for holding undivided first-differences
 *  - D2 (in):          scratch space for holding undivided second-differences
 *  - D3 (in):          scratch space for holding undivided third-differences
 *  - dx, dy, dz (in):  grid cell size
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void PQS3D_UPWIND_HJ_ENO3(
  PQS_REAL *phi_x,
  PQS_REAL *phi_y,
  PQS_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const PQS_REAL *vel_x,
  const PQS_REAL *vel_y,
  const PQS_REAL *vel_z,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  const int *klo_vel_gb,
  const int *khi_vel_gb,
  PQS_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  PQS_REAL *D2,
  const int *ilo_D2_gb,
  const int *ihi_D2_gb,
  const int *jlo_D2_gb,
  const int *jhi_D2_gb,
  const int *klo_D2_gb,
  const int *khi_D2_gb,
  PQS_REAL *D3,
  const int *ilo_D3_gb,
  const int *ihi_D3_gb,
  const int *jlo_D3_gb,
  const int *jhi_D3_gb,
  const int *klo_D3_gb,
  const int *khi_D3_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*! 
 * PQS3D_UPWIND_HJ_WENO5() computes the fifth-order Hamilton-Jacobi WENO
 * upwind approximation to the gradient of \f$ \phi \f$.
 * 
 * Arguments:
 *  - phi_* (out):      components of \f$ \nabla \phi \f$
 *  - phi (in):         \f$ \phi \f$
 *  - vel_* (in):       components of the velocity 
 *  - D1 (in):          scratch space for holding undivided first-differences
 *  - dx, dy, dz (in):  grid cell size
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 *
 * NOTES:
 *  - the fillbox is defined in terms of the index range for 
 *    cell-centered data
 */
void PQS3D_UPWIND_HJ_WENO5(
  PQS_REAL *phi_x,
  PQS_REAL *phi_y,
  PQS_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const PQS_REAL *vel_x,
  const PQS_REAL *vel_y,
  const PQS_REAL *vel_z,
  const int *ilo_vel_gb,
  const int *ihi_vel_gb,
  const int *jlo_vel_gb,
  const int *jhi_vel_gb,
  const int *klo_vel_gb,
  const int *khi_vel_gb,
  PQS_REAL *D1,
  const int *ilo_D1_gb,
  const int *ihi_D1_gb,
  const int *jlo_D1_gb,
  const int *jhi_D1_gb,
  const int *klo_D1_gb,
  const int *khi_D1_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*! 
 * PQS3D_CENTRAL_GRAD_ORDER2() computes the second-order, central, 
 * finite difference approximation to the gradient of \f$ \phi \f$ 
 * using the formula:
 * 
 *    \f[
 *
 *      \left( \frac{\partial \phi}{\partial x} \right)_i \approx
 *        \frac{ \phi_{i+1} - \phi_{i-1} }{ 2 dx }
 *
 *    \f]
 *
 * Arguments:
 *  - phi_* (out):      components of \f$ \nabla \phi \f$
 *  - phi (in):         \f$ \phi \f$
 *  - dx, dy, dz (in):  grid cell size
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 */
void PQS3D_CENTRAL_GRAD_ORDER2( 
  PQS_REAL *phi_x,
  PQS_REAL *phi_y,
  PQS_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*! 
 * PQS3D_CENTRAL_GRAD_ORDER4() computes the fourth-order, central,
 * finite difference approximation to the gradient of \f$ \phi \f$ 
 * using the formula:
 *
 *    \f[
 *
 *      \left( \frac{\partial \phi}{\partial x} \right)_i \approx
 *         \frac{ -\phi_{i+2} + 8 \phi_{i+1} - 8 \phi_{i-1} + \phi_{i-2} }
 *              { 12 dx }
 *
 *    \f]
 *
 * Arguments:
 *  - phi_* (out):      components of \f$ \nabla \phi \f$
 *  - phi (in):         \f$ \phi \f$
 *  - dx, dy, dz (in):  grid cell size
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 */
void PQS3D_CENTRAL_GRAD_ORDER4( 
  PQS_REAL *phi_x,
  PQS_REAL *phi_y,
  PQS_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*! 
 * PQS3D_LAPLACIAN_ORDER2() computes the second-order, central, 
 * finite difference approximation to the Laplacian of \f$ \phi \f$ 
 * using the formula:
 *
 *    \f[
 *
 *      \nabla^2 \phi \approx
 *         \frac{ \phi_{i+1,j,k} - 2 \phi_{i,j,k} + \phi_{i-1,j,k} }
 *              { dx^2 }
 *       + \frac{ \phi_{i,j+1,k} - 2 \phi_{i,j,k} + \phi_{i,j-1,k} }
 *              { dy^2 }
 *       + \frac{ \phi_{i,j,k+1} - 2 \phi_{i,j,k} + \phi_{i,j,k-1} }
 *              { dz^2 }
 *
 *    \f]
 *
 * Arguments:
 *  - laplacian_phi (out):  Laplacian of \f$ phi \f$
 *  - phi (in):             \f$ \phi \f$
 *  - dx (in):              grid cell size
 *  - *_gb (in):            index range for ghostbox
 *  - *_fb (in):            index range for fillbox
 *
 * Return value:            none
 */
void PQS3D_LAPLACIAN_ORDER2( 
  PQS_REAL *laplacian_phi,
  const int *ilo_laplacian_phi_gb,
  const int *ihi_laplacian_phi_gb,
  const int *jlo_laplacian_phi_gb,
  const int *jhi_laplacian_phi_gb,
  const int *klo_laplacian_phi_gb,
  const int *khi_laplacian_phi_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);


/*!
 * PQS3D_PHI_UPWIND_GRAD_F() computes the \f$ \phi \f$-upwind gradient of a
 * function, F, using the following "upwinding" scheme to compute
 * the normal:
 *
 *   if \f$ \phi > 0 \f$:  upwind direction is direction where 
 *                         \f$ \phi \f$ is smaller
 *
 *   if \f$ \phi < 0 \f$:  upwind direction is direction where 
 *                         \f$ \phi \f$ is larger
 *
 * Arguments:
 *  - F_* (out):        components of \f$ \phi \f$-upwinded \f$ \nabla F \f$
 *  - F_*_plus (in):    components of \f$ \nabla F \f$ in plus direction
 *  - F_*_minus (in):   components of \f$ \nabla F \f$ in minus direction
 *  - phi (in):         level set function
 *  - *_gb (in):        index range for ghostbox
 *  - *_fb (in):        index range for fillbox
 *
 * Return value:        none
 *
 * NOTES:
 *  - \f$ \phi \f$ is REQUIRED to have at least one ghost cell in each
 *    coordinate direction for upwinding
 *  - it is assumed that BOTH the plus AND minus derivatives have
 *    the same fillbox
 *
 */
void PQS3D_PHI_UPWIND_GRAD_F(
  PQS_REAL *F_x,
  PQS_REAL *F_y,
  PQS_REAL *F_z,
  const int *ilo_grad_F_gb,
  const int *ihi_grad_F_gb,
  const int *jlo_grad_F_gb,
  const int *jhi_grad_F_gb,
  const int *klo_grad_F_gb,
  const int *khi_grad_F_gb,
  PQS_REAL *F_x_plus,
  PQS_REAL *F_y_plus,
  PQS_REAL *F_z_plus,
  const int *ilo_grad_F_plus_gb,
  const int *ihi_grad_F_plus_gb,
  const int *jlo_grad_F_plus_gb,
  const int *jhi_grad_F_plus_gb,
  const int *klo_grad_F_plus_gb,
  const int *khi_grad_F_plus_gb,
  PQS_REAL *F_x_minus,
  PQS_REAL *F_y_minus,
  PQS_REAL *F_z_minus,
  const int *ilo_grad_F_minus_gb,
  const int *ihi_grad_F_minus_gb,
  const int *jlo_grad_F_minus_gb,
  const int *jhi_grad_F_minus_gb,
  const int *klo_grad_F_minus_gb,
  const int *khi_grad_F_minus_gb,
  const PQS_REAL *phi,
  const int *ilo_phi_gb,
  const int *ihi_phi_gb,
  const int *jlo_phi_gb,
  const int *jhi_phi_gb,
  const int *klo_phi_gb,
  const int *khi_phi_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb);


/*!
 * PQS3D_AVERAGE_GRAD_PHI() computes the average of the plus and minus
 * derivatives:
 *
 * \f[
 *
 *   \nabla \phi = (\nabla \phi_{plus} + \nabla \phi_{minus}) / 2
 *
 * \f]
 *
 * Arguments:
 *  - phi_* (out):       components of average \f$ \nabla \phi \f$
 *  - phi_*_plus (in):   components of \f$ \nabla \phi \f$ in plus direction
 *  - phi_*_minus (in):  components of \f$ \nabla \phi \f$ in minus direction
 *  - *_gb (in):         index range for ghostbox
 *  - *_fb (in):         index range for fillbox
 *
 * Return value:         none
 */
void PQS3D_AVERAGE_GRAD_PHI(
  PQS_REAL *phi_x,
  PQS_REAL *phi_y,
  PQS_REAL *phi_z,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  PQS_REAL *phi_x_plus,
  PQS_REAL *phi_y_plus,
  PQS_REAL *phi_z_plus,
  const int *ilo_grad_phi_plus_gb,
  const int *ihi_grad_phi_plus_gb,
  const int *jlo_grad_phi_plus_gb,
  const int *jhi_grad_phi_plus_gb,
  const int *klo_grad_phi_plus_gb,
  const int *khi_grad_phi_plus_gb,
  PQS_REAL *phi_x_minus,
  PQS_REAL *phi_y_minus,
  PQS_REAL *phi_z_minus,
  const int *ilo_grad_phi_minus_gb,
  const int *ihi_grad_phi_minus_gb,
  const int *jlo_grad_phi_minus_gb,
  const int *jhi_grad_phi_minus_gb,
  const int *klo_grad_phi_minus_gb,
  const int *khi_grad_phi_minus_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb);
  
/*!
 *
 *  PQS3D_GRADIENT_MAGNITUDE() computes magnitude of the gradient.
 *
 *  Arguments:
 *    phi_* (in):           components of grad(phi) 
 *    grad_phi_mag (out):   gradient magnitude
 *    *_gb (in):           index range for ghostbox
 *    *_fb (in):           index range for fillbox
 *
 */
void PQS3D_GRADIENT_MAGNITUDE(
  const PQS_REAL *phi_x,
  const PQS_REAL *phi_y,
  const PQS_REAL *phi_z,
  const PQS_REAL *grad_phi_mag,
  const int *ilo_grad_phi_gb,
  const int *ihi_grad_phi_gb,
  const int *jlo_grad_phi_gb,
  const int *jhi_grad_phi_gb,
  const int *klo_grad_phi_gb,
  const int *khi_grad_phi_gb,
  const int *ilo_grad_phi_fb,
  const int *ihi_grad_phi_fb,
  const int *jlo_grad_phi_fb,
  const int *jhi_grad_phi_fb,
  const int *klo_grad_phi_fb,
  const int *khi_grad_phi_fb);
  

/*! 
 * PQS3D_CENTRAL_HESSIAN() computes second-order, central, finite difference
 * approximations to the elements of the Hessian matrix of \f$ \phi \f$.
 *
 * Arguments:
 *   phi_* (out):  derivatives of phi
 *   phi (in):     function to compute derivatives for
 *   dx, dy, dz (in):  grid spacing
 *   *_gb (in):    index range for ghostbox
 *   *_fb (in):    index range for fillbox
 *
 */
void  PQS3D_CENTRAL_HESSIAN(
  PQS_REAL *phi_xx,
  PQS_REAL *phi_xy,
  PQS_REAL *phi_xz,
  PQS_REAL *phi_yy,
  PQS_REAL *phi_yz,
  PQS_REAL *phi_zz,
  const int *ilo_hessian_gb, 
  const int *ihi_hessian_gb,
  const int *jlo_hessian_gb,
  const int *jhi_hessian_gb,
  const int *klo_hessian_gb,
  const int *khi_hessian_gb,
  const PQS_REAL *phi,
  const int *ilo_gb, 
  const int *ihi_gb,
  const int *jlo_gb,
  const int *jhi_gb,
  const int *klo_gb,
  const int *khi_gb,
  const int *ilo_fb, 
  const int *ihi_fb,
  const int *jlo_fb,  
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const PQS_REAL *dx,
  const PQS_REAL *dy,
  const PQS_REAL *dz);

#ifdef __cplusplus
}
#endif

#endif
