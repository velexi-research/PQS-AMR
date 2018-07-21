/*! \file time_integration_3d.h
 *
 * \brief
 * Header files for numerical kernels for time integration.
 */

/*
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE.  This file is part of the XYZ package.  It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution.  No part of the XYZ
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

#ifndef INCLUDED_PQS_math_time_integration_3d_h
#define INCLUDED_PQS_math_time_integration_3d_h

#include "PQS/PQS_config.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief
 * @ref time_integration_3d.h provides support for time integration of
 * partial differential equations in three space dimensions. Support is
 * provided for
 *
 * - first-order Runge-Kutta method (i.e., forward Euler method) and
 *
 * - second-order and third-order total-variation diminishing (TVD)
 *   Runge-Kutta methods.
 *
 */


/* Link between C/C++ and Fortran function names
 *
 *      name in                             name in
 *      C/C++ code                          Fortran code
 *      ----------                          ------------
 */
#define PQS_MATH_3D_RK1_STEP                pqsmath3drk1step_
#define PQS_MATH_3D_TVD_RK2_STAGE1          pqsmath3dtvdrk2stage1_
#define PQS_MATH_3D_TVD_RK2_STAGE2          pqsmath3dtvdrk2stage2_
#define PQS_MATH_3D_TVD_RK3_STAGE1          pqsmath3dtvdrk3stage1_
#define PQS_MATH_3D_TVD_RK3_STAGE2          pqsmath3dtvdrk3stage2_
#define PQS_MATH_3D_TVD_RK3_STAGE3          pqsmath3dtvdrk3stage3_


/*!
 * Advance the solution 'u' through a single first-order Runge-Kutta
 * (i.e, forward Euler) step.
 *
 * Parameters
 * ----------
 * u_next: [output] u(t + dt)
 *
 * u_currentrent: u(t)
 *
 * rhs: right-hand side of time evolution equation
 *
 * dt: step size
 *
 * *_gb: index range for ghostbox
 *
 * *_fb: index range for fillbox
 *
 * Return value
 * ------------
 * None
 */
void PQS_MATH_3D_RK1_STEP(
  PQS_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const int *klo_u_next_gb,
  const int *khi_u_next_gb,
  const PQS_REAL *u_current,
  const int *ilo_u_current_gb,
  const int *ihi_u_current_gb,
  const int *jlo_u_current_gb,
  const int *jhi_u_current_gb,
  const int *klo_u_current_gb,
  const int *khi_u_current_gb,
  const PQS_REAL *rhs,
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const double *dt);

/*!
 * Advance the solution 'u' through the first stage of a second-order
 * TVD Runge-Kutta method step.
 *
 * Parameters
 * ----------
 * u_stage1: [output] u_approx(t + dt)
 *
 * u_current: u(t)
 *
 * rhs: right-hand side of time evolution equation
 *
 * dt: step size
 *
 * *_gb: index range for ghostbox
 *
 * *_fb: index range for fillbox
 *
 * Return value
 * ------------
 * None
 *
 * Notes
 * -----
 * - the first stage of TVD RK2 is identical to a single RK1 step
 */
void PQS_MATH_3D_TVD_RK2_STAGE1(
  PQS_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const int *klo_u_stage1_gb,
  const int *khi_u_stage1_gb,
  const PQS_REAL *u_current,
  const int *ilo_u_current_gb,
  const int *ihi_u_current_gb,
  const int *jlo_u_current_gb,
  const int *jhi_u_current_gb,
  const int *klo_u_current_gb,
  const int *khi_u_current_gb,
  const PQS_REAL *rhs,
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const double *dt);

/*!
 * Advance the solution 'u' through the second stage of a second-order
 * TVD Runge-Kutta method step.
 *
 * Parameters
 * ----------
 * u_next: [output]  u(t + dt)
 *
 * u_stage1: u_approx(t + dt)
 *
 * u_current: u(t)
 *
 * rhs: right-hand side of time evolution equation
 *
 * dt: step size
 *
 * *_gb: index range for ghostbox
 *
 * *_fb: index range for fillbox
 *
 * Return value
 * ------------
 * None
 */
void PQS_MATH_3D_TVD_RK2_STAGE2(
  PQS_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const int *klo_u_next_gb,
  const int *khi_u_next_gb,
  const PQS_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const int *klo_u_stage1_gb,
  const int *khi_u_stage1_gb,
  const PQS_REAL *u_current,
  const int *ilo_u_current_gb,
  const int *ihi_u_current_gb,
  const int *jlo_u_current_gb,
  const int *jhi_u_current_gb,
  const int *klo_u_current_gb,
  const int *khi_u_current_gb,
  const PQS_REAL *rhs,
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const double *dt);

/*!
 * Advance the solution 'u' through the first stage of a third-order
 * TVD Runge-Kutta method step.
 *
 * Parameters
 * ----------
 * u_stage1: [output] u_approx(t + dt)
 *
 * u_current: u(t)
 *
 * rhs: right-hand side of time evolution equation
 *
 * dt: step size
 *
 * *_gb: index range for ghostbox
 *
 * *_fb: index range for fillbox
 *
 * Return value
 * ------------
 * None
 *
 * Notes
 * -----
 * - the first stage of TVD RK3 is identical to a single RK1 step
 */
void PQS_MATH_3D_TVD_RK3_STAGE1(
  PQS_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const int *klo_u_stage1_gb,
  const int *khi_u_stage1_gb,
  const PQS_REAL *u_current,
  const int *ilo_u_current_gb,
  const int *ihi_u_current_gb,
  const int *jlo_u_current_gb,
  const int *jhi_u_current_gb,
  const int *klo_u_current_gb,
  const int *khi_u_current_gb,
  const PQS_REAL *rhs,
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const double *dt);

/*!
 * Advance the solution 'u' through the second stage of a third-order
 * TVD Runge-Kutta method step.
 *
 * Parameters
 * ----------
 * u_stage2: [output] u_approx(t + dt/2)
 *
 * u_stage1: u_approx(t + dt)
 *
 * u_current: u(t)
 *
 * rhs: right-hand side of time evolution equation
 *
 * dt: step size
 *
 * *_gb: index range for ghostbox
 *
 * *_fb: index range for fillbox
 *
 * Return value
 * ------------
 * None
 */
void PQS_MATH_3D_TVD_RK3_STAGE2(
  PQS_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const int *klo_u_stage2_gb,
  const int *khi_u_stage2_gb,
  const PQS_REAL *u_stage1,
  const int *ilo_u_stage1_gb,
  const int *ihi_u_stage1_gb,
  const int *jlo_u_stage1_gb,
  const int *jhi_u_stage1_gb,
  const int *klo_u_stage1_gb,
  const int *khi_u_stage1_gb,
  const PQS_REAL *u_current,
  const int *ilo_u_current_gb,
  const int *ihi_u_current_gb,
  const int *jlo_u_current_gb,
  const int *jhi_u_current_gb,
  const int *klo_u_current_gb,
  const int *khi_u_current_gb,
  const PQS_REAL *rhs,
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const double *dt);

/*!
 * Advance the solution 'u' through the third stage of a third-order
 * TVD Runge-Kutta method step.
 *
 * Parameters
 * ----------
 * u_next: [output] u(t + dt)
 *
 * u_stage2: u_approx(t + dt/2)
 *
 * u_current: u(t)
 *
 * rhs: right-hand side of time evolution equation
 *
 * dt: step size
 *
 * *_gb: index range for ghostbox
 *
 * *_fb: index range for fillbox
 *
 * Return value
 * ------------
 * None
 */
void PQS_MATH_3D_TVD_RK3_STAGE3(
  PQS_REAL *u_next,
  const int *ilo_u_next_gb,
  const int *ihi_u_next_gb,
  const int *jlo_u_next_gb,
  const int *jhi_u_next_gb,
  const int *klo_u_next_gb,
  const int *khi_u_next_gb,
  const PQS_REAL *u_stage2,
  const int *ilo_u_stage2_gb,
  const int *ihi_u_stage2_gb,
  const int *jlo_u_stage2_gb,
  const int *jhi_u_stage2_gb,
  const int *klo_u_stage2_gb,
  const int *khi_u_stage2_gb,
  const PQS_REAL *u_current,
  const int *ilo_u_current_gb,
  const int *ihi_u_current_gb,
  const int *jlo_u_current_gb,
  const int *jhi_u_current_gb,
  const int *klo_u_current_gb,
  const int *khi_u_current_gb,
  const PQS_REAL *rhs,
  const int *ilo_rhs_gb,
  const int *ihi_rhs_gb,
  const int *jlo_rhs_gb,
  const int *jhi_rhs_gb,
  const int *klo_rhs_gb,
  const int *khi_rhs_gb,
  const int *ilo_fb,
  const int *ihi_fb,
  const int *jlo_fb,
  const int *jhi_fb,
  const int *klo_fb,
  const int *khi_fb,
  const double *dt);

#ifdef __cplusplus
}
#endif

#endif
