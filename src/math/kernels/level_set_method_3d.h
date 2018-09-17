/*! \file level_set_method_3d.h
 *
 * \brief
 * Header file for numerical kernels for Level Set Method computations in
 * three space dimensions
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

#ifndef INCLUDED_PQS_math_level_set_method_3d_h
#define INCLUDED_PQS_math_level_set_method_3d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                          name in
 *      C/C++ code                       Fortran code
 *      ----------                       ------------
 */
#define LSM_3D_VOLUME_PHI_LESS_THAN_ZERO \
                                         lsm3dvolumephilessthanzero_
#define LSM_3D_VOLUME_PHI_GREATER_THAN_ZERO \
                                         lsm3dvolumephigreaterthanzero_
#define LSM_3D_COMPUTE_REINIT_EQN_SGN_PHI0_RHS \
                                         lsm3dcomputereiniteqnsgnphi0rhs_
#define LSM_3D_COMPUTE_REINIT_EQN_SGN_PHI_RHS \
                                         lsm3dcomputereiniteqnsgnphirhs_

/*!
 * See documentation in level_set_method_3d.h
 */
PQS_REAL LSM_3D_VOLUME_PHI_LESS_THAN_ZERO(
    const PQS_REAL* phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx,
    const double *eps);

/*!
 * See documentation in level_set_method_3d.h
 */
PQS_REAL LSM_3D_VOLUME_PHI_GREATER_THAN_ZERO(
    const PQS_REAL* phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx,
    const double *eps);

/*!
 * See documentation in level_set_method_3d.h
 */
void LSM_3D_COMPUTE_REINIT_EQN_SGN_PHI0_RHS(
    const PQS_REAL *rhs,
    const int *rhs_gb_lo,
    const int *rhs_gb_hi,
    const PQS_REAL *phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const PQS_REAL *phi_0,
    const int *phi_0_gb_lo,
    const int *phi_0_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx);

/*!
 * See documentation in level_set_method_3d.h
 */
void LSM_3D_COMPUTE_REINIT_EQN_SGN_PHI_RHS(
    const PQS_REAL *rhs,
    const int *rhs_gb_lo,
    const int *rhs_gb_hi,
    const PQS_REAL *phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx);

#ifdef __cplusplus
}
#endif

#endif // INCLUDED_PQS_math_level_set_method_3d_h
