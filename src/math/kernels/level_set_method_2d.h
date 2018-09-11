/*! \file level_set_method_2d.h
 *
 * \brief
 * Header file for numerical kernels for Level Set Method computations in
 * two space dimensions
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

#ifndef INCLUDED_PQS_math_level_set_method_2d_h
#define INCLUDED_PQS_math_level_set_method_2d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                          name in
 *      C/C++ code                       Fortran code
 *      ----------                       ------------
 */
#define LSM_2D_AREA_PHI_LESS_THAN_ZERO \
                                         lsm2dareaphilessthanzero_
#define LSM_2D_AREA_PHI_GREATER_THAN_ZERO \
                                         lsm2dareaphigreaterthanzero_

/*!
 *
 * Compute area of region where level set function is less than 0.
 *
 * Parameters
 * ----------
 * phi: level set function
 *
 * dx: array containing grid spacing in each coordinate direction
 *
 * eps: width of numerical smoothing to use for Heaviside function
 *
 * *_gb_lo: lower corner of index range for ghostbox
 *
 * *_gb_hi: upper corner of index range for ghostbox
 *
 * patch_box_lo: lower corner of index range for patch box
 *
 * patch_box_hi: upper corner of index range for patch box
 *
 * Return value
 * ------------
 * area of region where phi < 0
 *
 * Notes
 * -----
 * - When phi is a signed distance function, the numerical width of the
 *   Heaviside function as a function of spatial position is approximately
 *   equal to (2 * eps).
 */
PQS_REAL LSM_2D_AREA_PHI_LESS_THAN_ZERO(
    const PQS_REAL* phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx,
    const double *eps);

/*!
 *
 * Compute area of region where level set function is greater than 0.
 *
 * Parameters
 * ----------
 * phi: level set function
 *
 * dx: array containing grid spacing in each coordinate direction
 *
 * eps: width of numerical smoothing to use for Heaviside function
 *
 * *_gb_lo: lower corner of index range for ghostbox
 *
 * *_gb_hi: upper corner of index range for ghostbox
 *
 * patch_box_lo: lower corner of index range for patch box
 *
 * patch_box_hi: upper corner of index range for patch box
 *
 * Return value
 * ------------
 * area of region where phi > 0
 *
 * Notes
 * -----
 * - When phi is a signed distance function, the numerical width of the
 *   Heaviside function as a function of spatial position is approximately
 *   equal to (2 * eps).
 */
PQS_REAL LSM_2D_AREA_PHI_GREATER_THAN_ZERO(
    const PQS_REAL* phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx,
    const double *eps);

#ifdef __cplusplus
}
#endif

#endif // INCLUDED_PQS_math_level_set_method_2d_h
