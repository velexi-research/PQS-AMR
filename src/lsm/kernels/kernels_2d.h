/*! \file kernels_2d.h
 *
 * \brief
 * Header file numerical kernels for Level Set Method computations in
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

#ifndef INCLUDED_PQS_lsm_kernels_2d_h
#define INCLUDED_PQS_lsm_kernels_2d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                          name in
 *      C/C++ code                       Fortran code
 *      ----------                       ------------
 */
#define PQS_LSM_2D_AREA_PHI_LESS_THAN_ZERO \
                                         pqsLsm2dAreaPhiLessThanZero_
#define PQS_LSM_2D_AREA_PHI_GREATER_THAN_ZERO \
                                         pqsLsm2dAreaPhiGreaterThanZero_

/*!
 *
 * Compute area of region where signed distance function is less than 0.
 *
 * Parameters
 * ----------
 * phi: signed distance function
 *
 * dx: array containing grid spacing in each coordinate direction
 *
 * eps: width of numerical smoothing to use for Heaviside function
 *
 * *_gb_lo: lower corner of index range for ghostbox
 *
 * *_gb_hi: upper corner of index range for ghostbox
 *
 * ib_lo: lower corner of index range for interior box
 *
 * ib_hi: upper corner of index range for interior box
 *
 * Return value
 * ------------
 * area of region where phi < 0
 */
PQS_REAL PQS_LSM_2D_AREA_PHI_LESS_THAN_ZERO(
    const PQS_REAL* phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *ib_lo,
    const int * ib_hi,
    const double *dx,
    const double *eps);

/*!
 *
 * Compute area of region where signed distance function is greater than 0.
 *
 * Parameters
 * ----------
 * phi: signed distance function
 *
 * dx: array containing grid spacing in each coordinate direction
 *
 * eps: width of numerical smoothing to use for Heaviside function
 *
 * *_gb_lo: lower corner of index range for ghostbox
 *
 * *_gb_hi: upper corner of index range for ghostbox
 *
 * ib_lo: lower corner of index range for interior box
 *
 * ib_hi: upper corner of index range for interior box
 *
 * Return value
 * ------------
 * area of region where phi > 0
 */
PQS_REAL PQS_LSM_2D_AREA_PHI_GREATER_THAN_ZERO(
    const PQS_REAL* phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *ib_lo,
    const int * ib_hi,
    const double *dx,
    const double *eps);

#ifdef __cplusplus
}
#endif

#endif // INCLUDED_PQS_lsm_kernels_2d_h