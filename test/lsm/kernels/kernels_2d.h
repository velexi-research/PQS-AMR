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

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                          name in
 *      C/C++ code                       Fortran code
 *      ----------                       ------------
 */
#define SET_PHI_CIRCLE                   setphicircle_

/*!
 * Compute signed distance function for circle of specified radius
 * centered at the origin.
 *
 * Parameters
 * ----------
 * phi: [output] signed distance function
 *
 * x_lo: array containing lower corner of grid
 *
 * dx: array containing grid spacing in each coordinate direction
 *
 * r: radius of circle
 *
 * *_gb_lo: lower corner of index range for ghostbox
 *
 * *_gb_hi: upper corner of index range for ghostbox
 *
 * ib_lo: lower corner of index range for interior box
 *
 * ib_hi: upper corner of index range for interior box
 *
 */
void SET_PHI_CIRCLE(
    const PQS_REAL* phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *ib_lo,
    const int *ib_hi,
    const double *x_lo,
    const double *dx,
    const double *r);

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
