/*! \file utilities_2d.h
 *
 * \brief
 * Header file for numerical kernels math utility functions
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

#ifndef INCLUDED_PQS_math_utilities_2d_h
#define INCLUDED_PQS_math_utilities_2d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                          name in
 *      C/C++ code                       Fortran code
 *      ----------                       ------------
 */
#define PQS_MATH_2D_MAX_NORM_DIFF        pqsmath2dmaxnormdiff_
#define PQS_MATH_2D_MIN_UV               pqsmath2dminuv_

/*!
 *
 * Compute max norm of the difference between 'u' and 'v'.
 *
 * Parameters
 * ----------
 * u: first function in the expression |u - v|
 *
 * v: second function in the expression |u - v|
 *
 * patch_box_lo: lower corner of index range for patch box
 *
 * patch_box_hi: upper corner of index range for patch box
 *
 * *_gb_lo: lower corner of index range for ghost box
 *
 * *_gb_hi: upper corner of index range for ghost box
 *
 * Return value
 * ------------
 * max norm of (u - v)
 */
PQS_REAL PQS_MATH_2D_MAX_NORM_DIFF(
    const PQS_REAL* u,
    const int *u_gb_lo,
    const int *u_gb_hi,
    const PQS_REAL* v,
    const int *v_gb_lo,
    const int *v_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi);

/*!
 * Compute min(u, v).
 *
 * Parameters
 * ----------
 * min_uv: min(u, v)
 *
 * u: first function in the expression min(u, v)
 *
 * v: second function in the expression min(u, v)
 *
 * patch_box_lo: lower corner of index range for patch box
 *
 * patch_box_hi: upper corner of index range for patch box
 *
 * *_gb_lo: lower corner of index range for ghost box
 *
 * *_gb_hi: upper corner of index range for ghost box
 *
 */
void PQS_MATH_2D_MIN_UV(
    const PQS_REAL* min_uv,
    const int *min_uv_gb_lo,
    const int *min_uv_gb_hi,
    const PQS_REAL* u,
    const int *u_gb_lo,
    const int *u_gb_hi,
    const PQS_REAL* v,
    const int *v_gb_lo,
    const int *v_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi);

#ifdef __cplusplus
}
#endif

#endif // INCLUDED_PQS_math_utilities_2d_h
