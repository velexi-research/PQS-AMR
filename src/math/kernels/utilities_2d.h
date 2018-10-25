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
 * See documentation in utilities_2d.f.in
 */
PQS_REAL PQS_MATH_2D_MAX_NORM_DIFF(
    const PQS_REAL* u,
    const int *u_gb_lo,
    const int *u_gb_hi,
    const PQS_REAL* v,
    const int *v_gb_lo,
    const int *v_gb_hi,
    const PQS_REAL* control_volume,
    const int *control_volume_gb_lo,
    const int *control_volume_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi);

/*!
 * See documentation in utilities_2d.f.in
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
