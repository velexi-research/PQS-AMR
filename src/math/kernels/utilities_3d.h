/*! \file utilities_3d.h
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

#ifndef INCLUDED_PQS_math_utilities_3d_h
#define INCLUDED_PQS_math_utilities_3d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                          name in
 *      C/C++ code                       Fortran code
 *      ----------                       ------------
 */
#define PQS_MATH_3D_MAX_NORM_DIFF        pqsmath3dmaxnormdiff_

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
 * max norm of (u - v)
 */
PQS_REAL PQS_MATH_3D_MAX_NORM_DIFF(
    const PQS_REAL* u,
    const int *u_gb_lo,
    const int *u_gb_hi,
    const PQS_REAL* v,
    const int *v_gb_lo,
    const int *v_gb_hi,
    const int *ib_lo,
    const int *ib_hi);

#ifdef __cplusplus
}
#endif

#endif // INCLUDED_PQS_math_utilities_3d_h
