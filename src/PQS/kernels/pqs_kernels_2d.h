/*! \file kernels_2d.h
 *
 * \brief
 * Header files for numerical kernels that support PQS algorithm.
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

#ifndef INCLUDED_PQS_pqs_kernels_2d_h
#define INCLUDED_PQS_pqs_kernels_2d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                          name in
 *      C/C++ code                       Fortran code
 *      ----------                       ------------
 */
#define PQS_2D_COMPUTE_DN                pqs2dcomputedn_

/*!
 * TODO
 */
void PQS_2D_COMPUTE_DN(
   const int&, const double *, const double *, const double *,
   const int&, const int&,
   const int&, const int&,
   const int&,
   const int&,
   double *,
   const int&,
   const double *, const double *);

#ifdef __cplusplus
}
#endif

#endif
