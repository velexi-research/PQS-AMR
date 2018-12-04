/*! \file pore_kernels_3d.h
 *
 * \brief
 * Header files for numerical kernels to initialize fluid-fluid
 * interface for example application with pore space defined by a simple
 * packing of spheres.
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

#ifndef INCLUDED_pore_kernels_3d_h
#define INCLUDED_pore_kernels_3d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                    name in
 *      C/C++ code                 Fortran code
 *      ----------                 ------------
 */
#define INIT_PORE_SPACE_3D         init_pore_space_3d_

/*!
 * See documentation in pore_kernels_3d.f.in
 */
void INIT_PORE_SPACE_3D(
    const PQS_REAL *phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *fb_lo,
    const int *fb_hi,
    const double *x_lower,
    const double *dx,
    const double *radius,
    const double *center_1,
    const double *center_2,
    const double *center_3,
    const double *center_4,
    const double *center_5,
    const double *center_6,
    const double *center_7,
    const double *center_8);

#ifdef __cplusplus
}
#endif

#endif // INCLUDED_pore_kernels_3d_h
