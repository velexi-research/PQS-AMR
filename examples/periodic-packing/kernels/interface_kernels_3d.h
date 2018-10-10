/*! \file interface_kernels_3d.h
 *
 * \brief
 * Header files for numerical kernels to initialize fluid-fluid
 * interface for example application with pore space defined by a periodic
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

#ifndef INCLUDED_interface_kernels_3d_h
#define INCLUDED_interface_kernels_3d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                    name in
 *      C/C++ code                 Fortran code
 *      ----------                 ------------
 */
#define INIT_SPHERE                init_sphere_

/*!
 * See documentation in interface_kernels_3d.f.in
 */
void INIT_SPHERE(
    const PQS_REAL *phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *fb_lo,
    const int *fb_hi,
    const double *x_lower,
    const double *dx,
    const double *center,
    const double *radius);

#ifdef __cplusplus
}
#endif

#endif // INCLUDED_interface_kernels_3d_h
