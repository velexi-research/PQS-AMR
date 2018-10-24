/*! \file kernels_3d.h
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

#ifndef INCLUDED_PQS_pqs_kernels_3d_h
#define INCLUDED_PQS_pqs_kernels_3d_h

#ifdef __cplusplus
extern "C" {
#endif

/* Link between C/C++ and Fortran function names
 *
 *      name in                name in
 *      C/C++ code             Fortran code
 *      ----------             ------------
 */
#define PQS_3D_COMPRESSIBLE_MODEL_ZERO_CONTACT_ANGLE_RHS \
                               pqs3dcompressiblemodelzerocontactanglerhs_
#define PQS_3D_COMPRESSIBLE_MODEL_NONZERO_CONTACT_ANGLE_RHS \
                               pqs3dcompressiblemodelnonzerocontactanglerhs_
#define PQS_3D_CURVATURE_MODEL_ZERO_CONTACT_ANGLE_RHS \
                               pqs3dcurvaturemodelzerocontactanglerhs_
#define PQS_3D_CURVATURE_MODEL_NONZERO_CONTACT_ANGLE_RHS \
                               pqs3dcurvaturemodelnonzerocontactanglerhs_

#define PQS_3D_TAG_CELLS_FOR_REFINEMENT \
                               pqs3dtagcellsforrefinement_

/*!
 * See documentation in kernels_3d.f.in
 */
void PQS_3D_CURVATURE_MODEL_ZERO_CONTACT_ANGLE_RHS(
    const double *max_stable_dt,
    const PQS_REAL *rhs,
    const int *rhs_gb_lo,
    const int *rhs_gb_hi,
    const PQS_REAL *phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx,
    const double *pressure,
    const double *surface_tension);

/*!
 * See documentation in kernels_3d.f.in
 */
void PQS_3D_CURVATURE_MODEL_NONZERO_CONTACT_ANGLE_RHS(
    const double *max_stable_dt,
    const PQS_REAL *rhs,
    const int *rhs_gb_lo,
    const int *rhs_gb_hi,
    const PQS_REAL *phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const PQS_REAL *psi,
    const int *psi_gb_lo,
    const int *psi_gb_hi,
    const PQS_REAL *grad_psi_x,
    const PQS_REAL *grad_psi_y,
    const PQS_REAL *grad_psi_z,
    const int *grad_psi_gb_lo,
    const int *grad_psi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx,
    const double *pressure,
    const double *surface_tension,
    const double *contact_angle);

/*!
 * See documentation in kernels_3d.f.in
 */
void PQS_3D_COMPRESSIBLE_MODEL_ZERO_CONTACT_ANGLE_RHS(
    const double *max_stable_dt,
    const PQS_REAL *rhs,
    const int *rhs_gb_lo,
    const int *rhs_gb_hi,
    const PQS_REAL *phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx,
    const double *pressure,
    const double *volume,
    const double *target_volume,
    const double *bulk_modulus,
    const double *surface_tension);

/*!
 * See documentation in kernels_3d.f.in
 */
void PQS_3D_COMPRESSIBLE_MODEL_NONZERO_CONTACT_ANGLE_RHS(
    const double *max_stable_dt,
    const PQS_REAL *rhs,
    const int *rhs_gb_lo,
    const int *rhs_gb_hi,
    const PQS_REAL *phi,
    const int *phi_gb_lo,
    const int *phi_gb_hi,
    const PQS_REAL *psi,
    const int *psi_gb_lo,
    const int *psi_gb_hi,
    const PQS_REAL *grad_psi_x,
    const PQS_REAL *grad_psi_y,
    const PQS_REAL *grad_psi_z,
    const int *grad_psi_gb_lo,
    const int *grad_psi_gb_hi,
    const int *patch_box_lo,
    const int *patch_box_hi,
    const double *dx,
    const double *pressure,
    const double *volume,
    const double *target_volume,
    const double *bulk_modulus,
    const double *surface_tension,
    const double *contact_angle);

/*!
 * See documentation in kernels_3d.f.in
 */
void PQS_3D_TAG_CELLS_FOR_REFINEMENT(
     const int *tag,
     const int *tag_gb_lo,
     const int *tag_gb_hi,
     const PQS_REAL *phi,
     const int *phi_gb_lo,
     const int *phi_gb_hi,
     const PQS_REAL *psi,
     const int *psi_gb_lo,
     const int *psi_gb_hi,
     const int *patch_box_lo,
     const int *patch_box_hi,
     const double *refinement_cutoff);

#ifdef __cplusplus
}
#endif

#endif
