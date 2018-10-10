/*! \file Algorithms.cc
 *
 * \brief
 * Implementation file for Algorithms class.
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

// --- Headers, namespaces, and type declarations

// Standard library
#include <memory>
#include <ostream>

// SAMRAI
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities/error.h"
#include "PQS/utilities/macros.h"
#include "PQS/pqs/Algorithms.h"
#include "PQS/pqs/kernels/kernels_2d.h"
#include "PQS/pqs/kernels/kernels_3d.h"

// Class/type declarations


// --- Class implementation

namespace PQS {
namespace pqs {


// --- Public methods

Algorithms::Algorithms(
        const int lse_rhs_id,
        const int psi_id,
        const int grad_psi_id)
{
    // Check arguments
    if (lse_rhs_id < 0) {
        PQS_ERROR(this, "Algorithms", "'lse_rhs_id' must be non-negative");
    }

    if (psi_id < 0) {
        PQS_ERROR(this, "Algorithms", "'psi_id' must be non-negative");
    }

    // Set data members
    d_lse_rhs_id = lse_rhs_id;
    d_psi_id = psi_id;
    d_grad_psi_id = grad_psi_id;

} // Algorithms::Algorithms()

double Algorithms::computePrescribedCurvatureModelRHS(
        const shared_ptr<hier::Patch>& patch,
        const int phi_id,
        const int psi_id,
        const int grad_psi_id,
        const double target_curvature,
        const double surface_tension,
        const double contact_angle) const
{
    // --- Preparations

    // Get dimensionality of space
    const int dim = patch->getDim().getValue();

    // Compute pressure required to interface to have target curvature at
    // steady-state
    const double pressure = surface_tension * target_curvature;

    // Maximum stable time step
    double max_stable_dt;

    // Get geometry parameters
    shared_ptr<geom::CartesianPatchGeometry> patch_geom =
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                patch->getPatchGeometry());
    const double* dx = patch_geom->getDx();

    // --- Get pointers to data and index space ranges

    // RHS
    shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch->getPatchData(d_lse_rhs_id));

    hier::Box rhs_ghostbox = rhs_data->getGhostBox();
    const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
    const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(rhs_ghostbox_lo, rhs_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(rhs_ghostbox_hi, rhs_ghostbox_upper);

    PQS_REAL* rhs = rhs_data->getPointer();

    // phi
    shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch->getPatchData(phi_id));

    hier::Box phi_ghostbox = phi_data->getGhostBox();
    const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
    const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_lo, phi_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_hi, phi_ghostbox_upper);

    PQS_REAL* phi = phi_data->getPointer();

    // patch box
    hier::Box patch_box = patch->getBox();
    const hier::IntVector patch_box_lower = patch_box.lower();
    const hier::IntVector patch_box_upper = patch_box.upper();
    PQS_INT_VECT_TO_INT_ARRAY(patch_box_lo, patch_box_lower);
    PQS_INT_VECT_TO_INT_ARRAY(patch_box_hi, patch_box_upper);

    // --- Compute RHS and maximum stable time step

    if (contact_angle == 0.0) {
        if (dim == 2) {
            PQS_2D_CURVATURE_MODEL_ZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx,
                &pressure,
                &surface_tension);
        } else if (dim == 3) {
            PQS_3D_CURVATURE_MODEL_ZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx,
                &pressure,
                &surface_tension);
        }
    } else {
        // --- Get pointers to data and index space ranges for
        //     psi and grad(psi)

        // psi
        shared_ptr< pdat::CellData<PQS_REAL> > psi_data =
                SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData(d_psi_id));

        hier::Box psi_ghostbox = psi_data->getGhostBox();
        const hier::IntVector psi_ghostbox_lower = psi_ghostbox.lower();
        const hier::IntVector psi_ghostbox_upper = psi_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_lo, psi_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_hi, psi_ghostbox_upper);

        PQS_REAL* psi = psi_data->getPointer();

        // grad(psi)
        shared_ptr< pdat::CellData<PQS_REAL> > grad_psi_data =
                SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData(d_grad_psi_id));

        hier::Box grad_psi_ghostbox = grad_psi_data->getGhostBox();
        const hier::IntVector grad_psi_ghostbox_lower =
                grad_psi_ghostbox.lower();
        const hier::IntVector grad_psi_ghostbox_upper =
                grad_psi_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(grad_psi_ghostbox_lo,
                grad_psi_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(grad_psi_ghostbox_hi,
                grad_psi_ghostbox_upper);

        if (dim == 2) {
            // Get pointers to data arrays for grad(phi)
            PQS_REAL* grad_psi_x = grad_psi_data->getPointer(0);
            PQS_REAL* grad_psi_y = grad_psi_data->getPointer(1);

            PQS_2D_CURVATURE_MODEL_NONZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                psi,
                psi_ghostbox_lo, psi_ghostbox_hi,
                grad_psi_x, grad_psi_y,
                grad_psi_ghostbox_lo, grad_psi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx,
                &pressure,
                &surface_tension,
                &contact_angle);

        } else if (dim == 3) {
            // Get pointers to data arrays for grad(phi)
            PQS_REAL* grad_psi_x = grad_psi_data->getPointer(0);
            PQS_REAL* grad_psi_y = grad_psi_data->getPointer(1);
            PQS_REAL* grad_psi_z = grad_psi_data->getPointer(2);

            PQS_3D_CURVATURE_MODEL_NONZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                psi,
                psi_ghostbox_lo, psi_ghostbox_hi,
                grad_psi_x, grad_psi_y, grad_psi_z,
                grad_psi_ghostbox_lo, grad_psi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx,
                &pressure,
                &surface_tension,
                &contact_angle);
        }
    }

    return max_stable_dt;

} // Algorithms::computePrescribedCurvatureModelRHS()

double Algorithms::computeSlightlyCompressibleModelRHS(
        const shared_ptr<hier::Patch>& patch,
        const int phi_id,
        const int psi_id,
        const int grad_psi_id,
        const double target_curvature,
        const double surface_tension,
        const double bulk_modulus,
        const double volume,
        const double target_volume,
        const double contact_angle) const
{
    // --- Preparations

    // Get dimensionality of space
    const int dim = patch->getDim().getValue();

    // Compute pressure required to interface to have target curvature at
    // steady-state
    const double pressure = surface_tension * target_curvature;

    // Maximum stable time step
    double max_stable_dt;

    // Get geometry parameters
    shared_ptr<geom::CartesianPatchGeometry> patch_geom =
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                patch->getPatchGeometry());
    const double* dx = patch_geom->getDx();

    // --- Get pointers to data and index space ranges

    // RHS
    shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch->getPatchData(d_lse_rhs_id));

    hier::Box rhs_ghostbox = rhs_data->getGhostBox();
    const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
    const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(rhs_ghostbox_lo, rhs_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(rhs_ghostbox_hi, rhs_ghostbox_upper);

    PQS_REAL* rhs = rhs_data->getPointer();

    // phi
    shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch->getPatchData(phi_id));

    hier::Box phi_ghostbox = phi_data->getGhostBox();
    const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
    const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_lo, phi_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_hi, phi_ghostbox_upper);

    PQS_REAL* phi = phi_data->getPointer();

    // patch box
    hier::Box patch_box = patch->getBox();
    const hier::IntVector patch_box_lower = patch_box.lower();
    const hier::IntVector patch_box_upper = patch_box.upper();
    PQS_INT_VECT_TO_INT_ARRAY(patch_box_lo, patch_box_lower);
    PQS_INT_VECT_TO_INT_ARRAY(patch_box_hi, patch_box_upper);

    // --- Compute RHS and maximum stable time step

    if (contact_angle == 0.0) {
        if (dim == 2) {
            PQS_2D_COMPRESSIBLE_MODEL_ZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx,
                &pressure,
                &volume,
                &target_volume,
                &bulk_modulus,
                &surface_tension);
        } else if (dim == 3) {
            PQS_3D_COMPRESSIBLE_MODEL_ZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx,
                &pressure,
                &volume,
                &target_volume,
                &bulk_modulus,
                &surface_tension);
        }
    } else {
        // --- Get pointers to data and index space ranges for
        //     psi and grad(psi)

        // psi
        shared_ptr< pdat::CellData<PQS_REAL> > psi_data =
                SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData(d_psi_id));

        hier::Box psi_ghostbox = psi_data->getGhostBox();
        const hier::IntVector psi_ghostbox_lower = psi_ghostbox.lower();
        const hier::IntVector psi_ghostbox_upper = psi_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_lo, psi_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_hi, psi_ghostbox_upper);

        PQS_REAL* psi = psi_data->getPointer();

        // grad(psi)
        shared_ptr< pdat::CellData<PQS_REAL> > grad_psi_data =
                SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData(d_grad_psi_id));

        hier::Box grad_psi_ghostbox = grad_psi_data->getGhostBox();
        const hier::IntVector grad_psi_ghostbox_lower =
                grad_psi_ghostbox.lower();
        const hier::IntVector grad_psi_ghostbox_upper =
                grad_psi_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(grad_psi_ghostbox_lo,
                grad_psi_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(grad_psi_ghostbox_hi,
                grad_psi_ghostbox_upper);

        if (dim == 2) {
            // Get pointers to data arrays for grad(phi)
            PQS_REAL* grad_psi_x = grad_psi_data->getPointer(0);
            PQS_REAL* grad_psi_y = grad_psi_data->getPointer(1);

            PQS_2D_COMPRESSIBLE_MODEL_NONZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                psi,
                psi_ghostbox_lo, psi_ghostbox_hi,
                grad_psi_x, grad_psi_y,
                grad_psi_ghostbox_lo, grad_psi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx,
                &pressure,
                &volume,
                &target_volume,
                &bulk_modulus,
                &surface_tension,
                &contact_angle);

        } else if (dim == 3) {
            // Get pointers to data arrays for grad(phi)
            PQS_REAL* grad_psi_x = grad_psi_data->getPointer(0);
            PQS_REAL* grad_psi_y = grad_psi_data->getPointer(1);
            PQS_REAL* grad_psi_z = grad_psi_data->getPointer(2);

            PQS_3D_COMPRESSIBLE_MODEL_NONZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                psi,
                psi_ghostbox_lo, psi_ghostbox_hi,
                grad_psi_x, grad_psi_y, grad_psi_z,
                grad_psi_ghostbox_lo, grad_psi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx,
                &pressure,
                &volume,
                &target_volume,
                &bulk_modulus,
                &surface_tension,
                &contact_angle);
        }
    }

    return max_stable_dt;

} // Algorithms::computeSlightlyCompressibleModelRHS()

void Algorithms::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::Algorithms::printClassData..." << endl;
    os << "(Algorithms*) this = " << (Algorithms*) this << endl;
} // Algorithms::printClassData()

} // PQS::pqs namespace
} // PQS namespace
