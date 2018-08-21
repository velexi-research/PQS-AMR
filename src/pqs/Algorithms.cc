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
#include <sstream>
#include <stddef.h>
#include <string>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities/error.h"
#include "PQS/utilities/macros.h"
#include "PQS/pqs/Algorithms.h"
#include "PQS/pqs/kernels_2d.h"
#include "PQS/pqs/kernels_3d.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchData; } }
namespace SAMRAI { namespace hier { class PatchGeometry; } }


// --- Class implementation

namespace PQS {
namespace pqs {


// --- Public methods

Algorithms::Algorithms(
        const boost::shared_ptr<tbox::Database>& config_db,
        const int lse_rhs_id,
        const int psi_id,
        const int grad_psi_id)
{
    // Check arguments
    if (config_db == NULL) {
        PQS_ERROR(this, "Algorithms", "'config_db' must not be NULL");
    }
    verifyConfigurationDatabase(config_db);

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

    // Load configuration from config_db
    loadConfiguration(config_db);

} // Algorithms::Algorithms()

double Algorithms::computePrescribedCurvatureModelRHS(
        const boost::shared_ptr<hier::Patch> patch,
        const int phi_id) const
{
    patch->getPatchData(phi_id);
    return 1.0;
} // Algorithms::computePrescribedCurvatureModelRHS()

double Algorithms::computeSlightlyCompressibleModelRHS(
        const boost::shared_ptr<hier::Patch> patch,
        const int phi_id,
        const double volume) const
{
    // --- Preparations

    // Maximum stable time step
    double max_stable_dt;

    // Get geometry parameters
    boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom =
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch->getPatchGeometry());
    const double* dx = patch_geom->getDx();

    // ------- Get pointers to data and index space ranges

    // RHS
    boost::shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
            BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(d_lse_rhs_id));

    hier::Box rhs_ghostbox = rhs_data->getGhostBox();
    const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
    const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(rhs_ghostbox_lo, rhs_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(rhs_ghostbox_hi, rhs_ghostbox_upper);

    PQS_REAL* rhs = rhs_data->getPointer();

    // phi
    boost::shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(phi_id));

    hier::Box phi_ghostbox = phi_data->getGhostBox();
    const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
    const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_lo, phi_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_hi, phi_ghostbox_upper);

    PQS_REAL* phi = phi_data->getPointer();

    // fill box
    hier::Box fillbox = patch->getBox();
    const hier::IntVector fillbox_lower = fillbox.lower();
    const hier::IntVector fillbox_upper = fillbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(fillbox_lo, fillbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(fillbox_hi, fillbox_upper);

    // --- Compute RHS and maximum stable time step

    if (d_contact_angle == 0.0) {
        if (patch->getDim().getValue() == 2) {
            PQS_2D_COMPRESSIBLE_MODEL_ZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                fillbox_lo, fillbox_hi,
                dx,
                &volume,
                &d_target_volume,
                &d_reference_pressure,
                &d_bulk_modulus,
                &d_surface_tension);
        } else {
            PQS_3D_COMPRESSIBLE_MODEL_ZERO_CONTACT_ANGLE_RHS(
                &max_stable_dt,
                rhs,
                rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi,
                phi_ghostbox_lo, phi_ghostbox_hi,
                fillbox_lo, fillbox_hi,
                dx,
                &volume,
                &d_target_volume,
                &d_reference_pressure,
                &d_bulk_modulus,
                &d_surface_tension);
        }
    } else {
        // --- Get pointers to data and index space ranges for
        //     psi and grad(psi)

        // psi
        boost::shared_ptr< pdat::CellData<PQS_REAL> > psi_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                        patch->getPatchData(d_psi_id));

        hier::Box psi_ghostbox = psi_data->getGhostBox();
        const hier::IntVector psi_ghostbox_lower = psi_ghostbox.lower();
        const hier::IntVector psi_ghostbox_upper = psi_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_lo, psi_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_hi, psi_ghostbox_upper);

        PQS_REAL* psi = psi_data->getPointer();

        // grad(psi)
        boost::shared_ptr< pdat::CellData<PQS_REAL> > grad_psi_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
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

        if (patch->getDim().getValue() == 2) {
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
                fillbox_lo, fillbox_hi,
                dx,
                &volume,
                &d_target_volume,
                &d_reference_pressure,
                &d_bulk_modulus,
                &d_surface_tension,
                &d_contact_angle);

        } else {
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
                fillbox_lo, fillbox_hi,
                dx,
                &volume,
                &d_target_volume,
                &d_reference_pressure,
                &d_bulk_modulus,
                &d_surface_tension,
                &d_contact_angle);
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


// --- Private methods

void Algorithms::verifyConfigurationDatabase(
    const boost::shared_ptr<tbox::Database>& config_db) const
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'config_db' must not be NULL");
    }

    // --- Verify "SlightlyCompressibleModel" database

    if (!config_db->isDatabase("SlightlyCompressibleModel")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  string("'SlightlyCompressibleModel' database ") +
                  string("missing from 'config_db'"));
    }
    boost::shared_ptr<tbox::Database> scm_config_db =
            config_db->getDatabase("SlightlyCompressibleModel");

    if (!scm_config_db->isDouble("reference_pressure")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'reference_pressure' missing from 'PQS' database");
    }
    if (!scm_config_db->isDouble("target_volume")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'target_volume' missing from 'PQS' database");
    }
    if (!scm_config_db->isDouble("surface_tension")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'surface_tension' missing from 'PQS' database");
    }
    if (!scm_config_db->isDouble("bulk_modulus")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'bulk_modulus' missing from 'PQS' database");
    }

} // Algorithms::verifyConfigurationDatabase()

void Algorithms::loadConfiguration(
        const boost::shared_ptr<tbox::Database>& config_db)
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR(this, "loadConfiguration",
                  "'config_db' must not be NULL");
    }

    // --- Load configuration parameters

    if (!config_db->isDouble("contact_angle")) {
        d_contact_angle = 0.0;
    } else {
        d_contact_angle = config_db->getDouble("contact_angle");
        if (d_contact_angle < 0) {
            PQS_ERROR(this, "loadConfiguration",
                      "'contact_angle' must be non-negative.");
        }
    }

    // Slightly Compressible Model parameters

    // Prescribed Curvature parameters

} // Algorithms::loadConfiguration()

} // PQS::pqs namespace
} // PQS namespace
