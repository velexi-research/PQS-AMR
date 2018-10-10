/*! \file InterfaceInitModule.cc
 *
 * \brief
 * Implementation file for concrete subclass of pqs::InterfaceInitStrategy to
 * use in example application with pore space defined by a periodic packing of
 * circles (2D) or spheres (3D).
 */

/*
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE. This file is part of the XYZ package. It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution. No part of the XYZ
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

// --- Headers, namespaces, and type declarations

// Standard library
#include <stddef.h>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities/error.h"
#include "PQS/utilities/macros.h"

// PQS application
#include "InterfaceInitModule.h"
#include "kernels/interface_kernels_2d.h"
#include "kernels/interface_kernels_3d.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations


// --- Public Methods

InterfaceInitModule::InterfaceInitModule(
        const shared_ptr<tbox::Database>& config_db, const int dim,
        const vector<double>& X_lower, const vector<double>& X_upper)
{
    // Check arguments
    if (config_db == NULL) {
        PQS_ERROR(this, "InterfaceInitModule", "'config_db' must not be NULL");
    }
    verifyConfigurationDatabase(config_db);

    if ((dim != 2) and (dim != 3)) {
        PQS_ERROR(this, "InterfaceInitModule",
                  "'num_dimensions' must equal 2 or 3");
    }

    // Set problem dimension
    d_dim = dim;

    // Set X_lower and X_upper
    d_X_lower = new double[dim];
    d_X_upper = new double[dim];
    for (int i = 0; i < dim; i++) {
        d_X_lower[i] = X_lower[i];
        d_X_upper[i] = X_upper[i];
    }

    // Load configuration from config_db
    loadConfiguration(config_db);
} // IntefaceInitModule::IntefaceInitModule()

InterfaceInitModule::~InterfaceInitModule()
{
    // Clean up memory
    delete [] d_X_lower;
    delete [] d_X_upper;

} // IntefaceInitModule::~IntefaceInitModule()

void InterfaceInitModule::initializeInterface(
        hier::Patch& patch, int phi_id)
{
    // --- Preparations

    // Get geometry parameters
    shared_ptr<geom::CartesianPatchGeometry> patch_geom =
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                    patch.getPatchGeometry());
    const double* dx = patch_geom->getDx();
    const double* x_lower = patch_geom->getXLower();

    // ------- Get pointers to data and index space ranges

    // phi
    shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch.getPatchData(phi_id));

    hier::Box phi_ghostbox = phi_data->getGhostBox();
    const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
    const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_lo, phi_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_hi, phi_ghostbox_upper);

    PQS_REAL* phi = phi_data->getPointer();

    // fill box
    hier::Box fillbox = patch.getBox();
    const hier::IntVector fillbox_lower = fillbox.lower();
    const hier::IntVector fillbox_upper = fillbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(fillbox_lo, fillbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(fillbox_hi, fillbox_upper);

    // --- Initialize phi

    if (d_dim == 2) {
        const double center[2] = {
            0.5 * (d_X_lower[0] + d_X_upper[0]),
            d_X_lower[1]
        };

        INIT_CIRCLE(phi,
                    phi_ghostbox_lo, phi_ghostbox_hi,
                    fillbox_lo, fillbox_hi,
                    x_lower,
                    dx,
                    center,
                    &d_radius);
    } else if (d_dim == 3) {
        const double center[3] = {
            0.5 * (d_X_lower[0] + d_X_upper[0]),
            0.5 * (d_X_lower[1] + d_X_upper[1]),
            d_X_lower[2]
        };

        INIT_SPHERE(phi,
                    phi_ghostbox_lo, phi_ghostbox_hi,
                    fillbox_lo, fillbox_hi,
                    x_lower,
                    dx,
                    center,
                    &d_radius);
    }
} // InterfaceInitModule::initializeInterface()


// --- Private methods

void InterfaceInitModule::verifyConfigurationDatabase(
    const shared_ptr<tbox::Database>& config_db) const
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'config_db' must not be NULL");
    }

    // --- Verify configuration parameters

    if (!config_db->isDouble("radius")) {
        PQS_ERROR(this, "verifyConfigurationDatabase",
                  "'radius' missing from 'InterfaceInitModule' database");
    }

} // InterfaceInitModule::verifyConfigurationDatabase()

void InterfaceInitModule::loadConfiguration(
        const shared_ptr<tbox::Database>& config_db)
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR(this, "loadConfiguration",
                  "'config_db' must not be NULL");
    }

    // --- Load configuration parameters

    d_radius = config_db->getDouble("radius");

} // InterfaceInitModule::loadConfiguration()
