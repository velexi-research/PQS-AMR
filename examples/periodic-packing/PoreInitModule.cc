/*! \file PoreInitModule.cc
 *
 * \brief
 * Implementation file for concrete subclass of pqs::PoreInitStrategy to use
 * in example application with pore space defined by a periodic packing of
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
#include <memory>

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
#include "PoreInitModule.h"
#include "kernels/pore_kernels_2d.h"
#include "kernels/pore_kernels_3d.h"

// Namespaces
using namespace std;
using namespace SAMRAI;
using namespace PQS;

// Class/type declarations


// --- Public Methods

PoreInitModule::PoreInitModule(
        const shared_ptr<tbox::Database>& config_db, const int dim,
        const vector<double>& X_lower, const vector<double>& X_upper)
{
    // Check arguments
    if (config_db == NULL) {
        PQS_ERROR(this, "PoreInitModule", "'config_db' must not be NULL");
    }
    verifyConfigurationDatabase(config_db);

    if ((dim != 2) and (dim != 3)) {
        PQS_ERROR(this, "PoreInitModule",
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

PoreInitModule::~PoreInitModule()
{
    // Clean up memory
    delete [] d_X_lower;
    delete [] d_X_upper;

} // IntefaceInitModule::~IntefaceInitModule()

void PoreInitModule::initializePoreSpace(
        hier::Patch& patch, int psi_id)
{
    // --- Preparations

    // Get geometry parameters
    shared_ptr<geom::CartesianPatchGeometry> patch_geom =
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                    patch.getPatchGeometry());
    const double* dx = patch_geom->getDx();
    const double* x_lower = patch_geom->getXLower();

    // ------- Get pointers to data and index space ranges

    // psi
    shared_ptr< pdat::CellData<PQS_REAL> > psi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch.getPatchData(psi_id));

    hier::Box psi_ghostbox = psi_data->getGhostBox();
    const hier::IntVector psi_ghostbox_lower = psi_ghostbox.lower();
    const hier::IntVector psi_ghostbox_upper = psi_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_lo, psi_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_hi, psi_ghostbox_upper);

    PQS_REAL* psi = psi_data->getPointer();

    // fill box
    hier::Box fillbox = patch.getBox();
    const hier::IntVector fillbox_lower = fillbox.lower();
    const hier::IntVector fillbox_upper = fillbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(fillbox_lo, fillbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(fillbox_hi, fillbox_upper);

    // --- Initialize psi

    if (d_dim == 2) {
        const double center_1[2] = {d_X_lower[0], d_X_lower[1]};
        const double center_2[2] = {d_X_upper[0], d_X_lower[1]};
        const double center_3[2] = {d_X_lower[0], d_X_upper[1]};
        const double center_4[2] = {d_X_upper[0], d_X_upper[1]};

        INIT_CIRCLES_AT_CORNERS(psi,
                                psi_ghostbox_lo, psi_ghostbox_hi,
                                fillbox_lo, fillbox_hi,
                                x_lower,
                                dx,
                                &d_radius,
                                center_1,
                                center_2,
                                center_3,
                                center_4);
    } else if (d_dim == 3) {
        const double center_1[3] = {d_X_lower[0], d_X_lower[1], d_X_lower[2]};
        const double center_2[3] = {d_X_upper[0], d_X_lower[1], d_X_lower[2]};
        const double center_3[3] = {d_X_lower[0], d_X_upper[1], d_X_lower[2]};
        const double center_4[3] = {d_X_upper[0], d_X_upper[1], d_X_lower[2]};
        const double center_5[3] = {d_X_lower[0], d_X_lower[1], d_X_upper[2]};
        const double center_6[3] = {d_X_upper[0], d_X_lower[1], d_X_upper[2]};
        const double center_7[3] = {d_X_lower[0], d_X_upper[1], d_X_upper[2]};
        const double center_8[3] = {d_X_upper[0], d_X_upper[1], d_X_upper[2]};

        INIT_SPHERES_AT_CORNERS(psi,
                                psi_ghostbox_lo, psi_ghostbox_hi,
                                fillbox_lo, fillbox_hi,
                                x_lower,
                                dx,
                                &d_radius,
                                center_1,
                                center_2,
                                center_3,
                                center_4,
                                center_5,
                                center_6,
                                center_7,
                                center_8);
    }
} // PoreInitModule::initializeInterface()


// --- Private methods

void PoreInitModule::verifyConfigurationDatabase(
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
                  "'radius' missing from 'PoreInitModule' database");
    }

} // PoreInitModule::verifyConfigurationDatabase()

void PoreInitModule::loadConfiguration(
        const shared_ptr<tbox::Database>& config_db)
{
    // --- Check arguments

    if (config_db == NULL) {
        PQS_ERROR(this, "loadConfiguration",
                  "'config_db' must not be NULL");
    }

    // --- Load configuration parameters

    d_radius = config_db->getDouble("radius");

} // PoreInitModule::loadConfiguration()
