/*! \file fixture.cc
 *
 * \brief
 * Implementation of fixture for unit tests for PQS::pqs classes.
 */

/*
 * ---------------------------------------------------------------------
 * COPYRIGHT/LICENSE. This file is part of the FDK package. It is
 * subject to the license terms in the LICENSE file found in the
 * top-level directory of this distribution. No part of the FDK
 * package, including this file, may be copied, modified, propagated,
 * or distributed except according to the terms contained in the
 * LICENSE file.
 * ---------------------------------------------------------------------
 */

// --- Headers, namespaces, and type declarations

// Standard library

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep

// PQS test
#include "fixture.h"


// --- Fixtures

namespace pqsTests {

// Constructor (set up)
pqsTests::pqsTests() {

    // --- Initialize SAMRAI

    int argc = 0;
    char **argv = 0;
    tbox::SAMRAI_MPI::init(&argc, &argv);

    tbox::SAMRAIManager::initialize();
    tbox::SAMRAIManager::startup();
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

    // --- Simulation parameters

    // Problem dimension
    const tbox::Dimension dim(3);

    // Geometry
    boost::shared_ptr<tbox::MemoryDatabase> geometry_input_db =
        boost::shared_ptr<tbox::MemoryDatabase>(
            new tbox::MemoryDatabase("CartesianGridGeometry"));

    double x_lo[3] = {-1.0, -1.0, -1.0};
    double x_up[3] = {1.0, 1.0, 1.0};
    geometry_input_db->putDoubleArray("x_lo", x_lo, 3);
    geometry_input_db->putDoubleArray("x_up", x_up, 3);

    int box_lower[3] = {0, 0, 0};
    int box_upper[3] = {49, 49, 49};
    tbox::DatabaseBox domain_boxes(dim, box_lower, box_upper);
    geometry_input_db->putDatabaseBoxArray("domain_boxes", &domain_boxes, 1);

    // --- Initialize Geometry and PatchHierarchy

    // Geometry
    grid_geometry = boost::shared_ptr<geom::CartesianGridGeometry>(
        new geom::CartesianGridGeometry(
            dim, "CartesianGeometry", geometry_input_db));

    // PatchHierarchy
    patch_hierarchy = boost::shared_ptr<hier::PatchHierarchy>(
        new hier::PatchHierarchy("PatchHierarchy", grid_geometry));
}

pqsTests::~pqsTests() {
    // Shutdown SAMRAI
    tbox::SAMRAIManager::shutdown();
    tbox::SAMRAIManager::finalize();
    tbox::SAMRAI_MPI::finalize();
}

} // pqsTests namespace
