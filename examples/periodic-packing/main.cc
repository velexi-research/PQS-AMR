/*! \file pqs_amr.cc
 *
 * \brief
 * Example PQS AMR application with pore space defined by a periodic packing
 * of circles (2D) or spheres (3D).
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
#include <iosfwd>
#include <memory>
#include <string>
#include <stddef.h>
#include <vector>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/Solver.h"
#include "PQS/utilities/application.h"
#include "PQS/utilities/error.h"

// PQS Application
#include "InterfaceInitModule.h"
#include "PoreInitModule.h"

// Class/type declarations
namespace PQS { namespace pqs { class InterfaceInitStrategy; } }
namespace PQS { namespace pqs { class PoreInitStrategy; } }

// Namespaces
using namespace std;
using namespace PQS;
using namespace SAMRAI;


// --- Main program

int main(int argc, char *argv[])
{
    // --- Preparations

    // Initialize PQS simulation
    shared_ptr<tbox::Database> config_db;
    try {
        config_db = initialize_pqs(&argc, &argv);
    } catch (const PQS::exception& e) {
        PQS_ABORT(e.what());
    }

    // Get problem dimension
    shared_ptr<tbox::Database> samrai_config_db =
            config_db->getDatabase("SAMRAI");
    int dim = samrai_config_db->getDatabase("Geometry")->getInteger("dim");
    vector<double> x_lower =
        samrai_config_db->getDatabase("Geometry")->getDoubleVector("x_lo");
    vector<double> x_upper =
        samrai_config_db->getDatabase("Geometry")->getDoubleVector("x_up");

    // Get the base name for all name strings in application
    const string base_name = config_db->getString("base_name");

    // Get debug mode
    const bool enable_debug =
            config_db->getBoolWithDefault("enable_debug", false);

    // --- Construct fluid-fluid interface and pore interface
    //     initialization modules

    // ------ PoreInitModule (implements PoreInitStrategy)

    shared_ptr<tbox::Database> pore_init_module_config_db =
            config_db->getDatabase("PoreInitModule");
    shared_ptr<pqs::PoreInitStrategy> pore_init_strategy;
    try {
        pore_init_strategy = shared_ptr<pqs::PoreInitStrategy>(
                new PoreInitModule(pore_init_module_config_db, dim,
                                   x_lower, x_upper));
    } catch (const PQS::exception& e) {
        PQS_ABORT(e.what());
    }

    // ------ InterfaceInitModule (implements InterfaceInitStrategy)

    shared_ptr<tbox::Database> interface_init_module_config_db =
            config_db->getDatabase("InterfaceInitModule");
    shared_ptr<pqs::InterfaceInitStrategy> interface_init_strategy;
    try {
        interface_init_strategy =
            shared_ptr<pqs::InterfaceInitStrategy>(
                new InterfaceInitModule(interface_init_module_config_db,
                                        dim, x_lower, x_upper));
    } catch (const PQS::exception& e) {
        PQS_ABORT(e.what());
    }

    // ------ Construct PQS::Solver

    shared_ptr<pqs::Solver> pqs_solver;
    try {
        pqs_solver = shared_ptr<pqs::Solver>(
                new pqs::Solver(config_db, pore_init_strategy,
                                interface_init_strategy,
                                NULL,  // null PatchHierarchy
                                enable_debug));
    } catch (const PQS::exception& e) {
        PQS_ABORT(e.what());
    }

    // ------Set up VisIt data writer

    shared_ptr<appu::VisItDataWriter> visit_data_writer;
    bool enable_viz = config_db->getBoolWithDefault("enable_viz", false);
    if (enable_viz) {
        string visit_data_dir = base_name + ".visit";
        int viz_num_procs_per_file =
                config_db->getIntegerWithDefault("viz_num_procs_per_file",
                                                 1);
        visit_data_writer =
                shared_ptr<appu::VisItDataWriter>(
                        new appu::VisItDataWriter(tbox::Dimension(dim),
                                                  "VisIt Writer",
                                                  visit_data_dir,
                                                  viz_num_procs_per_file));

        // Register data to write to visualization file
        visit_data_writer->registerPlotQuantity(
                "phi", "SCALAR", pqs_solver->getInterfacePatchDataId());
        visit_data_writer->registerPlotQuantity(
                "psi", "SCALAR", pqs_solver->getPoreSpacePatchDataId());
    }

    // --- Run PQS simulation

    try {
        run_pqs(config_db, pqs_solver, visit_data_writer);
    } catch (const PQS::exception& e) {
        PQS_ABORT(e.what());
    }

    // --- Shutdown PQS simulation

    shutdown_pqs();

    return(0);
}
