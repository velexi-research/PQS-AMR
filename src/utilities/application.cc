/*! \file application.cc
 *
 * \brief
 * Implementation file for PQS application utility functions.
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

// Standard
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <string>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/math/LSMToolbox.h"
#include "PQS/pqs/Solver.h"
#include "PQS/utilities/error.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchHierarchy; } }

// Namespaces
using namespace std;
using namespace PQS;
using namespace SAMRAI;


namespace PQS {

// --- Function implementations

shared_ptr<tbox::Database> initialize_pqs(int *argc, char **argv[])
{
    // --- Initialize SAMRAI

    tbox::SAMRAI_MPI::init(argc, argv);
    tbox::SAMRAIManager::initialize();
    tbox::SAMRAIManager::startup();
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

    // --- Process command line

    string config_file;
    string restart_read_dir;
    int restore_num = 0;

    bool is_from_restart = false;

    if ( (*argc != 2) && (*argc != 4) ) {
        tbox::pout << "USAGE:  " << (*argv)[0] << " <config file> "
             << "\n"
             << "<restart dir> <restore number> [options]\n"
             << "  options:\n"
             << "  none at this time"
             << endl;
        PQS_ERROR_STATIC("pqs::application", "initialize_pqs",
                         "Invalid number of arguments passed to application.");
    } else {
        config_file = (*argv)[1];
        if (*argc == 4) {
            restart_read_dir = (*argv)[2];
            restore_num = atoi((*argv)[3]);
            is_from_restart = true;
        }
    }

    // --- Load parameters from configuration file

    shared_ptr<tbox::InputDatabase> config_db =
        shared_ptr<tbox::InputDatabase>(
            new tbox::InputDatabase("config_db"));
    tbox::InputManager::getManager()->parseInputFile(config_file, config_db);

    // Get the base name for all name strings in application
    const string base_name =
        config_db->getStringWithDefault("base_name", "PQS-OUTPUT");

    // --- Set is_from_restart in configuration database

    config_db->putBool("is_from_restart", is_from_restart);

    // --- Open the restart file (if restarting simulation)

    if (is_from_restart) {
        tbox::RestartManager::getManager()->openRestartFile(
                restart_read_dir, restore_num, mpi.getSize());
    }

    // --- Start logging

    const string log_file_name = base_name + ".log";
    const bool log_all_nodes = config_db->getBoolWithDefault("log_all_nodes",
                                                           false);
    if (log_all_nodes) {
        tbox::PIO::logAllNodes(log_file_name);
    } else {
        tbox::PIO::logOnlyNodeZero(log_file_name);
    }

    // Log the command-line args
    tbox::plog << "config_file= " << config_file << endl;
    tbox::plog << "restart_read_dir = " << restart_read_dir << endl;
    tbox::plog << "restore_num = " << restore_num << endl;

    return config_db;
}

void shutdown_pqs()
{
    // Shutdown SAMRAI
    tbox::SAMRAIManager::shutdown();
    tbox::SAMRAIManager::finalize();
    tbox::SAMRAI_MPI::finalize();
}

void run_pqs(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<pqs::Solver>& pqs_solver,
        const shared_ptr<appu::VisItDataWriter>& viz_data_writer)
{
    // --- Preparations

    // Get the base name for all name strings in application
    const string base_name = config_db->getString("base_name");

    // Get debugging parameters
    const bool enable_debug =
            config_db->getBoolWithDefault("enable_debug", false);

    // Set up restart parameters
    const bool is_from_restart =
        config_db->getBoolWithDefault("is_from_restart", false);

    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

    int restart_interval =
        config_db->getIntegerWithDefault("restart_interval", 0);

    const string restart_write_dir =
        config_db->getStringWithDefault("restart_dir",
                                        base_name + ".restart");

    const bool write_restart = (restart_interval > 0);

    // Set up visualization parameters
    const int viz_write_interval =
            config_db->getIntegerWithDefault("viz_write_interval", 0);

    const bool write_viz_files =
            (viz_data_writer != NULL) && (viz_write_interval > 0);

    if (write_viz_files && enable_debug) {
        viz_data_writer->registerPlotQuantity(
            "LSE RHS", "SCALAR", pqs_solver->getLSERHSPatchDataId());
    }

    // Emit contents of config database to log file.
    tbox::plog << "Configuration database..." << endl;
    config_db->printClassData(tbox::plog);

    // Set loop constants
    const shared_ptr<hier::PatchHierarchy> patch_hierarchy =
        pqs_solver->getPatchHierarchy();
    const double initial_curvature = pqs_solver->getInitialCurvature();
    const double final_curvature = pqs_solver->getFinalCurvature();
    const double curvature_increment = pqs_solver->getCurvatureIncrement();

    // Set up loop variables
    double curvature = pqs_solver->getCurvature();
    const int step_count = pqs_solver->getStepCount();

    // Initialize calculation
    if (!is_from_restart) {
        if (pqs_solver->initializeWithSlightlyCompressibleModel()) {
            // Emit status message
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << "  Equilibrating initial interface using "
                       << "Slightly Compressible Model ... " << endl;
            tbox::pout << "  Approximate curvature: " << curvature << endl;

            pqs_solver->equilibrateInterface(
                    curvature, pqs::SLIGHTLY_COMPRESSIBLE_MODEL);

        } else {
            // Emit status message
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << "  Skipping equilibration of initial " << endl;
            tbox::pout << "  interface using Slightly Compressible " << endl;
            tbox::pout << "  Model ... " << endl;
        }
    } else {
        // Simulation variables loaded from restart files, so fluid-fluid
        // interface already equilibrated. No additional equilibration
        // required.

        // Close restart file
        restart_manager->closeRestartFile();
    }

    // Write restart and viz files for initial conditions
    // (if this run is not from restart)
    if (!is_from_restart) {
        // Write restart files
        if (write_restart) {
            restart_manager->writeRestartFile(restart_write_dir,
                                              step_count);
        }

        // Write viz files
        if (write_viz_files) {
            viz_data_writer->writePlotData(patch_hierarchy, step_count,
                                           curvature);
        }
    }

    // --- Main curvature-stepping loop

    while (curvature <= final_curvature) {
        const int step_count = pqs_solver->getStepCount();

        // Emit status message
        tbox::pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
        if (curvature == initial_curvature) {
            tbox::pout << "  Initial curvature: " << curvature << endl;
        } else if (curvature == final_curvature) {
            tbox::pout << "  Final curvature: " << curvature << endl;
        } else {
            tbox::pout << "  Current curvature: " << curvature << endl;
        }
        tbox::pout << "  Step count: " << step_count << endl;

        // Update fluid-fluid interface
        pqs_solver->equilibrateInterface(curvature,
                                         pqs::PRESCRIBED_CURVATURE_MODEL);

        // Write restart file
        if (write_restart) {
            if (step_count % restart_interval == 0) {
                restart_manager->writeRestartFile(restart_write_dir,
                                                  step_count);
            }
        }

        // Write VisIt data
        if (write_viz_files) {
            if (step_count % viz_write_interval == 0) {
                viz_data_writer->writePlotData(patch_hierarchy, step_count,
                                               curvature);
            }
        }

        // Emit debugging messages
        if (enable_debug) {
            double solid_phase_volume = math::LSM::computeVolume(
                    pqs_solver->getPatchHierarchy(),
                    pqs_solver->getPoreSpacePatchDataId(),
                    -1, // compute volume for psi < 0
                    pqs_solver->getControlVolumePatchDataId());
            double non_wetting_phase_volume = math::LSM::computeVolume(
                    pqs_solver->getPatchHierarchy(),
                    pqs_solver->getInterfacePatchDataId(),
                    -1, // compute volume for phi < 0
                    pqs_solver->getControlVolumePatchDataId())
                - solid_phase_volume;

            tbox::pout << "-------------- DEBUGGING -----------------" << endl;
            tbox::pout << "Volume of non-wetting phase: "
                       << non_wetting_phase_volume << endl;
            tbox::pout << "------------------------------------------" << endl;
            tbox::pout << endl;
        }

        // Exit loop if we have reached the final curvature value
        if (curvature == final_curvature) {
           break;
        }

        // Compute next curvature target
        curvature += curvature_increment;
        if (final_curvature < curvature) {
            curvature = final_curvature;
        }
    }

    // Emit status message
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;

    // Write restart file for final time step
    if (write_restart) {
        if (step_count % restart_interval != 0) {
            restart_manager->writeRestartFile(restart_write_dir, step_count);
       }
    }

    // Write VisIt data for final time step
    if (write_viz_files) {
        if (step_count % viz_write_interval != 0) {
            viz_data_writer->writePlotData(patch_hierarchy, step_count,
                                           curvature);
        }
    }
}

}  // PQS namespace
