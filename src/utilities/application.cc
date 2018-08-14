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
#include <iosfwd>
#include <iostream>
#include <stdlib.h>
#include <string>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/hier/VariableDatabase.h"
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
#include "PQS/utilities/macros.h"
#include "PQS/pqs/Solver.h"

// Class/type declarations
namespace PQS { namespace pqs { class InterfaceInitStrategy; } }
namespace PQS { namespace pqs { class PoreInitStrategy; } }

// Namespaces
using namespace std;
using namespace PQS;
using namespace SAMRAI;


namespace PQS {

// --- Function implementations

boost::shared_ptr<tbox::Database> initialize_pqs(int argc, char *argv[])
{
    // --- Initialize SAMRAI

    tbox::SAMRAI_MPI::init(&argc, &argv);
    tbox::SAMRAIManager::initialize();
    tbox::SAMRAIManager::startup();
    const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

    // --- Process command line

    string config_file;
    string restart_read_dir;
    int restore_num = 0;

    bool is_from_restart = false;

    if ( (argc != 2) && (argc != 4) ) {
        tbox::pout << "USAGE:  " << argv[0] << " <config file> "
             << "\n"
             << "<restart dir> <restore number> [options]\n"
             << "  options:\n"
             << "  none at this time"
             << endl;
        tbox::SAMRAI_MPI::abort();
        PQS_ERROR_STATIC("pqs::application", "initialize_pqs",
                         "Invalid number of arguments passed to application.");
    } else {
        config_file = argv[1];
        if (argc == 4) {
            restart_read_dir = argv[2];
            restore_num = atoi(argv[3]);
            is_from_restart = true;
        }
    }

    // --- Load parameters from configuration file

    boost::shared_ptr<tbox::InputDatabase> config_db =
        boost::shared_ptr<tbox::InputDatabase>(
            new tbox::InputDatabase("config_db"));
    tbox::InputManager::getManager()->parseInputFile(config_file, config_db);

    // Get the base name for all name strings in application
    const string base_name =
        config_db->getStringWithDefault("base_name", "unnamed");

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
        const boost::shared_ptr<tbox::Database> config_db,
        const boost::shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const boost::shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy)
{
    // --- Preparations

    // Get the base name for all name strings in application
    const string base_name =
        config_db->getStringWithDefault("base_name", "unnamed");

    // Set up restart manager
    const bool is_from_restart =
        config_db->getBoolWithDefault("is_from_restart", false);

    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

    int restart_interval =
        config_db->getIntegerWithDefault("restart_interval", 0);

    const string restart_write_dir =
        config_db->getStringWithDefault("restart_dir",
                                        base_name + ".restart");

    const bool write_restart = (restart_interval > 0);

    // --- Create PQS Solver

    // Construct PQS::pqs::Solver object
    //solver = new pqs::Solver(config_db, pore_init_strategy,
    //                         interface_init_strategy);

    // Emit contents of config database and variable database to log file.
    tbox::plog << "Configuration database..." << endl;
    config_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

    // --- Set up visualization data writers

    bool use_visit = false;
    if (config_db->keyExists("use_visit")) {
        use_visit = config_db->getBool("use_visit");
    }

    // Set up viz write interval
    int viz_write_interval = -1;
    if (use_visit) {
        if (config_db->keyExists("viz_write_interval")) {
            viz_write_interval =
                config_db->getInteger("viz_write_interval");
        }
    }

    // Set up extra VisIt parameters
    int visit_number_procs_per_file = 1;
    if (use_visit) {
        if (config_db->keyExists("visit_number_procs_per_file")) {
            visit_number_procs_per_file =
                config_db->getInteger("visit_number_procs_per_file");
        }
    }

    /* TODO: activate
    boost::shared_ptr<appu::VisItDataWriter> visit_data_writer = 0;
    if (use_visit) {
        string visit_data_dir = base_name + ".visit";
        visit_data_writer = boost::shared_ptr<appu::VisItDataWriter>
           (new appu::VisItDataWriter(
                dim, "VisIt Writer", visit_data_dir,
                visit_number_procs_per_file));

        // Register data to write to visualization file
        // TODO
    }
    */

    // --- Initialize calculation

    // TODO: initialize simulation objects

    // Close restart file before starting main time-stepping loop
    restart_manager->closeRestartFile();

    // Set up loop variables
    int count = 0;
    double dt = 0; // TODO
    double current_time = 0; // TODO
    int cur_integrator_step = 0; // TODO

    // Output initial conditions (if this run is not from restart)
    if ( write_restart && (!is_from_restart) ) {
        restart_manager->writeRestartFile(restart_write_dir,
                                          cur_integrator_step);
    }

    // Write VisIt data for initial time step
    /* TODO: activate
    if ( use_visit && (!is_from_restart) ) {
        visit_data_writer->writePlotData(patch_hierarchy, cur_integrator_step,
                                         current_time);
    }
    */

    // --- Main time loop
    // TODO
    /*
    while ( !lsm_algorithm->endTimeReached() &&
            ((max_num_time_steps <= 0) || (count < max_num_time_steps)) ) {

        tbox::pout << "++++++++++++++++++++++++++++++++++++++++++"
            << endl;
        tbox::pout << "  Time step (in current run): " << count << endl;
        tbox::pout << "  Integrator time step: " << cur_integrator_step
            << endl;
        tbox::pout << "  Current time:  " << current_time << endl;

        // Compute next time step
        dt = lsm_algorithm->computeStableDt();
        double end_time = lsm_algorithm->getEndTime();
        if (end_time - current_time < dt) dt = end_time - current_time;
        tbox::pout << "  dt:  " << dt << endl;

        // Advance level set functions
        lsm_algorithm->advanceLevelSetFunctions(dt);

        // Add an extra line to output for aesthetic reasons
        tbox::pout << endl;

        // Output data for current time step if this is the
        // initial time step or if the next write interval has
        // been reached
        cur_integrator_step = lsm_algorithm->numIntegrationStepsTaken();

        // Write restart file
        if ( write_restart && (0==cur_integrator_step%restart_interval) ) {
            restart_manager->writeRestartFile(restart_write_dir,
                                          cur_integrator_step);
        }

        // Write VisIt data
        if ( use_visit && (0==cur_integrator_step%viz_write_interval) ) {
            visit_data_writer->writePlotData(patch_hierarchy,
                                             cur_integrator_step,
                                             lsm_algorithm->getCurrentTime());
      }

      // Update counter and current time
      count++;
      current_time = lsm_algorithm->getCurrentTime();

    }

    // Output information for final time step
    // (if it hasn't already been output)
    current_time = lsm_algorithm->getCurrentTime();
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;
    tbox::pout << "  Final time step (in current run): " << count << endl;
    tbox::pout << "  Final integrator time step: " << cur_integrator_step
        << endl;
    tbox::pout << "  Current time:  " << current_time << endl;
    tbox::pout << endl;
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++" << endl;

    // Write restart file for final time step
    if ( write_restart && (0!=cur_integrator_step%restart_interval) ) {
        restart_manager->writeRestartFile(restart_write_dir,
                                          cur_integrator_step);
    }

    // Write VisIt data for final time step
    if ( use_visit && (0!=cur_integrator_step%viz_write_interval) ) {
        visit_data_writer->writePlotData(patch_hierarchy, cur_integrator_step,
                                         lsm_algorithm->getCurrentTime());
    }
    */
}

}  // PQS namespace
