/*! \file TagAndInitModule.cc
 *
 * \brief
 * Implementation file for TagAndInitModule class.
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
#include <cstddef>
#include <sstream>
#include <string>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"
#include "PQS/pqs/TagAndInitModule.h"
#include "PQS/pqs/utilities.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class BoxContainer; } }
namespace SAMRAI { namespace hier { class PatchData; } }
namespace SAMRAI { namespace hier { class RefineOperator; } }
namespace SAMRAI { namespace hier { class Variable; } }
namespace SAMRAI { namespace tbox { class Database; } }


// --- Class implementation

namespace PQS {
namespace pqs {

// --- Implementation of public methods

// Constructor
TagAndInitModule::TagAndInitModule(
        const boost::shared_ptr<tbox::Database>& config_db,
        const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const boost::shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const boost::shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy,
        const int phi_id, const int psi_id):
    mesh::TagAndInitializeStrategy(d_object_name)
{
    // Check parameters
    if (config_db == NULL) {
        PQS_ERROR(this, "TagAndInitModule", "'config_db' must not be NULL");
    }
    if (pore_init_strategy == NULL) {
        PQS_ERROR(this, "TagAndInitModule",
                  "'pore_init_strategy' must not be NULL");
    }
    if (interface_init_strategy == NULL) {
        PQS_ERROR(this, "TagAndInitModule",
                  "'interface_init_strategy' must not be NULL");
    }
    if (phi_id < 0) {
        PQS_ERROR(this, "TagAndInitModule", "Invalid value for 'phi_id'");
    }
    if (psi_id < 0) {
        PQS_ERROR(this, "TagAndInitModule", "Invalid value for 'psi_id'");
    }

    // Set data members
    d_phi_id = phi_id;
    d_psi_id = psi_id;
    d_pore_init_strategy = pore_init_strategy;
    d_interface_init_strategy = interface_init_strategy;

    // Initialize communication objects
    initializeCommunicationObjects(patch_hierarchy->getGridGeometry());

} // TagAndInitModule::TagAndInitModule()

void TagAndInitModule::initializeLevelData(
        const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int patch_level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const boost::shared_ptr<hier::PatchLevel>& old_patch_level,
        const bool allocate_data)
{
    // Check parameters
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "'patch_hierarchy' may not be NULL");
    }
    if (patch_level_number < 0) {
        PQS_ERROR(this, "initializeLevelData",
                  "'patch_level_number' must non-negative");
    }
    if (patch_level_number > patch_hierarchy->getFinestLevelNumber()) {
        PQS_ERROR(this, "initializeLevelData",
                  "'patch_level_number' exceeds finest PatchLevel number");
    }
    if (old_patch_level != NULL) {
        if (patch_level_number != old_patch_level->getLevelNumber()) {
            PQS_ERROR(this, "initializeLevelData",
                      std::string("'patch_level_number' must equal ") +
                      std::string("PatchLevel number of 'old_patch_level'"));
        }
    }
    if (patch_hierarchy->getPatchLevel(patch_level_number) == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "PatchLevel at 'patch_level_number' must not be NULL");
    }

    // Get PatchLevel
    boost::shared_ptr<hier::PatchLevel> patch_level =
        patch_hierarchy->getPatchLevel(patch_level_number);

    // Allocate PatchData
    if (allocate_data) {
        patch_level->allocatePatchData(d_phi_id);
        patch_level->allocatePatchData(d_psi_id);
    } else {
        patch_level->setTime(init_data_time, d_phi_id);
        patch_level->setTime(init_data_time, d_psi_id);
    }

    if (initial_time) {
        // Initialize data on Patches on the PatchLevel.
        for (hier::PatchLevel::Iterator pi(patch_level->begin());
                pi!=patch_level->end(); pi++) {

            boost::shared_ptr<hier::Patch> patch = *pi;
            if (patch == NULL) {
                PQS_ERROR(this, "initializeLevelData",
                          "Null Patch pointer found when iterating over "
                          "PatchLevel.");
            }

            // Initialize level set functions for solid-pore and fluid-fluid
            // interfaces.
            d_pore_init_strategy->initializePoreSpace(*patch, d_psi_id);
            d_interface_init_strategy->initializeInterface(*patch, d_phi_id);
        }
    } else {
        // If appropriate, fill new PatchLevel with data from the old PatchLevel
        if ((patch_level_number > 0) || old_patch_level!=NULL) {
            boost::shared_ptr<xfer::RefineSchedule> sched =
                d_xfer_fill_new_level->createSchedule(
                patch_level, old_patch_level, patch_level_number-1,
                patch_hierarchy, NULL);

            sched->fillData(init_data_time);
        }
    }

    // deallocate data on old PatchLevel if it exists
    if (old_patch_level!=NULL) {
        old_patch_level->deallocatePatchData(d_phi_id);
        old_patch_level->deallocatePatchData(d_psi_id);
    }
} // TagAndInitModule::initializeLevelData()

void TagAndInitModule::resetHierarchyConfiguration(
        const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int coarsest_patch_level_number,
        const int finest_patch_level_number)
{
    // Check parameters
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  "'patch_hierarchy' may not be NULL");
    }
    if (coarsest_patch_level_number < 0) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  "'coarsest_patch_level_number' must be non-negative");
    }
    if (coarsest_patch_level_number > finest_patch_level_number) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  std::string("'coarsest_patch_level_number' must be less ") +
                  std::string("than or equal to 'finest_patch_level_number'"));
    }
    for (int ln = 0; ln <= finest_patch_level_number; ln++) {
        if (patch_hierarchy->getPatchLevel(ln) == NULL) {
            PQS_ERROR(this, "resetHierarchyConfiguration",
                      std::string("PatchLevel ") +
                      std::to_string(ln) +
                      std::string(" is NULL in PatchHierarchy"));
        }
    }

    // TODO: reset communication schedules
    // - reset communications schedules used to fill boundary data
    //   during time advance

    // recompute control volumes
    // LevelSetMethodToolbox::computeControlVolumes(
    // patch_hierarchy, d_control_volume_handle);

} // TagAndInitModule::resetHierarchyConfiguration()

void TagAndInitModule::tagCellsForRefinement(
        const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int patch_level_number,
        const int regrid_cycle,
        const double regrid_time,
        const int tag_id,
        const bool initial_time,
        const bool coarsest_sync_patch_level,
        const bool can_be_refined,
        const double regrid_start_time)
{
    // Check parameters
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "'patch_hierarchy' may not be NULL");
    }
    if (patch_level_number < 0) {
        PQS_ERROR(this, "initializeLevelData",
                  "'patch_level_number' must non-negative");
    }
    if (patch_level_number > patch_hierarchy->getFinestLevelNumber()) {
        PQS_ERROR(this, "initializeLevelData",
                  "'patch_level_number' exceeds finest PatchLevel number");
    }
    if (patch_hierarchy->getPatchLevel(patch_level_number) == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "PatchLevel at 'patch_level_number' must not be NULL");
    }

    // Get PatchLevel
    boost::shared_ptr<hier::PatchLevel> patch_level =
        patch_hierarchy->getPatchLevel(patch_level_number);

    // Initialize data on Patches on the PatchLevel.
    for (hier::PatchLevel::Iterator pi(patch_level->begin());
            pi!=patch_level->end(); pi++) {

        boost::shared_ptr<hier::Patch> patch = *pi;
        if (patch == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "Null Patch pointer found when iterating over "
                  "PatchLevel.");
        }

        boost::shared_ptr< pdat::CellData<int> > tag_data =
            BOOST_CAST<pdat::CellData<int>, hier::PatchData>(
                patch->getPatchData(tag_id));

        tag_data->fill(1);
    }
} // TagAndInitModule::tagCellsForRefinement()

bool TagAndInitModule::refineUserBoxInputOnly(int cycle, double time)
{
    // TODO
    return false;
} // TagAndInitModule::refineUserBoxInputOnly()

bool TagAndInitModule::getUserSuppliedRefineBoxes(
        hier::BoxContainer& refine_boxes,
        const int patch_level_number,
        const int cycle,
        const double time)
{
    // TODO
    return false;
} // TagAndInitModule::getUserSuppliedRefineBoxes()

void TagAndInitModule::resetRefineBoxes(
        const hier::BoxContainer& refine_boxes,
        const int patch_level_number)
{
} // TagAndInitModule::resetRefineBoxes()

void TagAndInitModule::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::TagAndInitModule::printClassData..." << endl;
    os << "(TagAndInitModule*) this = " << (TagAndInitModule*) this << endl;
    os << "d_pore_init_strategy = " << d_pore_init_strategy.get() << endl;
    os << "d_interface_init_strategy = " << d_interface_init_strategy.get()
       << endl;

    os << endl;
    d_pore_init_strategy->printClassData(os);
    d_interface_init_strategy->printClassData(os);
} // TagAndInitModule::printClassData()

// --- Implementation of private methods

void TagAndInitModule::loadConfiguration(
        const boost::shared_ptr<tbox::Database>& config_db)
{
} // TagAndInitModule::loadConfiguration()

void TagAndInitModule::initializeCommunicationObjects(
        const boost::shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{
    // Get pointer to VariableDatabase
    hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

    // Get variable associated with d_phi_id
    boost::shared_ptr< hier::Variable > phi_variable;
    if (!var_db->mapIndexToVariable(d_phi_id, phi_variable)) {
        PQS_ERROR(this, "initializeCommunicationObjects",
                  "'d_phi_id' not found in VariableDatabase");
    }

    // Lookup refine operations
    // TODO: review choice of refinement operator
    boost::shared_ptr<hier::RefineOperator> refine_op =
        grid_geometry->lookupRefineOperator(phi_variable, "CONSTANT_REFINE");

    // --- Set up communications objects for filling a new level
    //     (used during initialization of a PatchLevel)
    d_xfer_fill_new_level =
        boost::shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);

    d_xfer_fill_new_level->registerRefine(
        d_phi_id, d_phi_id, d_phi_id, refine_op);
    d_xfer_fill_new_level->registerRefine(
        d_psi_id, d_psi_id, d_psi_id, refine_op);

} // TagAndInitModule::initializeCommunicationObjects()

// Copy constructor
TagAndInitModule::TagAndInitModule(
        const TagAndInitModule& rhs):
    mesh::TagAndInitializeStrategy(d_object_name)
{
} // TagAndInitModule::TagAndInitModule()

} // PQS::pqs namespace
} // PQS namespace
