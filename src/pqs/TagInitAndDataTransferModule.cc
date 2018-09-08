/*! \file TagInitAndDataTransferModule.cc
 *
 * \brief
 * Implementation file for TagInitAndDataTransferModule class.
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
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// SAMRAI
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities/error.h"
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"
#include "PQS/pqs/Solver.h"
#include "PQS/pqs/TagInitAndDataTransferModule.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class RefineOperator; } }
namespace SAMRAI { namespace hier { class Variable; } }
namespace SAMRAI { namespace hier { class VariableContext; } }
namespace SAMRAI { namespace tbox { class Database; } }


// --- Class implementation

namespace PQS {
namespace pqs {

// --- Static data members

const string TagInitAndDataTransferModule::s_object_name =
    "PQS::pqs::TagInitAndDataTransferModule";

// --- Public methods

// Constructor
TagInitAndDataTransferModule::TagInitAndDataTransferModule(
        const std::shared_ptr<tbox::Database>& config_db,
        const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const std::shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const std::shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy,
        const int phi_pqs_id,
        const int phi_lsm_current_id,
        const int phi_lsm_next_id,
        const int psi_id,
        const int control_volume_id,
        const hier::IntVector& max_stencil_width):
    mesh::TagAndInitializeStrategy(s_object_name)
{
    // Check arguments
    if (config_db == NULL) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'config_db' must not be NULL");
    }
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'patch_hierarchy' must not be NULL");
    }
    if (pore_init_strategy == NULL) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'pore_init_strategy' must not be NULL");
    }
    if (interface_init_strategy == NULL) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'interface_init_strategy' must not be NULL");
    }
    if (phi_pqs_id < 0) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'phi_pqs_id' must be non-negative");
    }
    if (phi_lsm_current_id < 0) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'phi_lsm_current_id' must be non-negative");
    }
    if (phi_lsm_next_id < 0) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'phi_lsm_next_id' must be non-negative");
    }
    if (psi_id < 0) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'psi_id' must be non-negative");
    }
    if (control_volume_id < 0) {
        PQS_ERROR(this, "TagInitAndDataTransferModule",
                  "'control_volume_id' must be non-negative");
    }

    // Set data members
    d_phi_pqs_id = phi_pqs_id;
    d_phi_lsm_current_id = phi_lsm_current_id;
    d_phi_lsm_next_id = phi_lsm_next_id;
    d_psi_id = psi_id;
    d_control_volume_id = control_volume_id;
    d_patch_hierarchy = patch_hierarchy;
    d_pore_init_strategy = pore_init_strategy;
    d_interface_init_strategy = interface_init_strategy;

    // Set up data transfer objects
    setupDataTransferObjects(patch_hierarchy->getGridGeometry(),
                             max_stencil_width);

} // TagInitAndDataTransferModule::TagInitAndDataTransferModule()

TagInitAndDataTransferModule::~TagInitAndDataTransferModule()
{
    // Free memory allocated for simulation data
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels();
            level_num++) {

        std::shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        patch_level->deallocatePatchData(d_scratch_variables);
    }
} // TagInitAndDataTransferModule::~TagInitAndDataTransferModule()

void TagInitAndDataTransferModule::fillGhostCells(
        const int context) const
{
    // Check arguments
    if ( (context != LSM_CURRENT) && (context != LSM_NEXT) ) {
        PQS_ERROR(this, "fillGhostCells",
                  string("Invalid 'context': ") +
                  to_string(context) +
                  string(". Valid values: LSM_CURRENT (=1), ") +
                  string("LSM_NEXT (=2)"));
    }

    // Fill ghost cells
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels();
            level_num++) {

        if (context == LSM_CURRENT) {
            d_xfer_fill_bdry_schedule_lsm_current[level_num]->fillData(
                    0.0,    // not used
                    true);  // apply physical boundary conditions
        } else if (context == LSM_NEXT) {
            d_xfer_fill_bdry_schedule_lsm_next[level_num]->fillData(
                    0.0,    // not used
                    true);  // apply physical boundary conditions
        }
    }
} // TagInitAndDataTransferModule::fillGhostCells()

void TagInitAndDataTransferModule::initializeLevelData(
        const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int level_num,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const std::shared_ptr<hier::PatchLevel>& old_patch_level,
        const bool allocate_data)
{
    // Check arguments
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "'patch_hierarchy' may not be NULL");
    }
    if (level_num < 0) {
        PQS_ERROR(this, "initializeLevelData",
                  "'level_num' must non-negative");
    }
    if (level_num > patch_hierarchy->getFinestLevelNumber()) {
        PQS_ERROR(this, "initializeLevelData",
                  "'level_num' exceeds finest PatchLevel number");
    }
    if (old_patch_level != NULL) {
        if (level_num != old_patch_level->getLevelNumber()) {
            PQS_ERROR(this, "initializeLevelData",
                      string("'level_num' must equal ") +
                      string("PatchLevel number of 'old_patch_level'"));
        }
    }
    if (patch_hierarchy->getPatchLevel(level_num) == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "PatchLevel at 'level_num' must not be NULL");
    }

    // Get PatchLevel
    std::shared_ptr<hier::PatchLevel> patch_level =
        patch_hierarchy->getPatchLevel(level_num);

    // Allocate PatchData
    if (allocate_data) {
        patch_level->allocatePatchData(d_phi_pqs_id);
        patch_level->allocatePatchData(d_psi_id);
        patch_level->allocatePatchData(d_control_volume_id);
    } else {
        patch_level->setTime(init_data_time, d_phi_pqs_id);
        patch_level->setTime(init_data_time, d_psi_id);
        patch_level->setTime(init_data_time, d_control_volume_id);
    }

    if (initial_time) {
        // Initialize data on Patches on the PatchLevel.
        for (hier::PatchLevel::Iterator pi(patch_level->begin());
                pi!=patch_level->end(); pi++) {

            std::shared_ptr<hier::Patch> patch = *pi;
            if (patch == NULL) {
                PQS_ERROR(this, "initializeLevelData",
                          "Null Patch pointer found when iterating over "
                          "PatchLevel.");
            }

            // Initialize level set functions for solid-pore and fluid-fluid
            // interfaces.
            d_pore_init_strategy->initializePoreSpace(*patch, d_psi_id);
            d_interface_init_strategy->initializeInterface(*patch,
                                                           d_phi_pqs_id);
        }
    } else {
        // If appropriate, fill new PatchLevel with data from the old PatchLevel
        if ((level_num > 0) || (old_patch_level!=NULL)) {
            std::shared_ptr<xfer::RefineSchedule> sched =
                d_xfer_fill_new_level->createSchedule(
                patch_level, old_patch_level, level_num-1,
                patch_hierarchy, NULL);

            sched->fillData(init_data_time);
        }
    }

    // deallocate data on old PatchLevel if it exists
    if (old_patch_level!=NULL) {
        old_patch_level->deallocatePatchData(d_phi_pqs_id);
        old_patch_level->deallocatePatchData(d_psi_id);
        old_patch_level->deallocatePatchData(d_control_volume_id);
    }
} // TagInitAndDataTransferModule::initializeLevelData()

void TagInitAndDataTransferModule::resetHierarchyConfiguration(
        const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int coarsest_level_num,
        const int finest_level_num)
{
    // --- Check arguments

    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  "'patch_hierarchy' may not be NULL");
    }
    if (patch_hierarchy != d_patch_hierarchy) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  string("'patch_hierarchy' argument does not match ") +
                  string("'patch_hierarchy' argument passed to ") +
                  string("constructor"));
    }
    if (coarsest_level_num < 0) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  "'coarsest_level_num' must be non-negative");
    }
    if (coarsest_level_num > finest_level_num) {
        PQS_ERROR(this, "resetHierarchyConfiguration",
                  string("'coarsest_level_num' must be less ") +
                  string("than or equal to 'finest_level_num'"));
    }
    for (int level_num = 0;
            level_num <= finest_level_num;
            level_num++) {

        if (patch_hierarchy->getPatchLevel(level_num) == NULL) {
            PQS_ERROR(this, "resetHierarchyConfiguration",
                      string("PatchLevel ") +
                      to_string(level_num) +
                      string(" is NULL in PatchHierarchy"));
        }
    }

    // --- Reset data transfer schedules

    // data transfer schedules for filling ghost cells data
    d_xfer_fill_bdry_schedule_lsm_current.resize(
        patch_hierarchy->getNumberOfLevels());
    d_xfer_fill_bdry_schedule_lsm_next.resize(
        patch_hierarchy->getNumberOfLevels());

    for (int level_num = coarsest_level_num;
            level_num <= finest_level_num;
            level_num++) {

        std::shared_ptr<hier::PatchLevel> patch_level =
            patch_hierarchy->getPatchLevel(level_num);

        if (level_num == 0) {
            d_xfer_fill_bdry_schedule_lsm_current[level_num] =
                d_xfer_fill_bdry_lsm_current->createSchedule(
                    patch_level, NULL);  // TODO: change NULL to boundary
                                         // condition module

            d_xfer_fill_bdry_schedule_lsm_next[level_num] =
                d_xfer_fill_bdry_lsm_next->createSchedule(
                    patch_level, NULL);  // TODO: change NULL to boundary
                                         // condition module
        } else {
            d_xfer_fill_bdry_schedule_lsm_current[level_num] =
                d_xfer_fill_bdry_lsm_current->createSchedule(
                    patch_level, level_num-1,
                    patch_hierarchy, NULL);  // TODO: change NULL to boundary
                                             // condition module

            d_xfer_fill_bdry_schedule_lsm_next[level_num] =
                d_xfer_fill_bdry_lsm_next->createSchedule(
                    patch_level, level_num-1,
                    patch_hierarchy, NULL);  // TODO: change NULL to boundary
                                             // condition module
        }
    }

    // --- Recompute simulation parameters

    // recompute control volumes
    computeControlVolumes();

} // TagInitAndDataTransferModule::resetHierarchyConfiguration()

void TagInitAndDataTransferModule::tagCellsForRefinement(
        const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int level_num,
        const int regrid_cycle,
        const double regrid_time,
        const int tag_id,
        const bool initial_time,
        const bool coarsest_sync_patch_level,
        const bool can_be_refined,
        const double regrid_start_time)
{
    // Check arguments
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "'patch_hierarchy' may not be NULL");
    }
    if (level_num < 0) {
        PQS_ERROR(this, "initializeLevelData",
                  "'level_num' must non-negative");
    }
    if (level_num > patch_hierarchy->getFinestLevelNumber()) {
        PQS_ERROR(this, "initializeLevelData",
                  "'level_num' exceeds finest PatchLevel number");
    }
    if (patch_hierarchy->getPatchLevel(level_num) == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "PatchLevel at 'level_num' must not be NULL");
    }

    // Get PatchLevel
    std::shared_ptr<hier::PatchLevel> patch_level =
        patch_hierarchy->getPatchLevel(level_num);

    // Initialize data on Patches on the PatchLevel.
    for (hier::PatchLevel::Iterator pi(patch_level->begin());
            pi!=patch_level->end(); pi++) {

        std::shared_ptr<hier::Patch> patch = *pi;
        if (patch == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "Null Patch pointer found when iterating over "
                  "PatchLevel.");
        }

        // TODO: implement actual cell tagging algorithm
        std::shared_ptr< pdat::CellData<int> > tag_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<int> >(
                patch->getPatchData(tag_id));

        tag_data->fill(1);
    }
} // TagInitAndDataTransferModule::tagCellsForRefinement()

bool TagInitAndDataTransferModule::refineUserBoxInputOnly(
        int cycle, double time)
{
    // TODO
    return false;
} // TagInitAndDataTransferModule::refineUserBoxInputOnly()

bool TagInitAndDataTransferModule::getUserSuppliedRefineBoxes(
        hier::BoxContainer& refine_boxes,
        const int level_num,
        const int cycle,
        const double time)
{
    // TODO
    return false;
} // TagInitAndDataTransferModule::getUserSuppliedRefineBoxes()

void TagInitAndDataTransferModule::resetRefineBoxes(
        const hier::BoxContainer& refine_boxes,
        const int level_num)
{
} // TagInitAndDataTransferModule::resetRefineBoxes()

void TagInitAndDataTransferModule::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::TagInitAndDataTransferModule::printClassData..." << endl;
    os << "(TagInitAndDataTransferModule*) this = "
       << (TagInitAndDataTransferModule*) this << endl;
    os << "d_pore_init_strategy = " << d_pore_init_strategy.get() << endl;
    os << "d_interface_init_strategy = " << d_interface_init_strategy.get()
       << endl;

    os << endl;
    d_pore_init_strategy->printClassData(os);
    d_interface_init_strategy->printClassData(os);
} // TagInitAndDataTransferModule::printClassData()


// --- Private methods

void TagInitAndDataTransferModule::loadConfiguration(
        const std::shared_ptr<tbox::Database>& config_db)
{
} // TagInitAndDataTransferModule::loadConfiguration()

void TagInitAndDataTransferModule::setupDataTransferObjects(
        const std::shared_ptr<hier::BaseGridGeometry>& grid_geometry,
        const hier::IntVector& max_stencil_width)
{
    // --- Preparations

    // Get pointer to VariableDatabase
    hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

    hier::IntVector scratch_ghost_cell_width(max_stencil_width);

    // Initialize PatchData component selectors
    d_scratch_variables.clrAllFlags();

    // --- Create PatchData for data transfer scratch space

    // Get 'scratch' variable context
    std::shared_ptr<hier::VariableContext> scratch_context =
        var_db->getContext("scratch");

    // Get variable associated with d_phi_lsm_current_id
    std::shared_ptr< hier::Variable > phi_variable;
    if (!var_db->mapIndexToVariable(d_phi_lsm_current_id, phi_variable)) {
        PQS_ERROR(this, "setupDataTransferObjects",
                  "Error mapping 'd_phi_lsm_current_id' to Variable object.");
    }
    d_phi_scratch_id =
        var_db->registerVariableAndContext(phi_variable,
                                           scratch_context,
                                           scratch_ghost_cell_width);

    // Get variable associated with d_psi_id
    std::shared_ptr< hier::Variable > psi_variable;
    if (!var_db->mapIndexToVariable(d_psi_id, psi_variable)) {
        PQS_ERROR(this, "setupDataTransferObjects",
                  "Error mapping 'd_psi_id' to Variable object.");
    }
    d_psi_scratch_id =
        var_db->registerVariableAndContext(psi_variable,
                                           scratch_context,
                                           scratch_ghost_cell_width);


    // --- Get refinement operators for filling PatchData from coarser
    //     PatchLevels

    std::shared_ptr<hier::RefineOperator> phi_linear_refine_op =
        grid_geometry->lookupRefineOperator(phi_variable, "LINEAR_REFINE");

    std::shared_ptr<hier::RefineOperator> psi_linear_refine_op =
        grid_geometry->lookupRefineOperator(psi_variable, "LINEAR_REFINE");

    // --- Set up data transfer objects

    // Filling a new level (used during initialization of a PatchLevel)
    d_xfer_fill_new_level =
        std::shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);

    d_xfer_fill_new_level->registerRefine(
        d_phi_pqs_id, d_phi_pqs_id, d_phi_scratch_id, phi_linear_refine_op);
    d_xfer_fill_new_level->registerRefine(
        d_psi_id, d_psi_id, d_psi_scratch_id, psi_linear_refine_op);

    // Filling a ghost cells during time integration of level set functions
    d_xfer_fill_bdry_lsm_current =
        std::shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);
    d_xfer_fill_bdry_lsm_current->registerRefine(
        d_phi_lsm_current_id, d_phi_lsm_current_id, d_phi_scratch_id,
        phi_linear_refine_op);

    d_xfer_fill_bdry_lsm_next =
        std::shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);
    d_xfer_fill_bdry_lsm_next->registerRefine(
        d_phi_lsm_next_id, d_phi_lsm_next_id, d_phi_scratch_id,
        phi_linear_refine_op);

} // TagInitAndDataTransferModule::setupDataTransferObjects()

void TagInitAndDataTransferModule::computeControlVolumes() const
{
    // On every level, set control volume to volume of grid cell on level.
    const int finest_level_num = d_patch_hierarchy->getFinestLevelNumber();

    for (int level_num = finest_level_num; level_num >= 0; level_num--) {
        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(patch_level->begin());
                pi!=patch_level->end(); pi++) {

            shared_ptr<hier::Patch> patch = *pi;
            shared_ptr<geom::CartesianPatchGeometry> patch_geometry =
                SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                    patch->getPatchGeometry());
            int dim = d_patch_hierarchy->getDim().getValue();
            const double* dx = patch_geometry->getDx();

            double cell_volume = dx[0];
            for (int i = 1; i < dim; i++) {
                cell_volume *= dx[i];
            }

            shared_ptr< pdat::CellData<PQS_REAL> > control_volume_data =
                    SAMRAI_SHARED_PTR_CAST<pdat::CellData<PQS_REAL>>(
                            patch->getPatchData(d_control_volume_id));

            if (!control_volume_data) {
                PQS_ERROR(this, "computeControlVolumes",
                          string("'d_control_volume_id' does not refer to a ") +
                          string("pdat_CellVariable"));
            }

            control_volume_data->fillAll(cell_volume);
        }

        // On all but the finest level, set control volume to 0 for all
        // cells covered by finer cells.
        if (level_num < finest_level_num) {
            // Get boxes that describe index space of the next finer level
            // and coarsen them to describe corresponding index space at
            // this level.
            shared_ptr<hier::PatchLevel> next_finer_level =
                    d_patch_hierarchy->getPatchLevel(level_num + 1);
            hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
            hier::IntVector coarsen_ratio =
                    next_finer_level->getRatioToCoarserLevel();
            coarsened_boxes.coarsen(coarsen_ratio);

            // Set control volume to 0 wherever there is a nonempty intersection
            // with the next finer level.
            for (hier::PatchLevel::Iterator pi(patch_level->begin());
                    pi!=patch_level->end(); pi++) {

                shared_ptr<hier::Patch> patch = *pi;
                shared_ptr< pdat::CellData<PQS_REAL> > control_volume_data;

                for (hier::BoxContainer::iterator ib = coarsened_boxes.begin();
                        ib != coarsened_boxes.end(); ++ib) {

                    hier::Box coarse_box = *ib;
                    hier::Box intersection = coarse_box*(patch->getBox());
                    if (!intersection.empty()) {
                        control_volume_data =
                            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                                    patch->getPatchData(d_control_volume_id));

                        control_volume_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

    } // loop over PatchLevels
} // TagInitAndDataTransferModule::computeControlVolumes()

// Copy constructor
TagInitAndDataTransferModule::TagInitAndDataTransferModule(
        const TagInitAndDataTransferModule& rhs):
    mesh::TagAndInitializeStrategy(s_object_name)
{
} // TagInitAndDataTransferModule::TagInitAndDataTransferModule()

} // PQS::pqs namespace
} // PQS namespace
