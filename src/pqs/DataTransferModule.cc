/*! \file DataTransferModule.cc
 *
 * \brief
 * Implementation file for DataTransferModule class.
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
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/Solver.h"
#include "PQS/utilities/error.h"
#include "PQS/pqs/DataTransferModule.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class CoarsenOperator; } }
namespace SAMRAI { namespace hier { class PatchLevel; } }
namespace SAMRAI { namespace hier { class RefineOperator; } }
namespace SAMRAI { namespace hier { class Variable; } }
namespace SAMRAI { namespace hier { class VariableContext; } }
namespace SAMRAI { namespace tbox { class Database; } }


// --- Class implementation

namespace PQS {
namespace pqs {


    // --- Public methods

// Constructor
DataTransferModule::DataTransferModule(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int phi_pqs_id,
        const int phi_lsm_current_id,
        const int phi_lsm_next_id,
        const int psi_id,
        const int max_stencil_width)
{
    // Check arguments
    if (config_db == NULL) {
        PQS_ERROR(this, "DataTransferModule",
                  "'config_db' must not be NULL");
    }
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "DataTransferModule",
                  "'patch_hierarchy' must not be NULL");
    }
    if (phi_pqs_id < 0) {
        PQS_ERROR(this, "DataTransferModule",
                  "'phi_pqs_id' must be non-negative");
    }
    if (phi_lsm_current_id < 0) {
        PQS_ERROR(this, "DataTransferModule",
                  "'phi_lsm_current_id' must be non-negative");
    }
    if (phi_lsm_next_id < 0) {
        PQS_ERROR(this, "DataTransferModule",
                  "'phi_lsm_next_id' must be non-negative");
    }
    if (psi_id < 0) {
        PQS_ERROR(this, "DataTransferModule",
                  "'psi_id' must be non-negative");
    }

    // Set data members
    d_phi_pqs_id = phi_pqs_id;
    d_phi_lsm_current_id = phi_lsm_current_id;
    d_phi_lsm_next_id = phi_lsm_next_id;
    d_psi_id = psi_id;
    d_patch_hierarchy = patch_hierarchy;

    // Load configuration parameters
    loadConfiguration(config_db);

    // Set up data transfer objects
    setupDataTransferObjects(patch_hierarchy->getGridGeometry(),
                             max_stencil_width);

} // DataTransferModule::DataTransferModule()

void DataTransferModule::fillGhostCells(
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
            level_num < d_patch_hierarchy->getNumberOfLevels(); level_num++) {

        if (context == LSM_CURRENT) {
            d_xfer_fill_bdry_schedules_lsm_current[level_num]->fillData(
                    0.0,    // not used
                    true);  // apply physical boundary conditions
        } else if (context == LSM_NEXT) {
            d_xfer_fill_bdry_schedules_lsm_next[level_num]->fillData(
                    0.0,    // not used
                    true);  // apply physical boundary conditions
        }
    }
} // DataTransferModule::fillGhostCells()

void DataTransferModule::enforcePhiConsistency(
        const int context) const
{
    // Check arguments
    if ( (context != PQS) && (context != LSM_NEXT) ) {
        PQS_ERROR(this, "enforcePhiConsistency",
                  string("Invalid 'context': ") +
                  to_string(context) +
                  string(". Valid values: PQS (=0), ") +
                  string("LSM_NEXT (=2)"));
    }

    // Enforce consistency of phi across PatchLevels
    for (int level_num = d_patch_hierarchy->getNumberOfLevels()-1;
            level_num > 0; level_num--) {

        if (context == PQS) {
            d_xfer_enforce_phi_consistency_schedules_pqs[level_num]->
                coarsenData();
        } else if (context == LSM_NEXT) {
            d_xfer_enforce_phi_consistency_schedules_lsm_next[level_num]->
                coarsenData();
        }
    }
} // DataTransferModule::enforcePhiConsistency()

void DataTransferModule::enforcePsiConsistency() const
{
    // Enforce consistency of phi across PatchLevels
    for (int level_num = d_patch_hierarchy->getNumberOfLevels()-1;
            level_num > 0; level_num--) {

        d_xfer_enforce_psi_consistency_schedules[level_num]->coarsenData();
    }
} // DataTransferModule::enforcePsiConsistency()

void DataTransferModule::resetHierarchyConfiguration(
        const int coarsest_level_num,
        const int finest_level_num)
{
    // --- Check arguments

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

        if (d_patch_hierarchy->getPatchLevel(level_num) == NULL) {
            PQS_ERROR(this, "resetHierarchyConfiguration",
                      string("PatchLevel ") +
                      to_string(level_num) +
                      string(" is NULL in PatchHierarchy"));
        }
    }

    // --- Reset data transfer schedules

    // Data transfer schedules for filling ghost cells data
    d_xfer_fill_bdry_schedules_lsm_current.resize(
        d_patch_hierarchy->getNumberOfLevels());
    d_xfer_fill_bdry_schedules_lsm_next.resize(
        d_patch_hierarchy->getNumberOfLevels());

    for (int level_num = coarsest_level_num;
            level_num <= finest_level_num; level_num++) {

        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        if (level_num == 0) {
            d_xfer_fill_bdry_schedules_lsm_current[level_num] =
                d_xfer_fill_bdry_lsm_current->createSchedule(
                    patch_level, NULL);  // TODO: change NULL to boundary
                                         // condition module

            d_xfer_fill_bdry_schedules_lsm_next[level_num] =
                d_xfer_fill_bdry_lsm_next->createSchedule(
                    patch_level, NULL);  // TODO: change NULL to boundary
                                         // condition module
        } else {
            d_xfer_fill_bdry_schedules_lsm_current[level_num] =
                d_xfer_fill_bdry_lsm_current->createSchedule(
                    patch_level, level_num-1,
                    d_patch_hierarchy, NULL);  // TODO: change NULL to boundary
                                             // condition module

            d_xfer_fill_bdry_schedules_lsm_next[level_num] =
                d_xfer_fill_bdry_lsm_next->createSchedule(
                    patch_level, level_num-1,
                    d_patch_hierarchy, NULL);  // TODO: change NULL to boundary
                                             // condition module
        }
    }

    // Data transfer schedules for enforcing phi and psi consistency
    d_xfer_enforce_phi_consistency_schedules_pqs.resize(
        d_patch_hierarchy->getNumberOfLevels());
    d_xfer_enforce_phi_consistency_schedules_lsm_next.resize(
        d_patch_hierarchy->getNumberOfLevels());
    d_xfer_enforce_psi_consistency_schedules.resize(
        d_patch_hierarchy->getNumberOfLevels());

    for (int level_num = coarsest_level_num;
            level_num <= finest_level_num; level_num++) {

        // No consistency to enforce for coarsest level in PatchHierarchy
        if (level_num == 0) {
            continue;
        }

        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        shared_ptr<hier::PatchLevel> next_coarser_patch_level =
            d_patch_hierarchy->getPatchLevel(level_num-1);

        d_xfer_enforce_phi_consistency_schedules_pqs[level_num] =
            d_xfer_enforce_phi_consistency_pqs->createSchedule(
                next_coarser_patch_level, patch_level);

        d_xfer_enforce_phi_consistency_schedules_lsm_next[level_num] =
            d_xfer_enforce_phi_consistency_lsm_next->createSchedule(
                next_coarser_patch_level, patch_level);

        d_xfer_enforce_psi_consistency_schedules[level_num] =
            d_xfer_enforce_psi_consistency->createSchedule(
                next_coarser_patch_level, patch_level);
    }

} // DataTransferModule::resetHierarchyConfiguration()

// --- Private methods

void DataTransferModule::loadConfiguration(
        const shared_ptr<tbox::Database>& config_db)
{
    // --- Load parameters
    // TODO
} // DataTransferModule::loadConfiguration()

void DataTransferModule::setupDataTransferObjects(
        const shared_ptr<hier::BaseGridGeometry>& grid_geometry,
        const int max_stencil_width)
{
    // --- Preparations

    // Get dimensionality of PatchHierarchy
    tbox::Dimension dim = d_patch_hierarchy->getDim();

    // Get pointer to VariableDatabase
    hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

    hier::IntVector scratch_ghost_cell_width(dim, max_stencil_width);

    // --- Create PatchData for data transfer scratch space

    // Get 'scratch' variable context
    shared_ptr<hier::VariableContext> scratch_context =
        var_db->getContext("scratch");

    // Get variable associated with d_phi_lsm_current_id
    shared_ptr< hier::Variable > phi_variable;
    if (!var_db->mapIndexToVariable(d_phi_lsm_current_id, phi_variable)) {
        PQS_ERROR(this, "setupDataTransferObjects",
                  "Error mapping 'd_phi_lsm_current_id' to Variable object.");
    }
    d_phi_scratch_id =
        var_db->registerVariableAndContext(phi_variable,
                                           scratch_context,
                                           scratch_ghost_cell_width);

    // Get variable associated with d_psi_id
    shared_ptr< hier::Variable > psi_variable;
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

    shared_ptr<hier::RefineOperator> phi_linear_refine_op =
        grid_geometry->lookupRefineOperator(phi_variable, "LINEAR_REFINE");

    shared_ptr<hier::RefineOperator> psi_linear_refine_op =
        grid_geometry->lookupRefineOperator(psi_variable, "LINEAR_REFINE");

    // --- Get coarsening operators for enforcing consistency across
    //     PatchLevels
    //
    //     TODO: replace with custom coarsening operator based on
    //           interpolation

    shared_ptr<hier::CoarsenOperator> coarsen_op =
        grid_geometry->lookupCoarsenOperator(phi_variable,
                                             "CONSERVATIVE_COARSEN");

    // --- Set up data transfer objects

    // Filling a ghost cells during time integration of level set functions
    d_xfer_fill_bdry_lsm_current =
        shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);
    d_xfer_fill_bdry_lsm_current->registerRefine(
        d_phi_lsm_current_id, d_phi_lsm_current_id, d_phi_scratch_id,
        phi_linear_refine_op);

    d_xfer_fill_bdry_lsm_next =
        shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);
    d_xfer_fill_bdry_lsm_next->registerRefine(
        d_phi_lsm_next_id, d_phi_lsm_next_id, d_phi_scratch_id,
        phi_linear_refine_op);

    // Enforcing consistency of phi across PatchLevels
    d_xfer_enforce_phi_consistency_pqs =
        shared_ptr<xfer::CoarsenAlgorithm>(new xfer::CoarsenAlgorithm(dim));
    d_xfer_enforce_phi_consistency_pqs->registerCoarsen(
        d_phi_pqs_id, d_phi_pqs_id, coarsen_op);

    d_xfer_enforce_phi_consistency_lsm_next =
        shared_ptr<xfer::CoarsenAlgorithm>(new xfer::CoarsenAlgorithm(dim));
    d_xfer_enforce_phi_consistency_lsm_next->registerCoarsen(
        d_phi_lsm_next_id, d_phi_lsm_next_id, coarsen_op);

    // Enforcing consistency of psi across PatchLevels
    d_xfer_enforce_psi_consistency =
        shared_ptr<xfer::CoarsenAlgorithm>(new xfer::CoarsenAlgorithm(dim));
    d_xfer_enforce_psi_consistency->registerCoarsen(
        d_psi_id, d_psi_id, coarsen_op);

} // DataTransferModule::setupDataTransferObjects()

} // PQS::pqs namespace
} // PQS namespace
