/*! \file TagAndInitializeModule.cc
 *
 * \brief
 * Implementation file for TagAndInitializeModule class.
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

// SAMRAI
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/math/BoundaryConditions.h"
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"
#include "PQS/pqs/Solver.h"
#include "PQS/pqs/TagAndInitializeModule.h"
#include "PQS/pqs/kernels/kernels_2d.h"
#include "PQS/pqs/kernels/kernels_3d.h"
#include "PQS/utilities/error.h"
#include "PQS/utilities/macros.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class CoarsenOperator; } }
namespace SAMRAI { namespace hier { class RefineOperator; } }
namespace SAMRAI { namespace hier { class Variable; } }
namespace SAMRAI { namespace hier { class VariableContext; } }


// --- Class implementation

namespace PQS {
namespace pqs {

// --- Static data members

const string TagAndInitializeModule::s_object_name =
    "PQS::pqs::TagAndInitializeModule";

// --- Public methods

// Constructor
TagAndInitializeModule::TagAndInitializeModule(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        pqs::Solver* pqs_solver,
        const shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy,
        const int phi_pqs_id,
        const int phi_lsm_current_id,
        const int phi_lsm_next_id,
        const int psi_id,
        const int control_volume_id,
        const int max_stencil_width):
    mesh::TagAndInitializeStrategy(s_object_name)
{
    // Check arguments
    if (config_db == NULL) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'config_db' must not be NULL");
    }
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'patch_hierarchy' must not be NULL");
    }
    if (pqs_solver == NULL) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'pqs_solver' must not be NULL");
    }
    if (pore_init_strategy == NULL) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'pore_init_strategy' must not be NULL");
    }
    if (interface_init_strategy == NULL) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'interface_init_strategy' must not be NULL");
    }
    if (phi_pqs_id < 0) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'phi_pqs_id' must be non-negative");
    }
    if (phi_lsm_current_id < 0) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'phi_lsm_current_id' must be non-negative");
    }
    if (phi_lsm_next_id < 0) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'phi_lsm_next_id' must be non-negative");
    }
    if (psi_id < 0) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'psi_id' must be non-negative");
    }
    if (control_volume_id < 0) {
        PQS_ERROR(this, "TagAndInitializeModule",
                  "'control_volume_id' must be non-negative");
    }

    // Set data members
    d_max_stencil_width = max_stencil_width;

    d_phi_pqs_id = phi_pqs_id;
    d_phi_lsm_current_id = phi_lsm_current_id;
    d_phi_lsm_next_id = phi_lsm_next_id;
    d_psi_id = psi_id;
    d_control_volume_id = control_volume_id;
    d_patch_hierarchy = patch_hierarchy;
    d_pqs_solver = pqs_solver;
    d_pore_init_strategy = pore_init_strategy;
    d_interface_init_strategy = interface_init_strategy;

    // Load configuration parameters
    loadConfiguration(config_db);

    // Set up data transfer objects
    setupDataTransferObjects(patch_hierarchy->getGridGeometry(),
                             max_stencil_width);

} // TagAndInitializeModule::TagAndInitializeModule()

void TagAndInitializeModule::initializeLevelData(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int level_num,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const shared_ptr<hier::PatchLevel>& old_patch_level,
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
    shared_ptr<hier::PatchLevel> patch_level =
        patch_hierarchy->getPatchLevel(level_num);

    // Allocate PatchData
    if (allocate_data) {
        patch_level->allocatePatchData(d_phi_pqs_id);
        patch_level->allocatePatchData(d_psi_id);
        patch_level->allocatePatchData(d_control_volume_id);
    }

    // Set time for PatchData
    patch_level->setTime(init_data_time, d_phi_pqs_id);
    patch_level->setTime(init_data_time, d_psi_id);
    patch_level->setTime(init_data_time, d_control_volume_id);

    if (initial_time) {
        // Initialize data on Patches on the PatchLevel.
        for (hier::PatchLevel::Iterator pi(patch_level->begin());
                pi!=patch_level->end(); pi++) {

            shared_ptr<hier::Patch> patch = *pi;
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
        // Fill new PatchLevel with data from the old PatchLevel and coarser
        // levels in the PatchHierarchy
        shared_ptr<xfer::RefineSchedule> sched =
            d_xfer_fill_new_level->createSchedule(
            patch_level, old_patch_level, level_num-1,
            patch_hierarchy,
            this);

        sched->fillData(init_data_time, true);
    }

    // deallocate data on old PatchLevel if it exists
    if (old_patch_level!=NULL) {
        old_patch_level->deallocatePatchData(d_phi_pqs_id);
        old_patch_level->deallocatePatchData(d_psi_id);
        old_patch_level->deallocatePatchData(d_control_volume_id);
    }
} // TagAndInitializeModule::initializeLevelData()

void TagAndInitializeModule::resetHierarchyConfiguration(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
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

    // --- Recompute simulation parameters

    // recompute control volumes
    computeControlVolumes();

    // --- Call pqs::Solver::resetHierarchyConfiguration()

    d_pqs_solver->resetHierarchyConfiguration(
            coarsest_level_num, finest_level_num);

} // TagAndInitializeModule::resetHierarchyConfiguration()

void TagAndInitializeModule::tagCellsForRefinement(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
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

    // --- Preparations

    // Get problem dimension
    const int dim = d_patch_hierarchy->getDim().getValue();

    // Get PatchLevel
    const shared_ptr<hier::PatchLevel> patch_level =
        patch_hierarchy->getPatchLevel(level_num);

    // --- Set tag on Patches on the PatchLevel

    for (hier::PatchLevel::Iterator pi(patch_level->begin());
            pi!=patch_level->end(); pi++) {

        shared_ptr<hier::Patch> patch = *pi;
        if (patch == NULL) {
        PQS_ERROR(this, "initializeLevelData",
                  "Null Patch pointer found when iterating over "
                  "PatchLevel.");
        }

        // --- Compute refinement cutoff

        // Get geometry parameters
        shared_ptr<geom::CartesianPatchGeometry> patch_geom =
                SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                        patch->getPatchGeometry());
        const double* dx = patch_geom->getDx();

        // Compute maximum grid spacing
        double max_dx = dx[0];
        if (max_dx < dx[1]) {
            max_dx = dx[1];
        }
        if (dim == 3) {
            if (max_dx < dx[2]) {
                max_dx = dx[2];
            }
        }

        // Compute refinement_threshold
        const double refinement_cutoff =
                d_refinement_cutoff_multiplier * max_dx;

        // --- Get pointers to data and index space ranges

        // refinement tags
        shared_ptr< pdat::CellData<int> > tag_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<int> >(
                patch->getPatchData(tag_id));

        hier::Box tag_ghostbox = tag_data->getGhostBox();
        const hier::IntVector tag_ghostbox_lower = tag_ghostbox.lower();
        const hier::IntVector tag_ghostbox_upper = tag_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(tag_ghostbox_lo, tag_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(tag_ghostbox_hi, tag_ghostbox_upper);

        int* tag = tag_data->getPointer();

        // phi
        shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
                SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData(d_phi_pqs_id));

        hier::Box phi_ghostbox = phi_data->getGhostBox();
        const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
        const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_lo, phi_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_hi, phi_ghostbox_upper);

        PQS_REAL* phi = phi_data->getPointer();

        // psi
        shared_ptr< pdat::CellData<PQS_REAL> > psi_data =
                SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData(d_psi_id));

        hier::Box psi_ghostbox = psi_data->getGhostBox();
        const hier::IntVector psi_ghostbox_lower = psi_ghostbox.lower();
        const hier::IntVector psi_ghostbox_upper = psi_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_lo, psi_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(psi_ghostbox_hi, psi_ghostbox_upper);

        PQS_REAL* psi = psi_data->getPointer();

        // patch box
        hier::Box patch_box = patch->getBox();
        const hier::IntVector patch_box_lower = patch_box.lower();
        const hier::IntVector patch_box_upper = patch_box.upper();
        PQS_INT_VECT_TO_INT_ARRAY(patch_box_lo, patch_box_lower);
        PQS_INT_VECT_TO_INT_ARRAY(patch_box_hi, patch_box_upper);

        // --- Tag cells for refinement

        if (dim == 2) {
            PQS_2D_TAG_CELLS_FOR_REFINEMENT(
                tag, tag_ghostbox_lo, tag_ghostbox_hi,
                phi, phi_ghostbox_lo, phi_ghostbox_hi,
                psi, psi_ghostbox_lo, psi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                &refinement_cutoff);

        } else if (dim == 3) {
            PQS_3D_TAG_CELLS_FOR_REFINEMENT(
                tag, tag_ghostbox_lo, tag_ghostbox_hi,
                phi, phi_ghostbox_lo, phi_ghostbox_hi,
                psi, psi_ghostbox_lo, psi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                &refinement_cutoff);
        }
    }
} // TagAndInitializeModule::tagCellsForRefinement()

bool TagAndInitializeModule::refineUserBoxInputOnly(
        int cycle, double time)
{
    // TODO
    return false;
} // TagAndInitializeModule::refineUserBoxInputOnly()

bool TagAndInitializeModule::getUserSuppliedRefineBoxes(
        hier::BoxContainer& refine_boxes,
        const int level_num,
        const int cycle,
        const double time)
{
    // TODO
    return false;
} // TagAndInitializeModule::getUserSuppliedRefineBoxes()

void TagAndInitializeModule::resetRefineBoxes(
        const hier::BoxContainer& refine_boxes,
        const int level_num)
{
} // TagAndInitializeModule::resetRefineBoxes()

void TagAndInitializeModule::setPhysicalBoundaryConditions(
        hier::Patch& patch,
        const double fill_time,
        const hier::IntVector& ghost_width_to_fill)
{
    // Fill boundary data for level set functions by using linear
    // extrapolation to set values
    PQS::math::fillBdryDataLinearExtrapolation(patch, d_phi_scratch_id);
    PQS::math::fillBdryDataLinearExtrapolation(patch, d_psi_scratch_id);

} // TagAndInitializeModule::setPhysicalBoundaryConditions()

hier::IntVector TagAndInitializeModule::getRefineOpStencilWidth(
        const tbox::Dimension& dim) const {

    return hier::IntVector(dim, d_max_stencil_width);

} // TagAndInitializeModule::getRefineOpStencilWidth()

void TagAndInitializeModule::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::TagAndInitializeModule::printClassData..." << endl;
    os << "(TagAndInitializeModule*) this = "
       << (TagAndInitializeModule*) this << endl;
    os << "d_pore_init_strategy = " << d_pore_init_strategy.get() << endl;
    os << "d_interface_init_strategy = " << d_interface_init_strategy.get()
       << endl;

    os << endl;
    d_pore_init_strategy->printClassData(os);
    d_interface_init_strategy->printClassData(os);
} // TagAndInitializeModule::printClassData()


// --- Private methods

void TagAndInitializeModule::loadConfiguration(
        const shared_ptr<tbox::Database>& config_db)
{
    // --- Load parameters

    // AMR parameters
    if (config_db->keyExists("refinement_cutoff_multiplier")) {
        d_refinement_cutoff_multiplier =
                config_db->getInteger("refinement_cutoff_multiplier");
    } else {
        d_refinement_cutoff_multiplier = 5 * d_max_stencil_width;
    }

    if (d_refinement_cutoff_multiplier < d_max_stencil_width) {
        PQS_ERROR(this, "loadConfiguration",
                  string("'d_refinement_cutoff_multiplier' (=") +
                  to_string(d_refinement_cutoff_multiplier) +
                  string(") less than 'max_stencil_width' (=") +
                  to_string(d_max_stencil_width) +
                  string(")"));
    }
} // TagAndInitializeModule::loadConfiguration()

void TagAndInitializeModule::setupDataTransferObjects(
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

    // Filling a new level (used during initialization of a PatchLevel)
    d_xfer_fill_new_level =
        shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);

    d_xfer_fill_new_level->registerRefine(
        d_phi_pqs_id, d_phi_pqs_id, d_phi_scratch_id, phi_linear_refine_op);
    d_xfer_fill_new_level->registerRefine(
        d_psi_id, d_psi_id, d_psi_scratch_id, psi_linear_refine_op);

} // TagAndInitializeModule::setupDataTransferObjects()

void TagAndInitializeModule::computeControlVolumes() const
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
} // TagAndInitializeModule::computeControlVolumes()

// Copy constructor
TagAndInitializeModule::TagAndInitializeModule(
        const TagAndInitializeModule& rhs):
    mesh::TagAndInitializeStrategy(s_object_name)
{
} // TagAndInitializeModule::TagAndInitializeModule()

} // PQS::pqs namespace
} // PQS namespace
