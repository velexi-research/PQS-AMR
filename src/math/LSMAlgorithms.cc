/*! \file LSMAlgorithms.cc
 *
 * \brief
 * Implementation file for LSM::Algorithms class.
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
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// SAMRAI
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/math/LSMAlgorithms.h"
#include "PQS/math/kernels/level_set_method_2d.h"
#include "PQS/math/kernels/level_set_method_3d.h"
#include "PQS/math/TimeIntegration.h"
#include "PQS/math/Toolbox.h"
#include "PQS/utilities/error.h"
#include "PQS/utilities/macros.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class RefineOperator; } }
namespace SAMRAI { namespace hier { class Variable; } }
namespace SAMRAI { namespace hier { class VariableContext; } }

// --- Class implementation

namespace PQS {
namespace math {
namespace LSM {

// --- Public methods

Algorithms::Algorithms(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const shared_ptr<hier::IntVector>& max_stencil_width)
{
    // Check arguments
    if (patch_hierarchy == NULL) {
        PQS_ERROR(this, "Algorithms", "'patch_hierarchy' must not be NULL");
    }

    // Set data members
    d_patch_hierarchy = patch_hierarchy;
    d_max_stencil_width = max_stencil_width;

    // Set up simulation variables
    setupSimulationVariables();

    // Set up data transfer objects
    setupDataTransferObjects(d_patch_hierarchy->getGridGeometry());

    // Emit contents of variable database to log file.
    tbox::plog << "\nVariable database..." << endl;
    hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

} // Algorithms::Algorithms()

Algorithms::~Algorithms()
{
    // Free memory allocated for simulation data
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels(); level_num++) {
        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);
        patch_level->deallocatePatchData(d_lsm_algs_variables);
    }
} // Algorithms::~Algorithms()


void Algorithms::resetHierarchyConfiguration(
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

    // data transfer schedules for filling ghost cells data
    d_xfer_fill_bdry_schedule_current.resize(
        d_patch_hierarchy->getNumberOfLevels());
    d_xfer_fill_bdry_schedule_next.resize(
        d_patch_hierarchy->getNumberOfLevels());

    for (int level_num = coarsest_level_num;
            level_num <= finest_level_num;
            level_num++) {

        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        if (level_num == 0) {
            d_xfer_fill_bdry_schedule_current[level_num] =
                d_xfer_fill_bdry_current->createSchedule(
                    patch_level, NULL);  // TODO: change NULL to boundary
                                         // condition module

            d_xfer_fill_bdry_schedule_next[level_num] =
                d_xfer_fill_bdry_next->createSchedule(
                    patch_level, NULL);  // TODO: change NULL to boundary
                                         // condition module
        } else {
            d_xfer_fill_bdry_schedule_current[level_num] =
                d_xfer_fill_bdry_current->createSchedule(
                    patch_level, level_num-1,
                    d_patch_hierarchy, NULL);  // TODO: change NULL to boundary
                                               // condition module

            d_xfer_fill_bdry_schedule_next[level_num] =
                d_xfer_fill_bdry_next->createSchedule(
                    patch_level, level_num-1,
                    d_patch_hierarchy, NULL);  // TODO: change NULL to boundary
                                               // condition module
        }
    }
} // Algorithms::resetHierarchyConfiguration()

void Algorithms::reinitializeLevelSetFunction(
        const int phi_id,
        const int control_volume_id,
        const int time_integration_order,
        const REINIT_ALGORITHM_TYPE algorithm_type,
        const int max_time_steps,
        const double steady_state_condition,
        const double stop_distance)
{
    // --- Check arguments

    // Get dimensionality of problem
    const int dim = d_patch_hierarchy->getDim().getValue();
    if ((dim != 2) && (dim != 3)) {
        PQS_ERROR(this, "equilibrateInterface",
                  string("Invalid number of spatial dimensions (=") +
                  to_string(dim) +
                  string("). Valid values: 2, 3."));
    }

    // TODO: Constraints on max_time_steps, steady_state_condition,
    //       stop_distance

    // --- Preparations

    // Create SAMRAI::math::HierarchyCellDataOpsReal object
    shared_ptr< SAMRAI::math::HierarchyCellDataOpsReal<PQS_REAL> > math_ops =
            shared_ptr< SAMRAI::math::HierarchyCellDataOpsReal<PQS_REAL> >(
                    new SAMRAI::math::HierarchyCellDataOpsReal<PQS_REAL>(
                             d_patch_hierarchy));

    // Set phi_0_id
    int phi_0_id = -1;
    if (algorithm_type == REINIT_EQN_SGN_PHI0) {
        phi_0_id = phi_id;
    }

    // Compute t_max
    double t_max = std::numeric_limits<double>::max();
    if (stop_distance > 0) {
        t_max = stop_distance;
    }

    // Initialize loop variables
    double t = 0.0;
    int step = 0;

    double delta_phi = steady_state_condition + 1;

    // Allocate PatchData
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels();
            level_num++) {

        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        patch_level->allocatePatchData(d_lsm_algs_variables);
    }

    // Copy phi data from phi_id PatchData to lsm_algs PatchData
    math_ops->copyData(d_lsm_algs_current_id, phi_id);

    // Compute time step = min(dx) / sqrt(dim)
    double min_dx = std::numeric_limits<double>::max();
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels();
            level_num++) {

        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(patch_level->begin());
                pi!=patch_level->end(); pi++) {

            shared_ptr<hier::Patch> patch = *pi;

            shared_ptr<geom::CartesianPatchGeometry> patch_geom =
                    SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                            patch->getPatchGeometry());

            const double* dx = patch_geom->getDx();

            for (int i = 0; i < dim; i++) {
                if (min_dx > dx[i]) {
                    min_dx = dx[i];
                }
            }
        }
    }

    const double dt = min_dx / sqrt(dim);

    // --- Perform level set method computation

    while ( (step < max_time_steps) && (t < t_max) &&
            (delta_phi > steady_state_condition) ) {

        // --- Preparations

        // Use TVD Runge-Kutta integration in time to compute phi(t+dt)
        for (int rk_stage = 1; rk_stage <= time_integration_order; rk_stage++)
        {
            // --- Preparations

            // Configuration for computation of RHS of evolution equation
            int phi_reinit_id;
            if (time_integration_order == 1) {
                phi_reinit_id = d_lsm_algs_current_id;
            } else if (time_integration_order == 2) {
                if (rk_stage == 1) {
                    phi_reinit_id = d_lsm_algs_current_id;
                } else if (rk_stage == 2) {
                    phi_reinit_id = d_lsm_algs_next_id;
                }
            } else if (time_integration_order == 3) {
                if (rk_stage == 1) {
                    phi_reinit_id = d_lsm_algs_current_id;
                } else if (rk_stage == 2) {
                    phi_reinit_id = d_lsm_algs_next_id;
                } else if (rk_stage == 3) {
                    phi_reinit_id = d_lsm_algs_next_id;
                }
            }

            int ghost_cell_fill_context;
            if (rk_stage == 1) {
                ghost_cell_fill_context = CURRENT;
            } else {
                ghost_cell_fill_context = NEXT;
            }

            // --- Fill ghost cells

            fillGhostCells(ghost_cell_fill_context);

            // --- Compute RHS of level set evolution equation

            for (int level_num = 0;
                    level_num < d_patch_hierarchy->getNumberOfLevels();
                    level_num++) {

                shared_ptr<hier::PatchLevel> patch_level =
                    d_patch_hierarchy->getPatchLevel(level_num);

                // Compute RHS of level set evolution equation
                for (hier::PatchLevel::Iterator pi(patch_level->begin());
                        pi!=patch_level->end(); pi++) {

                    computeReinitEqnRHSOnPatch(*pi,
                                               d_lsm_algs_rhs_id,
                                               phi_reinit_id,
                                               phi_0_id);
                }
            }

            // --- Advance phi

            if (time_integration_order == 1) {
                math::TimeIntegration::RK1Step(
                        d_patch_hierarchy,
                        d_lsm_algs_next_id,
                        d_lsm_algs_current_id,
                        d_lsm_algs_rhs_id,
                        dt);

            } else if (time_integration_order == 2) {
                if (rk_stage == 1) {
                    math::TimeIntegration::TVDRK2Stage1(
                            d_patch_hierarchy,
                            d_lsm_algs_next_id,
                            d_lsm_algs_current_id,
                            d_lsm_algs_rhs_id,
                            dt);

                } else if (rk_stage == 2) {
                    math::TimeIntegration::TVDRK2Stage2(
                            d_patch_hierarchy,
                            d_lsm_algs_next_id,
                            d_lsm_algs_next_id,
                            d_lsm_algs_current_id,
                            d_lsm_algs_rhs_id,
                            dt);
                }
            } else if (time_integration_order == 3) {
                if (rk_stage == 1) {
                    math::TimeIntegration::TVDRK3Stage1(
                            d_patch_hierarchy,
                            d_lsm_algs_next_id,
                            d_lsm_algs_current_id,
                            d_lsm_algs_rhs_id,
                            dt);

                } else if (rk_stage == 2) {
                    math::TimeIntegration::TVDRK3Stage2(
                            d_patch_hierarchy,
                            d_lsm_algs_next_id,
                            d_lsm_algs_next_id,
                            d_lsm_algs_current_id,
                            d_lsm_algs_rhs_id,
                            dt);

                } else if (rk_stage == 3) {
                    math::TimeIntegration::TVDRK3Stage3(
                            d_patch_hierarchy,
                            d_lsm_algs_next_id,
                            d_lsm_algs_next_id,
                            d_lsm_algs_current_id,
                            d_lsm_algs_rhs_id,
                            dt);
                }
            }
        }

        // --- Update metrics used in stopping criteria

        // Compute max norm of change in phi
        delta_phi = math::computeMaxNormDiff(d_patch_hierarchy,
                                             d_lsm_algs_next_id,
                                             d_lsm_algs_current_id,
                                             control_volume_id);

        // --- Prepare for next iteration

        // Update time and step
        t += dt;
        step++;

        // Swap phi data from LSM current and LSM next contexts
        math_ops->swapData(d_lsm_algs_next_id, d_lsm_algs_current_id);
    }

    // --- Clean up

    // Copy phi data from lsm_algs PatchData to phi_id PatchData
    math_ops->copyData(phi_id, d_lsm_algs_current_id);

    // Deallocate PatchData
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels();
            level_num++) {

        shared_ptr<hier::PatchLevel> patch_level =
            d_patch_hierarchy->getPatchLevel(level_num);

        patch_level->deallocatePatchData(d_lsm_algs_variables);
    }
} // Algorithms::reinitializeLevelSetFunctions()

void Algorithms::computeExtensionField(
        const int S_id,
        const int phi_id,
        const int control_volume_id,
        const int max_time_steps,
        const double steady_state_condition,
        const double stop_distance)
{
} // Algorithms::computeExtensionField()

shared_ptr<hier::PatchHierarchy> Algorithms::getPatchHierarchy() const
{
    return d_patch_hierarchy;
} // Algorithms::getPatchHierarchy()

void Algorithms::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::math::LSM::Algorithms::printClassData..." << endl;
    os << "(LSM::Algorithms*) this = " << (LSM::Algorithms*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;

} // Algorithms::printClassData()


// --- Private methods

void Algorithms::setupSimulationVariables()
{
    // --- Preparations

    // Get dimensionality of problem
    tbox::Dimension dim = d_patch_hierarchy->getDim();

    // Create IntVector for zero ghost cell widths
    hier::IntVector zero_ghost_cell_width(dim, 0);

    // Initialize PatchData component selectors
    d_lsm_algs_variables.clrAllFlags();

    // --- Create PatchData for simulation variables

    // Get pointer to VariableDatabase
    hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

    // Get variable contexts
    shared_ptr<hier::VariableContext> current_context =
        var_db->getContext("lsm_algs_current");
    shared_ptr<hier::VariableContext> next_context =
        var_db->getContext("lsm_algs_next");

    // variables for level set method algorithms
    shared_ptr< pdat::CellVariable<PQS_REAL> > lsm_algs_variable;
    if (var_db->checkVariableExists("lsm_algs")) {
        lsm_algs_variable =
                SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
                        var_db->getVariable("lsm_algs"));
    } else {
        const int depth = 1;
        lsm_algs_variable = shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "lsm_algs", depth));
    }
    d_lsm_algs_current_id =
        var_db->registerVariableAndContext(lsm_algs_variable,
                                           current_context,
                                           *d_max_stencil_width);
    d_lsm_algs_next_id =
        var_db->registerVariableAndContext(lsm_algs_variable,
                                           next_context,
                                           *d_max_stencil_width);
    d_lsm_algs_variables.setFlag(d_lsm_algs_current_id);
    d_lsm_algs_variables.setFlag(d_lsm_algs_next_id);

    // RHS of level set evolution equation
    shared_ptr< pdat::CellVariable<PQS_REAL> > rhs_variable;
    if (var_db->checkVariableExists("rhs")) {
        rhs_variable =
                SAMRAI_SHARED_PTR_CAST< pdat::CellVariable<PQS_REAL> >(
        var_db->getVariable("rhs"));
    } else {
        const int depth = 1;
        rhs_variable = shared_ptr< pdat::CellVariable<PQS_REAL> >(
                new pdat::CellVariable<PQS_REAL>(dim, "rhs", depth));
    }
    d_lsm_algs_rhs_id =
        var_db->registerVariableAndContext(rhs_variable,
                                           current_context,
                                           zero_ghost_cell_width);
    d_lsm_algs_variables.setFlag(d_lsm_algs_rhs_id);

} // Algorithms::setupSimulationVariables()

void Algorithms::setupDataTransferObjects(
        const shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{
    // --- Preparations

    // Get pointer to VariableDatabase
    hier::VariableDatabase *var_db = hier::VariableDatabase::getDatabase();

    hier::IntVector scratch_ghost_cell_width(*d_max_stencil_width);

    // Initialize PatchData component selectors
    d_scratch_variables.clrAllFlags();

    // --- Create PatchData for data transfer scratch space

    // Get 'scratch' variable context
    shared_ptr<hier::VariableContext> scratch_context =
        var_db->getContext("lsm_algs_scratch");

    // Get variable associated with d_lsm_algs_current_id
    shared_ptr< hier::Variable > lsm_algs_variable;
    if (!var_db->mapIndexToVariable(d_lsm_algs_current_id, lsm_algs_variable)) {
        PQS_ERROR(this, "setupDataTransferObjects",
                  "Error mapping 'd_lsm_algs_current_id' to Variable object.");
    }
    d_lsm_algs_scratch_id =
        var_db->registerVariableAndContext(lsm_algs_variable,
                                           scratch_context,
                                           scratch_ghost_cell_width);

    // --- Get refinement operators for filling PatchData from coarser
    //     PatchLevels

    shared_ptr<hier::RefineOperator> linear_refine_op =
        grid_geometry->lookupRefineOperator(lsm_algs_variable, "LINEAR_REFINE");

    // --- Set up data transfer objects

    // Filling a ghost cells during time integration of level set functions
    d_xfer_fill_bdry_current =
        shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);
    d_xfer_fill_bdry_current->registerRefine(
        d_lsm_algs_current_id, d_lsm_algs_current_id, d_lsm_algs_scratch_id,
        linear_refine_op);

    d_xfer_fill_bdry_next =
        shared_ptr<xfer::RefineAlgorithm>(new xfer::RefineAlgorithm);
    d_xfer_fill_bdry_next->registerRefine(
        d_lsm_algs_next_id, d_lsm_algs_next_id, d_lsm_algs_scratch_id,
        linear_refine_op);

} // Algorithms::setupDataTransferObjects()

void Algorithms::fillGhostCells(const int context) const
{
    // Check arguments
    if ( (context != CURRENT) && (context != NEXT) ) {
        PQS_ERROR(this, "fillGhostCells",
                  string("Invalid 'context': ") +
                  to_string(context) +
                  string(". Valid values: CURRENT (=1), NEXT (=2)"));
    }

    // Fill ghost cells
    for (int level_num = 0;
            level_num < d_patch_hierarchy->getNumberOfLevels();
            level_num++) {

        if (context == CURRENT) {
            d_xfer_fill_bdry_schedule_current[level_num]->fillData(
                    0.0,    // not used
                    true);  // apply physical boundary conditions
        } else if (context == NEXT) {
            d_xfer_fill_bdry_schedule_next[level_num]->fillData(
                    0.0,    // not used
                    true);  // apply physical boundary conditions
        }
    }
} // Algorithms::fillGhostCells()

void Algorithms::computeReinitEqnRHSOnPatch(
        const shared_ptr<hier::Patch>& patch,
        const int rhs_id,
        const int phi_id,
        const int phi_0_id)
{
    // --- Preparations

    // Get dimensionality of space
    const int dim = patch->getDim().getValue();

    // Get geometry parameters
    shared_ptr<geom::CartesianPatchGeometry> patch_geom =
        SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry>(
                    patch->getPatchGeometry());

    const double* dx = patch_geom->getDx();

    // --- Get pointers to data and index space ranges

    // RHS
    shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch->getPatchData(rhs_id));

    hier::Box rhs_ghostbox = rhs_data->getGhostBox();
    const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
    const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(rhs_ghostbox_lo, rhs_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(rhs_ghostbox_hi, rhs_ghostbox_upper);

    PQS_REAL* rhs = rhs_data->getPointer();

    // phi
    shared_ptr< pdat::CellData<PQS_REAL> > phi_data =
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                    patch->getPatchData(phi_id));

    hier::Box phi_ghostbox = phi_data->getGhostBox();
    const hier::IntVector phi_ghostbox_lower = phi_ghostbox.lower();
    const hier::IntVector phi_ghostbox_upper = phi_ghostbox.upper();
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_lo, phi_ghostbox_lower);
    PQS_INT_VECT_TO_INT_ARRAY(phi_ghostbox_hi, phi_ghostbox_upper);

    PQS_REAL* phi = phi_data->getPointer();

    // patch box
    hier::Box patch_box = patch->getBox();
    const hier::IntVector patch_box_lower = patch_box.lower();
    const hier::IntVector patch_box_upper = patch_box.upper();
    PQS_INT_VECT_TO_INT_ARRAY(patch_box_lo, patch_box_lower);
    PQS_INT_VECT_TO_INT_ARRAY(patch_box_hi, patch_box_upper);

    // Compute RHS of reinitialization equation on Patch

    if (phi_0_id > 0) {
        // phi_0
        shared_ptr< pdat::CellData<PQS_REAL> > phi_0_data =
                SAMRAI_SHARED_PTR_CAST< pdat::CellData<PQS_REAL> >(
                        patch->getPatchData(phi_0_id));
        hier::Box phi_0_ghostbox = phi_0_data->getGhostBox();
        const hier::IntVector phi_0_ghostbox_lower = phi_0_ghostbox.lower();
        const hier::IntVector phi_0_ghostbox_upper = phi_0_ghostbox.upper();
        PQS_INT_VECT_TO_INT_ARRAY(phi_0_ghostbox_lo, phi_0_ghostbox_lower);
        PQS_INT_VECT_TO_INT_ARRAY(phi_0_ghostbox_hi, phi_0_ghostbox_upper);

        PQS_REAL* phi_0 = phi_0_data->getPointer();

        if (dim == 2) {
            LSM_2D_COMPUTE_REINIT_EQN_SGN_PHI0_RHS(
                rhs, rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi, phi_ghostbox_lo, phi_ghostbox_hi,
                phi_0, phi_0_ghostbox_lo, phi_0_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx);
        } else if (dim == 3) {
            LSM_3D_COMPUTE_REINIT_EQN_SGN_PHI0_RHS(
                rhs, rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi, phi_ghostbox_lo, phi_ghostbox_hi,
                phi_0, phi_0_ghostbox_lo, phi_0_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx);
        }
    } else {
        if (dim == 2) {
            LSM_2D_COMPUTE_REINIT_EQN_SGN_PHI_RHS(
                rhs, rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi, phi_ghostbox_lo, phi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx);
        } else if (dim == 3) {
            LSM_3D_COMPUTE_REINIT_EQN_SGN_PHI_RHS(
                rhs, rhs_ghostbox_lo, rhs_ghostbox_hi,
                phi, phi_ghostbox_lo, phi_ghostbox_hi,
                patch_box_lo, patch_box_hi,
                dx);
        }
    }
} // Algorithms::computeReinitEqnRHSOnPatch()

} // PQS::math::LSM namespace
} // PQS::math namespace
} // PQS namespace
