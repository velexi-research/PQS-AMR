/*! \file TagInitAndDataTransferModule.h
 *
 * \brief
 * Header file for TagInitAndDataTransferModule class.
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

#ifndef INCLUDED_PQS_pqs_TagInitAndDataTransferModule_h
#define INCLUDED_PQS_pqs_TagInitAndDataTransferModule_h

/*! \class PQS::pqs::TagInitAndDataTransferModule
 *
 * \brief
 * TODO: add description
 *
 * <h3> NOTES </h3>
 *
 *  - TODO
 *
 * <h3> USAGE </h3>
 *
 * TODO
 *
 * <h4> Method (A) </h4>
 *
 * <h3> User-specified parameters (input database field) </h3>
 *
 * <h4> Sample Input File </h4>
 *
 *  <pre>
 *  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *  LevelSetMethodAlgorithm{
 *
 *    LevelSetFunctionIntegrator {
 *      start_time  = 0.0
 *      end_time    = 0.5
 *
 *      cfl_number               = 0.5
 *      spatial_derivative_type  = "WENO"
 *      spatial_derivative_order = 5
 *      tvd_runge_kutta_order    = 3
 *
 *      reinitialization_interval = 0
 *      reinitialization_max_iters = 20
 *      reinitialization_stop_dist = 0.2
 *      orthogonalization_interval = 0
 *      orthogonalization_max_iters = 20
 *      orthogonalization_stop_dist = 0.2
 *
 *      lower_bc_phi_0 = 1, 1, 1
 *      upper_bc_phi_0 = 1, 1, 1
 *
 *      use_AMR = FALSE
 *      refinement_cutoff_value = 0.25
 *      tag_buffer = 2,2,2,2,2,2
 *
 *      verbose = false
 *
 *    } // end of LevelSetFunctionIntegrator database
 *
 *
 *    LevelSetMethodGriddingAlgorithm {
 *      max_levels = 4
 *
 *      ratio_to_coarser {
 *         level_1            = 2, 2
 *         level_2            = 2, 2
 *         level_3            = 2, 2
 *      }
 *
 *      largest_patch_size {
 *        level_0 = 50,50
 *        level_1 = 100,100
 *        // all finer levels will use same values as level_1...
 *      }
 *
 *      tagging_method = "GRADIENT_DETECTOR","REFINE_BOXES"
 *
 *      RefineBoxes {
 *        level_0 = [(15,0),(29,14)]
 *        level_1 = [(65,10),(114,40)]
 *      }
 *
 *      LoadBalancer {
 *        // load balancer input parameters
 *      }
 *
 *    } // end of LevelSetMethodGriddingAlgorithm database
 *
 *  } // end of LevelSetMethodAlgorithm database
 *
 *  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *  </pre>
 *
 */

// --- Headers, namespaces, and type declarations

// Standard
#include <ostream>
#include <vector>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

// PQS
#include "PQS/PQS_config.h"
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace SAMRAI { namespace hier { class BaseGridGeometry; } }


// --- PQS::pqs::TagInitAndDataTransferModule Class

namespace PQS {
namespace pqs {

class TagInitAndDataTransferModule:
    public mesh::TagAndInitializeStrategy
{
public:

    //! @{

    /*!
     ************************************************************************
     *
     * @name Constructor and destructor
     *
     ************************************************************************/

    /*!
     * This constructor for TagInitAndDataTransferModule creates a TODO
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     *
     * patch_hierarchy: pointer to PatchHierarchy that object manages
     *
     * pore_init_strategy: pointer to PoreInitStrategy object that
     *      implements algorithm for initializing level set function that
     *      defines solid-pore interface
     *
     * interface_init_strategy: pointer to PoreInitStrategy object that
     *      implements algorithm for initializing level set function that
     *      defines fluid-fluid interface
     *
     * phi_pqs_id: PatchData ID for the level set function for the
     *      fluid-fluid interface (at each PQS step)
     *
     * phi_lsm_current_id: PatchData ID for the level set function for the
     *      fluid-fluid interface (before each time step during evolution
     *      of the level set function)
     *
     * phi_lsm_next_id: PatchData ID for the level set function for the
     *      fluid-fluid interface (during and after each time step during
     *      evolution of the level set function)
     *
     * psi_id: PatchData ID for the level set function for the solid-pore
     *      interface
     *
     * max_stencil_width: maximum stencil width required for computations
     *
     */
    TagInitAndDataTransferModule(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy,
        const int phi_pqs_id,
        const int phi_lsm_current_id,
        const int phi_lsm_next_id,
        const int phi_id,
        const hier::IntVector& max_stencil_width);

    /*!
     * Default destructor frees memory allocated for data transfer scratch
     * space.
     */
    virtual ~TagInitAndDataTransferModule();

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Data transfer methods
     *
     ************************************************************************/

    /*!
     * Fill ghostcells for Patches in PatchLevel.
     *
     * Parameters
     * ----------
     * level_num: number of PatchLevel to fill ghost cells for
     *
     * context: variable context that ghost cell data should be filled for
     *
     * Return value
     * ------------
     * None
     */
    virtual void fillGhostCells(const int level_num, const int context) const;

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Methods inherited from TagAndInitializeStrategy class
     *
     ************************************************************************/

    /*!
     * Initialize data on a specified PatchLevel of PatchHierarchy.
     *
     * Parameters
     * ----------
     * patch_hierarchy: PatchHierarchy for AMR computation
     *
     * level_num: number of PatchLevel to initialize data for
     *
     * init_data_time: current simulation time
     *
     * can_be_refined: [unused] flag that indicates whether the PatchLevel can
     *      be further refined.
     *
     * initial_time: flag that indicates whether 'init_data_time' is equal to
     *      the initial simulation time. This information is useful if data
     *      values are initialized differently at the beginning of the
     *      simulation.
     *
     * old_patch_level: old PatchLevel in the PatchHierarchy at the specified
     *      'level_num'. If 'old_patch_level' is null, there was no
     *      PatchLevel in the PatchHierarchy prior to the call to
     *      'initializeLevelData()', so the data on the new PatchLevel is set
     *      by interpolating data from coarser PatchLevels in the
     *      PatchHierarchy. If 'old_patch_level' is not null, then the new
     *      Patchlevel is initialized by interpolating data from coarser
     *      PatchLevels and copying data from the old PatchLevel before it is
     *      destroyed.
     *
     * allocate_data: flag that indicates whether memory for the new PatchLevel
     *      at the specified 'level_num' must be allocated before
     *      data is initialized.
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - Application-specific pore initialization must be provided via a
     *   concrete implementation of the pure abstract 'initializePoreSpace()'
     *   method declared in the 'PQS::pqs::PoreInitStrategy' class.
     *
     * - Application-specific interface initialization must be provided via a
     *   concrete implementation of the pure abstract 'initializeInterface()'
     *   method declared in the 'PQS::pqs::InterfaceInitStrategy' class.
     *
     * - This method is invoked by GriddingAlgorithm objects when they insert
     *   new PatchLevels into PatchHierarchy objects.
     */
    virtual void initializeLevelData(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int level_num,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const shared_ptr<hier::PatchLevel>& old_patch_level =
            shared_ptr<hier::PatchLevel>(),
        const bool allocate_data = true);

    /*!
     * Reset any internal information that depends on the PatchHierarchy
     * or the data that resides on it.
     *
     * For the PQS algorithm, 'resetHierarchyConfiguration()' updates the
     * following:
     *
     * - communication schedules for PatchLevels in the PatchHierarchy
     *
     * - control volume values (TODO).
     *
     * Parameters
     * ----------
     * patch_hierarchy: PatchHierarchy for AMR computation
     *
     * coarsest_level_num: number of coarsest PatchLevel to reset
     *
     * finest_level_num: number of finest PatchLevel to reset
     *
     * Return value
     * ------------
     * None
     *
     */
    virtual void resetHierarchyConfiguration(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int coarsest_level_num,
        const int finest_level_num);

    /*!
     * Tag cells that should be refined on the specified PatchLevel.
     *
     * Parameters
     * ----------
     * patch_hierarchy: PatchHierarchy for AMR computation
     *
     * level_num: number of PatchLevel to tag cells on
     *
     * regrid_cycle: [unused] current cycle number
     *
     * regrid_time: current simulation time
     *
     * tag_id: PatchData ID containing the tag data
     *
     * initial_time: flag that indicates whether current call to
     *      'tagCellsForRefinement()' is at the initial simulation time.
     *
     * coarsest_sync_level_num: [unused] flag that indicates that
     *      the current PatchLevel is the coarsest PatchLevel in the
     *      PatchHierarchy in the current regridding process
     *
     * can_be_refined: [unused] flag that indicates whether the PatchLevel can
     *      be further refined.
     *
     * regrid_start_time: [unused] simulation time at the beginning of the
     *      regridding process.
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - Cells where refinement should occur are tagged by setting the integer
     *   PatchData associated with 'tag_id' to 1.
     */
    virtual void tagCellsForRefinement(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int level_num,
        const int regrid_cycle,
        const double regrid_time,
        const int tag_id,
        const bool initial_time,
        const bool coarsest_sync_patch_level,
        const bool can_be_refined = true,
        const double regrid_start_time = 0.);

    /*!
     * Return whether refinement is being performed using ONLY user-supplied
     * refine boxes. If any method is used that invokes tagging, this will
     * return false.
     *
     * Parameters
     * ----------
     * cycle: [unused] current cycle number
     *
     * time: [unused] current simulation time
     *
     * Return value
     * ------------
     * TODO
     */
    virtual bool refineUserBoxInputOnly(int cycle, double time);

    /*!
     * Set refine boxes to user supplied boxes for specified PatchLevel number.
     *
     * Parameters
     * ----------
     * refine_boxes: [output] refine boxes
     *
     * level_num: number of PatchLevel to set refine boxes for
     *
     * cycle: [unused] current cycle number
     *
     * time: [unused] current simulation time
     *
     * Return value
     * ------------
     * true the first time that getUserSuppliedRefineBoxes() is called; False
     * otherwise.
     *
     * Notes
     * -----
     * - The boolean return value specifies whether the boxes have been reset
     *   from the last time this method was called.
     */
    virtual bool getUserSuppliedRefineBoxes(
        hier::BoxContainer& refine_boxes,
        const int level_num,
        const int cycle,
        const double time);

    /*!
     * Reset the static refine boxes for the specified PatchLevel number in the
     * PatchHierarchy.
     *
     * Parameters
     * ----------
     * refine_boxes: [output] refine boxes
     *
     * level_num: number of PatchLevel to reset refine boxes for
     *
     * Notes
     * -----
     * - The PatchLevel number must be greater than or equal to zero.
     *
     * - Implementation is based on the resetRefineBoxes() method of the
     *   SAMRAI::mesh::StandardTagAndInitialize class.
     */
    virtual void resetRefineBoxes(
        const hier::BoxContainer& refine_boxes,
        const int level_num);

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Unused methods inherited from TagAndInitializeStrategy class
     *
     ************************************************************************/

    /*!
     * No-op because this class does not use error estimation to tag cells
     * for refinement.
     */
    virtual void preprocessErrorEstimation(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int level_num,
        const int cycle,
        const double regrid_time,
        const double regrid_start_time,
        const bool initial_time) {}

    /*!
     * Return false because grid management does not use time integration.
     *
     * Parameters
     * ----------
     * cycle: [unused] current cycle number
     *
     * time: [unused] current simulation time
     *
     * Return value
     * ------------
     * false
     */
    virtual bool usesTimeIntegration(int cycle, double time) {
        return false;
    }

    /*!
     * Return false because grid management never uses time integration.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * false
     */
    virtual bool everUsesTimeIntegration() const {
        return false;
    }

    /*!
     * Return true because all configurations of the coarsest PatchLevel
     * in the PatchHierarchy are acceptable for gridding.
     *
     * Parameters
     * ----------
     * boxes: [unused] boxes on coarsest PatchLevel of PatchHierarchy
     *
     * Return value
     * ------------
     * true
     */
    virtual bool coarsestLevelBoxesOK(const hier::BoxContainer& boxes) const {
        return true;
    }

    /*!
     * Return 1 because this class does not use error estimation to tag cells
     * for refinement.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * 1
     */
    virtual int getErrorCoarsenRatio() const {
        return 1;
    }

    /*!
     * No-op because this class does not use error estimation to tag cells
     * for refinement.
     */
    virtual void checkCoarsenRatios(
        const vector<hier::IntVector>& ratio_to_coarser) {}

    /*!
     * No-op because no special processing is required when swapping old and
     * new PatchLevels during regridding.
     */
    virtual void processHierarchyBeforeAddingNewLevel(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int level_num,
        const shared_ptr<hier::BoxLevel>& new_box_level) {}

    /*!
     * No-op because no special processing is required before a PatchLevel is
     * removed from the PatchHierarchy during regrid.
     */
    virtual void processLevelBeforeRemoval(
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int level_num,
        const shared_ptr<hier::PatchLevel>& old_patch_level =
            shared_ptr<hier::PatchLevel>()) {}

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Accessor methods for object parameters and state
     *
     ************************************************************************/

    /*!
     * Print the values of the data members the object.
     *
     * Parameters
     * ----------
     * os: output stream to write object information to
     *
     * Return value
     * ------------
     * None
     */
    virtual void printClassData(ostream& os) const;

    //! @}

protected:

    /****************************************************************
     *
     * Data Members
     *
     ****************************************************************/

    // --- Parameters

    // --- PatchData IDs

    // PatchData component selectors to organize variables by
    // data management cycle requirements
    hier::ComponentSelector d_scratch_variables;

    // fluid-fluid interface level set function
    int d_phi_pqs_id;
    int d_phi_lsm_current_id;
    int d_phi_lsm_next_id;
    int d_phi_scratch_id;

    // solid-pore interface level set function
    //
    // Note: the solid phase is defined by the region where psi < 0
    int d_psi_id;
    int d_psi_scratch_id;

    // --- Components

    // SAMRAI components
    shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;

    // Pore initialization
    shared_ptr<pqs::PoreInitStrategy> d_pore_init_strategy;

    // Interface initialization
    shared_ptr<pqs::InterfaceInitStrategy> d_interface_init_strategy;

    // Data transfer
    shared_ptr<xfer::RefineAlgorithm> d_xfer_fill_new_level;

    shared_ptr<xfer::RefineAlgorithm> d_xfer_fill_bdry_lsm_current;
    vector< shared_ptr<xfer::RefineSchedule> >
            d_xfer_fill_bdry_schedule_lsm_current;

    shared_ptr<xfer::RefineAlgorithm> d_xfer_fill_bdry_lsm_next;
    vector< shared_ptr<xfer::RefineSchedule> >
            d_xfer_fill_bdry_schedule_lsm_next;

    // --- Object name
    //
    // Note: used only to initialize TagAndInitializeStrategy base class
    const static string s_object_name;

private:

    /*
     * Load configuration parameters from specified database.
     */
    void loadConfiguration(const shared_ptr<tbox::Database>& config_db);

    /*
     * Set up data transfer objects and associated scratch space variables.
     *
     * Parameters
     * ----------
     * grid_geometry: grid geometry that defines refinement operators
     *
     * max_stencil_width: maximum stencil width required for computations
     */
    void setupDataTransferObjects(
        const shared_ptr<hier::BaseGridGeometry>& grid_geometry,
        const hier::IntVector& max_stencil_width);

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     */
    TagInitAndDataTransferModule(const TagInitAndDataTransferModule& rhs);

    /*
     * Private assignment operator to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object on right-hand side of assignment operator
     *
     * Return value
     * ------------
     * return object
     */
    const TagInitAndDataTransferModule& operator=(
            const TagInitAndDataTransferModule& rhs) {
        return *this;
    }

};  // PQS::pqs::TagInitAndDataTransferModule class

}  // PQS::pqs namespace
}  // PQS namespace

#endif
