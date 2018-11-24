/*! \file TagAndInitializeModule.h
 *
 * \brief
 * Header file for TagAndInitializeModule class.
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

#ifndef INCLUDED_PQS_pqs_TagAndInitializeModule_h
#define INCLUDED_PQS_pqs_TagAndInitializeModule_h

/*! \class PQS::pqs::TagAndInitializeModule
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
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

// PQS
#include "PQS/PQS_config.h"
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"
#include "PQS/pqs/Solver.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace SAMRAI { namespace hier { class BaseGridGeometry; } }


// --- PQS::pqs::TagAndInitializeModule Class

namespace PQS {
namespace pqs {

class TagAndInitializeModule:
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
     * This constructor for TagAndInitializeModule creates a TODO
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     *
     * patch_hierarchy: pointer to PatchHierarchy that object manages
     *
     * pqs_solver: pointer to pqs::Solver object
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
     * control_volume_id: PatchData ID for the control volume data
     *
     * max_stencil_width: maximum stencil width required for computations
     *
     */
    TagAndInitializeModule(
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
        const int max_stencil_width);

    /*!
     * Default destructor frees memory allocated for data transfer scratch
     * space.
     */
    virtual ~TagAndInitializeModule();

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
     * - control volume values.
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
     * For the PQS algorithm, we tag cells where the level set function
     * for either the fluid-fluid or pore-solid interface are below a
     *
     *     max{dx} * refinement_cutoff_multiplier
     *
     * where refinement_cutoff_multiplier can be specified as a configuration
     * parameter.
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

    // AMR parameters
    int d_refinement_cutoff_multiplier;  // cutoff value (in units of the
                                         // local grid spacing) to use when
                                         // tagging grid cells around the zero
                                         // level set of phi to refine.
                                         //
                                         // Note: should be set to be larger
                                         // that the maximum stencil width
                                         // required for the computation
                                         //
                                         // Default: 5 * maximum stencil width

    // --- PatchData IDs

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

    // control volume
    int d_control_volume_id;

    // --- Components

    // SAMRAI components
    shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;

    // PQS solver
    //
    // NOTE: pqs_solver is not a shared pointer because
    //       pqs::TagAndInitializeModule objects are created by
    //       pqs::Solver objects
    pqs::Solver* d_pqs_solver;

    // Pore initialization
    shared_ptr<pqs::PoreInitStrategy> d_pore_init_strategy;

    // Interface initialization
    shared_ptr<pqs::InterfaceInitStrategy> d_interface_init_strategy;

    // Data transfer: fill new level
    shared_ptr<xfer::RefineAlgorithm> d_xfer_fill_new_level;

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
        const int max_stencil_width);

    /*
     * Compute control volumes for grid cells.
     */
    void computeControlVolumes() const;

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     */
    TagAndInitializeModule(const TagAndInitializeModule& rhs);

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
    const TagAndInitializeModule& operator=(
            const TagAndInitializeModule& rhs) {
        return *this;
    }

};  // PQS::pqs::TagAndInitializeModule class

}  // PQS::pqs namespace
}  // PQS namespace

#endif // INCLUDED_PQS_pqs_TagAndInitializeModule_h
