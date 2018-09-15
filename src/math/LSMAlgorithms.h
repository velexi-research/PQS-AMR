/*! \file LSMAlgorithms.h
 *
 * \brief
 * Header file for LSMAlgorithms class.
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

#ifndef INCLUDED_PQS_math_LSMAlgorithms_h
#define INCLUDED_PQS_math_LSMAlgorithms_h

/*! \class PQS::math::LSM::Algorithms
 *
 * \brief
 * TODO: add description
 *
 * Notes
 * -----
 * - All level set method algorithms assume cell-centered PatchData.
 */

// --- Headers, namespaces, and type declarations

// Standard library
#include <ostream>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

// PQS
#include "PQS/PQS_config.h"
#include "PQS/pqs/Algorithms.h"
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations


// --- PQS::math::LSM::Algorithms Class

namespace PQS {
namespace math {
namespace LSM {

/*! \enum SPATIAL_DERIVATIVE_TYPE
 *
 * Enumerated type for the methods of computing spatial derivatives.
 *
 */
typedef enum {
    ENO1 = 1,
    ENO2 = 2,
    ENO3 = 3,
    WENO5 = 5,
    UNKNOWN = -1
} SPATIAL_DERIVATIVE_TYPE;

/*! \enum VARIABLE_CONTEXT
 *
 * Enumerated type for the contexts of simulation variables.
 *
 */
typedef enum {
    CURRENT = 1,
    NEXT = 2,
} VARIABLE_CONTEXT;


class Algorithms {

public:

    //! @{

    /*!
     ************************************************************************
     *
     * @name Constructor and destructor
     *
     ************************************************************************/

    /*!
     * This constructor for LSM::Algorithms creates a TODO
     *
     * - Create simulation variables.
     * - Register PatchData
     *
     * Parameters
     * ----------
     * patch_hierarchy: PatchHierarchy to use for simulation
     *
     * max_stencil_width: TODO
     */
    Algorithms(const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
               const shared_ptr<hier::IntVector>& max_stencil_width);

    /*!
     * Default destructor frees memory allocated for simulation.
     */
    virtual ~Algorithms();

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Algorithm state management
     *
     ************************************************************************/

    /*!
     * Reset any internal information that depends on the PatchHierarchy
     * or the data that resides on it.
     *
     * For LSAMAlgorithms, 'resetHierarchyConfiguration()' updates the
     * following:
     *
     * - communication schedules for PatchLevels in the PatchHierarchy.
     *
     * Parameters
     * ----------
     * coarsest_level_num: number of coarsest PatchLevel to reset
     *
     * finest_level_num: number of finest PatchLevel to reset
     *
     * Return value
     * ------------
     * None
     */
    virtual void resetHierarchyConfiguration(
            const int coarsest_level_num,
            const int finest_level_num);
    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Level Set Method algorithms
     *
     ************************************************************************/

    /*!
     * Reinitialize level set function to a signed distance function.
     *
     * Parameters
     * ----------
     * phi_id: PatchData id for phi
     *
     * control_volume_id: PatchData id for control volume
     *
     * max_time_steps: maximum number of time steps to take
     *
     * steady_state_threshold: phi is considered to have reached steady-state
     *      when max | (phi(t+dt) - phi(t)) / dt | < steady_state_threshold
     *
     * stop_distance: approximate distance out from interface that signed
     *      distance is computed. When stop_distance <= 0, it is ignored.
     *
     * Notes
     * -----
     * - When stop_distance is positive, it takes precedence over
     *   max_time_steps.
     */
    void reinitializeLevelSetFunction(
            const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int phi_id,
            const int control_volume_id,
            const int time_integration_order,
            const int max_time_steps = 20,
            const double steady_state_condition = 1e-4,
            const double stop_distance = -1);

    /*!
     * Compute extension field, S, that extrapolates values of S off of the
     * interface defined by the zero level set of phi.
     *
     * Parameters
     * ----------
     * patch_hierarchy: PatchHierarchy to perform computation on
     *
     * S_id: PatchData id for S
     *
     * phi_id: PatchData id for phi
     *
     * control_volume_id: PatchData id for control volume
     *
     * max_time_steps: maximum number of time steps to take
     *
     * steady_state_threshold: S is considered to have reached steady-state
     *      when *      max | (S(t+dt) - S(t)) / dt | < steady_state_threshold
     *
     * stop_distance: approximate distance out from interface that S is
     *      extended off of the zero level set of phi. When stop_distance <= 0,
     *      it is ignored.
     *
     * Notes
     * -----
     * - When stop_distance is positive, it takes precedence over
     *   max_time_steps.
     */
    void computeExtensionField(
            const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int S_id,
            const int phi_id,
            const int control_volume_id,
            const int max_time_steps,
            const double steady_state_condition,
            const double stop_distance);

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Accessor methods for object parameters and state
     *
     ************************************************************************/

    /*!
     * Get pointer to PatchHierarchy.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * pointer to PatchHierarchy
     */
    virtual shared_ptr<hier::PatchHierarchy> getPatchHierarchy() const;

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

    // ------ Data management

    // maximum stencil width (over all simulation variables and computations)
    shared_ptr<hier::IntVector> d_max_stencil_width;

    // ------ PatchData IDs

    // PatchData component selectors to organize variables by
    // data management cycle requirements
    hier::ComponentSelector d_lsm_algs_variables;
    hier::ComponentSelector d_scratch_variables;

    // grid functions that used when solving level set method evolution
    // equations
    int d_lsm_algs_current_id;  // depth = 1
    int d_lsm_algs_next_id;  // depth = 1
    int d_lsm_algs_scratch_id;  // depth = 1
    int d_lsm_algs_rhs_id;  // depth = 1

    // --- Components

    // SAMRAI components
    shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;

    // Data transfer
    shared_ptr<xfer::RefineAlgorithm> d_xfer_fill_bdry_current;
    vector< shared_ptr<xfer::RefineSchedule> >
            d_xfer_fill_bdry_schedule_current;

    shared_ptr<xfer::RefineAlgorithm> d_xfer_fill_bdry_next;
    vector< shared_ptr<xfer::RefineSchedule> >
            d_xfer_fill_bdry_schedule_next;

private:

    /*
     * Set up simulation variables.
     */
    void setupSimulationVariables();

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
            const shared_ptr<hier::BaseGridGeometry>& grid_geometry);

    /*!
     ************************************************************************
     *
     * @name Data transfer methods
     *
     ************************************************************************/

    /*
     * Fill ghostcells for Patches in PatchLevel.
     *
     * Parameters
     * ----------
     * context: variable context that ghost cell data should be filled for
     *
     * Return value
     * ------------
     * None
     */
    void fillGhostCells(const int context) const;

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     *
     */
    Algorithms(const Algorithms& rhs) {};

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
     *
     */
    const Algorithms& operator=(const Algorithms& rhs) {
        return *this;
    }

};  // PQS::math::LSM::Algorithms class

}  // PQS::math::LSM namespace
}  // PQS::math namespace
}  // PQS namespace

#endif // INCLUDED_PQS_math_LSMAlgorithms_h
