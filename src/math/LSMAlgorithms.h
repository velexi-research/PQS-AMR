/*! \file LSMAlgorithms.h
 *
 * \brief
 * Header file for LSM::Algorithms class.
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

// PQS
#include "PQS/PQS_config.h"
#include "PQS/pqs/Algorithms.h"
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"
#include "PQS/pqs/TagInitAndDataTransferModule.h"

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
     * @name Level Set Method Algorithms
     *
     ************************************************************************/

    /*!
     * Reinitialize level set function to a signed distance function.
     *
     * Parameters
     * ----------
     * phi_id: PatchData id for phi
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
     * phi_id: PatchData id for phi
     *
     * S_id: PatchData id for S
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
            const int phi_id,
            const int S_id,
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

    // --- State

    // --- SAMRAI parameters

    // ------ Data management

    // maximum stencil width (over all simulation variables and computations)
    shared_ptr<hier::IntVector> d_max_stencil_width;

    // ------ PatchData IDs

    // data management cycle requirements
    hier::ComponentSelector d_lsm_algorithm_variables;

    // grid functions that used when solving level set method evolution
    // equations
    int d_lsm_current_id;  // depth = 1
    int d_lsm_next_id;  // depth = 1
    int d_lsm_rhs_id;  // depth = 1

    // --- Components

    // SAMR grid
    shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;

private:

    /*
     * Set up simulation variables.
     */
    void setupSimulationVariables();

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
