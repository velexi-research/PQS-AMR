/*! \file DataInitStrategy.h
 *
 * \brief
 * Header file for DataInitStrategy class.
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

#ifndef INCLUDED_PQS_pqs_DataInitStrategy_h
#define INCLUDED_PQS_pqs_DataInitStrategy_h

/*! \class PQS::pqs::DataInitStrategy
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

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/Patch.h"

// PQS
#include "PQS/PQS_config.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations


// --- PQS::pqs::DataInitStrategy Class

namespace PQS {
namespace pqs {

class DataInitStrategy {

public:

    //! @{

    /*!
     ************************************************************************
     *
     * @name Constructor and destructor
     *
     ************************************************************************/

    /*!
     * Empty default constructor.
     */
    DataInitStrategy() {};

    /*!
     * Empty default destructor.
     */
    virtual ~DataInitStrategy() {};

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Methods for setting initial and boundary conditions
     *
     ************************************************************************/

    /*!
     * Initialize psi, the level set function that defines the solid-pore
     * interface.
     *
     * Parameters
     * ----------
     * patch: Patch on which to initialize psi
     *
     * psi_id: PatchData ID for psi
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - This is a pure abstract method that the user MUST override in
     *   order to build a PQS simulation.
     */
    virtual void initializePoreSpace(hier::Patch& patch, int psi_id) = 0;

    /*!
     * Initialize phi, the level set function that defines the fluid-fluid
     * interface.
     *
     * Parameters
     * ----------
     * patch: Patch on which to initialize phi
     *
     * phi_id: PatchData ID for phi
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - This is a pure abstract method that the user MUST override in
     *   order to build a PQS simulation.
     */
    virtual void initializeInterface(hier::Patch& patch, int phi_id) = 0;

    /*!
     * Set physical boundary conditions for phi, the level set function
     * that defines the fluid-fluid interace.
     *
     * Parameters
     * ----------
     * patch: Patch on which to set boundary conditions
     *
     * phi_id: PatchData ID for phi
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - No-op by default.
     */
    virtual void setBoundaryConditions(hier::Patch& patch, int phi_id) {};

    //! @}

    //! @{
    /*!
     ************************************************************************
     *
     * @name Accessor methods for object parameters and state
     *
     ************************************************************************/

    /*!
     * Print the values of the data members for object.
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

private:

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     *
     */
    DataInitStrategy(const DataInitStrategy& rhs){}

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
    const DataInitStrategy& operator=(const DataInitStrategy& rhs) {
        return *this;
    }

};  // PQS::pqs::DataInitStrategy class

}  // PQS::pqs namespace
}  // PQS namespace

#endif
