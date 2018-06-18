/*! \file Solver.h
 *
 * \brief
 * Header file for Solver class.
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

#ifndef INCLUDED_PQS_Solver_h
#define INCLUDED_PQS_Solver_h

/*! \class PQS::Solver
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

// Standard headers
#include <ostream>

// Boost headers
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"

// PQS headers
#include "PQS/PQS_config.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations


// --- PQS::Solver Class

namespace PQS {

class Solver {

public:

    //! @{

    //! @name Constructor and destructor

    /*!
     * This constructor for Solver creates a TODO
     */
    Solver(boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy);

    /*!
     * The destructor for Solver does nothing.
     */
    virtual ~Solver(){};

    //! @}

    //! @{
    /*!
     *******************************************************************
     *
     *  @name Accessor methods for solver state
     *
     *******************************************************************/

    /*!
     * printClassData() prints the values of the data members for
     * an instance of the LevelSetMethodAlgorithm class.
     *
     * Parameters
     * ----------
     * os: output stream to write object information
     *
     * Return value
     * ------------
     * None
     */
    virtual void printClassData(ostream& os) const;

    //! @}

    //! @{
    /*!
     ****************************************************************
     *
     * @ Accessor methods for solver parameters
     *
     ****************************************************************/

    // TODO
    //! @}

protected:

  /****************************************************************
   *
   * Data Members
   *
   ****************************************************************/

  // object name
  const string d_object_name = "PQS::Solver";

  // Grid management
  boost::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
  boost::shared_ptr<mesh::GriddingAlgorithm> d_gridding_alg;

  // TODO
  // Data management

private:

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: Solver object to copy
     *
     */
    Solver(const Solver& rhs){}

    /*
     * Private assignment operator to prevent use.
     *
     * Parameters
     * ----------
     * rhs: Solver object to copy
     *
     * Return value
     * ------------
     * return Solver object
     *
     */
    const Solver& operator=(const Solver& rhs) {
        return *this;
    }

};  // PQS::Solver

}  // PQS namespace

#endif
