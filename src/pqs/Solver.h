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

#ifndef INCLUDED_PQS_pqs_Solver_h
#define INCLUDED_PQS_pqs_Solver_h

/*! \class PQS::pqs::Solver
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

// Standard library
#include <ostream>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/Database.h"

// PQS
#include "PQS/PQS_config.h"
#include "PQS/pqs/GridManager.h"
#include "PQS/pqs/DataInitStrategy.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations


// --- PQS::pqs::Solver Class

namespace PQS {
namespace pqs {

class Solver {

public:

    //! @{

    /*!
     * This constructor for Solver creates a TODO
     */
    Solver(boost::shared_ptr<tbox::Database> config_db,
           boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
           boost::shared_ptr<pqs::DataInitStrategy> data_init_strategy);

    /*!
     * Empty default destructor.
     */
    virtual ~Solver() {};

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

    // TODO
    //! @}

protected:

    /****************************************************************
     *
     * Data Members
     *
     ****************************************************************/

    // Grid management
    boost::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
    boost::shared_ptr<mesh::GriddingAlgorithm> d_gridding_alg;
    boost::shared_ptr<pqs::GridManager> d_grid_manager;

    // TODO
    // Data management

private:

    /*
     * Load configuration parameters from specified database.
     */
    void loadConfiguration(const boost::shared_ptr<tbox::Database>& config_db);

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

};  // PQS::pqs::Solver class

}  // PQS::pqs namespace
}  // PQS namespace

#endif
