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
#include "PQS/pqs/DataInitStrategy.h"
#include "PQS/pqs/TagAndInitModule.h"

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
     ************************************************************************
     *
     * @name Constructor and destructor
     *
     ************************************************************************/

    /*!
     * This constructor for Solver creates a TODO
     *
     * - Create simulation variables.
     * - Register PatchData
     * - Create SAMRAI::mesh::GriddingAlgorithm
     */
    Solver(const boost::shared_ptr<tbox::Database>& config_db,
           const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
           const boost::shared_ptr<pqs::DataInitStrategy>& data_init_strategy);

    /*!
     * Empty default destructor.
     */
    virtual ~Solver() {};

    //! @}

    //! @{
    /*!
     ************************************************************************
     *
     * @name Interface dyanmics methods
     *
     ************************************************************************/

    /*!
     * Equilibrate fluid-fluid interface.
     *
     * Parameters
     * ----------
     * curvature: target curvature for steady-state interface
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - The steady-state interface is computed using the Slightly Compressible
     *   Model (Prodanovic and Bryant, 2006).
     */
    virtual void equilibrateInterface(const double curvature);

    /*!
     * Advance fluid-fluid interface.
     *
     * Parameters
     * ----------
     * delta_curvature: curvature increment to apply to current fluid-fluid
     *     interface
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - The new interface is computed using the Prescribed Curvature Model
     *   (Prodanovic and Bryant, 2006).
     */
    virtual void advanceInterface(const double delta_curvature);

    //! @}

    //! @{
    /*!
     ************************************************************************
     *
     * @name Accessor methods for object parameters and state
     *
     ************************************************************************/

    /*!
     * Get current curvature of fluid-fluid interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * current value of curvature
     *
     * Notes
     * -----
     * - After either equilibrateInterface() or advanceInterface() has been
     *   called, the returned curvature value is equal to the value of the
     *   fluid-fluid interface.
     *
     * - At the beginning of the simulation (i.e., before either
     *   equilibrateInterface() or advanceInterface() has been called), the
     *   returned curvature value is equal to the initial target value for
     *   the fluid-fluid interface.
     */
    virtual double getCurvature() const;

    /*!
     * Get number of curvature increments taken since the beginning of the
     * simulation.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * curvature
     */
    virtual int getStep() const;

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

    // PQS
    // TODO

    // --- PatchData IDs

    // PatchData component selectors to organize variables by
    // data management cycle requirements
    hier::ComponentSelector d_permanent_variables;
    hier::ComponentSelector d_intermediate_variables;

    // steady state fluid-fluid interface level set function before and
    // after increments of the interface curvature (which is related to
    // changes in the pressure difference across the fluid-fluid interface)
    int d_phi_pqs_current_id;  // depth = 1
    int d_phi_pqs_next_id;  // depth = 1

    // fluid-fluid interface level set function during evolution of the
    // interface towards steady-state
    int d_phi_lsm_current_id;  // depth = 1
    int d_phi_lsm_next_id;  // depth = 1

    // solid-pore interface level set function
    //
    // Note: the solid phase is defined by the region where psi < 0
    int d_psi_id;   // depth = 1
    int d_grad_psi_id;  // depth = number of spatial dimensions

    // velocities
    int d_normal_velocity_id;  // depth = 1
    int d_external_velocity_id;  // depth = number of spatial dimensions

    // geometry
    int d_curvature_id;  // depth = 1

    // AMR data
    int d_control_volume_id;  // depth = 1

    // TODO
    // int d_connectivity_label_id;  // depth = 1
    // int d_phi_binary;  // depth = 1

    // --- State

    // PQS
    double d_curvature;  // current value of prescribed curvature
    int d_num_steps;  // number of prescribed curvature steps taken

    // AMR
    int d_regrid_count;

    // --- Components

    // SAMR grid
    boost::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
    boost::shared_ptr<mesh::GriddingAlgorithm> d_gridding_alg;
    boost::shared_ptr<pqs::TagAndInitModule> d_tag_and_init_module;

    // Boundary conditions
    //boost::shared_ptr<BoundaryConditionModule> d_bc_module;
    //tbox::Array<hier::IntVector> d_lower_bc_phi;
    //tbox::Array<hier::IntVector> d_upper_bc_phi;
    //tbox::Array<hier::IntVector> d_lower_bc_psi;
    //tbox::Array<hier::IntVector> d_upper_bc_psi;

    // SAMRAI communication schedules
    // TODO

private:

    /*
     * Load configuration parameters from specified database.
     */
    void loadConfiguration(const boost::shared_ptr<tbox::Database>& config_db);

    /*
     * Set up simulation variables.
     */
    void setupSimulationVariables();

    /*
     * Set up grid management.
     */
    void setupGridManagement(
        const boost::shared_ptr<tbox::Database>& config_db,
        const boost::shared_ptr<pqs::DataInitStrategy>& data_init_strategy);

    /*
     * Initialize simulation data.
     */
    void initializeData();

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     *
     */
    Solver(const Solver& rhs) {};

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
    const Solver& operator=(const Solver& rhs) {
        return *this;
    }

};  // PQS::pqs::Solver class

}  // PQS::pqs namespace
}  // PQS namespace

#endif
