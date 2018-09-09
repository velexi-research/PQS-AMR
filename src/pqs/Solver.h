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
 *      time_integration_order    = 3
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


// --- PQS::pqs::Solver Class

namespace PQS {
namespace pqs {

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
    PQS = 0,
    LSM_CURRENT = 1,
    LSM_NEXT = 2,
    COMPUTED = 3,
    SAMR = 4
} VARIABLE_CONTEXT;

/*! \enum PQS_ALGORITHM_TYPE
 *
 * Enumerated type for PQS algorithms.
 *
 */
typedef enum {
    PRESCRIBED_CURVATURE_MODEL = 0,
    SLIGHTLY_COMPRESSIBLE_MODEL = 1
} PQS_ALGORITHM_TYPE;

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
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     *
     * pore_init_strategy: user implementation of PoreInitStrategy strategy
     *      interface
     *
     * interface_init_strategy: user implementation of
     *      InterfaceInitStrategy strategy interface
     *
     * patch_hierarchy: PatchHierarchy to use for simulation
     */
    Solver(const shared_ptr<tbox::Database>& config_db,
           const shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
           const shared_ptr<pqs::InterfaceInitStrategy>&
                   interface_init_strategy,
           const shared_ptr<hier::PatchHierarchy>& patch_hierarchy=NULL);

    /*!
     * Default destructor frees memory allocated for simulation.
     */
    virtual ~Solver();

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
     * curvature: target mean curvature for steady-state interface
     *
     * algorithm_type: algorithm to use to equilibrate fluid-fluid
     *      interface. Valid values: PRESCRIBED_CURVATURE_MODEL,
     *      SLIGHTLY_COMPRESSIBLE_MODEL.
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - When 'algorithm_type' is set to PRESCRIBED_CURVATURE_MODEL, the
     *   steady-state fluid-fluid interface is computed using the Prescribed
     *   Curvature Model (Prodanovic and Bryant, 2006).
     * - When 'algorithm_type' is set to SLIGHTLY_COMPRESSIBLE_MODEL, the
     *   steady-state fluid-fluid interface is computed using the Slightly
     *   Compressible Model (Prodanovic and Bryant, 2006).
     */
    virtual void equilibrateInterface(const double curvature,
                                      const PQS_ALGORITHM_TYPE algorithm_type);

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Accessor methods for object parameters and state
     *
     ************************************************************************/

    /*!
     * Get current mean curvature of fluid-fluid interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * current value of mean curvature
     *
     * Notes
     * -----
     * - After either equilibrateInterface() or advanceInterface() has been
     *   called, the returned mean curvature value is equal to the value of
     *   the fluid-fluid interface.
     *
     * - At the beginning of the simulation (i.e., before either
     *   equilibrateInterface() or advanceInterface() has been called), the
     *   returned mean curvature value is equal to the desired value for the
     *   initial curvature of the fluid-fluid interface.
     */
    virtual double getCurvature() const;

    /*!
     * Get initial mean curvature of fluid-fluid interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * initial value of mean curvature
     */
    virtual double getInitialCurvature() const;

    /*!
     * Get final mean curvature of fluid-fluid interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * final value of mean curvature
     */
    virtual double getFinalCurvature() const;

    /*!
     * Get mean curvature increment to use to advance fluid-fluid interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * mean curvature increment
     */
    virtual double getCurvatureStep() const;

    /*!
     * Get number of curvature steps taken since the beginning of the
     * simulation.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * number of steps taken since the beginning of the simulation
     */
    virtual int getStepCount() const;

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
     * Get PatchData ID for level set function that defines solid-pore
     * interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * PatchData ID for level set function that defines solid-pore interface
     */
    virtual int getPoreSpacePatchDataId() const;

    /*!
     * Get PatchData ID for level set function that defines fluid-fluid
     * interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * PatchData ID for level set function that defines fluid-fluid interface
     */
    virtual int getInterfacePatchDataId() const;

    /*!
     * Get PatchData ID for control volume.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * PatchData ID for control volume
     */
    virtual int getControlVolumePatchDataId() const;

    /*!
     * Get PatchData ID for RHS of level set evolution equation.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * PatchData ID for RHS of level set evolution equation
     *
     * Notes
     * -----
     * - Intended to be used primarily for debugging purposes.
     */
    virtual int getLSERHSPatchDataId() const;

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
    double d_initial_curvature;
    double d_final_curvature;
    double d_curvature_step;

    // Level set method parameters
    double d_lsm_t_max;
    int d_lsm_max_iterations;
    double d_lsm_min_delta_phi;
    double d_lsm_min_delta_saturation;

    int d_lsm_spatial_derivative_type;  // ENO1, ENO2, ENO3, or WENO5

    // Numerical method parameters
    int d_time_integration_order;  // for TVD Runge-Kutta time integration

    // --- State

    // PQS
    double d_curvature;  // current value of prescribed mean curvature
    int d_step_count;  // number of prescribed curvature steps taken

    // AMR
    int d_regrid_count;

    // --- SAMRAI parameters

    // ------ Data management

    // maximum stencil width (over all simulation variables and computations)
    shared_ptr<hier::IntVector> d_max_stencil_width;

    // ------ PatchData IDs

    // PatchData component selectors to organize variables by
    // data management cycle requirements
    hier::ComponentSelector d_permanent_variables;
    hier::ComponentSelector d_intermediate_variables;

    // level set function that defines steady-state fluid-fluid interface
    // before increment of the interface curvature (which is related to
    // changes in the pressure difference across the fluid-fluid interface)
    //
    // Note: the non-wetting phase is defined by the region where phi < 0
    int d_phi_pqs_id;  // depth = 1

    // level set function that defines fluid-fluid interface during evolution
    // of the interface towards steady-state
    //
    // Note: d_phi_lsm_next_id is also used to hold intermediate
    //       Runge-Kutta steps
    int d_phi_lsm_current_id;  // depth = 1
    int d_phi_lsm_next_id;  // depth = 1

    // level set function that defines solid-pore interface
    //
    // Note: the solid phase is defined by the region where psi < 0
    int d_psi_id;   // depth = 1
    int d_grad_psi_id;  // depth = number of spatial dimensions

    // level set method
    int d_lse_rhs_id;  // depth = 1

    // AMR data
    int d_control_volume_id;  // depth = 1

    // TODO
    // int d_curvature_id;  // depth = 1
    // int d_connectivity_label_id;  // depth = 1
    // int d_phi_binary;  // depth = 1

    // velocities
    //int d_normal_velocity_id;  // depth = 1
    //int d_vector_velocity_id;  // depth = number of spatial dimensions

    // --- Components

    // PQS algorithms
    shared_ptr<pqs::Algorithms> d_pqs_algorithms;

    // SAMR grid
    shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
    shared_ptr<mesh::GriddingAlgorithm> d_gridding_algorithm;
    shared_ptr<pqs::TagInitAndDataTransferModule>
            d_tag_init_and_data_xfer_module;

    // Boundary conditions
    //shared_ptr<BoundaryConditionModule> d_bc_module;
    //tbox::Array<hier::IntVector> d_lower_bc_phi;
    //tbox::Array<hier::IntVector> d_upper_bc_phi;
    //tbox::Array<hier::IntVector> d_lower_bc_psi;
    //tbox::Array<hier::IntVector> d_upper_bc_psi;

private:

    /*
     * Verify that configuration database from specified database.
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     *
     * verify_patch_hierarchy: flag indicating whether or not to check
     *      PatchHierarchy and Geometry parameters
     *
     * Return value
     * ------------
     * true
     *
     * Notes
     * -----
     * - If 'config_db' is not valid, this method throws an exception
     *   containing error information.
     */
    void verifyConfigurationDatabase(
        const shared_ptr<tbox::Database>& config_db,
        const bool verify_patch_hierarchy=false) const;

    /*
     * Load configuration parameters from specified database.
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     */
    void loadConfiguration(const shared_ptr<tbox::Database>& config_db);

    /*
     * Create PatchHierarchy.
     */
    void createPatchHierarchy(
        const shared_ptr<tbox::Database>& config_db);

    /*
     * Set up simulation variables.
     */
    void setupSimulationVariables();

    /*
     * Set up grid management.
     *
     * Parameters
     * ----------
     * config_db: Database containing configuration parameters
     *
     * pore_init_strategy: PoreInitStrategy object to use for initialization
     *     of level set function defining solid-pore interface
     *
     * interface_init_strategy: InterfaceInitStrategy object to use for
     *     initialization of level set function defining fluid-fluid interface
     *
     * Return value
     * ------------
     * None
     */
    void setupGridManagement(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
        const shared_ptr<pqs::InterfaceInitStrategy>&
            interface_init_strategy);

    /*
     * Initialize simulation.
     */
    void initializeSimulation();

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

#endif // INCLUDED_PQS_pqs_Solver_h
