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
#include "PQS/math/LSMAlgorithms.h"
#include "PQS/pqs/Algorithms.h"
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace PQS { namespace pqs { class TagAndInitializeModule; } }
namespace PQS { namespace pqs { class DataTransferModule; } }


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
     * Constants
     *
     ************************************************************************/

    static constexpr double DEFAULT_LSM_STEADY_STATE_CONDITION = 0.01;
    static const int DEFAULT_LSM_MAX_ITERATIONS = 100;
    static constexpr double DEFAULT_LSM_STOP_TIME = 0.0;
    static constexpr double DEFAULT_LSM_SATURATION_STEADY_STATE_CONDITION = 0.0;

    static const int DEFAULT_REINITIALIZATION_INTERVAL = 5;
    static const int DEFAULT_TAG_BUFFER = 2;
    static const int DEFAULT_REGRID_INTERVAL = 5;

    //! @}

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
     *
     * enable_debug: flag indicating whether debug mode should be enabled
     */
    Solver(const shared_ptr<tbox::Database>& config_db,
           const shared_ptr<pqs::PoreInitStrategy>& pore_init_strategy,
           const shared_ptr<pqs::InterfaceInitStrategy>&
                   interface_init_strategy,
           const shared_ptr<hier::PatchHierarchy>& patch_hierarchy=NULL,
           const bool enable_debug=false);

    /*!
     * Default destructor frees memory allocated for simulation.
     */
    virtual ~Solver();

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Solver management methods
     *
     ************************************************************************/


    /*!
     * Reset any internal information that depends on the PatchHierarchy
     * or the data that resides on it.
     *
     * For the Solver class, 'resetHierarchyConfiguration()' updates the
     * following:
     *
     * - call resetHierarchyConfiguration() for LSM::Algorithms object.
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
     *
     */
    virtual void resetHierarchyConfiguration(
            const int coarsest_level_num,
            const int finest_level_num);

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
     * steady_state_condition: phi is considered to have reached steady-state
     *      when
     *
     *          max | (phi(t+dt) - phi(t)) / dt | < steady_state_condition.
     *
     *      When 'steady_state_condition' is set to a non-positive number,
     *      the 'lsm_steady_state_condition' value from the configuration
     *      database is used.
     *
     * max_num_iterations: maximum number of time steps to evolve level set
     *      function. When 'max_num_iterations' is set to zero, the number of
     *      time steps taken is not used as a stopping criteria.  When
     *      'max_num_iterations' is set to a negative number, the
     *      'lsm_max_num_iterations' value from the configuration database is
     *      used.
     *
     * stop_time: time at which evolution of the level set function is
     *      stopped (even if the fluid-fluid interface has not yet reached
     *      steady-state). When 'stop_time' is set to zero, the time is not
     *      used as a stopping criteria.  When 'stop_time' is set to a
     *      negative number, the 'lsm_stop_time' value from the configuration
     *      database is used.
     *
     * saturation_steady_state_condition: phi is considered to have reached
     *      steady-state when the saturation
     *
     *          max | (saturation(t+dt) - saturationphi(t)) / dt | <
     *              saturation_steady_state_condition.
     *
     *      When 'saturation_steady_state_condition' is set zero, the
     *      saturation is not used as a stopping criteria.  When
     *      'saturation_steady_state_condition' is set to a negative number,
     *      the 'lsm_saturation_steady_state_condition' value from the
     *      configuration database is used.
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
     *
     * - When 'algorithm_type' is set to SLIGHTLY_COMPRESSIBLE_MODEL, the
     *   steady-state fluid-fluid interface is computed using the Slightly
     *   Compressible Model (Prodanovic and Bryant, 2006).
     */
    virtual void equilibrateInterface(
            const double curvature,
            const PQS_ALGORITHM_TYPE algorithm_type,
            double steady_state_condition = -1,
            int max_num_iterations = -1,
            double stop_time = -1,
            double saturation_steady_state_condition = -1);

    /*!
     * Reinitialize fluid-fluid interface to be a signed distance function.
     *
     * Parameters
     * ----------
     * algorithm_type: algorithm to use to reinitialize interface.
     *      Valid values: math::LSM::REINIT_EQN_SGN_PHI0,
     *      math::LSM::REINIT_EQN_SGN_PHI.
     *
     * steady_state_condition: phi is considered to have reached steady-state
     *      when
     *
     *          max | (phi(t+dt) - phi(t)) / dt | < steady_state_condition
     *
     *      When 'steady_state_condition' is set to a non-positive number,
     *      the 'reinitialization_steady_state_condition' value from the
     *      configuration database is used.
     *
     * stop_distance: approximate distance out from interface that signed
     *      distance is computed. When 'stop_distance' is set to zero, the
     *      time is not used as a stopping criteria.  When 'stop_distance'
     *      is set to a negative number, the 'reinitialization_stop_distance'
     *      value from the configuration database is used.
     *
     * max_num_iterations: maximum number of reinitialization steps to take.
     *      When 'max_num_iterations' is set zero, the number of steps is not
     *      used as a stopping criteria.  When 'max_num_iterations' is set to
     *      a negative number, the 'reinitialization_max_num_iterations' value
     *      from the configuration database is used.
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - TODO: phi0 vs phi
     */
    virtual void reinitializeInterface(
            const math::LSM::REINIT_ALGORITHM_TYPE algorithm_type =
                    math::LSM::REINIT_EQN_SGN_PHI0,
            double steady_state_condition = -1,
            double stop_distance = -1,
            int max_num_iterations = -1);

    //! @}

    //! @{
    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Accessor methods for object parameters and state
     *
     ************************************************************************/

    // --- Solver Parameters

    /*!
     * Return whether or not to use Slightly Compressible Model to
     * equilibrate initial fluid-fluid interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * true if initial fluid-fluid interface should be equilibrated using
     * Slightly Compressible Model; false otherwise.
     */
    virtual bool initializeWithSlightlyCompressibleModel() const;

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
    virtual double getCurvatureIncrement() const;

    /*!
     * Get contact angle of wetting phase with solid phase.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * contact angle (in degrees)
     */
    virtual double getContactAngle() const;

    /*!
     * Get surface tension of fluid-fluid interface.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * surface tension
     */
    virtual double getSurfaceTension() const;

    /*!
     * Get bulk modulus of the non-wetting phase.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * bulk modulus of non-wetting phase in Slightly Compressible Model
     */
    virtual double getBulkModulus() const;

    /*!
     * Get target volume for the non-wetting phase.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * target volue of non-wetting phase in Slightly Compressible Model
     */
    virtual double getTargetVolume() const;

    /*!
     * Get maximum stencil width.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * maximum stencil width
     */
    virtual int getMaxStencilWidth() const;

    // --- Solver State

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

    // Physical parameters
    bool d_initialize_with_slightly_compressible_model;

    double d_contact_angle;  // units: degrees. default: 0
    double d_surface_tension;

    double d_bulk_modulus;  // only used for Slightly Compressible Model
    double d_target_volume;  // only used for Slightly Compressible Model

    // Level set method parameters
    double d_lsm_steady_state_condition;
    double d_lsm_stop_time;
    int d_lsm_max_num_iterations;
    double d_lsm_saturation_steady_state_condition;

    // Reinitialization parameters
    int d_reinitialization_interval;
    math::LSM::REINIT_ALGORITHM_TYPE d_reinitialization_algorithm_type;
    double d_reinitialization_steady_state_condition;
    double d_reinitialization_stop_distance;
    int d_reinitialization_max_num_iterations;

    // Numerical method parameters
    int d_lsm_spatial_derivative_type;  // ENO1, ENO2, ENO3, or WENO5
    int d_time_integration_order;  // for TVD Runge-Kutta time integration

    // AMR parameters
    int d_tag_buffer;  // number of cells to buffer tagged cells with
                       // during regridding operations
    int d_regrid_interval;  // number of equilibration steps to take
                            // between regridding operations

    // Debugging parameters
    bool d_enable_debug;

    // --- State

    // PQS
    double d_curvature;  // current value of prescribed mean curvature
    int d_step_count;  // number of prescribed curvature steps taken

    // --- SAMRAI parameters

    // ------ Data management

    // maximum stencil width (over all simulation variables and computations)
    int d_max_stencil_width;

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

    // PQS data transfer module
    shared_ptr<pqs::DataTransferModule> d_data_xfer_module;

    // LSM algorithms
    shared_ptr<PQS::math::LSM::Algorithms> d_lsm_algorithms;

    // SAMR grid
    shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
    shared_ptr<mesh::GriddingAlgorithm> d_gridding_algorithm;
    shared_ptr<pqs::TagAndInitializeModule> d_tag_and_init_module;

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
