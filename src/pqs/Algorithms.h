/*! \file Algorithms.h
 *
 * \brief
 * Header file for Algorithms class.
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

#ifndef INCLUDED_PQS_pqs_Algorithms_h
#define INCLUDED_PQS_pqs_Algorithms_h

/*! \class PQS::pqs::Algorithms
 *
 * \brief
 * TODO:
 * - provide interface velocity models as static methods for computing RHS of
 *   level set evolution equations
 *
 * <h3> NOTES </h3>
 *
 *  - TODO
 *
 * <h3> USAGE </h3>
 *
 * Configuration Parameters
 * ------------------------
 *
 * ### Prescribed Curvature Model
 * ------------------------------
 * surface_tension: surface tension of fluid-fluid interface
 *
 * ### Slightly Compressible Model
 * -------------------------------
 * bulk_modulus: dimensionless bulk modulus
 * target_volume: target volume of non-wetting phase (phi < 0)
 * surface_tension: surface tension of fluid-fluid interface
 *
 */

// --- Headers, namespaces, and type declarations

// Standard
#include <ostream>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

// PQS
#include "PQS/PQS_config.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace SAMRAI { namespace hier { class BaseGridGeometry; } }


// --- PQS::pqs::Algorithms Class

namespace PQS {
namespace pqs {

class Algorithms
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
     * This constructor for pqs::Algorithms creates a TODO
     *
     * - Create simulation variables.
     * - Register PatchData
     * - Create SAMRAI::mesh::GriddingAlgorithm
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     *
     * lse_rhs_id: PatchData ID for the right-hand side of level set evolution
     *     equation
     *
     * psi_id: PatchData ID for the level set function for the solid-pore
     *     interface
     *
     * grad_psi_id: PatchData ID for the gradient of the level set
     *     function for the solid-pore interface
     */
    Algorithms(const shared_ptr<tbox::Database>& config_db,
               const int lse_rhs_id,
               const int psi_id,
               const int grad_psi_id = -1);

    /*!
     * Empty default destructor.
     */
    virtual ~Algorithms() {};

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Interface velocity models
     *
     ************************************************************************/

    /*!
     * Compute RHS of level set equation for "Prescribed Curvature Model".
     *
     * Parameters
     * ----------
     * patch: Patch to compute RHS of level set equation on
     *
     * phi_id: PatchData ID for level set function that defines fluid-fluid
     *     interface
     *
     * psi_id: PatchData ID for level set function that defines solid-pore
     *     interface
     *
     * grad_psi_id: PatchData ID for gradient of level set function that
     *     defines solid-pore interface
     *
     * target_curvature: target curvature for fluid-fluid interface
     *
     * Return value
     * ------------
     * maximum stable timestep on Patch
     */
    double computePrescribedCurvatureModelRHS(
        const shared_ptr<hier::Patch>& patch,
        const int phi_id,
        const int psi_id,
        const int grad_psi_id,
        const double target_curvature) const;

    /*!
     * Compute RHS of level set equation for "Slightly Compressible Model".
     *
     * Parameters
     * ----------
     * patch: Patch to compute RHS of level set equation on
     *
     * phi_id: PatchData ID for level set function that defines fluid-fluid
     *     interface
     *
     * psi_id: PatchData ID for level set function that defines solid-pore
     *     interface
     *
     * grad_psi_id: PatchData ID for gradient of level set function that
     *     defines solid-pore interface
     *
     * target_curvature: target curvature for fluid-fluid interface
     *
     * volume: current volume of non-wetting phase
     *
     * Return value
     * ------------
     * maximum stable timestep on Patch
     */
    double computeSlightlyCompressibleModelRHS(
        const shared_ptr<hier::Patch>& patch,
        const int phi_id,
        const int psi_id,
        const int grad_psi_id,
        const double target_curvature,
        const double volume) const;

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

    // ------ PQS

    // Physical parameters
    double d_contact_angle;  // default: 0
    double d_surface_tension;

    // Slightly Compressible Model
    double d_scm_bulk_modulus;
    double d_scm_target_volume;

    // --- SAMRAI parameters

    // ------ PatchData IDs

    // right-hand side of level set evolution equation
    int d_lse_rhs_id;   // depth = 1

    // level set function that defines solid-pore interface
    //
    // Note: the solid phase is defined by the region where psi < 0
    int d_psi_id;   // depth = 1
    int d_grad_psi_id;  // depth = number of spatial dimensions

private:

    /*
     * Verify that configuration database from specified database.
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
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
            const shared_ptr<tbox::Database>& config_db) const;

    /*
     * Load configuration parameters from specified database.
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     */
    void loadConfiguration(const shared_ptr<tbox::Database>& config_db);

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     */
    Algorithms(const Algorithms& rhs);

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
    const Algorithms& operator=(const Algorithms& rhs) {
        return *this;
    }

};  // PQS::pqs::Algorithms class

}  // PQS::pqs namespace
}  // PQS namespace

#endif // INCLUDED_PQS_pqs_Algorithms_h
