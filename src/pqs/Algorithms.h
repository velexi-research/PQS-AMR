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
 * "P_reference": reference pressure
 * "V_target": target volume of phase where phi < 0
 * "surface_tensions": surface tension
 * "bulk_modulus": dimensionless bulk modulus
 *
 */

// --- Headers, namespaces, and type declarations

// Standard
#include <ostream>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

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
     * psi_id: PatchData ID of the level set function for the solid-pore
     *      interface
     *
     * grad_psi_id: PatchData ID of the gradient of the level set function
     *      for the solid-pore interface
     */
    Algorithms(const boost::shared_ptr<tbox::Database>& config_db,
               const int psi_id,
               const int grad_psi_id);

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
     * phi_id: PatchData ID of level set function that defines fluid-fluid
     *      interface
     *
     * psi_id: PatchData ID of level set function that defines solid-pore
     *      interface
     *
     * grad_psi_id: PatchData ID of gradient of level set function that
     *      defines solid-pore interface
     *
     * Return value
     * ------------
     * maximum stable timestep on Patch
     */
    double computePrescribedCurvatureModelRHS(
        const boost::shared_ptr<hier::Patch> patch,
        const int phi_id) const;

    /*!
     * Compute RHS of level set equation for "Slightly Compressible Model".
     *
     * Parameters
     * ----------
     * patch: Patch to compute RHS of level set equation on
     *
     * phi_id: PatchData ID of level set function that defines fluid-fluid
     *      interface
     *
     * psi_id: PatchData ID of level set function that defines solid-pore
     *      interface
     *
     * grad_psi_id: PatchData ID of gradient of level set function that
     *      defines solid-pore interface
     *
     * Return value
     * ------------
     * maximum stable timestep on Patch
     */
    double computeSlightlyCompressibleModelRHS(
        const boost::shared_ptr<hier::Patch> patch,
        const int phi_id) const;

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

    // Slightly Compressible Model
    double d_P_reference;
    double d_V_target;
    double d_surface_tension;
    double d_bulk_modulus;

    // --- SAMRAI parameters

    // ------ PatchData IDs

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
            const boost::shared_ptr<tbox::Database>& config_db) const;

    /*
     * Load configuration parameters from specified database.
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     */
    void loadConfiguration(const boost::shared_ptr<tbox::Database>& config_db);

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

#endif
