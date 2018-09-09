/*! \file PoreInitStrategy.h
 *
 * \brief
 * Header file for PoreInitStrategy class.
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

#ifndef INCLUDED_PQS_pqs_PoreInitStrategy_h
#define INCLUDED_PQS_pqs_PoreInitStrategy_h

/*! \class PQS::pqs::PoreInitStrategy
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


// --- PQS::pqs::PoreInitStrategy Class

namespace PQS {
namespace pqs {

class PoreInitStrategy {

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
    PoreInitStrategy() {};

    /*!
     * Empty default destructor.
     */
    virtual ~PoreInitStrategy() {};

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
     * Set physical boundary conditions for phi, the level set function
     * that defines the solid-pore interace.
     *
     * Parameters
     * ----------
     * patch: Patch on which to set boundary conditions
     *
     * psi_id: PatchData ID for psi
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - No-op by default.
     */
    virtual void setBoundaryConditions(hier::Patch& patch, int psi_id) {};

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
    PoreInitStrategy(const PoreInitStrategy& rhs){}

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
    const PoreInitStrategy& operator=(const PoreInitStrategy& rhs) {
        return *this;
    }

};  // PQS::pqs::PoreInitStrategy class

}  // PQS::pqs namespace
}  // PQS namespace

#endif INCLUDED_PQS_pqs_PoreInitStrategy_h
