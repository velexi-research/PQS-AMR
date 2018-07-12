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
 * TODO
 *
 */

// --- Headers, namespaces, and type declarations

// Standard
#include <ostream>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"

// PQS
#include "PQS/PQS_config.h"
#include "PQS/pqs/InterfaceInitStrategy.h"
#include "PQS/pqs/PoreInitStrategy.h"

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
     * Empty default constructor.
     */
    Algorithms() {};

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
     * TODO
     *
     * Parameters
     * ----------
     */
    static void computePrescribedCurvatureModelRHS();

    /*!
     * TODO
     *
     * Parameters
     * ----------
     */
    static void computeSlightlyCompressibleModelRHS();

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

    // PQS
    // TODO

private:

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
