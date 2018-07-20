/*! \file TimeIntegration.h
 *
 * \brief
 * Header file for TimeIntegration class.
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

#ifndef INCLUDED_PQS_math_TimeIntegration_h
#define INCLUDED_PQS_math_TimeIntegration_h

/*! \class PQS::math::TimeIntegration
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

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep

// PQS
#include "PQS/PQS_config.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchHierarchy; } }


// --- PQS::math::TimeIntegration Class

namespace PQS {
namespace math {

class TimeIntegration {

public:

    //! @{

    /*!
     ************************************************************************
     *
     * @name Forward Euler methods
     *
     ************************************************************************/

    /*!
     * Advance solution through a single forward Euler time step.
     *
     * Parameters
     * ----------
     * hierarchy: pointer to PatchHierarchy containing
     *
     * u_next_id: PatchData ID for u(t+dt)
     *
     * u_cur_id: PatchData ID for u(t)
     *
     * rhs_id: PatchData ID for rhs(t)
     *
     * dt: time step to advance u(t) by
     *
     * Return value
     * ------------
     * None
     *
     */
    static void forwardEulerStep(
            boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_current_id,
            const int rhs_id,
            const double dt);

    //! @}

private:

    /*
     * Private default constructor to prevent use.
     */
    TimeIntegration(){}

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     *
     */
    TimeIntegration(const TimeIntegration& rhs){}

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
    const TimeIntegration& operator=(const TimeIntegration& rhs) {
        return *this;
    }

};  // PQS::math::TimeIntegration class

}  // PQS::math namespace
}  // PQS namespace

#endif
