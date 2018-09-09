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
     * @name Forward Euler method
     *
     ************************************************************************/

    /*!
     * Advance solution 'u' through a single first-order Runge-Kutta
     * (i.e., forward Euler) step.
     *
     * Parameters
     * ----------
     * patch_hierarchy: pointer to PatchHierarchy containing
     *
     * u_next_id: PatchData ID for u(t + dt)
     *
     * u_current_id: PatchData ID for u(t)
     *
     * rhs_id: PatchData ID for rhs(t)
     *
     * dt: time step to advance u(t) by
     *
     * Return value
     * ------------
     * None
     */
    static void RK1Step(
            shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_current_id,
            const int rhs_id,
            const double dt);
    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Second-order TVD Runge-Kutta method
     *
     ************************************************************************/

    /*!
     * Advance solution 'u' through a first stage of second-order
     * TVD Runge-Kutta step.
     *
     * Parameters
     * ----------
     * patch_hierarchy: pointer to PatchHierarchy containing
     *
     * u_stage1_id: PatchData ID for u_approx(t + dt)
     *
     * u_current_id: PatchData ID for u(t)
     *
     * rhs_id: PatchData ID for rhs(t)
     *
     * dt: time step to advance u(t) by
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - the first stage of TVD RK2 is identical to a single RK1 step
     */
    static void TVDRK2Stage1(
            shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_stage1_id,
            const int u_current_id,
            const int rhs_id,
            const double dt);

    /*!
     * Advance solution 'u' through a second stage of second-order
     * TVD Runge-Kutta step.
     *
     * Parameters
     * ----------
     * patch_hierarchy: pointer to PatchHierarchy containing
     *
     * u_next_id: PatchData ID for u(t + dt)
     *
     * u_stage1_id: PatchData ID for u_approx(t + dt)
     *
     * u_current_id: PatchData ID for u(t)
     *
     * rhs_id: PatchData ID for rhs(t)
     *
     * dt: time step to advance u(t) by
     *
     * Return value
     * ------------
     * None
     */
    static void TVDRK2Stage2(
            shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_stage1_id,
            const int u_current_id,
            const int rhs_id,
            const double dt);
    //! @}

    /*!
     ************************************************************************
     *
     * @name Third-order TVD Runge-Kutta method
     *
     ************************************************************************/

    /*!
     * Advance solution 'u' through a first stage of third-order
     * TVD Runge-Kutta step.
     *
     * Parameters
     * ----------
     * patch_hierarchy: pointer to PatchHierarchy containing
     *
     * u_stage1_id: PatchData ID for u_approx(t + dt)
     *
     * u_current_id: PatchData ID for u(t)
     *
     * rhs_id: PatchData ID for rhs(t)
     *
     * dt: time step to advance u(t) by
     *
     * Return value
     * ------------
     * None
     *
     * Notes
     * -----
     * - the first stage of TVD RK2 is identical to a single RK1 step
     */
    static void TVDRK3Stage1(
            shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_stage1_id,
            const int u_current_id,
            const int rhs_id,
            const double dt);

    /*!
     * Advance solution 'u' through a second stage of third-order
     * TVD Runge-Kutta step.
     *
     * Parameters
     * ----------
     * patch_hierarchy: pointer to PatchHierarchy containing
     *
     * u_stage2_id: PatchData ID for u(t + dt/2)
     *
     * u_stage1_id: PatchData ID for u_approx(t + dt)
     *
     * u_current_id: PatchData ID for u(t)
     *
     * rhs_id: PatchData ID for rhs(t)
     *
     * dt: time step to advance u(t) by
     *
     * Return value
     * ------------
     * None
     */
    static void TVDRK3Stage2(
            shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_stage2_id,
            const int u_stage1_id,
            const int u_current_id,
            const int rhs_id,
            const double dt);

    /*!
     * Advance solution 'u' through a third stage of third-order
     * TVD Runge-Kutta step.
     *
     * Parameters
     * ----------
     * patch_hierarchy: pointer to PatchHierarchy containing
     *
     * u_next_id: PatchData ID for u(t + dt)
     *
     * u_stage2_id: PatchData ID for u(t + dt/2)
     *
     * u_current_id: PatchData ID for u(t)
     *
     * rhs_id: PatchData ID for rhs(t)
     *
     * dt: time step to advance u(t) by
     *
     * Return value
     * ------------
     * None
     */
    static void TVDRK3Stage3(
            shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_stage2_id,
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

#endif // INCLUDED_PQS_math_TimeIntegration_h
