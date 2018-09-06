/*! \file Toolbox.h
 *
 * \brief
 * Header file for Toolbox class.
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

#ifndef INCLUDED_PQS_lsm_Toolbox_h
#define INCLUDED_PQS_lsm_Toolbox_h

/*! \class PQS::pqs::Toolbox
 *
 * \brief
 * TODO: add description
 */

// --- Headers, namespaces, and type declarations

// Standard headers
#include <ostream>

// SAMRAI headers
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"

// PQS headers
#include "PQS/PQS_config.h"
#include "PQS/lsm/ReinitializationAlgorithm.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations


// --- PQS::lsm::Toolbox Class

namespace PQS {
namespace lsm {

class Toolbox {

public:

    //! @{

    /*!
     * This constructor for Toolbox creates a TODO
     */
    Toolbox();

    /*!
     * Empty default destructor.
     */
    virtual ~Toolbox() {};

    //! @}

    //! @{
    /*!
     ************************************************************************
     *
     * @name Utility functions
     *
     ************************************************************************/

    /*!
     * Compute the volume of one of two following regions:
     * region with phi < 0 or region with phi > 0.
     *
     * Parameters
     * ----------
     * patch_hierarchy: PatchHierarchy on which to compute the volume
     * 
     * phi_id: PatchData id for phi
     *
     * control_volume_id: PatchData id for control volume
     *
     * region_indicator: integer indicating which region to integrate over
     *      - region_indicator >0:  integration region = {x | phi(x) > 0}
     *      - region_indicator <=0: integration region = {x | phi(x) <= 0}
     *
     * Return value
     * ------------
     * volume of the specified domain that is enclosed by the zero level set
     *
     * Notes
     * -----
     * - When phi is a signed distance function, the smoothed Heaviside
     *   function used compute the integral has a width approximately equal
     *   to the three grid cells widths.
     */
    static PQS_REAL computeVolume(
            const shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int phi_id,
            const int control_volume_id,
            const int region_indicator);

    //! @}

    //! @{
    /*!
     ************************************************************************
     *
     * @name Accessor methods for object parameters and state
     *
     ************************************************************************/

    /*!
     * printClassData() prints the values of the data members for
     * an instance of the LevelSetMethodAlgorithm class.
     *
     * Parameters
     * ----------
     * os: output stream to write object information
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

    bool d_use_reinitialization;
    bool d_use_reinitialization_stop_tol;
    bool d_use_reinitialization_stop_dist;
    bool d_use_reinitialization_max_iters;

    // --- PQS PatchData IDs

    // gradient of level set function
    // TODO: is this needed?
    int d_grad_phi_id;  // depth = number of spatial dimensions

    // right-hand side of level set evolution equation
    int d_lse_rhs_id;  // depth = 1

    // --- State

    // Reinitialization
    int d_reinitialization_count;

    // --- Components

    shared_ptr<lsm::ReinitializationAlgorithm> d_reinitialization_alg;

private:

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     *
     */
    Toolbox(const Toolbox& rhs){}

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
    const Toolbox& operator=(const Toolbox& rhs) {
        return *this;
    }

};  // PQS::lsm::Toolbox class

}  // PQS::lsm namespace
}  // PQS namespace

#endif
