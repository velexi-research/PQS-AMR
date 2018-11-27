/*! \file DataTransferModule.h
 *
 * \brief
 * Header file for DataTransferModule class.
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

#ifndef INCLUDED_PQS_pqs_DataTransferModule_h
#define INCLUDED_PQS_pqs_DataTransferModule_h

/*! \class PQS::pqs::DataTransferModule
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
#include <vector>

// SAMRAI
#include "SAMRAI/SAMRAI_config.h"  // IWYU pragma: keep
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"

// PQS
#include "PQS/PQS_config.h"
#include "PQS/pqs/Solver.h"

// Namespaces
using namespace std;
using namespace SAMRAI;

// Class/type declarations
namespace SAMRAI { namespace hier { class BaseGridGeometry; } }


// --- PQS::pqs::DataTransferModule Class

namespace PQS {
namespace pqs {

class DataTransferModule
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
     * This constructor for DataTransferModule TODO
     *
     * Parameters
     * ----------
     * config_db: database containing configuration parameters
     *
     * patch_hierarchy: pointer to PatchHierarchy that object manages
     *
     * phi_pqs_id: PatchData ID for the level set function for the
     *      fluid-fluid interface (at each PQS step)
     *
     * phi_lsm_current_id: PatchData ID for the level set function for the
     *      fluid-fluid interface (before each time step during evolution
     *      of the level set function)
     *
     * phi_lsm_next_id: PatchData ID for the level set function for the
     *      fluid-fluid interface (during and after each time step during
     *      evolution of the level set function)
     *
     * psi_id: PatchData ID for the level set function for the solid-pore
     *      interface
     *
     * max_stencil_width: maximum stencil width required for computations
     *
     */
    DataTransferModule(
        const shared_ptr<tbox::Database>& config_db,
        const shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
        const int phi_pqs_id,
        const int phi_lsm_current_id,
        const int phi_lsm_next_id,
        const int psi_id,
        const int max_stencil_width);

    /*!
     * Empty default destructor.
     */
    virtual ~DataTransferModule() {};

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Data transfer methods
     *
     ************************************************************************/

    /*!
     * Fill ghostcells for Patches in PatchLevel.
     *
     * Parameters
     * ----------
     * context: variable context that ghost cell data should be filled for
     *
     * Return value
     * ------------
     * None
     */
    virtual void fillGhostCells(const int context) const;

    /*!
     * Enforce consistency of phi across PatchLevels.
     *
     * Parameters
     * ----------
     * context: variable context that data should be made consistent for
     *
     * Return value
     * ------------
     * None
     */
    virtual void enforcePhiConsistency(const int context) const;

    /*!
     * Enforce consistency of psi across PatchLevels.
     *
     * Parameters
     * ----------
     * None
     *
     * Return value
     * ------------
     * None
     */
    virtual void enforcePsiConsistency() const;

    //! @}

    //! @{

    /*!
     ************************************************************************
     *
     * @name Module state and parameter methods
     *
     ************************************************************************/

    /*!
     * Reset any internal information that depends on the PatchHierarchy
     * or the data that resides on it:
     *
     * - communication schedules for PatchLevels in the PatchHierarchy.
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
     */
    virtual void resetHierarchyConfiguration(
            const int coarsest_level_num,
            const int finest_level_num);

    //! @}

protected:

    /****************************************************************
     *
     * Data Members
     *
     ****************************************************************/

    // --- PatchData IDs

    // PatchData component selectors to organize variables by
    // data management cycle requirements

    // fluid-fluid interface level set function
    int d_phi_pqs_id;
    int d_phi_lsm_current_id;
    int d_phi_lsm_next_id;
    int d_phi_scratch_id;

    // solid-pore interface level set function
    //
    // Note: the solid phase is defined by the region where psi < 0
    int d_psi_id;
    int d_psi_scratch_id;

    // --- Components

    // SAMRAI components
    shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;

    // Data transfer: fill boundary data for phi when solving level set
    // evolution equation
    shared_ptr<xfer::RefineAlgorithm> d_xfer_fill_bdry_lsm_current;
    vector< shared_ptr<xfer::RefineSchedule> >
            d_xfer_fill_bdry_schedules_lsm_current;

    shared_ptr<xfer::RefineAlgorithm> d_xfer_fill_bdry_lsm_next;
    vector< shared_ptr<xfer::RefineSchedule> >
            d_xfer_fill_bdry_schedules_lsm_next;

    // Data transfer: enforce consistency of phi when solving level set
    // evolution equation
    shared_ptr<xfer::CoarsenAlgorithm>
            d_xfer_enforce_phi_consistency_lsm_next;
    vector< shared_ptr<xfer::CoarsenSchedule> >
            d_xfer_enforce_phi_consistency_schedules_lsm_next;

    // Data transfer: enforce consistency of phi when initializing simulation
    shared_ptr<xfer::CoarsenAlgorithm>
            d_xfer_enforce_phi_consistency_pqs;
    vector< shared_ptr<xfer::CoarsenSchedule> >
            d_xfer_enforce_phi_consistency_schedules_pqs;

    // Data transfer: enforce consistency of psi
    shared_ptr<xfer::CoarsenAlgorithm> d_xfer_enforce_psi_consistency;
    vector< shared_ptr<xfer::CoarsenSchedule> >
            d_xfer_enforce_psi_consistency_schedules;

private:

    /*
     * Set up data transfer objects and associated scratch space variables.
     *
     * Parameters
     * ----------
     * grid_geometry: grid geometry that defines refinement operators
     *
     * max_stencil_width: maximum stencil width required for computations
     */
    void setupDataTransferObjects(
        const shared_ptr<hier::BaseGridGeometry>& grid_geometry,
        const int max_stencil_width);

    /*
     * Private copy constructor to prevent use.
     *
     * Parameters
     * ----------
     * rhs: object to copy
     */
    DataTransferModule(const DataTransferModule& rhs){}

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
    const DataTransferModule& operator=(const DataTransferModule& rhs) {
        return *this;
    }

};  // PQS::pqs::DataTransferModule class

}  // PQS::pqs namespace
}  // PQS namespace

#endif // INCLUDED_PQS_pqs_DataTransferModule_h
