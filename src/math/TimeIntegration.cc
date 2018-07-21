/*! \file TimeIntegration.cc
 *
 * \brief
 * Implementation file for TimeIntegration class.
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

// --- Headers, namespaces, and type declarations

// Standard library
#include <stddef.h>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/Utilities.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/utilities.h"
#include "PQS/math/TimeIntegration.h"
#include "PQS/math/time_integration_2d.h"
#include "PQS/math/time_integration_3d.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchData; } }


// --- Class implementation

namespace PQS {
namespace math {

// --- Implementation of public methods

void TimeIntegration::RK1Step(
            boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_current_id,
            const int rhs_id,
            const double dt)
{
    // --- Check arguments

    if (patch_hierarchy == NULL) {
        PQS_ERROR_STATIC("TimeIntegration", "RK1Step",
                         "'patch_hierarchy' must not be NULL");
    }
    if (u_next_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "RK1Step",
                         "'u_next_id' must be non-negative");
    }
    if (u_current_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "RK1Step",
                         "'u_current_id' must be non-negative");
    }
    if (rhs_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "RK1Step",
                         "'rhs_id' must be non-negative");
    }
    if (dt <= 0) {
        PQS_ERROR_STATIC("TimeIntegration", "RK1Step",
                         "'dt' must be positive");
    }

    // --- Preparations

    int dim = patch_hierarchy->getDim().getValue();

    // --- Advance solution by time step dt

    // loop over PatchHierarchy and take Runge-Kutta step
    // by calling Fortran routines
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    for (int level_num=0 ; level_num < num_levels; level_num++) {
        boost::shared_ptr<hier::PatchLevel> level =
            patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(level->begin());
                pi!=level->end(); pi++) {

            boost::shared_ptr<hier::Patch> patch = *pi;

            // --- Get pointers to data and index space ranges

            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_next_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_next_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_current_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_current_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(rhs_id));

            // ghost box
            hier::Box u_next_ghostbox = u_next_data->getGhostBox();
            const hier::IntVector u_next_ghostbox_lower =
                u_next_ghostbox.lower();
            const hier::IntVector u_next_ghostbox_upper =
                u_next_ghostbox.upper();

            hier::Box u_current_ghostbox = u_current_data->getGhostBox();
            const hier::IntVector u_current_ghostbox_lower =
                u_current_ghostbox.lower();
            const hier::IntVector u_current_ghostbox_upper =
                u_current_ghostbox.upper();

            hier::Box rhs_ghostbox = rhs_data->getGhostBox();
            const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
            const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

            // fill box
            hier::Box fillbox = rhs_data->getBox();
            const hier::IntVector fillbox_lower = fillbox.lower();
            const hier::IntVector fillbox_upper = fillbox.upper();

            // pointers to data
            PQS_REAL* u_next = u_next_data->getPointer();
            PQS_REAL* u_current = u_current_data->getPointer();
            PQS_REAL* rhs = rhs_data->getPointer();


            if (dim == 2) {
                PQS_MATH_2D_RK1_STEP(
                    u_next,
                    &u_next_ghostbox_lower[0],
                    &u_next_ghostbox_upper[0],
                    &u_next_ghostbox_lower[1],
                    &u_next_ghostbox_upper[1],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &dt);
            } else if (dim == 3) {
                PQS_MATH_3D_RK1_STEP(
                    u_next,
                    &u_next_ghostbox_lower[0],
                    &u_next_ghostbox_upper[0],
                    &u_next_ghostbox_lower[1],
                    &u_next_ghostbox_upper[1],
                    &u_next_ghostbox_lower[2],
                    &u_next_ghostbox_upper[2],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    &u_current_ghostbox_lower[2],
                    &u_current_ghostbox_upper[2],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &rhs_ghostbox_lower[2],
                    &rhs_ghostbox_upper[2],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &fillbox_lower[2],
                    &fillbox_upper[2],
                    &dt);
            }

        }  // loop over Patches
    }  // loop over PatchLevels
}  // TimeIntegration::RK1Step()

void TimeIntegration::TVDRK2Stage1(
            boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_stage1_id,
            const int u_current_id,
            const int rhs_id,
            const double dt)
{
    // --- Check arguments

    if (patch_hierarchy == NULL) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage1",
                         "'patch_hierarchy' must not be NULL");
    }
    if (u_stage1_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage1",
                         "'u_stage1_id' must be non-negative");
    }
    if (u_current_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage1",
                         "'u_current_id' must be non-negative");
    }
    if (rhs_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage1",
                         "'rhs_id' must be non-negative");
    }
    if (dt <= 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage1",
                         "'dt' must be positive");
    }

    // --- Preparations

    int dim = patch_hierarchy->getDim().getValue();

    // --- Advance solution by time step dt

    // loop over PatchHierarchy and take Runge-Kutta step
    // by calling Fortran routines
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    for (int level_num=0 ; level_num < num_levels; level_num++) {
        boost::shared_ptr<hier::PatchLevel> level =
            patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(level->begin());
                pi!=level->end(); pi++) {

            boost::shared_ptr<hier::Patch> patch = *pi;

            // --- Get pointers to data and index space ranges

            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_stage1_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_stage1_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_current_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_current_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(rhs_id));

            // ghost box
            hier::Box u_stage1_ghostbox = u_stage1_data->getGhostBox();
            const hier::IntVector u_stage1_ghostbox_lower =
                u_stage1_ghostbox.lower();
            const hier::IntVector u_stage1_ghostbox_upper =
                u_stage1_ghostbox.upper();

            hier::Box u_current_ghostbox = u_current_data->getGhostBox();
            const hier::IntVector u_current_ghostbox_lower =
                u_current_ghostbox.lower();
            const hier::IntVector u_current_ghostbox_upper =
                u_current_ghostbox.upper();

            hier::Box rhs_ghostbox = rhs_data->getGhostBox();
            const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
            const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

            // fill box
            hier::Box fillbox = rhs_data->getBox();
            const hier::IntVector fillbox_lower = fillbox.lower();
            const hier::IntVector fillbox_upper = fillbox.upper();

            // pointers to data
            PQS_REAL* u_stage1 = u_stage1_data->getPointer();
            PQS_REAL* u_current = u_current_data->getPointer();
            PQS_REAL* rhs = rhs_data->getPointer();


            if (dim == 2) {
                PQS_MATH_2D_TVD_RK2_STAGE1(
                    u_stage1,
                    &u_stage1_ghostbox_lower[0],
                    &u_stage1_ghostbox_upper[0],
                    &u_stage1_ghostbox_lower[1],
                    &u_stage1_ghostbox_upper[1],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &dt);
            } else if (dim == 3) {
                PQS_MATH_3D_TVD_RK2_STAGE1(
                    u_stage1,
                    &u_stage1_ghostbox_lower[0],
                    &u_stage1_ghostbox_upper[0],
                    &u_stage1_ghostbox_lower[1],
                    &u_stage1_ghostbox_upper[1],
                    &u_stage1_ghostbox_lower[2],
                    &u_stage1_ghostbox_upper[2],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    &u_current_ghostbox_lower[2],
                    &u_current_ghostbox_upper[2],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &rhs_ghostbox_lower[2],
                    &rhs_ghostbox_upper[2],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &fillbox_lower[2],
                    &fillbox_upper[2],
                    &dt);
            }

        }  // loop over Patches
    }  // loop over PatchLevels
}  // TimeIntegration::TVDRK2Stage1(

void TimeIntegration::TVDRK2Stage2(
            boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_stage1_id,
            const int u_current_id,
            const int rhs_id,
            const double dt)
{
    // --- Check arguments

    if (patch_hierarchy == NULL) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage2",
                         "'patch_hierarchy' must not be NULL");
    }
    if (u_next_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage2",
                         "'u_next_id' must be non-negative");
    }
    if (u_stage1_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage2",
                         "'u_stage1_id' must be non-negative");
    }
    if (u_current_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage2",
                         "'u_current_id' must be non-negative");
    }
    if (rhs_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage2",
                         "'rhs_id' must be non-negative");
    }
    if (dt <= 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK2Stage2",
                         "'dt' must be positive");
    }

    // --- Preparations

    int dim = patch_hierarchy->getDim().getValue();

    // --- Advance solution by time step dt

    // loop over PatchHierarchy and take Runge-Kutta step
    // by calling Fortran routines
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    for (int level_num=0 ; level_num < num_levels; level_num++) {
        boost::shared_ptr<hier::PatchLevel> level =
            patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(level->begin());
                pi!=level->end(); pi++) {

            boost::shared_ptr<hier::Patch> patch = *pi;

            // --- Get pointers to data and index space ranges

            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_next_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_next_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_stage1_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_stage1_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_current_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_current_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(rhs_id));

            // ghost box
            hier::Box u_next_ghostbox = u_next_data->getGhostBox();
            const hier::IntVector u_next_ghostbox_lower =
                u_next_ghostbox.lower();
            const hier::IntVector u_next_ghostbox_upper =
                u_next_ghostbox.upper();

            hier::Box u_stage1_ghostbox = u_stage1_data->getGhostBox();
            const hier::IntVector u_stage1_ghostbox_lower =
                u_stage1_ghostbox.lower();
            const hier::IntVector u_stage1_ghostbox_upper =
                u_stage1_ghostbox.upper();

            hier::Box u_current_ghostbox = u_current_data->getGhostBox();
            const hier::IntVector u_current_ghostbox_lower =
                u_current_ghostbox.lower();
            const hier::IntVector u_current_ghostbox_upper =
                u_current_ghostbox.upper();

            hier::Box rhs_ghostbox = rhs_data->getGhostBox();
            const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
            const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

            // fill box
            hier::Box fillbox = rhs_data->getBox();
            const hier::IntVector fillbox_lower = fillbox.lower();
            const hier::IntVector fillbox_upper = fillbox.upper();

            // pointers to data
            PQS_REAL* u_next = u_next_data->getPointer();
            PQS_REAL* u_stage1 = u_stage1_data->getPointer();
            PQS_REAL* u_current = u_current_data->getPointer();
            PQS_REAL* rhs = rhs_data->getPointer();


            if (dim == 2) {
                PQS_MATH_2D_TVD_RK2_STAGE2(
                    u_next,
                    &u_next_ghostbox_lower[0],
                    &u_next_ghostbox_upper[0],
                    &u_next_ghostbox_lower[1],
                    &u_next_ghostbox_upper[1],
                    u_stage1,
                    &u_stage1_ghostbox_lower[0],
                    &u_stage1_ghostbox_upper[0],
                    &u_stage1_ghostbox_lower[1],
                    &u_stage1_ghostbox_upper[1],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &dt);
            } else if (dim == 3) {
                PQS_MATH_3D_TVD_RK2_STAGE2(
                    u_next,
                    &u_next_ghostbox_lower[0],
                    &u_next_ghostbox_upper[0],
                    &u_next_ghostbox_lower[1],
                    &u_next_ghostbox_upper[1],
                    &u_next_ghostbox_lower[2],
                    &u_next_ghostbox_upper[2],
                    u_stage1,
                    &u_stage1_ghostbox_lower[0],
                    &u_stage1_ghostbox_upper[0],
                    &u_stage1_ghostbox_lower[1],
                    &u_stage1_ghostbox_upper[1],
                    &u_stage1_ghostbox_lower[2],
                    &u_stage1_ghostbox_upper[2],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    &u_current_ghostbox_lower[2],
                    &u_current_ghostbox_upper[2],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &rhs_ghostbox_lower[2],
                    &rhs_ghostbox_upper[2],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &fillbox_lower[2],
                    &fillbox_upper[2],
                    &dt);
            }

        }  // loop over Patches
    }  // loop over PatchLevels
}  // TimeIntegration::TVDRK2Stage2(

void TimeIntegration::TVDRK3Stage1(
            boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_stage1_id,
            const int u_current_id,
            const int rhs_id,
            const double dt)
{
    // --- Check arguments

    if (patch_hierarchy == NULL) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage1",
                         "'patch_hierarchy' must not be NULL");
    }
    if (u_stage1_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage1",
                         "'u_stage1_id' must be non-negative");
    }
    if (u_current_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage1",
                         "'u_current_id' must be non-negative");
    }
    if (rhs_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage1",
                         "'rhs_id' must be non-negative");
    }
    if (dt <= 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage1",
                         "'dt' must be positive");
    }

    // --- Preparations

    int dim = patch_hierarchy->getDim().getValue();

    // --- Advance solution by time step dt

    // loop over PatchHierarchy and take Runge-Kutta step
    // by calling Fortran routines
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    for (int level_num=0 ; level_num < num_levels; level_num++) {
        boost::shared_ptr<hier::PatchLevel> level =
            patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(level->begin());
                pi!=level->end(); pi++) {

            boost::shared_ptr<hier::Patch> patch = *pi;

            // --- Get pointers to data and index space ranges

            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_stage1_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_stage1_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_current_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_current_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(rhs_id));

            // ghost box
            hier::Box u_stage1_ghostbox = u_stage1_data->getGhostBox();
            const hier::IntVector u_stage1_ghostbox_lower =
                u_stage1_ghostbox.lower();
            const hier::IntVector u_stage1_ghostbox_upper =
                u_stage1_ghostbox.upper();

            hier::Box u_current_ghostbox = u_current_data->getGhostBox();
            const hier::IntVector u_current_ghostbox_lower =
                u_current_ghostbox.lower();
            const hier::IntVector u_current_ghostbox_upper =
                u_current_ghostbox.upper();

            hier::Box rhs_ghostbox = rhs_data->getGhostBox();
            const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
            const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

            // fill box
            hier::Box fillbox = rhs_data->getBox();
            const hier::IntVector fillbox_lower = fillbox.lower();
            const hier::IntVector fillbox_upper = fillbox.upper();

            // pointers to data
            PQS_REAL* u_stage1 = u_stage1_data->getPointer();
            PQS_REAL* u_current = u_current_data->getPointer();
            PQS_REAL* rhs = rhs_data->getPointer();


            if (dim == 2) {
                PQS_MATH_2D_TVD_RK3_STAGE1(
                    u_stage1,
                    &u_stage1_ghostbox_lower[0],
                    &u_stage1_ghostbox_upper[0],
                    &u_stage1_ghostbox_lower[1],
                    &u_stage1_ghostbox_upper[1],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &dt);
            } else if (dim == 3) {
                PQS_MATH_3D_TVD_RK3_STAGE1(
                    u_stage1,
                    &u_stage1_ghostbox_lower[0],
                    &u_stage1_ghostbox_upper[0],
                    &u_stage1_ghostbox_lower[1],
                    &u_stage1_ghostbox_upper[1],
                    &u_stage1_ghostbox_lower[2],
                    &u_stage1_ghostbox_upper[2],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    &u_current_ghostbox_lower[2],
                    &u_current_ghostbox_upper[2],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &rhs_ghostbox_lower[2],
                    &rhs_ghostbox_upper[2],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &fillbox_lower[2],
                    &fillbox_upper[2],
                    &dt);
            }

        }  // loop over Patches
    }  // loop over PatchLevels
}  // TimeIntegration::TVDRK3Stage1(

void TimeIntegration::TVDRK3Stage2(
            boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_stage1_id,
            const int u_current_id,
            const int rhs_id,
            const double dt)
{
    // --- Check arguments

    if (patch_hierarchy == NULL) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage2",
                         "'patch_hierarchy' must not be NULL");
    }
    if (u_next_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage2",
                         "'u_next_id' must be non-negative");
    }
    if (u_stage1_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage2",
                         "'u_stage1_id' must be non-negative");
    }
    if (u_current_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage2",
                         "'u_current_id' must be non-negative");
    }
    if (rhs_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage2",
                         "'rhs_id' must be non-negative");
    }
    if (dt <= 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage2",
                         "'dt' must be positive");
    }

    // --- Preparations

    int dim = patch_hierarchy->getDim().getValue();

    // --- Advance solution by time step dt

    // loop over PatchHierarchy and take Runge-Kutta step
    // by calling Fortran routines
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    for (int level_num=0 ; level_num < num_levels; level_num++) {
        boost::shared_ptr<hier::PatchLevel> level =
            patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(level->begin());
                pi!=level->end(); pi++) {

            boost::shared_ptr<hier::Patch> patch = *pi;

            // --- Get pointers to data and index space ranges

            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_next_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_next_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_stage1_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_stage1_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_current_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_current_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(rhs_id));

            // ghost box
            hier::Box u_next_ghostbox = u_next_data->getGhostBox();
            const hier::IntVector u_next_ghostbox_lower =
                u_next_ghostbox.lower();
            const hier::IntVector u_next_ghostbox_upper =
                u_next_ghostbox.upper();

            hier::Box u_stage1_ghostbox = u_stage1_data->getGhostBox();
            const hier::IntVector u_stage1_ghostbox_lower =
                u_stage1_ghostbox.lower();
            const hier::IntVector u_stage1_ghostbox_upper =
                u_stage1_ghostbox.upper();

            hier::Box u_current_ghostbox = u_current_data->getGhostBox();
            const hier::IntVector u_current_ghostbox_lower =
                u_current_ghostbox.lower();
            const hier::IntVector u_current_ghostbox_upper =
                u_current_ghostbox.upper();

            hier::Box rhs_ghostbox = rhs_data->getGhostBox();
            const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
            const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

            // fill box
            hier::Box fillbox = rhs_data->getBox();
            const hier::IntVector fillbox_lower = fillbox.lower();
            const hier::IntVector fillbox_upper = fillbox.upper();

            // pointers to data
            PQS_REAL* u_next = u_next_data->getPointer();
            PQS_REAL* u_stage1 = u_stage1_data->getPointer();
            PQS_REAL* u_current = u_current_data->getPointer();
            PQS_REAL* rhs = rhs_data->getPointer();


            if (dim == 2) {
                PQS_MATH_2D_TVD_RK3_STAGE2(
                    u_next,
                    &u_next_ghostbox_lower[0],
                    &u_next_ghostbox_upper[0],
                    &u_next_ghostbox_lower[1],
                    &u_next_ghostbox_upper[1],
                    u_stage1,
                    &u_stage1_ghostbox_lower[0],
                    &u_stage1_ghostbox_upper[0],
                    &u_stage1_ghostbox_lower[1],
                    &u_stage1_ghostbox_upper[1],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &dt);
            } else if (dim == 3) {
                PQS_MATH_3D_TVD_RK3_STAGE2(
                    u_next,
                    &u_next_ghostbox_lower[0],
                    &u_next_ghostbox_upper[0],
                    &u_next_ghostbox_lower[1],
                    &u_next_ghostbox_upper[1],
                    &u_next_ghostbox_lower[2],
                    &u_next_ghostbox_upper[2],
                    u_stage1,
                    &u_stage1_ghostbox_lower[0],
                    &u_stage1_ghostbox_upper[0],
                    &u_stage1_ghostbox_lower[1],
                    &u_stage1_ghostbox_upper[1],
                    &u_stage1_ghostbox_lower[2],
                    &u_stage1_ghostbox_upper[2],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    &u_current_ghostbox_lower[2],
                    &u_current_ghostbox_upper[2],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &rhs_ghostbox_lower[2],
                    &rhs_ghostbox_upper[2],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &fillbox_lower[2],
                    &fillbox_upper[2],
                    &dt);
            }

        }  // loop over Patches
    }  // loop over PatchLevels
}  // TimeIntegration::TVDRK3Stage2(

void TimeIntegration::TVDRK3Stage3(
            boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
            const int u_next_id,
            const int u_stage2_id,
            const int u_current_id,
            const int rhs_id,
            const double dt)
{
    // --- Check arguments

    if (patch_hierarchy == NULL) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage3",
                         "'patch_hierarchy' must not be NULL");
    }
    if (u_next_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage3",
                         "'u_next_id' must be non-negative");
    }
    if (u_stage2_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage3",
                         "'u_stage2_id' must be non-negative");
    }
    if (u_current_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage3",
                         "'u_current_id' must be non-negative");
    }
    if (rhs_id < 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage3",
                         "'rhs_id' must be non-negative");
    }
    if (dt <= 0) {
        PQS_ERROR_STATIC("TimeIntegration", "TVDRK3Stage3",
                         "'dt' must be positive");
    }

    // --- Preparations

    int dim = patch_hierarchy->getDim().getValue();

    // --- Advance solution by time step dt

    // loop over PatchHierarchy and take Runge-Kutta step
    // by calling Fortran routines
    const int num_levels = patch_hierarchy->getNumberOfLevels();
    for (int level_num=0 ; level_num < num_levels; level_num++) {
        boost::shared_ptr<hier::PatchLevel> level =
            patch_hierarchy->getPatchLevel(level_num);

        for (hier::PatchLevel::Iterator pi(level->begin());
                pi!=level->end(); pi++) {

            boost::shared_ptr<hier::Patch> patch = *pi;

            // --- Get pointers to data and index space ranges

            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_next_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_next_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_stage2_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_stage2_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > u_current_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(u_current_id));
            boost::shared_ptr< pdat::CellData<PQS_REAL> > rhs_data =
                BOOST_CAST<pdat::CellData<PQS_REAL>, hier::PatchData>(
                    patch->getPatchData(rhs_id));

            // ghost box
            hier::Box u_next_ghostbox = u_next_data->getGhostBox();
            const hier::IntVector u_next_ghostbox_lower =
                u_next_ghostbox.lower();
            const hier::IntVector u_next_ghostbox_upper =
                u_next_ghostbox.upper();

            hier::Box u_stage2_ghostbox = u_stage2_data->getGhostBox();
            const hier::IntVector u_stage2_ghostbox_lower =
                u_stage2_ghostbox.lower();
            const hier::IntVector u_stage2_ghostbox_upper =
                u_stage2_ghostbox.upper();

            hier::Box u_current_ghostbox = u_current_data->getGhostBox();
            const hier::IntVector u_current_ghostbox_lower =
                u_current_ghostbox.lower();
            const hier::IntVector u_current_ghostbox_upper =
                u_current_ghostbox.upper();

            hier::Box rhs_ghostbox = rhs_data->getGhostBox();
            const hier::IntVector rhs_ghostbox_lower = rhs_ghostbox.lower();
            const hier::IntVector rhs_ghostbox_upper = rhs_ghostbox.upper();

            // fill box
            hier::Box fillbox = rhs_data->getBox();
            const hier::IntVector fillbox_lower = fillbox.lower();
            const hier::IntVector fillbox_upper = fillbox.upper();

            // pointers to data
            PQS_REAL* u_next = u_next_data->getPointer();
            PQS_REAL* u_stage2 = u_stage2_data->getPointer();
            PQS_REAL* u_current = u_current_data->getPointer();
            PQS_REAL* rhs = rhs_data->getPointer();


            if (dim == 2) {
                PQS_MATH_2D_TVD_RK3_STAGE3(
                    u_next,
                    &u_next_ghostbox_lower[0],
                    &u_next_ghostbox_upper[0],
                    &u_next_ghostbox_lower[1],
                    &u_next_ghostbox_upper[1],
                    u_stage2,
                    &u_stage2_ghostbox_lower[0],
                    &u_stage2_ghostbox_upper[0],
                    &u_stage2_ghostbox_lower[1],
                    &u_stage2_ghostbox_upper[1],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &dt);
            } else if (dim == 3) {
                PQS_MATH_3D_TVD_RK3_STAGE3(
                    u_next,
                    &u_next_ghostbox_lower[0],
                    &u_next_ghostbox_upper[0],
                    &u_next_ghostbox_lower[1],
                    &u_next_ghostbox_upper[1],
                    &u_next_ghostbox_lower[2],
                    &u_next_ghostbox_upper[2],
                    u_stage2,
                    &u_stage2_ghostbox_lower[0],
                    &u_stage2_ghostbox_upper[0],
                    &u_stage2_ghostbox_lower[1],
                    &u_stage2_ghostbox_upper[1],
                    &u_stage2_ghostbox_lower[2],
                    &u_stage2_ghostbox_upper[2],
                    u_current,
                    &u_current_ghostbox_lower[0],
                    &u_current_ghostbox_upper[0],
                    &u_current_ghostbox_lower[1],
                    &u_current_ghostbox_upper[1],
                    &u_current_ghostbox_lower[2],
                    &u_current_ghostbox_upper[2],
                    rhs,
                    &rhs_ghostbox_lower[0],
                    &rhs_ghostbox_upper[0],
                    &rhs_ghostbox_lower[1],
                    &rhs_ghostbox_upper[1],
                    &rhs_ghostbox_lower[2],
                    &rhs_ghostbox_upper[2],
                    &fillbox_lower[0],
                    &fillbox_upper[0],
                    &fillbox_lower[1],
                    &fillbox_upper[1],
                    &fillbox_lower[2],
                    &fillbox_upper[2],
                    &dt);
            }

        }  // loop over Patches
    }  // loop over PatchLevels
}  // TimeIntegration::TVDRK3Stage3(

} // PQS::math namespace
} // PQS namespace
