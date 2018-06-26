/*! \file GridManager.cc
 *
 * \brief
 * Implementation file for GridManager class.
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
#include <cstddef>
#include <sstream>

// Boost
#include <boost/smart_ptr/shared_ptr.hpp>

// SAMRAI
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"

// PQS
#include "PQS/PQS_config.h"  // IWYU pragma: keep
#include "PQS/pqs/GridManager.h"
#include "PQS/pqs/DataInitStrategy.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace tbox { class Database; } }


// --- Implementation for PQS::pqs::GridManager methods

namespace PQS {
namespace pqs {

// --- Implementation of public PQS::pqs::GridManager methods

// Constructor
GridManager::GridManager(
    boost::shared_ptr<tbox::Database> config_db,
    boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy):
    mesh::TagAndInitializeStrategy(d_object_name)
{
    // Check parameters
    if (patch_hierarchy == NULL) {
        // TODO
    }

    d_patch_hierarchy = patch_hierarchy;
}

void GridManager::initializeLevelData(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const int patch_level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const boost::shared_ptr<hier::PatchLevel>& old_patch_level,
    const bool allocate_data)
{
    // TODO
}

void GridManager::resetHierarchyConfiguration(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const int coarsest_patch_level_number,
    const int finest_patch_level_number)
{
    // TODO
}

void GridManager::tagCellsForRefinement(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const int patch_level_number,
    const int regrid_cycle,
    const double regrid_time,
    const int tag_id,
    const bool initial_time,
    const bool coarsest_sync_patch_level,
    const bool can_be_refined,
    const double regrid_start_time)
{
    // TODO
}

void GridManager::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::GridManager::printClassData..." << endl;
    os << "(GridManager*) this = " << (GridManager*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;
    os << "d_data_init_strategy = " << d_data_init_strategy.get() << endl;

    os << endl;
    d_data_init_strategy->printClassData(os);
}

bool GridManager::refineUserBoxInputOnly(int cycle, double time)
{
    // TODO
    return false;
}

bool GridManager::getUserSuppliedRefineBoxes(
    hier::BoxContainer& refine_boxes,
    const int patch_level_number,
    const int cycle,
    const double time)
{
    // TODO
    return false;
}

void GridManager::resetRefineBoxes(
    const hier::BoxContainer& refine_boxes,
    const int patch_level_number)
{
}

// --- Implementation of private PQS::pqs::GridManager methods

// Copy constructor
GridManager::GridManager(const GridManager& rhs):
    mesh::TagAndInitializeStrategy(d_object_name)
{}

} // PQS::pqs namespace
} // PQS namespace
