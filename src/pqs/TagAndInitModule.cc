/*! \file TagAndInitModule.cc
 *
 * \brief
 * Implementation file for TagAndInitModule class.
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
#include "PQS/pqs/TagAndInitModule.h"
#include "PQS/pqs/DataInitStrategy.h"

// Class/type declarations
namespace SAMRAI { namespace hier { class PatchHierarchy; } }
namespace SAMRAI { namespace tbox { class Database; } }


// --- Implementation for PQS::pqs::TagAndInitModule methods

namespace PQS {
namespace pqs {

// --- Implementation of public PQS::pqs::TagAndInitModule methods

// Constructor
TagAndInitModule::TagAndInitModule(
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

void TagAndInitModule::initializeLevelData(
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

void TagAndInitModule::resetHierarchyConfiguration(
    const boost::shared_ptr<hier::PatchHierarchy>& patch_hierarchy,
    const int coarsest_patch_level_number,
    const int finest_patch_level_number)
{
    // TODO
}

void TagAndInitModule::tagCellsForRefinement(
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

bool TagAndInitModule::refineUserBoxInputOnly(int cycle, double time)
{
    // TODO
    return false;
}

bool TagAndInitModule::getUserSuppliedRefineBoxes(
    hier::BoxContainer& refine_boxes,
    const int patch_level_number,
    const int cycle,
    const double time)
{
    // TODO
    return false;
}

void TagAndInitModule::resetRefineBoxes(
    const hier::BoxContainer& refine_boxes,
    const int patch_level_number)
{
}

void TagAndInitModule::printClassData(ostream& os) const
{
    os << endl;
    os << "PQS::pqs::TagAndInitModule::printClassData..." << endl;
    os << "(TagAndInitModule*) this = " << (TagAndInitModule*) this << endl;
    os << "d_patch_hierarchy = " << d_patch_hierarchy.get() << endl;
    os << "d_data_init_strategy = " << d_data_init_strategy.get() << endl;

    os << endl;
    d_data_init_strategy->printClassData(os);
}

// --- Implementation of private PQS::pqs::TagAndInitModule methods

// Copy constructor
TagAndInitModule::TagAndInitModule(const TagAndInitModule& rhs):
    mesh::TagAndInitializeStrategy(d_object_name)
{}

} // PQS::pqs namespace
} // PQS namespace
