/**
 *  \file FilteredListGroupnameContainer.h
 *  \brief Store a list of Classnames filtered based on another list
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPCORE_FILTERED_LIST_GROUPNAME_CONTAINER_H
#define IMPCORE_FILTERED_LIST_GROUPNAME_CONTAINER_H

#include "core_exports.h"
#include "GroupnameContainer.h"
#include "internal/core_version_info.h"

IMPCORE_BEGIN_NAMESPACE

//! Store a list of Classnames filtered by other lists
/** This class stores a list of Classnames and a list of
    GroupnameContainers with the invariant that none of the
    GroupnameContainers contain any of the Classnames stored.

    \note Currently the filter is only applied upon addition
    of a Value to the container. So adding more sets to the
    filter afterwards won't remove objects. Nor will changing
    the filtering sets.

    \note The indexes can change when particles are inserted
    as the list is maintained in sorted order.

    \verbinclude simple_examples/groupname_filtered_container.py
 */
class IMPCOREEXPORT FilteredListGroupnameContainer
  : public GroupnameContainer
{
  std::vector<Value> data_;
public:
  //! cannot pass a Groupnames on construction
  FilteredListGroupnameContainer();

  virtual ~FilteredListGroupnameContainer();

  IMP_GROUPNAME_CONTAINER(internal::core_version_info);

  //! Add vt if none of the referenced containers already contains it
  void add_classname(Value vt);

  //! remove all objects from the container
  void clear_classnames() {
    data_.clear();
  }

  IMP_CONTAINER(GroupnameContainer, groupname_container,
                GroupnameContainerIndex);
};


IMPCORE_END_NAMESPACE

#endif  /* IMPCORE_FILTERED_LIST_GROUPNAME_CONTAINER_H */
