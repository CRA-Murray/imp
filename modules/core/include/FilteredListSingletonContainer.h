/**
 *  \file FilteredListSingletonContainer.h
 *  \brief Store a list of Particles filtered based on another list
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPCORE_FILTERED_LIST_SINGLETON_CONTAINER_H
#define IMPCORE_FILTERED_LIST_SINGLETON_CONTAINER_H

#include "config.h"
#include "internal/version_info.h"
#include <IMP/SingletonContainer.h>
#include <IMP/container_macros.h>

IMPCORE_BEGIN_NAMESPACE

//! Store a list of Particles filtered by other lists
/** This class stores a list of Particles and a list of
    SingletonContainers with the invariant that none of the
    SingletonContainers contain any of the Particles stored.

    \note Currently the filter is only applied upon addition
    of a Particle* to the container. So adding more sets to the
    filter afterwards won't remove objects. Nor will changing
    the filtering sets.

    \note The indexes can change when particles are inserted
    as the list is maintained in sorted order.

    \see ListSingletonContainer
 */
class IMPCOREEXPORT FilteredListSingletonContainer
  : public SingletonContainer
{
  std::vector<Particle*> data_;
public:
  //! cannot pass a Singletons on construction
  FilteredListSingletonContainer();

  virtual ~FilteredListSingletonContainer();

  IMP_SINGLETON_CONTAINER(internal::version_info);

  //! Add vt if none of the referenced containers already contains it
  void add_particle(Particle* vt);

  //! remove all objects from the container
  void clear_particles() {
    data_.clear();
  }

  IMP_LIST(public, SingletonContainer, singleton_container,
           SingletonContainer*);
};


IMPCORE_END_NAMESPACE

#endif  /* IMPCORE_FILTERED_LIST_SINGLETON_CONTAINER_H */
