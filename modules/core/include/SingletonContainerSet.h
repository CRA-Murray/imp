/**
 *  \file SingletonContainerSet.h
 *  \brief Store a set of SingletonContainers
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-9 Sali Lab. All rights reserved.
 */

#ifndef IMPCORE_SINGLETON_CONTAINER_SET_H
#define IMPCORE_SINGLETON_CONTAINER_SET_H

#include "config.h"
#include <IMP/SingletonContainer.h>
#include <IMP/container_macros.h>

IMPCORE_BEGIN_NAMESPACE

//! Stores a set of SingletonContainers
/**
 */
class IMPCOREEXPORT SingletonContainerSet
  : public SingletonContainer
{
  unsigned int rev_;
public:
  //! Construct and empty set
  SingletonContainerSet(std::string name="SingletonContainerSet %1%");

  IMP_SINGLETON_CONTAINER(SingletonContainerSet, get_module_version_info());
 /** @name Methods to control the nested container

     This container merges a set of nested containers. To add
     or remove nested containers, use the methods below.
  */
  /**@{*/
  IMP_LIST(public, SingletonContainer, singleton_container,
           SingletonContainer*, SingletonContainers);
  /**@}*/
};


IMPCORE_END_NAMESPACE

#endif  /* IMPCORE_SINGLETON_CONTAINER_SET_H */
