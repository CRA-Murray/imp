/**
 *  \file SingletonContainerSet.h
 *  \brief Store a set of SingletonContainers
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPCONTAINER_SINGLETON_CONTAINER_SET_H
#define IMPCONTAINER_SINGLETON_CONTAINER_SET_H

#include "container_config.h"
#include <IMP/SingletonContainer.h>
#include <IMP/container_macros.h>
#include <IMP/internal/container_helpers.h>
#include <IMP/scoped.h>

IMPCONTAINER_BEGIN_NAMESPACE

//! Stores a set of SingletonContainers
/** The input sets must be disjoint. This can change if there is
    demand for it.

    \usesconstraint
*/
class IMPCONTAINEREXPORT SingletonContainerSet
  : public SingletonContainer
{
  IMP_CONTAINER_DEPENDENCIES(SingletonContainerSet,
                             {
                               ret.insert(ret.end(),
                                          back_->singleton_containers_begin(),
                                          back_->singleton_containers_end());
                             });
  static SingletonContainerSet* get_set(SingletonContainer* c) {
    return dynamic_cast<SingletonContainerSet*>(c);
  }
 public:
  //! Construct and empty set
  SingletonContainerSet(Model *m,
                        std::string name="SingletonContainerSet %1%");

  SingletonContainerSet(const SingletonContainersTemp &pc,
                        std::string name="SingletonContainerSet %1%");

  bool get_contains_particle(Particle*) const;
  void apply(const SingletonModifier *sm) const;
  void apply(const SingletonDerivativeModifier *sm,
             DerivativeAccumulator &da) const;
  double evaluate(const SingletonScore *s,
                  DerivativeAccumulator *da) const;
  double evaluate_if_good(const SingletonScore *s,
                          DerivativeAccumulator *da,
                          double max) const;
 template <class SM>
  void template_apply(const SM *sm,
                      DerivativeAccumulator &da) const {
   for (unsigned int i=0; i< get_number_of_singleton_containers(); ++i) {
     get_singleton_container(i)->apply(sm, da);
   }
 }
  template <class SM>
  void template_apply(const SM *sm) const {
    for (unsigned int i=0; i< get_number_of_singleton_containers(); ++i) {
      get_singleton_container(i)->apply(sm);
    }
  }
  template <class SS>
  double template_evaluate(const SS *s,
                           DerivativeAccumulator *da) const {
    double ret=0;
    for (unsigned int i=0; i< get_number_of_singleton_containers(); ++i) {
      ret+=get_singleton_container(i)->evaluate(s, da);
    }
    return ret;
  }
  template <class SS>
    double template_evaluate_if_good(const SS *s,
                                 DerivativeAccumulator *da, double max) const {
    double ret=0;
    for (unsigned int i=0; i< get_number_of_singleton_containers(); ++i) {
      double cur=get_singleton_container(i)->evaluate_if_good(s, da, max);
      ret+=cur;
      max-=cur;
      if (max < 0) break;
    }
    return ret;
  }
  bool get_is_changed() const;
  ParticlesTemp get_all_possible_particles() const;
  IMP_OBJECT(SingletonContainerSet);

  /** @name Methods to control the nested container

      This container merges a set of nested containers. To add
      or remove nested containers, use the methods below.
  */
  /**@{*/
  IMP_LIST_ACTION(public, SingletonContainer, SingletonContainers,
                  singleton_container, singleton_containers,
                  SingletonContainer*, SingletonContainers,
                  {
                    obj->set_was_used(true);
                    set_is_changed(true);
                    get_model()->reset_dependencies();
                  },{},
                  );
  /**@}*/
#ifndef IMP_DOXYGEN
  ParticleIndexes get_indexes() const;
  ParticleIndexes get_all_possible_indexes() const;
  ContainersTemp get_input_containers() const;
  void do_before_evaluate();
#endif
};


IMPCONTAINER_END_NAMESPACE

#endif  /* IMPCONTAINER_SINGLETON_CONTAINER_SET_H */
