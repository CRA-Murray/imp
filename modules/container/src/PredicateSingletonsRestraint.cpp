/**
 *  \file PredicateSingletonsRestraint.cpp   \brief Container for singleton.
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include "IMP/container/PredicateSingletonsRestraint.h"
#include <IMP/singleton_macros.h>

IMPCONTAINER_BEGIN_NAMESPACE

PredicateSingletonsRestraint
::PredicateSingletonsRestraint(SingletonPredicate *pred,
                               SingletonContainerAdaptor input,
                               std::string name):
  Restraint(input->get_model(), name),
  predicate_(pred), input_(input), updated_(false),
  error_on_unknown_(true){}

double
PredicateSingletonsRestraint
::unprotected_evaluate(DerivativeAccumulator *da) const {
  update_lists_if_necessary();
  double ret=0;
  for (unsigned int i=0; i< restraints_.size(); ++i) {
    ret+=restraints_[i]->unprotected_evaluate(da);
  }
  return ret;
}

ParticlesTemp PredicateSingletonsRestraint
::get_input_particles() const {
  // not correct, but correct is complicated
  return input_->get_all_possible_particles();
}
ContainersTemp PredicateSingletonsRestraint
::get_input_containers() const {
  // List containers don't do anything interesting
  return ContainersTemp(1, input_);
}

Restraints PredicateSingletonsRestraint
::do_create_current_decomposition() const {
  Restraints ret;
  for (unsigned int i=0; i< restraints_.size(); ++i) {
    Pointer<Restraint> r=restraints_[i]->create_current_decomposition();
    if (r) {
      RestraintSet *rs= dynamic_cast<RestraintSet*>(r.get());
      if (rs) {
        ret+=rs->get_restraints();
      } else {
        ret.push_back(r);
      }
    }
  }
  return ret;
}

bool PredicateSingletonsRestraint
::assign_pair(ParticleIndex index) const {
  int bin=predicate_->get_value_index(get_model(), index);
  Map::const_iterator it= containers_.find(bin);
  if (it == containers_.end()) {
    if (unknown_container_) {
      unknown_container_->add_particle(index);
      return true;
    } else if (error_on_unknown_) {
      IMP_THROW("Invalid predicate value of " << bin
                << " encounted for " << index,
                ValueException);
      return true;
    } else {
      return false;
    }
  } else {
    it->second->add_particle(index);
    return true;
  }
}
void PredicateSingletonsRestraint
::update_lists_if_necessary() const {
  if (updated_ && !input_->get_is_changed()) return;
  updated_=true;
  if (unknown_container_) {
    unknown_container_->clear_particles();
  }
  for (Map::const_iterator it= containers_.begin();
       it != containers_.end(); ++it) {
    it->second->clear_particles();
  }
  int dropped=0;
  IMP_FOREACH_SINGLETON_INDEX(input_, {
      bool added=assign_pair(_1);
      if (!added) ++dropped;
    });
  IMP_IF_CHECK(USAGE_AND_INTERNAL) {
    unsigned int total=dropped;
    for (Map::const_iterator it= containers_.begin();
         it != containers_.end(); ++it) {
      total+=it->second->get_number();
    }
    total+= unknown_container_->get_number();
    IMP_INTERNAL_CHECK(input_->get_number()==total,
                       "Wrong number of particles "
                       << total << "!=" << input_->get_number());
  }
}

void PredicateSingletonsRestraint::do_show(std::ostream &) const {
}
IMPCONTAINER_END_NAMESPACE
