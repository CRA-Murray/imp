/**
 *  \file MinimumSingletonRestraint.cpp
 *  \brief The minimum score over a container.
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2011 IMP Inventors. All rights reserved.
 *
 */

#include "IMP/container/MinimumSingletonRestraint.h"
#include "IMP/algebra/internal/MinimalSet.h"
#include <IMP/internal/container_helpers.h>
#include <IMP/core/SingletonRestraint.h>

IMPCONTAINER_BEGIN_NAMESPACE

MinimumSingletonRestraint
::MinimumSingletonRestraint(SingletonScore *f,
                           SingletonContainer *c,
                           unsigned int n,
                           std::string name):
  Restraint(name),
  f_(f), c_(c), n_(n){
}

namespace {
  typedef algebra::internal::MinimalSet<double,
          Particle*, std::less<double> > SingletonMinimumMS;
  template <class C, class F>
  SingletonMinimumMS find_minimal_set_SingletonMinimum(C* c, F *f,
                                                         unsigned int n) {
    IMP_LOG(VERBOSE, "Finding Minimum " << n << " of "
            << c->get_number() << std::endl);
    SingletonMinimumMS bestn(n);
    IMP_FOREACH_SINGLETON(c, {
        double score= f->evaluate(_1, NULL);
        IMP_LOG(VERBOSE, "Found " << score << " for "
                << _1 << std::endl);
        bestn.insert(score, _1);
      });
    return bestn;
  }
}

double MinimumSingletonRestraint
::unprotected_evaluate(DerivativeAccumulator *da) const {
  IMP_OBJECT_LOG;
  SingletonMinimumMS bestn
    = find_minimal_set_SingletonMinimum(c_.get(),
                                         f_.get(), n_);

  double score=0;
  for (unsigned int i=0; i< bestn.size(); ++i) {
    if (da) {
      f_->evaluate(bestn[i].second, da);
    }
    score+= bestn[i].first;
  }
  IMP_LOG(VERBOSE, "Total score is " << score << std::endl);
  return score;
}

double MinimumSingletonRestraint
::unprotected_evaluate_if_good(DerivativeAccumulator *da,
                               double max) const {
  IMP_OBJECT_LOG;
  SingletonMinimumMS bestn
    = find_minimal_set_SingletonMinimum(c_.get(),
                                         f_.get(), n_);

  double score=0;
  for (unsigned int i=0; i< bestn.size(); ++i) {
    if (da) {
      f_->evaluate(bestn[i].second, da);
    }
    score+= bestn[i].first;
    if (score > max) break;
  }
  IMP_LOG(VERBOSE, "Total score is " << score << std::endl);
  return score;
}

Restraints MinimumSingletonRestraint
::do_create_current_decomposition() const {
  IMP_OBJECT_LOG;
  SingletonMinimumMS bestn
    = find_minimal_set_SingletonMinimum(c_.get(),
                                         f_.get(), n_);
  Restraints ret;
  for (unsigned int i=0; i< bestn.size(); ++i) {
    ret.push_back(new core::SingletonRestraint(f_, bestn[i].second));
  }
  return ret;
}


void MinimumSingletonRestraint::do_show(std::ostream &out) const {
  out << "container " << *c_ << std::endl;
  out << "score " << *f_ << std::endl;
}


ParticlesTemp MinimumSingletonRestraint::get_input_particles() const
{
  return IMP::internal::get_input_particles(f_.get(),
                                            c_->get_contained_particles());
}

ContainersTemp MinimumSingletonRestraint::get_input_containers() const
{
  return ContainersTemp(1, c_);
}


IMPCONTAINER_END_NAMESPACE
