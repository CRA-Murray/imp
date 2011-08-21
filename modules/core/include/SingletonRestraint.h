/**
 *  \file SingletonRestraint.h
 *  \brief Apply a SingletonScore to a Singleton.
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2011 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPCORE_SINGLETON_RESTRAINT_H
#define IMPCORE_SINGLETON_RESTRAINT_H

#include "core_config.h"

#include <IMP/Restraint.h>
#include <IMP/Pointer.h>
#include <IMP/SingletonScore.h>
#include "internal/singleton_helpers.h"
#include <IMP/internal/container_helpers.h>

#include <iostream>

IMPCORE_BEGIN_NAMESPACE

//! Applies a SingletonScore to a Singleton.
/** This restraint stores a Singleton.
    \see SingletonRestraint
 */
class IMPCOREEXPORT SingletonRestraint :
  public SingletonScoreRestraint
{
  IMP::OwnerPointer<SingletonScore> ss_;
  ParticleIndex v_;
public:
  //! Create the restraint.
  /** This function takes the function to apply to the
      stored Singleton and the Singleton.
   */
  SingletonRestraint(SingletonScore *ss,
                     Particle* vt,
                     std::string name="SingletonRestraint %1%");

  SingletonScore* get_score() const {
    return ss_;
  }
  Particle* get_argument() const {
    return IMP::internal::get_particle(get_model(), v_);
  }

  IMP_RESTRAINT(SingletonRestraint);

  double unprotected_evaluate_if_good(DerivativeAccumulator *da,
                                      double max) const;

  Restraints get_current_decomposition() const;
};

IMPCORE_END_NAMESPACE

#endif  /* IMPCORE_SINGLETON_RESTRAINT_H */
