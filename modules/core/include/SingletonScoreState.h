/**
 *  \file SingletonScoreState.h
 *  \brief Use a SingletonModifier applied to a Particles to
 *  maintain an invariant
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 */

#ifndef IMPCORE_SINGLETON_SCORE_STATE_H
#define IMPCORE_SINGLETON_SCORE_STATE_H

#include "config.h"
#include "internal/version_info.h"
#include <IMP/SingletonModifier.h>
#include <IMP/ScoreState.h>

IMP_BEGIN_NAMESPACE
// for swig
class SingletonModifier;
IMP_END_NAMESPACE

IMPCORE_BEGIN_NAMESPACE
//! Apply a SingletonFunction to a Singleton
/** \ingroup restraint
    \see SingletonsScoreState
 */
class IMPCOREEXPORT SingletonScoreState : public ScoreState
{
  Pointer<SingletonModifier> f_;
  Pointer<SingletonModifier> af_;
  Particle* v_;
public:
  /** \param[in] v The Particle to modify
      \param[in] before The SingletonModifier to apply to all elements
       before evaluate.
      \param[in] after The SingletonModifier to apply to all elements
       after evaluate.
   */
  SingletonScoreState(SingletonModifier *before,
                       SingletonModifier *after, Particle* v);

  //! Apply this modifier to all the elements after an evaluate
  void set_after_evaluate_modifier(SingletonModifier* f) {
    af_=f;
  }

  //! Apply this modifier to all the elements before an evaluate
  void set_before_evaluate_modifier(SingletonModifier* f) {
    f_=f;
  }

  virtual ~SingletonScoreState();

  IMP_SCORE_STATE(internal::version_info)
};


IMPCORE_END_NAMESPACE

#endif  /* IMPCORE_SINGLETON_SCORE_STATE_H */
