/**
 *  \file PairsScoreState.h
 *  \brief Use a PairModifier applied to a ParticlePairs to
 *  maintain an invariant
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-9 Sali Lab. All rights reserved.
 */

#ifndef IMPCORE_PAIRS_SCORE_STATE_H
#define IMPCORE_PAIRS_SCORE_STATE_H

#include "config.h"
#include "internal/version_info.h"
#include <IMP/PairContainer.h>
#include <IMP/PairModifier.h>
#include <IMP/ScoreState.h>

IMP_BEGIN_NAMESPACE
// for swig
class PairContainer;
class PairModifier;
IMP_END_NAMESPACE

IMPCORE_BEGIN_NAMESPACE
//! Apply a PairFunction to a PairContainer to maintain an invariant
/** \ingroup restraint
    The score state is passed up to two PairModifiers, one to
    apply before evaluation and the other after. The one after
    should take a DerivativeAccumulator as its last argument for
    PairModifier::apply() and will only be called if
    the score was computed with derivatives.

    An example showing a how to use such a score state to maintain a cover
    of the atoms of a protein by a sphere per residue.
    \htmlinclude cover_particles.py.html
    \see PairScoreState
 */
class IMPCOREEXPORT PairsScoreState : public ScoreState
{
  Pointer<PairModifier> f_;
  Pointer<PairModifier> af_;
  Pointer<PairContainer> c_;
public:
  /** \param[in] c The Container to hold the elements to process
      \param[in] before The PairModifier to apply to all elements
       before evaluate.
      \param[in] after The PairModifier to apply to all elements
       after evaluate.
   */
  PairsScoreState(PairContainer *c, PairModifier *before,
                       PairModifier *after);

  //! Apply this modifier to all the elements after an evaluate
  void set_after_evaluate_modifier(PairModifier* f) {
    af_=f;
  }

  //! Apply this modifier to all the elements before an evaluate
  void set_before_evaluate_modifier(PairModifier* f) {
    f_=f;
  }

  IMP_SCORE_STATE(PairsScoreState, internal::version_info)
};


IMPCORE_END_NAMESPACE

#endif  /* IMPCORE_PAIRS_SCORE_STATE_H */
