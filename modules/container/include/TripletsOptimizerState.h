/**
 *  \file TripletsOptimizerState.h
 *  \brief Use a TripletModifier applied to a ParticleTripletsTemp to
 *  maintain an invariant
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPCONTAINER_TRIPLETS_OPTIMIZER_STATE_H
#define IMPCONTAINER_TRIPLETS_OPTIMIZER_STATE_H

#include "container_config.h"
#include <IMP/TripletContainer.h>
#include <IMP/TripletModifier.h>
#include <IMP/OptimizerState.h>
#include <IMP/optimizer_state_macros.h>

IMP_BEGIN_NAMESPACE
// for swig
class TripletContainer;
class TripletModifier;
IMP_END_NAMESPACE

IMPCONTAINER_BEGIN_NAMESPACE
//! Apply a TripletFunction to a TripletContainer to maintain an invariant
/** \ingroup restraint
    \see TripletOptimizerState
 */
class IMPCONTAINEREXPORT TripletsOptimizerState : public OptimizerState
{
  IMP::OwnerPointer<TripletModifier> f_;
  IMP::OwnerPointer<TripletContainer> c_;
public:
  /** \param[in] c The Container to hold the elements to process
      \param[in] gf The TripletModifier to apply to all elements.
      \param[in] name The name to use for this Object
   */
  TripletsOptimizerState(TripletContainerAdaptor c, TripletModifier *gf,
                           std::string name="TripletsOptimizerState %1%");

  IMP_OPTIMIZER_STATE(TripletsOptimizerState);
};

IMP_OBJECTS(TripletsOptimizerState,TripletsOptimizerStates);


IMPCONTAINER_END_NAMESPACE

#endif  /* IMPCONTAINER_TRIPLETS_OPTIMIZER_STATE_H */
