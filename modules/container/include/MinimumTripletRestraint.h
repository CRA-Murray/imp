/**
 *  \file MinimumTripletRestraint.h
 *  \brief Score based on the minimum score over a set of Triplets
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2011 IMP Inventors. All rights reserved.
 */

#ifndef IMPCONTAINER_MINIMUM_TRIPLET_RESTRAINT_H
#define IMPCONTAINER_MINIMUM_TRIPLET_RESTRAINT_H

#include "container_config.h"
#include <IMP/Restraint.h>
#include <IMP/TripletScore.h>
#include <IMP/TripletContainer.h>

IMPCONTAINER_BEGIN_NAMESPACE

//! Score based on the min or max TripletScore over a set
/** The score is evaluated for each of the VALUETYPE in the container
    and the value of the min or max n scores is used. That is,
    if n is 1, the value of the restraint is the value of the min or max
    score over the container.
 */
class IMPCONTAINEREXPORT MinimumTripletRestraint
: public Restraint
{
  IMP::OwnerPointer<TripletScore> f_;
  IMP::OwnerPointer<TripletContainer> c_;
  unsigned int n_;
  double unprotected_evaluate_if_good(DerivativeAccumulator *da,
                                      double max) const;
public:
  /** n is the number of LCMinimum scores to use.
   */
  MinimumTripletRestraint(TripletScore *f,
                                 TripletContainer *c,
                                 unsigned int n=1,
                                 std::string name
                                 ="MinimumTripletRestraint %1%");

  IMP_RESTRAINT(MinimumTripletRestraint);

  //! Set the number of lowest scores to use.
  void set_n(unsigned int n) { n_=n;}
#ifndef IMP_DOXYGEN
  Restraints do_create_current_decomposition() const;
#endif
};

IMP_OBJECTS(MinimumTripletRestraint,MinimumTripletRestraints);


IMPCONTAINER_END_NAMESPACE

#endif  /* IMPCONTAINER_MINIMUM_TRIPLET_RESTRAINT_H */
