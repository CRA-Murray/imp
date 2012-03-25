/**
 *  \file IMP/declare_ScoringFunction.h
 *  \brief Storage of a model, its restraints,
 *                         constraints and particles.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPKERNEL_DECLARE_SCORING_FUNCTION_H
#define IMPKERNEL_DECLARE_SCORING_FUNCTION_H

#include "kernel_config.h"
#include "base_types.h"
#include "dependency_graph.h"
#include "declare_Restraint.h"
#include <IMP/base/tracking.h>
#include <IMP/base/Pointer.h>

#include <limits>


IMP_BEGIN_NAMESPACE
class Model;

/** A ScoringFunction represents a scoring function on the model.
    The Model has a default scoring function
    (Model::get_model_scoring_function()), but it can be useful to use
    others in different contexts during a samping process.

    \note This may eventually become an abstract base class to support
    certain optimizations.
*/
class IMPEXPORT ScoringFunction:
#if defined(IMP_DOXYGEN) || defined(SWIG)
    public base::Object
#else
    public base::TrackedObject<ScoringFunction, Model>
#endif
 {
   typedef  base::TrackedObject<ScoringFunction, Model> Tracked;

  friend class Model;
  // kept alive in model
  ScoreStatesTemp ss_;
  double last_score_;
  bool last_was_good_;
  inline void ensure_dependencies();
  void update_score_states(const DependencyGraph &dg);
public:
  typedef std::pair<double, bool> ScoreIsGoodPair;
protected:
  virtual ScoreIsGoodPair do_evaluate_if_good(bool derivatives,
                                              const ScoreStatesTemp &ss)=0;
  virtual ScoreIsGoodPair do_evaluate(bool derivatives,
                                      const ScoreStatesTemp &ss)=0;
  virtual ScoreIsGoodPair do_evaluate_if_below(bool derivatives,
                                               double max,
                                               const ScoreStatesTemp &ss)=0;
  /** Return any extra score states that should be included in the list
      generated by the model for this ScoringFunction in addition to
      the ones implied by the restraints.

      By default it just calls create_restraints() and gets the score states
      needed for those restraints.*/
  virtual ScoreStatesTemp
    get_required_score_states(const DependencyGraph &dg) const;

 public:
  ScoringFunction(Model *m, std::string name);
  IMP_OBJECT_INLINE(ScoringFunction, out << create_restraints(),);
  inline Model *get_model() const;
  inline double evaluate_if_good(bool derivatives);
  inline double evaluate(bool derivatives);
  inline double evaluate_if_below(bool derivatives, double max);
  /** Return true if the last evaluate satisfied all the restraint
      thresholds.*/
  bool get_had_good_score() const {
    return last_was_good_;
  }
  double get_last_score() const {
    return last_score_;
  }
  /** Return a set of restraints equivalent to this scoring function.
   */
  virtual Restraints create_restraints() const=0;
  /** Return the score states needed to evaluate this ScoringFunction.*/
  const ScoreStatesTemp& get_score_states();
};


/** Return a list of ScoringFunction objects where each is as simple
    as possible and evaluating the sum (and anding the good score bits)
    is exactly like evaluating the one ScoringFunction.*/
IMPEXPORT ScoringFunctions create_decomposition(ScoringFunction *sf);


/** Return a list of ScoringFunction objects where each is as simple
    as possible and evaluating the sum (and anding the good score bits)
    is exactly like evaluating the one ScoringFunction.

    These scoring functions are optimized for incremental evaluation.
*/
IMPEXPORT ScoringFunctions
create_incremental_decomposition(ScoringFunction *sf);

/** This class is to provide a backward compatible interface for things
    that take ScoringFunctions as arguments, but used to take
    RestraintsTemp or a RestraintSet. */
class IMPEXPORT ScoringFunctionInput
#ifndef SWIG
    : public base::OwnerPointer<ScoringFunction>
#endif
{
  typedef base::OwnerPointer<ScoringFunction> P;
 public:
  ScoringFunctionInput(ScoringFunction *sf): P(sf){}
  ScoringFunctionInput(Model *sf);
  ScoringFunctionInput(const RestraintsTemp &sf);
  ScoringFunctionInput(RestraintSet *sf);
};

//! Print the hierarchy of restraints
IMPEXPORT void show_restraint_hierarchy(ScoringFunctionInput rs,
                                     std::ostream &out=std::cout);


IMP_END_NAMESPACE

#endif  /* IMPKERNEL_DECLARE_SCORING_FUNCTION_H */
