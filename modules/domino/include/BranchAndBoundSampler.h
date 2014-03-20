/**
 *  \file IMP/domino/BranchAndBoundSampler.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2014 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPDOMINO_BRANCH_AND_BOUND_SAMPLER_H
#define IMPDOMINO_BRANCH_AND_BOUND_SAMPLER_H

#include <IMP/domino/domino_config.h>
//#include "Evaluator.h"
#include "DiscreteSampler.h"
#include <IMP/Sampler.h>
#include <IMP/macros.h>
#include <IMP/base/Pointer.h>

IMPDOMINO_BEGIN_NAMESPACE

//! Sample best solutions using BranchAndBound
/** Find all good configurations of the model using branch and bound.
    Searches are truncated when the score is worse than the the thresholds
    in the Sampler or when two particles with the same kernel::ParticlesState
    are assigned the same state.
 */
class IMPDOMINOEXPORT BranchAndBoundSampler : public DiscreteSampler {
 public:
  BranchAndBoundSampler(kernel::Model *m,
                        std::string name = "BranchAndBoundSampler %1%");
  BranchAndBoundSampler(kernel::Model *m, ParticleStatesTable *pst,
                        std::string name = "BranchAndBoundSampler %1%");
  Assignments do_get_sample_assignments(const IMP::domino::Subset &known) const
      IMP_OVERRIDE;
  IMP_OBJECT_METHODS(BranchAndBoundSampler);
};

IMP_OBJECTS(BranchAndBoundSampler, BranchAndBoundSamplers);

IMPDOMINO_END_NAMESPACE

#endif /* IMPDOMINO_BRANCH_AND_BOUND_SAMPLER_H */
