/**
 *  \file Sampler.cpp
 *  \brief Storage of a model, its restraints,
 *                         constraints and particles.
 *
 *  Copyright 2007-2014 IMP Inventors. All rights reserved.
 *
 */

#include "IMP/kernel/Sampler.h"
#include <limits>

IMPKERNEL_BEGIN_NAMESPACE

Sampler::Sampler(kernel::Model *m, std::string nm) : Object(nm), model_(m) {
  sf_ = m->create_model_scoring_function();
}

ConfigurationSet *Sampler::create_sample() const {
  IMP_OBJECT_LOG;
  set_was_used(true);
  // IMP_LOG_TERSE( "Sampling " << num_opt << " particles."<<std::endl);
  /*if (num_opt == 0) {
    IMP_WARN("There are no particles to optimize."<<std::endl);
    return nullptr;
    }*/
  return do_sample();
}

void Sampler::set_scoring_function(ScoringFunctionAdaptor sf) { sf_ = sf; }

ConfigurationSet *Sampler::get_sample() const {
  IMPKERNEL_DEPRECATED_METHOD_DEF(2.1, "Use create_sample instead");
  return create_sample();
}

Sampler::~Sampler() {}

IMPKERNEL_END_NAMESPACE
