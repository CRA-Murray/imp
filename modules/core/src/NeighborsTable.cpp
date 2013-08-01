/**
 *  \file BallMover.cpp  \brief A modifier which perturbs a discrete variable.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/core/NeighborsTable.h>
#include <IMP/kernel/container_macros.h>

IMPCORE_BEGIN_NAMESPACE

void NeighborsTable::do_before_evaluate() {
  if (updated_ && !input_->get_is_changed()) return;
  updated_ = true;

  data_.clear();

  IMP_CONTAINER_FOREACH(PairContainer, input_, {
      data_[_1[0]].push_back(_1[1]);
      data_[_1[1]].push_back(_1[0]);
  });
}

NeighborsTable::NeighborsTable(kernel::PairContainer *input,
                               std::string name):
  ScoreState(input->get_model(), name),
  input_(input), updated_(false) {

}

IMPCORE_END_NAMESPACE
