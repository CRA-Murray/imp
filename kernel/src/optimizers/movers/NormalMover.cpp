/**
 *  \file NormalMover.cpp
 *  \brief A modifier which perturbs a point with a gaussian.
 *
 *  Copyright 2007-8 Sali Lab. All rights reserved.
 *
 */

#include <boost/random/normal_distribution.hpp>
#include "IMP/optimizers/movers/NormalMover.h"
#include "IMP/random.h"

namespace IMP
{

NormalMover::NormalMover(const Particles &pis,
                         const FloatKeys &vars,
                         Float max)
{
  IMP_assert(max != 0, "Must have some width");
  add_particles(pis);
  for (unsigned int i = 0; i < vars.size(); ++i) {
    add_key(vars[i]);
  }
  stddev_ = max;
}

void NormalMover::generate_move(float scale)
{
  std::vector<Float> center(number_of_float_keys());
  boost::uniform_01<RandomNumberGenerator> u01(random_number_generator);
  boost::normal_distribution<double> mrng(0, stddev_);
  for (unsigned int i = 0; i < number_of_particles(); ++i) {
    for (unsigned int j = 0; j < number_of_float_keys(); ++j) {
      float c = get_float(i, j);
      float r = mrng(u01);
      IMP_assert(!std::isnan(r), "Bad random");
      IMP_assert(!std::isnan(c), "Bad stored");
      propose_value(i, j, c + r);
    }
  }
}

} // namespace IMP
