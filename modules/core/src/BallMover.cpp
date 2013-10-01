/**
 *  \file BallMover.cpp  \brief A modifier which perturbs a discrete variable.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/core/BallMover.h>

#include <IMP/base/random.h>
#include <IMP/core/XYZ.h>
#include <IMP/algebra/vector_generators.h>
#include <boost/random/uniform_real.hpp>

IMPCORE_BEGIN_NAMESPACE

namespace {
std::string get_ball_mover_name(kernel::Model *m, kernel::ParticleIndex pi) {
  return "BallMover-" + m->get_particle(pi)->get_name();
}
}

void BallMover::initialize(kernel::ParticleIndexes pis, FloatKeys keys,
                           double radius) {
  pis_ = pis;
  keys_ = keys;
  set_radius(radius);
  originals_.resize(pis.size(), algebra::get_zero_vector_kd(keys.size()));
}

BallMover::BallMover(kernel::Model *m, kernel::ParticleIndex pi,
                     const FloatKeys &keys, double radius)
    : MonteCarloMover(m, get_ball_mover_name(m, pi)) {
  initialize(kernel::ParticleIndexes(1, pi), keys, radius);
}

BallMover::BallMover(kernel::Model *m, kernel::ParticleIndex pi, double radius)
    : MonteCarloMover(m, get_ball_mover_name(m, pi)) {
  initialize(kernel::ParticleIndexes(1, pi), XYZ::get_xyz_keys(), radius);
}

// backwards compat
BallMover::BallMover(const kernel::ParticlesTemp &sc, const FloatKeys &vars,
                     double max)
    : MonteCarloMover(sc[0]->get_model(), "BallMover%1%") {
  initialize(kernel::get_indexes(sc), vars, max);
}

// backwards compat
BallMover::BallMover(const kernel::ParticlesTemp &sc, double max)
    : MonteCarloMover(sc[0]->get_model(), "XYZBallMover%1%") {
  initialize(kernel::get_indexes(sc), XYZ::get_xyz_keys(), max);
}

MonteCarloMoverResult BallMover::do_propose() {
  IMP_OBJECT_LOG;
  algebra::SphereKD ball(algebra::get_zero_vector_kd(keys_.size()), radius_);
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    for (unsigned int j = 0; j < keys_.size(); ++j) {
      originals_[i][j] = get_model()->get_attribute(keys_[j], pis_[i]);
    }
    algebra::VectorKD nv =
        originals_[i] + IMP::algebra::get_random_vector_in(ball);
    for (unsigned int j = 0; j < keys_.size(); ++j) {
      IMP_USAGE_CHECK(
          get_model()->get_is_optimized(keys_[j], pis_[i]),
          "BallMover can't move non-optimized attribute. "
              << "particle: " << get_model()->get_particle_name(pis_[i])
              << "attribute: " << keys_[j]);
      get_model()->set_attribute(keys_[j], pis_[i], nv[j]);
    }
  }
  return MonteCarloMoverResult(pis_, 1.0);
}

void BallMover::do_reject() {
  IMP_OBJECT_LOG;
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    for (unsigned int j = 0; j < keys_.size(); ++j) {
      get_model()->set_attribute(keys_[j], pis_[i], originals_[i][j]);
    }
  }
}

kernel::ModelObjectsTemp BallMover::do_get_inputs() const {
  kernel::ModelObjectsTemp ret(pis_.size());
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    ret[i] = get_model()->get_particle(pis_[i]);
  }
  return ret;
}

IMPCORE_END_NAMESPACE
