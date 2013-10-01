/**
 *  \file CloseBipartitePairContainer.cpp   \brief A list of
 *kernel::ParticlePairs.
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2013 IMP Inventors. Close rights reserved.
 *
 */

#include "IMP/core/internal/CoreCloseBipartitePairContainer.h"
#include <IMP/core/BoxSweepClosePairsFinder.h>
#include <IMP/core/GridClosePairsFinder.h>
#include <IMP/container/ListPairContainer.h>
#include <IMP/core/internal/close_pairs_helpers.h>
#include <IMP/base/warning_macros.h>
#include <IMP/PairModifier.h>
#include <algorithm>

IMPCORE_BEGIN_INTERNAL_NAMESPACE

IMP_LIST_IMPL(CoreCloseBipartitePairContainer, PairFilter, pair_filter,
              PairFilter *, PairFilters);

CoreCloseBipartitePairContainer::CoreCloseBipartitePairContainer(
    SingletonContainer *a, SingletonContainer *b, double distance, double slack,
    std::string name)
    : P(a->get_model(), name) {
  std::ostringstream oss;
  oss << "BCPC " << get_name() << " hierarchy " << this;
  ObjectKey key = ObjectKey(oss.str());
  initialize(a, b, base::get_invalid_index<kernel::ParticleIndexTag>(),
             base::get_invalid_index<kernel::ParticleIndexTag>(), distance,
             slack, key);
}

CoreCloseBipartitePairContainer::CoreCloseBipartitePairContainer(
    kernel::SingletonContainer *a, kernel::SingletonContainer *b,
    kernel::ParticleIndex cover_a, kernel::ParticleIndex cover_b,
    kernel::ObjectKey key, double distance, double slack, std::string name)
    : P(a->get_model(), name) {
  initialize(a, b, cover_a, cover_b, distance, slack, key);
}

void CoreCloseBipartitePairContainer::initialize(kernel::SingletonContainer *a,
                                                 kernel::SingletonContainer *b,
                                                 kernel::ParticleIndex cover_a,
                                                 kernel::ParticleIndex cover_b,
                                                 double distance, double slack,
                                                 kernel::ObjectKey key) {
  IMP_OBJECT_LOG;
  slack_ = slack;
  distance_ = distance;
  key_ = key;
  sc_[0] = a;
  sc_[1] = b;
  were_close_ = false;
  reset_ = false;
  covers_[0] = cover_a;
  covers_[1] = cover_b;
  for (unsigned int i = 0; i < 2; ++i) {
    internal::initialize_particles(sc_[i], key_, xyzrs_[i], rbs_[i],
                                   constituents_, rbs_backup_sphere_[i],
                                   rbs_backup_rot_[i], xyzrs_backup_[i]);
  }
}

ModelObjectsTemp CoreCloseBipartitePairContainer::do_get_inputs() const {
  kernel::ModelObjectsTemp ret;
  ret += internal::get_inputs(get_model(), sc_[0], access_pair_filters());
  ret += internal::get_inputs(get_model(), sc_[1], access_pair_filters());
  if (covers_[0] != base::get_invalid_index<kernel::ParticleIndexTag>()) {
    ret.push_back(get_model()->get_particle(covers_[0]));
    ret.push_back(get_model()->get_particle(covers_[1]));
  }
  return ret;
}

void CoreCloseBipartitePairContainer::do_before_evaluate() {
  IMP_OBJECT_LOG;
  IMP_IF_LOG(VERBOSE) {
    algebra::Sphere3Ds coords[2];
    for (unsigned int i = 0; i < 2; ++i) {
      for (unsigned int j = 0; j < xyzrs_[i].size(); ++j) {
        coords[i].push_back(get_model()->get_sphere(xyzrs_[i][j]));
      }
    }
    Floats distances;
    for (unsigned int i = 0; i < coords[0].size(); ++i) {
      for (unsigned int j = 0; j < coords[1].size(); ++j) {
        distances.push_back(algebra::get_distance(coords[0][i], coords[1][j]));
      }
    }
    IMP_LOG_VERBOSE(xyzrs_[0] << " " << coords[0] << " " << xyzrs_backup_[0]
                              << std::endl);
    IMP_LOG_VERBOSE(xyzrs_[1] << " " << coords[1] << " " << xyzrs_backup_[1]
                              << std::endl);
    IMP_LOG_VERBOSE(distances << std::endl);
  }
  if (covers_[0] == base::get_invalid_index<kernel::ParticleIndexTag>() ||
      algebra::get_distance(get_model()->get_sphere(covers_[0]),
                            get_model()->get_sphere(covers_[1])) < distance_ ||
      reset_) {
    if (!reset_ && were_close_ &&
        !internal::get_if_moved(get_model(), slack_, xyzrs_[0], rbs_[0],
                                constituents_, rbs_backup_sphere_[0],
                                rbs_backup_rot_[0], xyzrs_backup_[0]) &&
        !internal::get_if_moved(get_model(), slack_, xyzrs_[1], rbs_[1],
                                constituents_, rbs_backup_sphere_[1],
                                rbs_backup_rot_[1], xyzrs_backup_[1])) {
      IMP_LOG_TERSE("Nothing to update" << std::endl);
      // all ok
    } else {
      // rebuild
      IMP_LOG_TERSE("Recomputing bipartite close pairs list." << std::endl);
      internal::reset_moved(get_model(), xyzrs_[0], rbs_[0], constituents_,
                            rbs_backup_sphere_[0], rbs_backup_rot_[0],
                            xyzrs_backup_[0]);
      internal::reset_moved(get_model(), xyzrs_[1], rbs_[1], constituents_,
                            rbs_backup_sphere_[1], rbs_backup_rot_[1],
                            xyzrs_backup_[1]);
      ParticleIndexPairs pips;
      internal::fill_list(get_model(), access_pair_filters(), key_,
                          2 * slack_ + distance_, xyzrs_, rbs_, constituents_,
                          pips);
      reset_ = false;

      swap(pips);

      IMP_LOG_VERBOSE("List is " << get_access() << std::endl);
      IMP_IF_CHECK(base::USAGE_AND_INTERNAL) {
        kernel::ParticleIndexes sc0p = sc_[0]->get_indexes();
        kernel::ParticleIndexes sc1p = sc_[1]->get_indexes();
        kernel::ParticleIndexPairs unfound;
        for (unsigned int i = 0; i < sc0p.size(); ++i) {
          XYZR d0(get_model(), sc0p[i]);
          for (unsigned int j = 0; j < sc1p.size(); ++j) {
            XYZR d1(get_model(), sc1p[j]);
            double dist = get_distance(d0, d1);
            if (dist < .9 * (distance_ + 2 * slack_)) {
              ParticleIndexPair pip(sc0p[i], sc1p[j]);
              bool filtered = false;
              IMP_CHECK_VARIABLE(filtered);
              for (unsigned int i = 0; i < get_number_of_pair_filters(); ++i) {
                if (get_pair_filter(i)->get_value_index(get_model(), pip)) {
                  filtered = true;
                  break;
                }
              }
              if (!filtered &&
                  std::find(get_access().begin(), get_access().end(), pip) ==
                      get_access().end()) {
                unfound.push_back(pip);
              }
            }
          }
        }
        IMP_INTERNAL_CHECK(unfound.empty(),
                           "Missing particle pairs: " << unfound);
      }
    }
    were_close_ = true;
  } else {
    IMP_LOG_TERSE("Covers are well separated." << std::endl);
    kernel::ParticleIndexPairs none;
    swap(none);
  }
  IMP_LOG_VERBOSE("List is " << get_access() << std::endl);
  IMP_IF_CHECK(base::USAGE_AND_INTERNAL) {
    kernel::ParticleIndexes sc0p = sc_[0]->get_indexes();
    kernel::ParticleIndexes sc1p = sc_[1]->get_indexes();
    kernel::ParticleIndexPairs unfound;
    for (unsigned int i = 0; i < sc0p.size(); ++i) {
      XYZR d0(get_model(), sc0p[i]);
      for (unsigned int j = 0; j < sc1p.size(); ++j) {
        XYZR d1(get_model(), sc1p[j]);
        double dist = get_distance(d0, d1);
        if (dist < .9 * distance_) {
          kernel::ParticleIndexPair pip(sc0p[i], sc1p[j]);
          bool filtered = false;
          IMP_CHECK_VARIABLE(filtered);
          for (unsigned int i = 0; i < get_number_of_pair_filters(); ++i) {
            if (get_pair_filter(i)->get_value_index(get_model(), pip)) {
              filtered = true;
              break;
            }
          }
          if (!filtered && std::find(get_access().begin(), get_access().end(),
                                     pip) == get_access().end()) {
            unfound.push_back(pip);
          }
        }
      }
    }
    IMP_INTERNAL_CHECK(unfound.empty(), "Missing particle pairs: " << unfound);
  }
}

ParticleIndexes CoreCloseBipartitePairContainer::get_all_possible_indexes()
    const {
  kernel::ParticleIndexes ret = sc_[0]->get_indexes();
  ret += sc_[1]->get_indexes();
  return ret;
}

ParticleIndexPairs CoreCloseBipartitePairContainer::get_range_indexes() const {
  kernel::ParticleIndexes pis = sc_[0]->get_range_indexes();
  kernel::ParticleIndexes pjs = sc_[1]->get_range_indexes();
  kernel::ParticleIndexPairs ret;
  ret.reserve(pis.size() * pjs.size());
  for (unsigned int i = 0; i < pis.size(); ++i) {
    for (unsigned int j = 0; j < pjs.size(); ++j) {
      ret.push_back(kernel::ParticleIndexPair(pis[i], pjs[j]));
    }
  }
  return ret;
}

IMPCORE_END_INTERNAL_NAMESPACE
