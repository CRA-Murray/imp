/**
 *  \file BondPairContainer.cpp
 *  \brief A fake container that returns true if a pair of particles are bonded
 *
 *  Copyright 2007-9 Sali Lab. All rights reserved.
 *
 */

#include "IMP/atom/BondPairContainer.h"

IMPATOM_BEGIN_NAMESPACE

BondPairContainer
::BondPairContainer(SingletonContainer *sc): sc_(sc){
}

bool BondPairContainer
::get_contains_particle_pair(ParticlePair pp) const {
  if (!Bonded::is_instance_of(pp.first)
      || ! Bonded::is_instance_of(pp.second)) {
    return false;
  }

  Bonded ba(pp.first);
  Bonded bb(pp.second);
  Bond bd=get_bond(ba, bb);
  return sc_->get_contains_particle(bd);
}

unsigned int BondPairContainer
::get_number_of_particle_pairs() const {
  return sc_->get_number_of_particles();
}

ParticlePair BondPairContainer
::get_particle_pair(unsigned int i) const {
  Bond bd(sc_->get_particle(i));
  return ParticlePair(bd.get_bonded(0), bd.get_bonded(1));
}


void BondPairContainer::show(std::ostream &out) const {
  out << "BondPairContainer" << std::endl;
}

IMPATOM_END_NAMESPACE
