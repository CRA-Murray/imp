/**
 *  \file ListQuadContainer.cpp   \brief A list of ParticleQuadsTemp.
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2011 IMP Inventors. All rights reserved.
 *
 */

#include "IMP/core/internal/CoreListQuadContainer.h"
#include "IMP/QuadModifier.h"
#include "IMP/QuadScore.h"
#include <algorithm>


IMPCORE_BEGIN_INTERNAL_NAMESPACE


IMP_ACTIVE_CONTAINER_DEF(CoreListQuadContainer,);


CoreListQuadContainer
::CoreListQuadContainer(Model *m, std::string name):
  internal::ListLikeQuadContainer(m, name){
  initialize_active_container(m);
}


CoreListQuadContainer
::CoreListQuadContainer(Model *m, const char *name):
  internal::ListLikeQuadContainer(m, name){
  initialize_active_container(m);
}


void CoreListQuadContainer::do_show(std::ostream &out) const {
  IMP_CHECK_OBJECT(this);
  out << get_number_of_particle_quads()
      << " Quads." << std::endl;
}



void CoreListQuadContainer
::remove_particle_quads(const ParticleQuadsTemp &c) {
  if (c.empty()) return;
  ParticleIndexQuads cp= IMP::internal::get_index(c);
  remove_from_list(cp);
  IMP_IF_CHECK(USAGE) {
    for (unsigned int i=0; i< c.size(); ++i) {
      IMP_USAGE_CHECK(IMP::internal::is_valid(c[i]),
                    "Passed Quad cannot be nullptr (or None)");
    }
  }
}

ParticleIndexQuads
CoreListQuadContainer::get_all_possible_indexes() const {
    return get_indexes();
  }

void CoreListQuadContainer::do_before_evaluate() {
  internal::ListLikeQuadContainer::do_before_evaluate();
}

void CoreListQuadContainer::do_after_evaluate() {
  internal::ListLikeQuadContainer::do_after_evaluate();
}

ParticlesTemp CoreListQuadContainer::get_state_input_particles() const {
  return ParticlesTemp();
}

ContainersTemp CoreListQuadContainer::get_state_input_containers() const {
  return ContainersTemp();
}

IMPCORE_END_INTERNAL_NAMESPACE
