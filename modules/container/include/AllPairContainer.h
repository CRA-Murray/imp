/**
 *  \file IMP/container/AllPairContainer.h
 *  \brief Return all pairs from a SingletonContainer
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2014 IMP Inventors. All rights reserved.
 */

#ifndef IMPCONTAINER_ALL_PAIR_CONTAINER_H
#define IMPCONTAINER_ALL_PAIR_CONTAINER_H

#include <IMP/container/container_config.h>
#include <IMP/generic.h>
#include <IMP/PairContainer.h>
#include <IMP/SingletonContainer.h>
#include <IMP/container/ListPairContainer.h>
#include <IMP/base/Pointer.h>
#include <IMP/pair_macros.h>
IMPCONTAINER_BEGIN_NAMESPACE

//! Return all unordered pairs of particles taken from the SingletonContainer
/*  See also AllBipartitePairContainer, ClosePairContainer,
    CloseBipartitePairContainer for variants on the functionality provided.
 */
class IMPCONTAINEREXPORT AllPairContainer : public PairContainer {
  IMP::base::PointerMember<SingletonContainer> c_;

 protected:
  virtual std::size_t do_get_contents_hash() const IMP_OVERRIDE {
    return c_->get_contents_hash();
  }

 public:
  template <class F>
  void apply_generic(F* f) const {
    validate_readable();
    const kernel::ParticleIndexes& pis = c_->get_contents();
    for (unsigned int i = 0; i < pis.size(); ++i) {
      for (unsigned int j = 0; j < i; ++j) {
        f->apply_index(get_model(), kernel::ParticleIndexPair(pis[i], pis[j]));
      }
    }
  }
  //! Get the individual particles from the passed SingletonContainer
  AllPairContainer(SingletonContainerAdaptor c,
                   std::string name = "AllPairContainer%1%");
  virtual kernel::ParticleIndexPairs get_indexes() const IMP_OVERRIDE;
  virtual kernel::ParticleIndexPairs get_range_indexes() const IMP_OVERRIDE;
  virtual kernel::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
  virtual kernel::ParticleIndexes get_all_possible_indexes() const IMP_OVERRIDE;
  IMP_PAIR_CONTAINER_METHODS(AllPairContainer);
  IMP_OBJECT_METHODS(AllPairContainer);
};

IMP_OBJECTS(AllPairContainer, AllPairContainers);

IMPCONTAINER_END_NAMESPACE

#endif /* IMPCONTAINER_ALL_PAIR_CONTAINER_H */
