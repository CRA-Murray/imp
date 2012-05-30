/**
 *  \file declare_QuadContainer.h    \brief A container for Quads.
 *
 *  This file is generated by a script (core/tools/make-container).
 *  Do not edit directly.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 */

#ifndef IMPKERNEL_DECLARE_QUAD_CONTAINER_H
#define IMPKERNEL_DECLARE_QUAD_CONTAINER_H

#include "kernel_config.h"
#include "internal/IndexingIterator.h"
#include "declare_Particle.h"
#include "container_base.h"
#include "internal/container_helpers.h"
#include "DerivativeAccumulator.h"
#include "internal/OwnerPointer.h"
#include "ParticleTuple.h"
#include <IMP/base/ref_counted_macros.h>
#include <IMP/base/check_macros.h>
#include <IMP/base/Pointer.h>
#include <IMP/base/InputAdaptor.h>


IMP_BEGIN_NAMESPACE
class QuadModifier;
class QuadDerivativeModifier;
class QuadScore;

//! A shared container for Quads
/** Stores a searchable shared collection of Quads.
    \ingroup restraints

    \implementationwithoutexample{QuadContainer, IMP_QUAD_CONTAINER}
 */
class IMPEXPORT QuadContainer : public Container
{
 protected:
  QuadContainer(Model *m,
                     std::string name="QuadContainer %1%");
public:
  typedef ParticleQuad ContainedType;
  typedef ParticleQuadsTemp ContainedTypes;
  typedef ParticleIndexQuads ContainedIndexTypes;
  typedef ParticleIndexQuad ContainedIndexType;
  /** \note This function may be linear. Be aware of the complexity
      bounds of your particular container.
   */
  virtual bool get_contains_particle_quad(const ParticleQuad& v) const =0;

  ParticleQuadsTemp get_particle_quads() const {
    return IMP::internal::get_particle(get_model(),
                                       get_indexes());
  }
#ifndef IMP_DOXGEN
  //! return the number of Quads in the container
  /** \note this isn't always constant time
   */
  virtual unsigned int get_number_of_particle_quads() const {
    return get_number();
  }

  virtual ParticleQuad get_particle_quad(unsigned int i) const {
    return get(i);
  }

#endif

  //! Apply a SingletonModifier to the contents
  virtual void apply(const QuadModifier *sm) const=0;
  //! Apply a SingletonModifier to the contents
  virtual void apply(const QuadDerivativeModifier *sm,
                     DerivativeAccumulator &da) const=0;

  //! Evaluate a score on the contents
  virtual double evaluate(const QuadScore *s,
                          DerivativeAccumulator *da) const=0;

  //! Evaluate a score on the contents
  virtual double evaluate_if_good(const QuadScore *s,
                                  DerivativeAccumulator *da,
                                  double max) const=0;

#ifndef IMP_DOXYGEN
  ParticleQuad get(unsigned int i) const {
    return IMP::internal::get_particle(get_model(),
                                       get_indexes()[i]);
  }
  ParticleQuadsTemp get() const {
    return IMP::internal::get_particle(get_model(),
                                       get_indexes());
  }
  bool get_contains(const ParticleQuad& v) const {
    return get_contains_particle_quad(v);
  }
  virtual bool get_contains_index(ParticleIndexQuad v) const {
    return get_contains_particle_quad(IMP::internal
                                     ::get_particle(get_model(),
                                                    v));
  }
  unsigned int get_number() const {return get_indexes().size();}
  virtual ParticleIndexQuads get_indexes() const=0;
  virtual ParticleIndexQuads get_all_possible_indexes() const=0;
#ifndef SWIG
  virtual bool get_provides_access() const {return false;}
  virtual const ParticleIndexQuads& get_access() const {
    IMP_THROW("Object not implemented properly.", base::IndexException);
  }
#endif
#endif

  IMP_REF_COUNTED_NONTRIVIAL_DESTRUCTOR(QuadContainer);
};


/** This class allows either a list or a container to be
    accepted as input.
*/
class IMPEXPORT QuadContainerAdaptor:
#if !defined(SWIG) && !defined(IMP_DOXYGEN)
public base::Pointer<QuadContainer>
#else
public base::InputAdaptor
#endif
{
  typedef base::Pointer<QuadContainer> P;
 public:
  QuadContainerAdaptor(){}
  QuadContainerAdaptor(QuadContainer *c);
  template <class C>
  QuadContainerAdaptor(base::internal::PointerBase<C> c): P(c){}
  QuadContainerAdaptor(const ParticleQuadsTemp &t,
                          std::string name="QuadContainerAdaptor%1%");
};


IMP_END_NAMESPACE

#endif  /* IMPKERNEL_DECLARE_QUAD_CONTAINER_H */
