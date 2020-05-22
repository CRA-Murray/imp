/**
 *  \file IMP/em/CoarseL2Norm.h
 *  \brief Perform coarse fitting between two density objects.
 *
 *  Copyright 2007-2020 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPEM_COARSE_L2_NORM_H
#define IMPEM_COARSE_L2_NORM_H

#include <IMP/em/em_config.h>
#include <IMP/em/KernelParameters.h>
#include "DensityMap.h"
#include "def.h"
#include <vector>
#include <IMP/core/rigid_bodies.h>
#include <IMP/Refiner.h>
#include <IMP/em/MRCReaderWriter.h>

IMPEM_BEGIN_NAMESPACE

class IMPEMEXPORT CoarseL2Norm : public IMP::Object {

 public:
 CoarseL2Norm() : Object("CoarseL2Norm%1%") {}

  static FloatPair logabssumexp(double x,  double y,
				double sx, double sy);

  static FloatPair logabssumprodexp(double x,  double y,
				    double wx, double wy);

  static double logsumexp(double x, double y);

  static std::vector<double> get_value(Particle *p,
				       const algebra::Vector3D &pt,
				       double mass_ii, 
				       IMP::em::KernelParameters kps_);

  static std::vector<double> get_value_no_deriv(Particle *p,
						const algebra::Vector3D &pt,
						double mass_ii,
						IMP::em::KernelParameters kps_);
  
  static std::pair<double, algebra::Vector3Ds> calc_score_and_derivative(DensityMap *em,
									 const IMP::ParticlesTemp &ps,
									 double resolution,
									 const algebra::Vector3Ds &dv);
  IMP_OBJECT_METHODS(CoarseL2Norm);

};
IMPEM_END_NAMESPACE
#endif 
