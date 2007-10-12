/**
 *  \file Steepest_Descent.h   Simple steepest descent optimizer.
 *
 *  Copyright 2007 Sali Lab. All rights reserved.
 *
 */

#ifndef __IMP_STEEPEST_DESCENT_H
#define __IMP_STEEPEST_DESCENT_H

#include "../IMP_config.h"
#include "Optimizer.h"

namespace imp
{

/** A simple steepest descent optimizer */
class IMPDLLEXPORT Steepest_Descent : public Optimizer
{
public:
  Steepest_Descent();
  virtual ~Steepest_Descent();

  virtual Float optimize(Model& model, int max_steps, Float threshold=0.0);
  virtual std::string version(void) {
    return "0.5.0";
  }
  virtual std::string last_modified_by(void) {
    return "Bret Peterson";
  }

protected:
  Float step_size;
};

} // namespace imp

#endif  /* __IMP_STEEPEST_DESCENT_H */
