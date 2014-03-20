/**
 *  \file IMP/benchmark/flags.h
 *  \brief Various general useful macros for IMP.
 *
 *  Copyright 2007-2014 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBENCHMARK_FLAGS_H
#define IMPBENCHMARK_FLAGS_H

#include <IMP/benchmark/benchmark_config.h>
#include <IMP/base/Flag.h>
#include <boost/cstdint.hpp>

IMPBENCHMARK_BEGIN_INTERNAL_NAMESPACE

extern IMPBENCHMARKEXPORT base::AdvancedFlag<int> run_only;

#if IMP_BASE_HAS_GPERFTOOLS
extern IMPBENCHMARKEXPORT base::AdvancedFlag<bool> cpu_profile_benchmarks;
#endif

#if IMP_BASE_HAS_TCMALLOC_HEAPPROFILER
extern IMPBENCHMARKEXPORT base::AdvancedFlag<bool> heap_profile_benchmarks;
#endif

#if IMP_BASE_HAS_TCMALLOC_HEAPCHECKER
extern IMPBENCHMARKEXPORT base::AdvancedFlag<bool> heap_check_benchmarks;
#endif

IMPBENCHMARK_END_INTERNAL_NAMESPACE

#endif /* IMPBENCHMARK_FLAGS_H */
