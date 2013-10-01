/**
 *  \file exception.cpp   \brief Check handling.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include "IMP/base/statistics.h"
#include "IMP/base/internal/static.h"
#include <boost/format.hpp>
#include <boost/foreach.hpp>

IMPBASE_BEGIN_NAMESPACE
void clear_statistics() { internal::timings.clear(); }

void show_timings(TextOutput out) {
  out.get_stream() << (boost::format("%-60s%10s%8s") % "Operation," %
                       "seconds," % "calls,") << std::endl;
  typedef std::pair<std::string, internal::Timing> VT;
  BOOST_FOREACH(VT tp, internal::timings) {
    std::string name = tp.first;
    if (name.size() > 60) {
      name = std::string(name.begin(), name.begin() + 60);
    }
    out.get_stream() << (boost::format("%-61s%10f,%8d") % (name + ",") %
                         tp.second.total_time % tp.second.calls) << std::endl;
  }
}

IMPBASE_END_NAMESPACE
