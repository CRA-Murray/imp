/**
 *  \file ConjugateGradients.cpp  \brief Simple conjugate gradients optimizer.
 *
 *  Copyright 2007-2010 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino2/DominoSampler.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/domino2/utility.h>
#include <IMP/domino2/internal/tree_inference.h>
#include <IMP/internal/graph_utility.h>
#include <IMP/file.h>
#include <boost/graph/connected_components.hpp>
#if BOOST_VERSION > 103900
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#include <boost/vector_property_map.hpp>
#endif


IMPDOMINO2_BEGIN_NAMESPACE


DominoSampler::DominoSampler(Model *m, ParticleStatesTable* pst,
                             std::string name):
  DiscreteSampler(m, pst, name){
}

DominoSampler::DominoSampler(Model *m, std::string name):
  DiscreteSampler(m, new ParticleStatesTable(), name){
}



namespace {
  SubsetStatesList
  get_solutions_from_tree(const SubsetGraph &jt,
                          const Subset &known_particles,
                          SubsetStatesTable *sst,
                          const SubsetFilterTables &sfts) {
    const SubsetStatesList pd
      = internal::get_best_conformations(jt, 0,
                                         known_particles,
                                         sfts, sst);
    return pd;
  }

  bool get_is_tree(const SubsetGraph &g) {
    // check connected components too
    if  (boost::num_edges(g)+1 != boost::num_vertices(g)) return false;
    else {
      boost::vector_property_map<int> comp(boost::num_vertices(g));
      int cc= boost::connected_components(g, comp);
      return cc==1;
    }
  }
}

SubsetStatesList DominoSampler
::do_get_sample_states(const Subset &known_particles) const {
  IMP_LOG(TERSE, "Sampling with " << known_particles.size()
          << " particles as " << known_particles << std::endl);
  Pointer<RestraintSet> rs= get_model()->get_root_restraint_set();
  OptimizeContainers co(rs, get_particle_states_table());
  OptimizeRestraints ro(rs, get_particle_states_table()->get_particles());
  ParticlesTemp pt(known_particles.begin(), known_particles.end());
  Pointer<SubsetGraphTable> sgt;
  if (sgt_) {
    sgt= sgt_;
  } else {
    sgt= new JunctionTreeTable(rs);
  }
  sgt->set_was_used(true);
  SubsetGraph jt=sgt->get_subset_graph(get_particle_states_table());
  IMP_IF_LOG(VERBOSE) {
    IMP_LOG(VERBOSE, "Subset graph is ");
    //std::ostringstream oss;
    IMP::internal::show_as_graphviz(jt, std::cout);
    //oss << std::endl;
    //IMP_LOG(TERSE, oss.str() << std::endl);
  }
  IMP::internal::OwnerPointer<SubsetEvaluatorTable> set
    = get_subset_evaluator_table_to_use();
  SubsetFilterTables sfts= get_subset_filter_tables_to_use(set);
  IMP::internal::OwnerPointer<SubsetStatesTable> sst
    = DiscreteSampler::get_subset_states_table_to_use(sfts);

  SubsetStatesList final_solutions;
  if (get_is_tree(jt)) {
    final_solutions
      = get_solutions_from_tree(jt, known_particles,
                                sst,
                                sfts);
  } else {
    // loopy case
    IMP_NOT_IMPLEMENTED;
    /* - For each subset, keep a list of possibly valid conformations
       - iterate for each node, n
          - for each neighbor, m
              - for each state in n, check if their is a compatible one
              in my (which passes the filters), discard if not.
       - when nothing changes on one run, stop
     */
  }
  return final_solutions;
}

void DominoSampler::set_subset_graph_table(SubsetGraphTable *s) {
  sgt_=s;
}

void DominoSampler::do_show(std::ostream &out) const {
}

IMPDOMINO2_END_NAMESPACE
