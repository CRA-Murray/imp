/**
 *  \file Model.cpp \brief Storage of a model, its restraints,
 *                         constraints and particles.
 *
 *  Copyright 2007-2012 IMP Inventors. All rights reserved.
 *
 */

#include "IMP/dependency_graph.h"
#include "IMP/Model.h"
#include "IMP/RestraintSet.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/visitors.hpp>
#include <IMP/internal/graph_utility.h>
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/dynamic_bitset.hpp>
#include <IMP/base/file.h>
//#include <boost/graph/lookup_edge.hpp>
#include <IMP/compatibility/vector_property_map.h>
#include <boost/graph/reverse_graph.hpp>

IMP_BEGIN_NAMESPACE

template <class Graph, class Type, class Types>
class DirectCollectVisitor: public boost::default_dfs_visitor {
  typename boost::property_map<Graph,
                      boost::vertex_name_t>::const_type vm_;
  Types &vals_;
public:
  const Types &get_collected() {
    std::sort(vals_.begin(), vals_.end());
    vals_.erase(std::unique(vals_.begin(), vals_.end()), vals_.end());
    return vals_;
  }
  DirectCollectVisitor(const Graph &g, Types &vals): vals_(vals)
    {
      vm_=boost::get(boost::vertex_name, g);
    }
  template <class TG>
  void discover_vertex(typename boost::graph_traits<Graph>::vertex_descriptor u,
                       const TG&) {
    base::Object *o= vm_[u];
    //std::cout << "Visiting " << o->get_name() << std::endl;
    Type *p=dynamic_cast<Type*>(o);
    if (p) {
      //IMP_LOG(VERBOSE, "Found vertex " << o->get_name() << std::endl);
      vals_.push_back(p);
    } else {
      //IMP_LOG(VERBOSE, "Visited vertex " << o->get_name() << std::endl);
    }
  }
};

namespace {
  template <class ResultType, class Type, bool REVERSE>
  ResultType get_dependent(const base::ObjectsTemp &p,
                         const base::ObjectsTemp &all,
                         const DependencyGraph &dg,
                         const DependencyGraphVertexIndex &index) {
  IMP_FUNCTION_LOG;
  // find p in graph, ick
  DependencyGraphConstVertexName dpm= boost::get(boost::vertex_name, dg);
  std::pair<DependencyGraphTraits::vertex_iterator,
    DependencyGraphTraits::vertex_iterator> be
    = boost::vertices(dg);
  boost::vector_property_map<int> color(boost::num_vertices(dg));
  for (unsigned int i=0; i< all.size(); ++i) {
    IMP_USAGE_CHECK(index.find(all[i]) != index.end(),
                    "Blocker node not found in index");
    DependencyGraphVertex blocked=
      index.find(all[i])->second;
    IMP_INTERNAL_CHECK(color[blocked]==boost::color_traits<int>::white(),
                       "Vertex does not start white");
    color[blocked]= boost::color_traits<int>::black();
  }
  ResultType pt;
  DirectCollectVisitor<DependencyGraph, Type, ResultType> cv(dg, pt);
  for (unsigned int i=0; i< p.size(); ++i) {
    IMP_USAGE_CHECK(index.find(p[i]) != index.end(),
                    "Object " << p[i] << " not found in dependency graph");
    DependencyGraphVertex start= index.find(p[i])->second;
    if (REVERSE) {
      boost::depth_first_visit(boost::make_reverse_graph(dg), start, cv, color);
    } else {
      boost::depth_first_visit(dg, start, cv, color);
    }
  }
  return cv.get_collected();
}


}







ParticlesTemp get_dependent_particles(base::Object *p,
                                      const base::ObjectsTemp &all,
                                      const DependencyGraph &dg,
                       const DependencyGraphVertexIndex &index) {
  return get_dependent<ParticlesTemp, Particle, false>(base::ObjectsTemp(1,p),
                                                       all, dg,index);
}



RestraintsTemp get_dependent_restraints(base::Object *p,
                                        const base::ObjectsTemp &all,
                                        const DependencyGraph &dg,
                       const DependencyGraphVertexIndex &index) {
  return get_dependent<RestraintsTemp, Restraint, false>(base::ObjectsTemp(1,p),
                                                         all, dg, index);
}
ScoreStatesTemp get_dependent_score_states(base::Object *p,
                                           const base::ObjectsTemp &all,
                                           const DependencyGraph &dg,
                       const DependencyGraphVertexIndex &index) {
  return get_dependent<ScoreStatesTemp,
      ScoreState, false>(base::ObjectsTemp(1,p),all, dg, index);
}





ParticlesTemp get_required_particles(base::Object *p,
                                     const base::ObjectsTemp &all,
                                     const DependencyGraph &dg,
                       const DependencyGraphVertexIndex &index) {
  return get_dependent<ParticlesTemp, Particle, true>(base::ObjectsTemp(1,p),
                                                      all, dg, index);
}


ScoreStatesTemp get_required_score_states(base::Object *p,
                                          const base::ObjectsTemp &all,
                                          const DependencyGraph &dg,
                                 const DependencyGraphVertexIndex &index ) {
  return get_dependent<ScoreStatesTemp, ScoreState,
      true>(base::ObjectsTemp(1,p),all, dg,index);
}




namespace {


  template <class C>
  C filter(C c) {
    std::sort(c.begin(), c.end());
    c.erase(std::unique(c.begin(), c.end()), c.end());
    IMP_INTERNAL_CHECK(c.empty() || c[0],
                       "nullptr returned for dependencies.");
    return c;
  }

  bool get_has_edge(const DependencyGraph &graph,
                    DependencyGraphTraits::vertex_descriptor va,
                    DependencyGraphTraits::vertex_descriptor vb) {
    /*std::pair<MDependencyGraphTraits::out_edge_iterator,
      MDependencyGraphTraits::out_edge_iterator> edges
      = boost::out_edges(va, graph);
    for (; edges.first != edges.second;++edges.first) {
      if (boost::target(*edges.first, graph) == vb) return true;
    }
    return false;*/
    return boost::edge(va, vb, graph).second;
  }

  void add_edge(DependencyGraph &graph,
                DependencyGraphTraits::vertex_descriptor va,
                DependencyGraphTraits::vertex_descriptor vb) {
    if (get_has_edge(graph, va, vb)) return;
    IMP_INTERNAL_CHECK(va != vb, "Can't dependend on itself " << va);
    IMP_INTERNAL_CHECK(!get_has_edge(graph, va, vb),
                       "Already has edge between " << va << " and " << vb);
    boost::add_edge(va, vb, graph);
    IMP_INTERNAL_CHECK(get_has_edge(graph, va, vb),
                       "No has edge between " << va << " and " << vb);

  }

  DependencyGraphTraits::vertex_descriptor get_vertex(DependencyGraph &dg,
                                         DependencyGraphVertexIndex &dgi,
                                         base::Object *o) {
    DependencyGraphVertexIndex::const_iterator it=dgi.find(o);
    if (it==dgi.end()) {
      boost::property_map<DependencyGraph, boost::vertex_name_t>::type vm
        = boost::get(boost::vertex_name, dg);
      DependencyGraphTraits::vertex_descriptor v= boost::add_vertex(dg);
      vm[v]=o;
      dgi[o]=v;
      return v;
    } else {
      return it->second;
    }
  }

template <class It>
void add_out_edges(DependencyGraphTraits::vertex_descriptor rv,
                   It b, It e,
                   DependencyGraph &dg,
                   DependencyGraphVertexIndex &dgi) {
  for (It c=b; c!= e; ++c) {
    DependencyGraphTraits::vertex_descriptor cv= get_vertex(dg, dgi, *c);
    if (!get_has_edge(dg, rv, cv)) {
      add_edge(dg, cv, rv);
    }
  }
}

  template <class It>
  void build_inputs_graph(It b, It e,
                          DependencyGraph &dg,
                          DependencyGraphVertexIndex &dgi) {
    for (It c= b; c != e; ++c) {
      DependencyGraphTraits::vertex_descriptor rv= dgi.find(*c)->second;
      base::Object *o= *c;
      if (dynamic_cast<RestraintSet*>(o)) {
        RestraintSet *rs=dynamic_cast<RestraintSet*>(o);
        add_out_edges(rv, rs->restraints_begin(),
                      rs->restraints_end(), dg, dgi);
      } else {
        /*IMP_LOG(VERBOSE, "Processing inputs for \""
          << (*c)->get_name() << "\" ");*/
        {
          ContainersTemp ct= filter((*c)->get_input_containers());
          /*if (!ct.empty()) {
            IMP_LOG(VERBOSE, ", containers are "
            << ct);
            }*/
          add_out_edges(rv, ct.begin(), ct.end(), dg, dgi);
        }
        {
          ParticlesTemp pt= filter((*c)->get_input_particles());
          /*if (!pt.empty()) {
            IMP_LOG(VERBOSE, ", particles are " << pt);
            }*/
          add_out_edges(rv, pt.begin(), pt.end(), dg, dgi);
        }
      }
      //IMP_LOG(VERBOSE, std::endl);
    }
  }

  template <class It>
  void build_outputs_graph(It b, It e,
                           DependencyGraph &dg,
                           DependencyGraphVertexIndex &dgi) {
    for (It c= b; c != e; ++c) {
      DependencyGraphTraits::vertex_descriptor rv= dgi.find(*c)->second;
      /*IMP_LOG(VERBOSE, "Processing outputs for \""
        << (*c)->get_name()  << "\"");*/
      {
        ContainersTemp ct= filter((*c)->get_output_containers());
        /*IMP_IF_LOG(VERBOSE) {
          if (!ct.empty()) {
            IMP_LOG(VERBOSE, ", containers are "
                    << ct);
          }
          }*/
        for (unsigned int j=0; j < ct.size(); ++j) {
          DependencyGraphTraits::vertex_descriptor cv
            = get_vertex(dg, dgi, ct[j]);
          if (!get_has_edge(dg, cv, rv)) {
            add_edge(dg, rv, cv);
          }
        }
      }
      {
        ParticlesTemp pt= filter((*c)->get_output_particles());
        /*if (!pt.empty()) {
          IMP_LOG(VERBOSE, ", particles are "
                  << pt);
                  }*/
        for (unsigned int j=0; j < pt.size(); ++j) {
          DependencyGraphTraits::vertex_descriptor cv
            = get_vertex(dg, dgi, pt[j]);
          if (!get_has_edge(dg, cv, rv)) {
            add_edge(dg, rv, cv);
          }
        }
      }
      //IMP_LOG(VERBOSE, std::endl);
    }
  }
}
DependencyGraph
get_dependency_graph(Model *m) {
  ScoreStatesTemp ss(m->score_states_begin(),
                     m->score_states_end());
  RestraintsTemp rs=m->get_known_restraints();
  DependencyGraphVertexIndex index;
  DependencyGraph ret(ss.size()+rs.size());
  DependencyGraphVertexName
    vm = boost::get(boost::vertex_name, ret);
  for (unsigned int i=0; i< ss.size(); ++i) {
    vm[i]= ss[i];
    index[ss[i]]=i;
  }
  for (unsigned int i=0; i< rs.size(); ++i) {
    vm[i+ss.size()]= rs[i];
    index[rs[i]]=i+ss.size();
  }
  // Very important to do outputs first
  build_outputs_graph(ss.begin(), ss.end(), ret, index);
  build_inputs_graph(ss.begin(), ss.end(), ret, index);
  build_inputs_graph(rs.begin(), rs.end(), ret, index);
  base::Vector<std::pair<base::Object*, base::Object*> > extra;
  for (unsigned int i=0; i< extra.size(); ++i) {
    int va= index[extra[i].first];
    int vb= index[extra[i].second];
    boost::add_edge(va, vb, ret);
  }
  return ret;
}

namespace {
  template <class P>
  bool get_range_is_empty(const P &p) {
    return p.first==p.second;
  }
}

namespace {
  struct Connections {
    Ints in, out;
    Connections(int v, const DependencyGraph &g) {
      typedef boost::graph_traits<DependencyGraph> T;
      {
        typedef T::in_edge_iterator IIT;
        for (std::pair<IIT,IIT> be= boost::in_edges(v, g);
             be.first != be.second; ++be.first) {
          in.push_back(boost::source(*be.first, g));
        }
        std::sort(in.begin(), in.end());
      }
      {
        typedef T::out_edge_iterator OIT;
        for (std::pair<OIT,OIT> be= boost::out_edges(v, g);
             be.first != be.second; ++be.first) {
          out.push_back(boost::target(*be.first, g));
        }
        std::sort(out.begin(), out.end());
      }
      IMP_INTERNAL_CHECK(*this == *this, "Not equal");
    }
    int compare(const Connections &o) const {
      if (in.size() < o.in.size()) return -1;
      else if (in.size() > o.in.size()) return 1;
      else if (out.size() < o.out.size()) return -1;
      else if (out.size() > o.out.size()) return 1;
      else {
        for (unsigned int i=0; i< in.size(); ++i) {
          if (in[i] < o.in[i]) return -1;
          else if (in[i] > o.in[i]) return 1;
        }
        for (unsigned int i=0; i< out.size(); ++i) {
          if (out[i] < o.out[i]) return -1;
          else if (out[i] > o.out[i]) return 1;
        }
        return 0;
      }
    }
    IMP_HASHABLE_INLINE(Connections,
                        return boost::hash_range(in.begin(),
                                                 in.end())
                        + boost::hash_range(out.begin(),
                                            out.end()));
    IMP_COMPARISONS(Connections);
  };
  //IMP_VALUES(Connections, ConnectionsList);
  inline std::size_t hash_value(const Connections &t) {
    return t.__hash__();
  }
}




DependencyGraph
get_pruned_dependency_graph(Model *m) {
  IMP_FUNCTION_LOG;
  DependencyGraph full= get_dependency_graph(m);
  typedef boost::graph_traits<DependencyGraph> T;
  bool changed=true;
  while (changed) {
    changed=false;
    IMP_LOG(VERBOSE, "Searching for vertices to prune" << std::endl);
    compatibility::set<Connections> connections;
    for (unsigned int i=0; i< boost::num_vertices(full); ++i) {
      Connections c(i, full);
      if (connections.find(c) != connections.end()) {
        boost::property_map<DependencyGraph, boost::vertex_name_t>::type
          vm = boost::get(boost::vertex_name, full);
        IMP_LOG(VERBOSE, "Removing object " << vm[i]->get_name() << std::endl);
        for (unsigned int j=0; j< c.in.size(); ++j) {
          for (unsigned int k=0; k< c.out.size(); ++k) {
            //if (!boost::lookup_edge(c.in[j], c.out[k], full).second) {
              // why am I doing this anyway?
              //boost::add_edge(c.in[j], c.out[k], full);
            //}
          }
        }
        boost::clear_vertex(i, full);
        boost::remove_vertex(i, full);
        changed=true;
      } else {
        connections.insert(c);
      }
    }
  }
  return full;
}




struct cycle_detector : public boost::default_dfs_visitor {
  base::Vector<DependencyGraphVertex> cycle_;
  template <class DGEdge>
  void tree_edge(DGEdge e, const DependencyGraph&g) {
    DependencyGraphVertex t= boost::target(e, g);
    //MDGVertex s= boost::source(e, g);
    cycle_.push_back(t);
  }
  template <class DependencyGraphVertex>
  void finish_vertex(DependencyGraphVertex v, const DependencyGraph&) {
    IMP_USAGE_CHECK(cycle_.back()==v, "They don't match");
    cycle_.pop_back();
  }
  template <class ED>
  void back_edge(ED e, const DependencyGraph&g) {
    DependencyGraphVertex t= boost::target(e, g);
    //MDGVertex s= boost::source(e, g);
    base::Vector<DependencyGraphVertex>::iterator it
        = std::find(cycle_.begin(), cycle_.end(), t);
    //std::cout << s << " " << cycle_.back() << std::endl;
    if (it != cycle_.end()) {
      cycle_.erase(cycle_.begin(), it);
      cycle_.push_back(t);
      throw cycle_;
    } else {
      //std::cout << "non-loop " << s << " " << t << std::endl;
    }
  }
};

namespace {

  base::Vector<DependencyGraphVertex> get_cycle(const DependencyGraph &g) {
    cycle_detector vis;
    try {
      boost::vector_property_map<int> color(boost::num_vertices(g));
      boost::depth_first_search(g, boost::visitor(vis).color_map(color));
    } catch (base::Vector<DependencyGraphVertex> cycle) {
      //std::cerr << "Caught cycle " << cycle << std::endl;
      return cycle;
    }
    //std::cerr << "No cycle found" << std::endl;
    return base::Vector<DependencyGraphVertex>();
  }

#pragma GCC diagnostic ignored "-Wunused-parameter"
  void order_score_states(const DependencyGraph &dg,
                          ScoreStatesTemp &out) {
    base::Vector<DependencyGraphVertex> sorted;
    DependencyGraphConstVertexName om= boost::get(boost::vertex_name, dg);
    ScoreStatesTemp ret;
    try {
      boost::topological_sort(dg, std::back_inserter(sorted));
    } catch (...) {
      base::TextOutput out=base::create_temporary_file();
      base::internal::show_as_graphviz(dg, out);
      base::Vector<DependencyGraphVertex> cycle= get_cycle(dg);
      //std::ostringstream oss;
      std::cerr << "[";
      for (unsigned int i=0; i< cycle.size(); ++i) {
        std::cerr << om[cycle[i]]->get_name() << " -- ";
      }
      std::cerr << "]";
      IMP_THROW("Topological sort failed, probably due to loops in "
                << " dependency graph. See \"" << out.get_name() << "\"",
                base::ValueException);
    }
    for (int i=sorted.size()-1; i > -1; --i) {
      base::Object *o= om[sorted[i]];
      ScoreState *s=dynamic_cast<ScoreState*>(o);
      if (s) {
        out.push_back(s);
      }
    }
    //out= ScoreStatesTemp(out.rbegin(), out.rend());
  }
}

ScoreStatesTemp get_ordered_score_states(const DependencyGraph &dg) {
  IMP_FUNCTION_LOG;
  ScoreStatesTemp ret;
  order_score_states(dg, ret);
  return ret;
}

namespace {
struct LessOrder {
  bool operator()(const ScoreState*a, const ScoreState*b) const {
    return a->order_ < b->order_;
  }
};
}

IMPEXPORT
ScoreStatesTemp get_required_score_states(const RestraintsTemp &irs,
                                          const DependencyGraph &dg,
                             const DependencyGraphVertexIndex &index) {
  ScoreStatesTemp sst
      =  get_dependent<ScoreStatesTemp, ScoreState, true>(irs,
                                                      ObjectsTemp(), dg, index);
  return get_ordered_score_states(sst);
}

IMP_END_NAMESPACE
