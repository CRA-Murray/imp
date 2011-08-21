import IMP
import IMP.atom
import IMP.container

def create_representation():
    m= IMP.Model()
    all=IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
    all.set_name("the universe")
    def create_protein(name, ds):
        h=IMP.atom.create_protein(m, name, 10, ds)
        leaves= IMP.atom.get_leaves(h)
        all.add_child(h)
        r=IMP.atom.create_connectivity_restraint([IMP.atom.Selection(c)\
                                                  for c in h.get_children()],
                                                 1)
        if r:
            m.add_restraint(r)
    def create_protein_from_pdbs(name, files):
        def create_from_pdb(file):
            sls=IMP.SetLogState(IMP.NONE)
            t=IMP.atom.read_pdb( IMP.get_example_path("data/"+file), m,
                                 IMP.atom.ATOMPDBSelector())
            del sls
            #IMP.atom.show_molecular_hierarchy(t)
            c=IMP.atom.Chain(IMP.atom.get_by_type(t, IMP.atom.CHAIN_TYPE)[0])
            if c.get_number_of_children()==0:
                IMP.atom.show_molecular_hierarchy(t)
            # there is no reason to use all atoms, just approximate the pdb shape instead
            s=IMP.atom.create_simplified_along_backbone(c,
                                                        10.0/2.0)
            IMP.atom.destroy(t)
            # make the simplified structure rigid
            rb=IMP.atom.create_rigid_body(s)
            rb.set_coordinates_are_optimized(True)
            return s
        if len(files) >1:
            p= IMP.Particle(m)
            h= IMP.atom.Hierarchy.setup_particle(p)
            h.set_name(name)
            for i, f in enumerate(files):
                c=create_from_pdb(f)
                h.add_child(c)
                c.set_name(name+" chain "+str(i))
            r=IMP.atom.create_connectivity_restraint([IMP.atom.Selection(c)\
                                                      for c in h.get_children()],
                                                     1)
            if r:
                m.add_restraint(r)
        else:
            h= create_from_pdb(files[0])
            h.set_name(name)
        all.add_child(h)
    create_protein("Nup85", 570)
    create_protein("Nup84", 460)
    create_protein("Nup145C", 442)
    create_protein("Nup120", [0, 500, 761])
    create_protein("Nup133", [0, 450, 778, 1160])
    create_protein_from_pdbs("Seh1", ["seh1.pdb"])
    create_protein_from_pdbs("Sec13", ["sec13.pdb"])
    return (m, all)

def create_restraints(m, all):
    def add_connectivity_restraint(s):
        r= IMP.atom.create_connectivity_restraint(s, 1)
        m.add_restraint(r)
    def add_distance_restraint(s0, s1):
        r=IMP.atom.create_distance_restraint(s0,s1, 0, 1)
        m.add_restraint(r)
    evr=IMP.atom.create_excluded_volume_restraint([all])
    m.add_restraint(evr)
    S= IMP.atom.Selection
    s0=S(hierarchy=all, molecule="Nup145C", residue_indexes=[(0,423)])
    s1=S(hierarchy=all, molecule="Nup84")
    s2=S(hierarchy=all, molecule="Sec13")
    add_connectivity_restraint([s0,s1,s2])
    add_distance_restraint(S(hierarchy=all, molecule="Nup145C", residue_indexes=[(0,423)]),
                           S(hierarchy=all, molecule="Nup85"))
    add_distance_restraint(S(hierarchy=all, molecule="Nup145C", residue_indexes=[(0,423)]),
                           S(hierarchy=all, molecule="Nup120",
                             residue_indexes= [(500, 762)]))
    add_distance_restraint(S(hierarchy=all, molecule="Nup84"),
                           S(hierarchy=all, molecule="Nup133",
                             residue_indexes=[(778, 1160)]))
    add_distance_restraint(S(hierarchy=all, molecule="Nup85"),
                           S(hierarchy=all, molecule="Seh1"))
    add_distance_restraint(S(hierarchy=all, molecule="Nup145C",
                             residue_indexes=[(0,423)]),
                           S(hierarchy=all, molecule="Sec13"))

# now do the actual work
(m,all)= create_representation()
create_restraints(m, all)

# we can get the full dependency graph for the whole model with all the restraints
# but it is pretty complex
dg= IMP.get_dependency_graph([m.get_root_restraint_set()])
IMP.show_graphviz(dg);

# better thing to do is to get the "pruned" graph
pdg= IMP.get_pruned_dependency_graph([m.get_root_restraint_set()])
try:
    # these all open new windows which must be closed to continue
    # also, the graph is no where near as nice as displayed by
    # IMP.show_graphviz below

    #import matplotlib
    # the engine to be used must be selected before pyplot is imported
    #matplotlib.use("macosx")
    #import matplotlib.pyplot as plt

    # the method below requires the altgraph python package
    #xg=IMP.get_networkx_graph(pdg)

    #import networkx
    #networkx.draw_spectral(xg)
    #plt.show()
    pass
except:
    try:
        IMP.show_altgraph(pdg)
    except:
        print 'Need networkx and matplotlib or altgraph to display graphs interactively'

IMP.show_graphviz(pdg)
