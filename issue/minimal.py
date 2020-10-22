

import IMP
import IMP.saxs
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.saxs
import sys
pdb = 'pnkp.pdb'
saxs_data = 'smooth.dat'


model = IMP.Model()
system = IMP.pmi.topology.System(model)
state = system.create_state()
dof_obj = IMP.pmi.dof.DegreesOfFreedom(model)


pnkp = state.create_molecule("PNKP", "LGSRGRLWLQSPTGGPPPIFLPSDGQALVLGRGPLTQVTDRKCSRNQVELIADPESRTVAVKQLGVNPSTVGVHELKPGLSGSLSLGDVLYLVNGLYPLTLRWEELSTSGSQPDAPPDTPGDPEEGEDTEPQKKRVRKSSLGWESLKKLLVFTASGVKPQGKVAAFDLDGTLITTRSGKVFPTSPSDWRILYPEIPKKLQELAAEGYKLVIFTNQMGIGRGKLPAEVFKGKVEAVLEKLGVPFQVLVATHAGLNRKPVSGMWDHLQEQANEGIPISVEDSVFVGDAAGRLANWAPGRKKKDFSCADRLFALNVGLPFATPEEFFLKWPAARFELPAFDPRTISSAGPLYLPESSSLLSPNPEVVVAVGFPGAGKSTFIQEHLVSAGYVHVNRDTLGSWQRCVSSCQAALRQGKRVVIDNTNPDVPSRARYIQCAKDAGVPCRCFNFCATIEQARHNNRFREMTDPSHAPVSDMVMFSYRKQFEPPTLAEGFLEILEIPFRLQEHLDPALQRLYRQFSEG", "D")

pnkp.add_structure(pdb, 'D', soft_check=True)

pnkp.add_representation(resolutions=0, color='medium purple')

#build molecules?
root_hierarchy = system.build()

dof_obj.create_rigid_body(pnkp)

###################### RESTRAINTS #####################
output_objects = [] # keep a list of functions that need to be reported

# SAXS restraint
sr = IMP.pmi.restraints.saxs.SAXSRestraint(root_hierarchy, saxs_data, maxq=0.3)
output_objects.append(sr)
IMP.pmi.tools.add_restraint_to_model(model, sr.rs)

dof_obj.optimize_flexible_beads(100)

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(model,
                                    root_hier=root_hierarchy,                          # pass the root hierarchy
                                    #rmf_output_objects=crs,                     # will display like XLs
                                    monte_carlo_sample_objects=dof_obj.get_movers(),  # pass MC movers
                                    global_output_directory='output/',
                                    output_objects=output_objects,
                                    num_sample_rounds=1,
                                    monte_carlo_steps=50,
                                    number_of_best_scoring_models=0,      # set >0 to store best PDB files (but this is slow to do online)
                                    number_of_frames=31,
                                    nframes_write_coordinates=2)                   # increase number of frames to get better results!
rex.execute_macro()


