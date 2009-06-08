import unittest
import sys,os
import IMP
import IMP.em
import IMP.test
import IMP.core
import IMP.atom
import IMP.helper

class ProteinRigidFittingTest(IMP.test.TestCase):
    """Class to test EM correlation restraint"""

    def load_density_map(self):
        self.scene = IMP.em.DensityMap()
        mrw = IMP.em.MRCReaderWriter()
        self.scene.Read(self.get_input_file_name("1z5s_10.mrc"), mrw)
        self.scene.get_header_writable().set_resolution(10.)
        self.scene.get_header_writable().set_spacing(2.0)
        self.scene.set_origin(34.0,8.0,-92.0)

    def load_protein(self,pdb_filename):
        self.m = IMP.Model()
        self.mp= IMP.atom.read_pdb(self.open_input_file(pdb_filename),
                              self.m, IMP.atom.CAlphaSelector())#IMP.atom.NonWaterSelector())
        self.mps = IMP.atom.Hierarchys()
        self.mps.append(self.mp)
        self.radius_key = IMP.FloatKey("radius")
        self.weight_key = IMP.FloatKey("weight")
        #add radius and weight attributes
        self.particles = IMP.Particles()
        for p in IMP.core.get_leaves(self.mp):
            self.particles.append(p)
        for p in self.particles:
            p.add_attribute(self.radius_key, 1.5)
            p.add_attribute(self.weight_key, 1.0)

    def setUp(self):
        """Build test model and optimizer"""
        IMP.test.TestCase.setUp(self)
        IMP.set_log_level(IMP.SILENT)
        self.imp_model = IMP.Model()
        self.load_density_map()
        self.load_protein("1z5s_A.pdb")

    def test_local_fitting(self):
        """Check that the local rigid fitting functionality works"""
        #create a rigid body
        fr = IMP.em.FittingSolutions()
        print "start : "
        print fr.get_number_of_solutions()
        IMP.em.local_rigid_fitting_grid_search(
            self.particles,self.imp_model,
            self.radius_key, self.weight_key,
            self.scene,
            fr,1,1,30)

        #test that if you apply the transformation on the original configuration you get the same result
        # (in rmsd and score)

        #second, test that the optimization gets you close.
        score = self.imp_model.evaluate(False)
        print fr.get_number_of_solutions()
        self.assert_(fr.get_number_of_solutions() >= 0, "not enough solutions returned")
        self.assert_(fr.get_score(0) < fr.get_score(1), "solutions are not sorted")
        for i in xrange(fr.get_number_of_solutions()):
            self.assert_(fr.get_score(i) < 1.0, "wrong CC values")

if __name__ == '__main__':
    unittest.main()
