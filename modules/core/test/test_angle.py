import IMP
import IMP.test
import IMP.core
import io
import math


class Tests(IMP.test.TestCase):

    """Tests for angle restraints"""

    def _setup_particles(self, system_angle, scored_angle):
        """Build a system of 3 xyz particles that form (by construction) the
           angle `system_angle`, and add an angle restraint with mean
           `scored_angle`."""
        model = IMP.Model()
        particles = []
        particles.append(self.create_point_particle(model, -1.0, 0.0, 0.0))
        particles.append(self.create_point_particle(model, 0.0, 0.0, 0.0))
        particles.append(self.create_point_particle(model,
                                                    -math.cos(system_angle),
                                                    math.sin(system_angle),
                                                    0.0))
        k = IMP.core.Harmonic.get_k_from_standard_deviation(0.1)
        rs = [IMP.core.DistanceRestraint(IMP.core.Harmonic(1.0, k),
                                         particles[0], particles[1]),
              IMP.core.DistanceRestraint(IMP.core.Harmonic(1.0, k),
                                         particles[1], particles[2]),
              IMP.core.AngleRestraint(IMP.core.Harmonic(scored_angle, k),
                                      particles[0], particles[1], particles[2])]
        sf = IMP.core.RestraintsScoringFunction(rs)
        return model, rs[2], sf, particles

    def test_deriv(self):
        """Check derivatives of angle restraints"""
        angles = [0.25 * math.pi, 0.3 * math.pi, 0.6 * math.pi, 0.75 * math.pi]
        for system_angle in angles:
            for score_angle in angles:
                model, rsr, sf, ps = self._setup_particles(system_angle,
                                                           score_angle)
                self.assertXYZDerivativesInTolerance(
                    sf, IMP.core.XYZ(ps[0]),
                    0.3, 5.0)
                self.assertXYZDerivativesInTolerance(
                    sf, IMP.core.XYZ(ps[1]),
                    0.3, 5.0)
                self.assertXYZDerivativesInTolerance(
                    sf, IMP.core.XYZ(ps[2]),
                    0.3, 5.0)

    def test_score(self):
        """Check score of angle restraints"""
        angles = [0.25 * math.pi, 0.3 * math.pi, 0.6 * math.pi, 0.75 * math.pi]
        for i in range(len(angles)):
            # Score of model with the same angle as the scoring function's mean
            # should be zero:
            model, rsr, sf, ps = self._setup_particles(angles[i], angles[i])
            self.assertLess(sf.evaluate(False), 1e-6)
            # When the angle is different, score should be far from zero:
            model, rsr, sf, ps = self._setup_particles(angles[i],
                                                       angles[-i - 1])
            self.assertGreater(sf.evaluate(False), 10.0)
            # Optimizing should reduce the score to zero:
            opt = IMP.core.ConjugateGradients(model)
            opt.set_scoring_function(sf)
            self.assertLess(opt.optimize(50), 1e-6)
            self.assertLess(sf.evaluate(False), 1e-6)

    def test_show(self):
        """Check AngleRestraint::show() method"""
        model, rsr, sf, ps = self._setup_particles(math.pi / 2.0, math.pi / 2.0)
        s = io.BytesIO()
        rsr.show(s)
        # no reason to check the show value
        #self.assertEqual(s.getvalue().split('\n')[0], "angle restraint:")

if __name__ == '__main__':
    IMP.test.main()
