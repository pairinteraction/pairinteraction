import numpy as np
import unittest

import @LIBNAME@.picomplex as pi

class TestState(unittest.TestCase):

    def test_instantiation(self):
        s = pi.StateOne("Sr3", 79, 1, 2, 0)
        self.assertEqual(type(s), pi.StateOne)

    def test_properties(self):
        s = pi.StateOne("Sr3", 79, 1, 2, 0)
        self.assertEqual(s.getSpecies(), "Sr3")
        self.assertEqual(s.species, "Sr3")
        self.assertEqual(s.element, "Sr")
        self.assertEqual(s.s, 1)
        self.assertEqual(s.getN(), 79)
        self.assertEqual(s.n, 79)
        self.assertEqual(s.getL(), 1)
        self.assertEqual(s.l, 1)
        self.assertEqual(s.getJ(), 2)
        self.assertEqual(s.j, 2)
        self.assertEqual(s.getM(), 0)
        self.assertEqual(s.m, 0)

class TestPairState(unittest.TestCase):

    def test_combination(self):
        s1 = pi.StateOne("Sr3", 80, 2, 1, 0)
        s2 = pi.StateOne("Sr3", 79, 1, 2, 0)
        s = pi.StateTwo(s1, s2)
        self.assertEqual(type(s), pi.StateTwo)

    def test_properties(self):
        s1 = pi.StateOne("Sr3", 80, 2, 1, 0)
        s2 = pi.StateOne("Sr3", 79, 1, 2, 0)
        s = pi.StateTwo(s1, s2)

        self.assertEqual(s.first(), s1)
        self.assertEqual(s.second(), s2)

        self.assertEqual(s.getSpecies(), ("Sr3", "Sr3"))
        self.assertEqual(s.species[0], "Sr3")
        self.assertEqual(s.species[1], "Sr3")
        self.assertEqual(s.element[0], "Sr")
        self.assertEqual(s.element[1], "Sr")
        self.assertEqual(s.s[0], 1)
        self.assertEqual(s.s[1], 1)
        np.testing.assert_equal(s.getN(), [80, 79])
        self.assertEqual(s.n[0], 80)
        self.assertEqual(s.n[1], 79)
        np.testing.assert_equal(s.getL(), [2, 1])
        self.assertEqual(s.l[0], 2)
        self.assertEqual(s.l[1], 1)
        np.testing.assert_equal(s.getJ(), [1, 2])
        self.assertEqual(s.j[0], 1)
        self.assertEqual(s.j[1], 2)
        np.testing.assert_equal(s.getM(), [0, 0])
        self.assertEqual(s.m[0], 0)
        self.assertEqual(s.m[1], 0)

if __name__ == '__main__':
    unittest.main()
