import unittest

import pairinteraction.picomplex as pi
au2GHz = 6579683.920711

class TestState(unittest.TestCase):

    def test_instantiation(self):
        s = pi.StateOne("Rb", 79, 1, 1.5, 0.5)
        self.assertEqual(type(s), pi.StateOne)

    def test_properties(self):
        s = pi.StateOne("Rb", 79, 1, 1.5, 0.5)
        self.assertEqual(s.element, "Rb")
        self.assertEqual(s.n, 79)
        self.assertEqual(s.l, 1)
        self.assertEqual(s.j, 1.5)
        self.assertEqual(s.m, 0.5)

class TestPairState(unittest.TestCase):

    def test_combination(self):
        s1 = pi.StateOne("Rb", 80, 0, 0.5, 0.5)
        s2 = pi.StateOne("Rb", 79, 1, 1.5, 0.5)
        s = pi.StateTwo(s1, s2)
        self.assertEqual(type(s), pi.StateTwo)

    def test_properties(self):
        s1 = pi.StateOne("Rb", 80, 0, 0.5, 0.5)
        s2 = pi.StateOne("Rb", 79, 1, 1.5, 0.5)
        s = pi.StateTwo(s1, s2)

        self.assertEqual(s.first(), s1)
        self.assertEqual(s.second(), s2)
        
        self.assertEqual(s.element[0], "Rb")
        self.assertEqual(s.element[1], "Rb")
        self.assertEqual(s.n[0], 80)
        self.assertEqual(s.n[1], 79)
        self.assertEqual(s.l[0], 0)
        self.assertEqual(s.l[1], 1)
        self.assertEqual(s.j[0], 0.5)
        self.assertEqual(s.j[1], 1.5)
        self.assertEqual(s.m[0], 0.5)
        self.assertEqual(s.m[1], 0.5)
        
if __name__ == '__main__':
    unittest.main()
