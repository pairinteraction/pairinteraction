import unittest

import numpy as np

from pairinteraction import pireal as pi


class TestState(unittest.TestCase):
    def setUp(self):
        self.s = pi.StateOne("Sr3", 79, 1, 2, 0)
        self.s_artificial = pi.StateOne("G")

    def test_instantiation(self):
        self.assertEqual(type(self.s), pi.StateOne)

    def test_properties(self):
        self.assertEqual(self.s.getSpecies(), "Sr3")
        self.assertEqual(self.s.getElement(), "Sr")
        self.assertEqual(self.s.getS(), 1)
        self.assertEqual(self.s.getN(), 79)
        self.assertEqual(self.s.getL(), 1)
        self.assertEqual(self.s.getJ(), 2)
        self.assertEqual(self.s.getM(), 0)
        self.assertEqual(self.s_artificial.getLabel(), "G")
        self.assertEqual(self.s_artificial.isArtificial(), True)

    def test_comparison(self):
        self.assertTrue(self.s == pi.StateOne("Sr3", 79, 1, 2, 0))
        self.assertFalse(self.s != pi.StateOne("Sr3", 79, 1, 2, 0))
        self.assertTrue(self.s != pi.StateOne("Sr3", 79, 1, 2, 1))
        self.assertFalse(self.s == pi.StateOne("Sr3", 79, 1, 2, 1))
        self.assertTrue(self.s ^ pi.StateOne("Sr3", 79, 1, 2, 0))
        self.assertFalse(self.s ^ pi.StateOne("Sr3", 79, 2, 2, 0))
        self.assertTrue(self.s ^ pi.StateOne("Sr3", pi.ARB, 1, 2, 0))
        self.assertFalse(self.s ^ pi.StateOne("Sr3", pi.ARB, 2, 2, 0))
        self.assertFalse(self.s == self.s_artificial)

    def test_output(self):
        self.assertEqual(str(self.s), "|Sr3, 79 P_2, mj=0>")

    def test_reflection(self):
        self.assertEqual(pi.StateOne("Sr3", 79, 1, 2, 1).getReflected(), pi.StateOne("Sr3", 79, 1, 2, -1))


class TestPairState(unittest.TestCase):
    def setUp(self):
        self.s1 = pi.StateOne("Sr3", 80, 2, 1, 0)
        self.s2 = pi.StateOne("Sr3", 79, 1, 2, 0)
        self.s = pi.StateTwo(self.s1, self.s2)
        self.s_artificial = pi.StateTwo(["G", "G"])

    def test_combination(self):
        self.assertEqual(type(self.s), pi.StateTwo)

    def test_properties(self):
        self.assertEqual(self.s.getFirstState(), self.s1)
        self.assertEqual(self.s.getSecondState(), self.s2)
        self.assertEqual(self.s.getSpecies(), ["Sr3", "Sr3"])
        self.assertEqual(self.s.getElement(), ["Sr", "Sr"])
        np.testing.assert_equal(self.s.getS(), [1, 1])
        np.testing.assert_equal(self.s.getN(), [80, 79])
        np.testing.assert_equal(self.s.getL(), [2, 1])
        np.testing.assert_equal(self.s.getJ(), [1, 2])
        np.testing.assert_equal(self.s.getM(), [0, 0])
        self.assertEqual(self.s_artificial.getLabel(), ["G", "G"])
        self.assertEqual(self.s_artificial.isArtificial(), [True, True])

        for i in range(2):
            self.assertEqual(self.s.getSpecies(i), ["Sr3", "Sr3"][i])
            self.assertEqual(self.s.getElement(i), ["Sr", "Sr"][i])
            np.testing.assert_equal(self.s.getS(i), [1, 1][i])
            np.testing.assert_equal(self.s.getN(i), [80, 79][i])
            np.testing.assert_equal(self.s.getL(i), [2, 1][i])
            np.testing.assert_equal(self.s.getJ(i), [1, 2][i])
            np.testing.assert_equal(self.s.getM(i), [0, 0][i])
            self.assertEqual(self.s_artificial.getLabel(i), ["G", "G"][i])
            self.assertEqual(self.s_artificial.isArtificial(i), [True, True][i])

    def test_comparison(self):
        self.assertTrue(self.s == pi.StateTwo(["Sr3", "Sr3"], [80, 79], [2, 1], [1, 2], [0, 0]))
        self.assertFalse(self.s != pi.StateTwo(["Sr3", "Sr3"], [80, 79], [2, 1], [1, 2], [0, 0]))
        self.assertTrue(self.s != pi.StateTwo(["Sr3", "Sr3"], [80, 79], [2, 1], [1, 2], [0, 1]))
        self.assertFalse(self.s == pi.StateTwo(["Sr3", "Sr3"], [80, 79], [2, 1], [1, 2], [0, 1]))
        self.assertTrue(self.s ^ pi.StateTwo(["Sr3", "Sr3"], [80, 79], [2, 1], [1, 2], [0, 0]))
        self.assertFalse(self.s ^ pi.StateTwo(["Sr3", "Sr3"], [80, 59], [2, 1], [1, 2], [0, 0]))
        self.assertTrue(self.s ^ pi.StateTwo(["Sr3", "Sr3"], [pi.ARB, 79], [2, 1], [1, 2], [0, 0]))
        self.assertFalse(self.s ^ pi.StateTwo(["Sr3", "Sr3"], [pi.ARB, 59], [2, 1], [1, 2], [0, 0]))
        self.assertFalse(self.s == self.s_artificial)

    def test_output(self):
        self.assertEqual(str(self.s), "|Sr3, 80 D_1, mj=0>|Sr3, 79 P_2, mj=0>")

    def test_reflection(self):
        self.assertEqual(
            pi.StateTwo(["Sr3", "Sr3"], [80, 79], [2, 1], [1, 2], [-1, 2]).getReflected(),
            pi.StateTwo(["Sr3", "Sr3"], [80, 79], [2, 1], [1, 2], [1, -2]),
        )


if __name__ == "__main__":
    unittest.main()
