import unittest
import sqlite3
from @LIBNAME@ import picomplex as pi

class TestQuantumDefect(unittest.TestCase):

    def test_comparison(self):
        qd = pi.QuantumDefect("Rb", 78, 1, 0.5)

        conn = sqlite3.connect("@LIBNAME@/databases/quantum_defects.db")
        stmt = conn.cursor()
        stmt.execute("select ac,Z,a1,a2,a3,a4,rc from model_potential"
                  + " where ( (element = 'Rb') and (L = 1) );")
        ac, Z, a1, a2, a3, a4, rc = stmt.fetchone()
        
        self.assertEqual(qd.ac, ac)
        self.assertEqual(qd.Z , Z )
        self.assertEqual(qd.a1, a1)
        self.assertEqual(qd.a2, a2)
        self.assertEqual(qd.a3, a3)
        self.assertEqual(qd.a4, a4)
        self.assertEqual(qd.rc, rc)

if __name__ == '__main__':
    unittest.main()
