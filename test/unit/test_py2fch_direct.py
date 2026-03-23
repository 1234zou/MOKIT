import unittest
import os

try:
    from pyscf import gto, scf, mcscf
    from mokit.lib.py2fch_direct import pchk2fch
    pyscf_available = True
except ImportError:
    pyscf_available = False

@unittest.skipUnless(pyscf_available, "pyscf not available")
class TestPchk2fchDirect(unittest.TestCase):
    def setUp(self):
        # Create a simple molecule for testing
        self.mol = gto.Mole()
        self.mol.atom = '''
        O 0.0 0.0 0.0
        H 0.757 0.586 0.0
        H -0.757 0.586 0.0
        '''
        self.mol.basis = 'sto-3g'
        self.mol.spin = 0
        self.mol.charge = 0
        self.mol.output = '/dev/null'
        self.mol.build()

    def test_rhf_pchk_to_fch(self):
        """Test pchk2fch with RHF checkpoint file."""
        # 1. Run RHF calculation and save checkpoint
        pchk_path = 'test_rhf.pchk'
        mf = scf.RHF(self.mol)
        mf.chkfile = pchk_path
        mf.kernel()

        # 2. Convert pchk to fch
        fch_path = 'test_rhf.fch'
        pchk2fch(pchk_path, fch_path)

        # 3. Check if fch file is created
        self.assertTrue(os.path.exists(fch_path))
        
        # Cleanup
        if os.path.exists(pchk_path):
            os.remove(pchk_path)
        if os.path.exists(fch_path):
            os.remove(fch_path)

    def test_uhf_pchk_to_fch(self):
        """Test pchk2fch with UHF checkpoint file."""
        # 1. Run UHF calculation and save checkpoint
        # Use a radical (CH3) to have an odd number of electrons (9)
        mol_uhf = gto.Mole()
        mol_uhf.atom = '''
        C 0.0 0.0 0.0
        H 1.0 0.0 0.0
        H 0.0 1.0 0.0
        H 0.0 0.0 1.0
        '''
        mol_uhf.basis = 'sto-3g'
        mol_uhf.spin = 1
        mol_uhf.charge = 0
        mol_uhf.output = '/dev/null'
        mol_uhf.build()
        
        pchk_path = 'test_uhf.pchk'
        mf = scf.UHF(mol_uhf)
        mf.chkfile = pchk_path
        mf.kernel()

        # 2. Convert pchk to fch
        fch_path = 'test_uhf.fch'
        pchk2fch(pchk_path, fch_path)

        # 3. Check if fch file is created
        self.assertTrue(os.path.exists(fch_path))
        
        # Cleanup
        if os.path.exists(pchk_path):
            os.remove(pchk_path)
        if os.path.exists(fch_path):
            os.remove(fch_path)

    def test_casscf_pchk_to_fch(self):
        """Test pchk2fch with CASSCF checkpoint file."""
        # 1. Run CASSCF calculation and save checkpoint
        self.mol.build()
        mf = scf.RHF(self.mol).run()
        mc = mcscf.CASSCF(mf, 2, 2) 
        mc.natorb = True 
        
        pchk_path = 'test_casscf.pchk'
        mc.chkfile = pchk_path
        mc.kernel()

        # 2. Convert pchk to fch
        fch_path = 'test_casscf.fch'
        pchk2fch(pchk_path, fch_path)

        # 3. Check if fch file is created
        self.assertTrue(os.path.exists(fch_path))
        
        # Cleanup
        if os.path.exists(pchk_path):
            os.remove(pchk_path)
        if os.path.exists(fch_path):
            os.remove(fch_path)

if __name__ == '__main__':
    unittest.main()
