import argparse
from mokit.lib.gaussian import uno, permute_orb
from mokit.lib.py2fch_direct import pchk2fch

descriptions = {
    'uno': 'uno: Generate UHF natural orbitals (UNOs) from a Gaussian .fch(k) file',
    'permute_orb': 'permute_orb: Swap/exchange two orbitals in a given .fch(k) file',
    'pchk2fch': 'pchk2fch: Convert a PySCF chkfile to a .fch file'
}

def uno_cli():
    parser = argparse.ArgumentParser(description=descriptions['uno'])
    parser.add_argument('filename', help='fch filename')
    args = parser.parse_args()

    uno(args.filename)

def permute_orb_cli():
    parser = argparse.ArgumentParser(description=descriptions['permute_orb'])
    parser.add_argument('filename', help='fch filename')
    parser.add_argument('orb1', type=int, help='index of the first orbital')
    parser.add_argument('orb2', type=int, help='index of the second orbital')
    args = parser.parse_args()

    permute_orb(args.filename, args.orb1, args.orb2)

def pchk2fch_cli():
    parser = argparse.ArgumentParser(description=descriptions['pchk2fch'])
    parser.add_argument('pchk', help='PySCF chkfile name')
    parser.add_argument('fchname', nargs='?', 
        help='.fch filename. If not provided, it will be derived from the pchk filename.')
    parser.add_argument('--density', action='store_true', help='save density matrix', default=False)
    args = parser.parse_args()

    pchk2fch(args.pchk, fchname=args.fchname, density=args.density)
