import argparse
from mokit.lib.gaussian import uno, permute_orb

descriptions = {
    'uno': 'uno: Generate UHF natural orbitals (UNOs) from a Gaussian .fch(k) file',
    'permute_orb': 'permute_orb: Swap/exchange two orbitals in a given .fch(k) file'
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
