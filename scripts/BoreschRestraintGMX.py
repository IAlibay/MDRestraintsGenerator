"""
Basic script for finding a Boresch restraint from an MD trajectory

Usage:
    python BoreschRestraintGMX.py --top file.gro --traj file.xtc --l1 2611 \
                                  --l2 2609 --l3 2607 \
                                  --host_selection "protein and name CA"
"""
import MDAnalysis as mda
from MDRestraintsGenerator.restraints import FindBoreschRestraint
import argparse


if __name__ == "__main__":

    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--top', default='file.tpr',
                            help=('path to input structure topology file '
                                  '(e.g. TPR, PARM7, etc...)'))
        parser.add_argument('--traj', default='file.xtc',
                            help=('path to input trajectory file'
                                  '(e.g. XTC, NC, TRJ')
        parser.add_argument('--l1', default=None,
                            help=('zero formated index of ligand bond '
                                  'forming atom'))
        parser.add_argument('--l2', default=None,
                            help=('zero formatted index of ligand angle '
                                  'forming atom')
        parser.add_argument('--l3', default=None,
                            help=('zero formatted index of the ligand '
                                  'dihedral forming atom'))
        parser.add_argument('--host_selection', default="protein and name CA",
                            help='host atom selection string')
        args = parser.parse_args()
        return args

    args = parse_args()

    ligand_atoms = [int(args.l1), int(args.l2), int(args.l3)]

    if None in ligand_atoms:
        errmsg = "Missing ligand atoms"
        raise IOError(errmsg)

    u = mda.Universe(args.top, args.traj)

    # Create the boresch finder analysis object
    find = FindBoreschRestraint(u, l_atoms=ligand_atoms,
                                p_selection=args.host_selection)

    # Run the restraint analysis
    find.run()
