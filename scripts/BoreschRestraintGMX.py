"""
Basic script for finding a Boresch restraint from an MD trajectory

Usage:
    python BoreschRestraintGMX.py --top file.gro --traj file.xtc --l1 2611 \
                                  --l2 2609 --l3 2607 \
                                  --host_selection "protein and name CA"
"""
import MDAnalysis as mda
from MDRestraintsGenerator import search
from MDRestraintsGenerator.restraints import FindBoreschRestraint
import argparse
import warnings


if __name__ == "__main__":

    def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument('--top', default='file.tpr',
                            help=('path to input structure topology file '
                                  '(e.g. TPR, PARM7, etc...)'))
        parser.add_argument('--traj', default='file.xtc',
                            help=('path to input trajectory file'
                                  '(e.g. XTC, NC, TRJ'))
        parser.add_argument('--l1', default=None,
                            help=('zero formated index of ligand bond '
                                  'forming atom'))
        parser.add_argument('--l2', default=None,
                            help=('zero formatted index of ligand angle '
                                  'forming atom'))
        parser.add_argument('--l3', default=None,
                            help=('zero formatted index of the ligand '
                                  'dihedral forming atom'))
        parser.add_argument('--ligand_selection', default="resname LIG",
                            help='ligand selection string')
        parser.add_argument('--host_selection', default="protein and name CA",
                            help='host atom selection string')
        parser.add_argument('--temperature', type=float, default=298.15,
                            help='simulation temperature')
        parser.add_argument('--force_constant', type=float, default=10.0,
                            help='restraint force constant')
        parser.add_argument('--outpath', default='./',
                            help='output path for writing files')
        args = parser.parse_args()
        return args

    args = parse_args()

    u = mda.Universe(args.top, args.traj)

    if None in [args.l1, args.l2, args.l3]:
        wmsg = "Ligand atoms not passed, will search for suitable atoms instead"
        warnings.warn(wmsg)
        # by default we will exclude H* named atoms
        l_sel = f"{args.ligand_selection} and not name H*"
        # We align based on the host selection
        ligand_atoms = search.find_ligand_atoms(u, l_selection=l_sel,
                                                p_align=args.host_selection)
    else:
        ligand_atoms = [[int(args.l1), int(args.l2), int(args.l3)]]

    # find protein atoms
    atom_set = []
    for l_atoms in ligand_atoms:
        psearch = search.FindHostAtoms(u, l_atoms[0],
                                       p_selection=args.host_selection)
        psearch.run()
        atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])

    # Create the boresch finder analysis object
    boresch = FindBoreschRestraint(u, atom_set)

    # Run the restraint analysis
    boresch.run()

    # Plot out the statistics
    boresch.restraint.plot(path=args.outpath)

    # Write out the intermolecular section to a topology
    boresch.restraint.write(path=args.outpath,
                            force_constant=args.force_constant)

    dG_off = boresch.restraint.standard_state(force_constant=args.force_constant,
                                              temperature=args.temperature)

    print(f"dG_off: {dG_off}, dG_on: {-dG_off}")
