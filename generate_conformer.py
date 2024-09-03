import argparse

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
from rdkit.Chem.rdForceFieldHelpers import (UFFOptimizeMoleculeConfs, MMFFOptimizeMoleculeConfs)

RDLogger.DisableLog("rdApp.*")


def get_rdkit_mol(sdf):
    mol = Chem.SDMolSupplier(sdf)[0]
    if mol is None:
        raise ValueError(f"Failed to read molecule from file: {sdf}")
    return mol


def generate_conformers(mol, num_confs, n_cpu):
    mol = Chem.AddHs(mol, addCoords=True)
    mol_conf = Chem.Mol(mol)
    EmbedMultipleConfs(mol_conf, numConfs=num_confs, numThreads=n_cpu)
    return mol, mol_conf


def optimize_conformers_UFF(mol_conf, max_iter=500, n_cpu=0):
    mol_conf_UFF = Chem.Mol(mol_conf)
    success = UFFOptimizeMoleculeConfs(mol_conf_UFF, maxIters=max_iter, numThreads=n_cpu)
    return (mol_conf_UFF, success)


def optimize_conformers_MMFF(mol_conf, max_iter=500, n_cpu=0):
    mol_conf_MMFF = Chem.Mol(mol_conf)
    success = MMFFOptimizeMoleculeConfs(mol_conf_MMFF, maxIters=max_iter, numThreads=n_cpu)
    return (mol_conf_MMFF, success)


def align_and_merge_conformers(mol, mol_conf, mol_UFF, mol_MMFF):
    combined_mol = mol
    mol_UFF, success_UFF = mol_UFF
    mol_MMFF, success_MMFF = mol_MMFF

    for conf_id in range(mol_conf.GetNumConformers()):
        conformer = mol_conf.GetConformer(conf_id)
        combined_mol.AddConformer(conformer, assignId=True)

    for conf_id, (success, energy) in zip( range(mol_UFF.GetNumConformers()), success_UFF ) :
        if success == 0:
            conformer = mol_UFF.GetConformer(conf_id)
            combined_mol.AddConformer(conformer, assignId=True)

    for conf_id, (success, energy) in zip( range(mol_MMFF.GetNumConformers()), success_MMFF ):
        if success == 0:
            conformer = mol_MMFF.GetConformer(conf_id)
            combined_mol.AddConformer(conformer, assignId=True)

    AllChem.AlignMolConformers(combined_mol)

    return combined_mol


def calculate_rmsds(mol):
    return [ AllChem.GetConformerRMS(mol, 0, i, prealigned=True) for i in range(mol.GetNumConformers()) ]


def write_conformers(mol, rmsds, out_file):
    writer = Chem.SDWriter(out_file)
    for i in range(mol.GetNumConformers()):
        conf_mol = Chem.Mol(mol)
        conf_mol.RemoveAllConformers()
        conf_mol.AddConformer(mol.GetConformer(i), assignId=True)
        conf_mol.SetProp("RMSD", str(rmsds[i]))
        writer.write(conf_mol)
    writer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate molecular conformations using the given molecule (SDF)")
    parser.add_argument("--sdf", required=True, type=str, help="Input molecule SDF file")
    parser.add_argument("--out", default=r"", type=str, help="Output SDF file name")
    parser.add_argument("--num", default=100, type=int, help="Number of conformations to generate")
    parser.add_argument("--max_iter", default=200, type=int, help="Maximum number of iterations for force field optimization")
    parser.add_argument("--n_cpu", default=0, type=int, help="Maximum number of CPU cores to use")
    parser.add_argument("--verbose", action="store_true", help="Verbosity (True or False)")
    args = parser.parse_args()

    mol = get_rdkit_mol(args.sdf)
    if args.out == "":
        args.out = args.sdf.rsplit(".", 1)[0] + "_conf.sdf"

    if args.verbose:
        smiles = Chem.MolToSmiles(mol)
        print(f"[INFO] Input SDF File: {args.sdf}")
        print(f"[INFO] Output SDF File: {args.out}")
        print(f"[INFO] Number of Conformers: {args.num}")
        print(f"[INFO] SMILES: {smiles}")

    mol, mol_conf = generate_conformers (mol, num_confs=args.num, n_cpu=args.n_cpu )
    mol_UFF = optimize_conformers_UFF( mol_conf, max_iter=args.max_iter, n_cpu=args.n_cpu )
    mol_MMFF = optimize_conformers_MMFF( mol_conf, max_iter=args.max_iter, n_cpu=args.n_cpu )

    combined_mol = align_and_merge_conformers( mol, mol_conf, mol_UFF, mol_MMFF )

    rmsds = calculate_rmsds( combined_mol )
    write_conformers( combined_mol, rmsds, args.out )

    if args.verbose:
        print(f"[INFO] Molecule Conformer Generation DONE: ", args.out)
