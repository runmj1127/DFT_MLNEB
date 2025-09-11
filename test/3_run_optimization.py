
# 3_run_optimization.py (Auto-generated)
import copy, sys, traceback
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso

# Part 1: Global Settings
UNOPTIMIZED_INITIAL = 'initial_unoptimized.traj'
UNOPTIMIZED_FINAL = 'final_unoptimized.traj'
OPTIMIZED_INITIAL = 'initial.traj'
OPTIMIZED_FINAL = 'final.traj'
OPTIMIZE_FMAX = 0.1
N_CORES = 32

# Part 2: Quantum Espresso Calculator

pseudopotentials = {
            'H': 'H.pbe-rrkjus.UPF',
            'O': 'O.pbe-rrkjus.UPF'
}
input_data = {
    'control': {
        'calculation': 'relax'
    },
    'system': {
        'ibrav': 0,
        'nat': 3,
        'ntyp': 2,
        'ecutwfc': 30.0,
        'ecutrho': 240.0
    },
    'electrons': {
        'mixing_beta': 0.7
    }
}
kpts = (1, 1, 1)

ase_calculator = Espresso(
    label='qe_calc_opt',
    nprocs=N_CORES,
    executable='pw.x',
    pseudopotentials=pseudopotentials,
    pseudo_dir='./pseudo/',
    input_data=input_data,
    kpts=kpts)

# Part 3: Endpoint Optimization Logic
def optimize_endpoints():
    try:
        print("\n>>> Optimizing Initial structure...")
        initial_atoms = read(UNOPTIMIZED_INITIAL)
        initial_atoms.calc = copy.deepcopy(ase_calculator)
        optimizer_initial = BFGS(initial_atoms, trajectory=OPTIMIZED_INITIAL)
        optimizer_initial.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Initial structure optimization complete! Saved to '{OPTIMIZED_INITIAL}'.")

        print("\n>>> Optimizing Final structure...")
        final_atoms = read(UNOPTIMIZED_FINAL)
        final_atoms.calc = copy.deepcopy(ase_calculator)
        optimizer_final = BFGS(final_atoms, trajectory=OPTIMIZED_FINAL)
        optimizer_final.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Final structure optimization complete! Saved to '{OPTIMIZED_FINAL}'.")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    optimize_endpoints()
