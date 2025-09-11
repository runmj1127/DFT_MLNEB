
# 3_run_optimization.py (Auto-generated)
import copy, sys, traceback
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso, EspressoProfile

# Part 1: Global Settings
UNOPTIMIZED_INITIAL = 'initial_unoptimized.traj'
UNOPTIMIZED_FINAL = 'final_unoptimized.traj'
OPTIMIZED_INITIAL = 'initial.traj'
OPTIMIZED_FINAL = 'final.traj'
OPTIMIZE_FMAX = 0.1
N_CORES = 12

# Part 2: Quantum Espresso Calculator
pseudopotentials = {
    'O': 'O.pbe-rrkjus.UPF',
    'Cu': 'Cu.pbe-n-van_ak.UPF',
    'N': 'N.pbe-rrkjus.UPF',
    'H': 'H.pbe-rrkjus.UPF',
}
qe_input_data = {
    'control': {
        'verbosity': "low"
    },
    'system': {
        'degauss': 0.01,
        'ecutrho': 225.0,
        'ecutwfc': 25.0,
        'ibrav': 0,
        'nat': 134,
        'nspin': 2,
        'ntyp': 4,
        'occupations': "smearing",
        'smearing': "gaussian",
        'starting_magnetization(1)': 0.0,
        'starting_magnetization(2)': 0.2
    },
    'electrons': {
        'conv_thr': 1.0e-4,
        'diagonalization': "david",
        'electron_maxstep': 200,
        'mixing_beta': 0.1,
        'startingpot': "atomic",
        'startingwfc': "atomic+random"
    }
}
command = f'mpirun -np {N_CORES} pw.x'
profile = EspressoProfile(command=command, pseudo_dir='./)
ase_calculator = Espresso(
    label='qe_calc_opt', pseudopotentials=pseudopotentials,
    input_data=qe_input_data, kpts=(1, 1, 1), profile=profile)

# Part 3: Endpoint Optimization Logic
def optimize_endpoints():
    try:
        print("\n>>> Optimizing Initial structure...")
        initial_atoms = read(UNOPTIMIZED_INITIAL)
        initial_atoms.set_calculator(copy.deepcopy(ase_calculator))
        optimizer_initial = BFGS(initial_atoms, trajectory=OPTIMIZED_INITIAL)
        optimizer_initial.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Initial structure optimization complete! Saved to '{OPTIMIZED_INITIAL}'.")

        print("\n>>> Optimizing Final structure...")
        final_atoms = read(UNOPTIMIZED_FINAL)
        final_atoms.set_calculator(copy.deepcopy(ase_calculator))
        optimizer_final = BFGS(final_atoms, trajectory=OPTIMIZED_FINAL)
        optimizer_final.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Final structure optimization complete! Saved to '{OPTIMIZED_FINAL}'.")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    optimize_endpoints()
