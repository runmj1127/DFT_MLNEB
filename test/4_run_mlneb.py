
# 4_run_mlneb.py (Auto-generated)
import copy, sys, traceback
from ase.io import read
from ase.calculators.espresso import Espresso
from catlearn.optimize.mlneb import MLNEB

# Part 1: Global Settings
OPTIMIZED_INITIAL = 'initial.traj'
OPTIMIZED_FINAL = 'final.traj'
N_IMAGES = 5
NEB_FMAX = 0.1
TRAJECTORY_FILE = 'mlneb_final.traj'
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
    label='qe_calc_neb',
    nprocs=N_CORES,
    executable='pw.x',
    pseudopotentials=pseudopotentials,
    pseudo_dir='./pseudo/',
    input_data=input_data,
    kpts=kpts)

# Part 3: ML-NEB Run
def run_main_mlneb():
    print("\n>>> ML-NEB calculation starting.")
    initial_atoms = read(OPTIMIZED_INITIAL)
    final_images = read(OPTIMIZED_FINAL)
    initial_atoms.calc = copy.deepcopy(ase_calculator)
    final_images.calc = copy.deepcopy(ase_calculator)
    mlneb = MLNEB(start=initial_atoms, end=final_images,
                  ase_calc=copy.deepcopy(ase_calculator), n_images=N_IMAGES, k=0.1)
    try:
        mlneb.run(fmax=NEB_FMAX, trajectory=TRAJECTORY_FILE)
        print("\n>>> ML-NEB calculation completed successfully!")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    run_main_mlneb()
