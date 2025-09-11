
# 4_run_mlneb.py (Auto-generated)
import copy, sys, traceback
from ase.io import read
from ase.calculators.espresso import Espresso, EspressoProfile
from catlearn.optimize.mlneb import MLNEB

# Part 1: Global Settings
OPTIMIZED_INITIAL = 'initial.traj'
OPTIMIZED_FINAL = 'final.traj'
N_IMAGES = 5
NEB_FMAX = 0.1
TRAJECTORY_FILE = 'mlneb_final.traj'
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
    label='qe_calc_neb', pseudopotentials=pseudopotentials,
    input_data=qe_input_data, kpts=(1, 1, 1), profile=profile)

# Part 3: ML-NEB Run
def run_main_mlneb():
    print("\n>>> ML-NEB calculation starting.")
    mlneb = MLNEB(start=OPTIMIZED_INITIAL, end=OPTIMIZED_FINAL,
                  ase_calc=copy.deepcopy(ase_calculator), n_images=N_IMAGES, k=0.1)
    try:
        mlneb.run(fmax=NEB_FMAX, trajectory=TRAJECTORY_FILE)
        print("\n>>> ML-NEB calculation completed successfully!")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    run_main_mlneb()
