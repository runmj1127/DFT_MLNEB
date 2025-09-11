# 4_run_mlneb.py (최종 수정본)
import copy, sys, traceback
from ase.io import read
# EspressoProfile은 더 이상 사용하지 않으므로 import에서 제거합니다.
from ase.calculators.espresso import Espresso
from catlearn.optimize.mlneb import MLNEB

# Part 1: Global Settings
OPTIMIZED_INITIAL = 'initial.traj'
OPTIMIZED_FINAL = 'final.traj'
N_IMAGES = 5
NEB_FMAX = 0.1
TRAJECTORY_FILE = 'mlneb_final.traj'
N_CORES = 32 # 사용할 코어 수

# Part 2: Quantum Espresso Calculator (수정된 부분)
pseudopotentials = {
    'O': 'O.pbe-rrkjus.UPF',
    'Cu': 'Cu.pbe-n-van_ak.UPF',
    'N': 'N.pbe-rrkjus.UPF',
    'H': 'H.pbe-rrkjus.UPF',
}
input_data = {
    'control': {
        'verbosity': "'low'" # 문자열 값은 따옴표로 감싸야 합니다.
    },
    'system': {
        'degauss': 0.01,
        'ecutrho': 225.0,
        'ecutwfc': 25.0,
        'ibrav': 0,
        'nspin': 2,
        'occupations': "'smearing'",
        'smearing': "'gaussian'",
        'starting_magnetization(1)': 0.0,
        'starting_magnetization(2)': 0.2
    },
    'electrons': {
        'conv_thr': 1.0e-4,
        'diagonalization': "'david'",
        'electron_maxstep': 200,
        'mixing_beta': 0.1,
        'startingpot': "'atomic'",
        'startingwfc': "'atomic+random'"
    }
}
# EspressoProfile을 사용하지 않고, Espresso에 직접 command와 pseudo_dir를 전달합니다.
command = f'mpirun -np {N_CORES} pw.x -in PREFIX.pwi > PREFIX.pwo'
ase_calculator = Espresso(
    label='qe_calc_neb',
    command=command,
    pseudopotentials=pseudopotentials,
    pseudo_dir='./pseudo/', # 오타 수정 및 경로 완성
    input_data=input_data,
    kpts=(1, 1, 1))

# Part 3: ML-NEB Run (사용자께서 수정한 올바른 버전)
def run_main_mlneb():
    print("\\n>>> ML-NEB calculation starting.")
    
    # 1. 최적화된 초기/최종 구조 파일을 직접 읽어옵니다.
    initial_atoms = read(OPTIMIZED_INITIAL)
    final_images = read(OPTIMIZED_FINAL)
    
    # 2. 읽어온 구조에 계산기를 명시적으로 연결해줍니다.
    initial_atoms.calc = copy.deepcopy(ase_calculator)
    final_images.calc = copy.deepcopy(ase_calculator)

    # 3. 파일 경로 대신, 계산기가 연결된 Atoms 객체를 직접 전달합니다.
    mlneb = MLNEB(start=initial_atoms, end=final_images,
                  ase_calc=copy.deepcopy(ase_calculator), n_images=N_IMAGES, k=0.1)
    try:
        mlneb.run(fmax=NEB_FMAX, trajectory=TRAJECTORY_FILE)
        print("\\n>>> ML-NEB calculation completed successfully!")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    run_main_mlneb()
