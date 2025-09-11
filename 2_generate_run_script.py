# 2_generate_run_scripts.py (NameError 수정된 최종 버전)

import sys

def parse_qe_input(filename='espresso.neb.in'):
    """
    espresso.neb.in 파일을 파싱하여 주요 계산 파라미터를 추출합니다.
    """
    settings = {
        'control': {}, 'system': {}, 'electrons': {},
        'pseudos': {}, 'kpoints_str': '(1, 1, 1)'
    }
    current_block = None
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                stripped_line = line.strip()
                if not stripped_line or stripped_line.startswith(('!', '#')): continue
                line_upper = stripped_line.upper()
                if line_upper.startswith(('&', 'ATOMIC_SPECIES', 'K_POINTS', 'BEGIN_POSITIONS')): current_block = None
                if line_upper.startswith('&'): current_block = line_upper.strip('&')
                elif 'ATOMIC_SPECIES' in line_upper: current_block = 'SPECIES'
                elif 'K_POINTS' in line_upper:
                    k_type = stripped_line.split()[-1]
                    if 'gamma' in k_type.lower(): settings['kpoints_str'] = '(1, 1, 1)'
                elif current_block and line_upper.startswith('/'): current_block = None
                if current_block in ['CONTROL', 'SYSTEM', 'ELECTRONS']:
                    line_no_comment = stripped_line.split('!')[0]
                    line_no_comma = line_no_comment.split(',')[0]
                    if '=' in line_no_comma:
                        key, value = [p.strip() for p in line_no_comma.split('=')]
                        settings[current_block.lower()][key.lower()] = value
                elif current_block == 'SPECIES':
                    parts = stripped_line.split()
                    if len(parts) >= 3 and parts[0].isalpha():
                        symbol, pseudo_file = parts[0], parts[2]
                        settings['pseudos'][symbol] = pseudo_file
    except FileNotFoundError:
        print(f"오류: 입력 파일 '{filename}'을 찾을 수 없습니다."); sys.exit(1)
    return settings

def create_run_scripts(settings, opt_filename='3_run_optimization.py', neb_filename='4_run_mlneb.py'):
    """
    파싱된 설정을 바탕으로 최적화 및 NEB 실행 스크립트 두 개를 생성합니다.
    """
    def format_dict_items(d, indent=8):
        prefix = ' ' * indent
        return ',\n'.join([f"{prefix}'{k}': {v}" for k, v in d.items()])
    def format_pseudos(d, indent=4):
        prefix = ' ' * indent
        return '\n'.join([f"{prefix}'{k}': '{v}'," for k, v in d.items()])

    pseudos_str = format_pseudos(settings['pseudos'])
    pseudo_dir_str = "'./"
    if 'pseudo_dir' in settings['control']: del settings['control']['pseudo_dir']

        # ibrav = 0 과 충돌하는 격자 상수 관련 키들을 제거합니다.
    lattice_keys_to_remove = ['a', 'b', 'c', 'cosab', 'cosac', 'cosbc']
    for key in lattice_keys_to_remove:
        if key in settings['system']:
            del settings['system'][key]
            
    # ASE 정책에 따라 ibrav=0으로 강제 설정
    settings['system']['ibrav'] = 0

    control_str = format_dict_items(settings['control'])
    system_str = format_dict_items(settings['system'])
    electrons_str = format_dict_items(settings['electrons'])
    kpts_str = settings['kpoints_str']

    # --- 템플릿 1: 최적화 스크립트 ---
    optimization_template = f"""
# {opt_filename} (Auto-generated)
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
pseudopotentials = {{
{pseudos_str}
}}
qe_input_data = {{
    'control': {{
{control_str}
    }},
    'system': {{
{system_str}
    }},
    'electrons': {{
{electrons_str}
    }}
}}
command = f'mpirun -np {{N_CORES}} pw.x'
profile = EspressoProfile(command=command, pseudo_dir={pseudo_dir_str})
ase_calculator = Espresso(
    label='qe_calc_opt', pseudopotentials=pseudopotentials,
    input_data=qe_input_data, kpts={kpts_str}, profile=profile)

# Part 3: Endpoint Optimization Logic
def optimize_endpoints():
    try:
        print("\\n>>> Optimizing Initial structure...")
        initial_atoms = read(UNOPTIMIZED_INITIAL)
        initial_atoms.set_calculator(copy.deepcopy(ase_calculator))
        optimizer_initial = BFGS(initial_atoms, trajectory=OPTIMIZED_INITIAL)
        optimizer_initial.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Initial structure optimization complete! Saved to '{{OPTIMIZED_INITIAL}}'.")

        print("\\n>>> Optimizing Final structure...")
        final_atoms = read(UNOPTIMIZED_FINAL)
        final_atoms.set_calculator(copy.deepcopy(ase_calculator))
        optimizer_final = BFGS(final_atoms, trajectory=OPTIMIZED_FINAL)
        optimizer_final.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Final structure optimization complete! Saved to '{{OPTIMIZED_FINAL}}'.")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    optimize_endpoints()
"""

    # --- 템플릿 2: ML-NEB 스크립트 ---
    mlneb_template = f"""
# {neb_filename} (Auto-generated)
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
pseudopotentials = {{
{pseudos_str}
}}
qe_input_data = {{
    'control': {{
{control_str}
    }},
    'system': {{
{system_str}
    }},
    'electrons': {{
{electrons_str}
    }}
}}
command = f'mpirun -np {{N_CORES}} pw.x'
profile = EspressoProfile(command=command, pseudo_dir={pseudo_dir_str})
ase_calculator = Espresso(
    label='qe_calc_neb', pseudopotentials=pseudopotentials,
    input_data=qe_input_data, kpts={kpts_str}, profile=profile)

# Part 3: ML-NEB Run
def run_main_mlneb():
    print("\\n>>> ML-NEB calculation starting.")
    mlneb = MLNEB(start=OPTIMIZED_INITIAL, end=OPTIMIZED_FINAL,
                  ase_calc=copy.deepcopy(ase_calculator), n_images=N_IMAGES, k=0.1)
    try:
        mlneb.run(fmax=NEB_FMAX, trajectory=TRAJECTORY_FILE)
        print("\\n>>> ML-NEB calculation completed successfully!")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    run_main_mlneb()
"""

    with open(opt_filename, 'w', encoding='utf-8') as f:
        f.write(optimization_template)
    print(f"'{opt_filename}' 파일이 생성되었습니다.")
    
    with open(neb_filename, 'w', encoding='utf-8') as f:
        f.write(mlneb_template)
    print(f"'{neb_filename}' 파일이 생성되었습니다.")

if __name__ == "__main__":
    parsed_settings = parse_qe_input()

    create_run_scripts(parsed_settings)
