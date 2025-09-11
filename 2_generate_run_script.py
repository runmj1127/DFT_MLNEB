# 2_generate_run_scripts.py (모든 오류 수정된 최종 버전)

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
            lines = f.readlines()
            in_species_block = False
            in_kpoints_block = False

            for line in lines:
                stripped_line = line.strip()
                if not stripped_line or stripped_line.startswith(('!', '#')):
                    continue

                line_upper = stripped_line.upper()

                if line_upper.startswith('&'):
                    current_block = line_upper.strip('&')
                    in_species_block = False
                    in_kpoints_block = False
                elif 'ATOMIC_SPECIES' in line_upper:
                    current_block = 'SPECIES'
                    in_species_block = True
                    in_kpoints_block = False
                elif 'K_POINTS' in line_upper:
                    current_block = 'K_POINTS'
                    in_kpoints_block = True
                    in_species_block = False
                elif line_upper.startswith('/'):
                    current_block = None
                    in_species_block = False
                    in_kpoints_block = False
                
                if current_block in ['CONTROL', 'SYSTEM', 'ELECTRONS']:
                    line_no_comment = stripped_line.split('!')[0]
                    if '=' in line_no_comment:
                        key, value = [p.strip() for p in line_no_comment.split('=')]
                        settings[current_block.lower()][key.lower().strip(',')] = value
                elif in_species_block and 'ATOMIC_SPECIES' not in line_upper:
                    parts = stripped_line.split()
                    if len(parts) >= 3:
                        symbol, pseudo_file = parts[0], parts[2]
                        settings['pseudos'][symbol] = pseudo_file
                elif in_kpoints_block and 'K_POINTS' not in line_upper:
                    k_values = stripped_line.split()[:3]
                    if all(k.isdigit() for k in k_values):
                         settings['kpoints_str'] = f"({', '.join(k_values)})"
                    in_kpoints_block = False

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

    def format_pseudos(d, indent=12):
        prefix = ' ' * indent
        return ',\n'.join([f"{prefix}'{k}': '{v}'" for k, v in d.items()])

    pseudos_str = format_pseudos(settings['pseudos'])
    
    # 입력 파일에서 pseudo_dir을 읽어오거나 기본값 설정
    pseudo_dir_str = settings['control'].get('pseudo_dir', "'./pseudo/'")
    if 'pseudo_dir' in settings['control']:
        del settings['control']['pseudo_dir']

    # ASE 정책에 따라 ibrav=0으로 강제 설정하고 관련 키 제거
    lattice_keys_to_remove = ['a', 'b', 'c', 'cosab', 'cosac', 'cosbc']
    for key in lattice_keys_to_remove:
        if key in settings['system']:
            del settings['system'][key]
    settings['system']['ibrav'] = 0

    control_str = format_dict_items(settings['control'])
    system_str = format_dict_items(settings['system'])
    electrons_str = format_dict_items(settings['electrons'])
    kpts_str = settings['kpoints_str']

    # --- 공통 설정 텍스트 생성 ---
    common_settings_template = f"""
# Part 2: Quantum Espresso Calculator
N_CORES = 32
command = f'mpirun -np {{N_CORES}} pw.x -in PREFIX.pwi > PREFIX.pwo'
pseudopotentials = {{
{pseudos_str}
}}
input_data = {{
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
"""

    # --- 템플릿 1: 최적화 스크립트 ---
    optimization_template = f"""
# {opt_filename} (Auto-generated)
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
{common_settings_template}
ase_calculator = Espresso(
    label='qe_calc_opt', command=command, pseudopotentials=pseudopotentials,
    pseudo_dir={pseudo_dir_str}, input_data=input_data, kpts={kpts_str})

# Part 3: Endpoint Optimization Logic
def optimize_endpoints():
    try:
        print("\\n>>> Optimizing Initial structure...")
        initial_atoms = read(UNOPTIMIZED_INITIAL)
        initial_atoms.calc = copy.deepcopy(ase_calculator)
        optimizer_initial = BFGS(initial_atoms, trajectory=OPTIMIZED_INITIAL)
        optimizer_initial.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Initial structure optimization complete! Saved to '{{OPTIMIZED_INITIAL}}'.")

        print("\\n>>> Optimizing Final structure...")
        final_atoms = read(UNOPTIMIZED_FINAL)
        final_atoms.calc = copy.deepcopy(ase_calculator)
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
from ase.calculators.espresso import Espresso
from catlearn.optimize.mlneb import MLNEB

# Part 1: Global Settings
OPTIMIZED_INITIAL = 'initial.traj'
OPTIMIZED_FINAL = 'final.traj'
N_IMAGES = 5
NEB_FMAX = 0.1
TRAJECTORY_FILE = 'mlneb_final.traj'
{common_settings_template}
ase_calculator = Espresso(
    label='qe_calc_neb', command=command, pseudopotentials=pseudopotentials,
    pseudo_dir={pseudo_dir_str}, input_data=input_data, kpts={kpts_str})

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
