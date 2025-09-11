# 2_generate_run_scripts.py (pseudo_dir 고정된 최종 버전)

import sys
import re

def parse_qe_input(filename='espresso.neb.in'):
    """
    espresso.neb.in 파일을 안정적으로 파싱하여 ASE와 호환되는 설정을 추출합니다.
    """
    settings = {
        'control': {}, 'system': {}, 'electrons': {},
        'pseudos': {}, 'kpoints_str': '(1, 1, 1)'
    }
    
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()

    # --- Namelist 파싱 ---
    for section in ['CONTROL', 'SYSTEM', 'ELECTRONS']:
        match = re.search(f'&{section}(.*?)/', content, re.IGNORECASE | re.DOTALL)
        if match:
            for line in match.group(1).splitlines():
                line = line.split('!')[0].strip()
                if '=' in line:
                    parts = [p.strip().strip(',') for p in line.split('=')]
                    key, value = parts[0].lower(), parts[1]
                    try:
                        float(value)
                    except ValueError:
                        if not (value.startswith("'") and value.endswith("'")):
                            value = f"'{value}'"
                    settings[section.lower()][key] = value

    # --- ATOMIC_SPECIES 파싱 ---
    match = re.search(r'ATOMIC_SPECIES.*?((?:\n\s*\w+\s+[\d.]+\s+\S+)+)', content, re.IGNORECASE | re.DOTALL)
    if match:
        for line in match.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) >= 3:
                settings['pseudos'][parts[0]] = parts[2]

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

    # --- ASE와 호환되도록 데이터 정리 ---
    # 1. pseudo_dir은 Espresso()에 직접 전달할 것이므로 input_data에서 제거
    if 'pseudo_dir' in settings['control']:
        del settings['control']['pseudo_dir']
    
    # 2. ASE 정책에 따라 ibrav=0으로 강제 설정하고 충돌 가능성 있는 키 제거
    lattice_keys_to_remove = ['a', 'b', 'c', 'cosab', 'cosac', 'cosbc']
    for key in lattice_keys_to_remove:
        if key in settings['system']:
            del settings['system'][key]
    settings['system']['ibrav'] = 0
    
    # 3. 호환성을 위해 smearing 방식 수정
    if 'smearing' in settings['system']:
        settings['system']['smearing'] = "'methfessel-paxton'"

    # --- 템플릿에 들어갈 문자열 준비 ---
    pseudos_str = format_pseudos(settings['pseudos'])
    control_str = format_dict_items(settings['control'])
    system_str = format_dict_items(settings['system'])
    electrons_str = format_dict_items(settings['electrons'])
    kpts_str = settings.get('kpoints_str', '(1, 1, 1)')


    # --- 공통 설정 템플릿 ---
    common_template = f"""
# Part 2: Quantum Espresso Calculator
N_CORES = 32
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
kpts = {kpts_str}

ase_calculator = Espresso(
    pseudopotentials=pseudopotentials,
    input_data=input_data,
    kpts=kpts,
    pseudo_dir='./',
    nprocs=N_CORES,
    executable='pw.x')
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
{common_template}
ase_calculator.set_label('qe_calc_opt')

# Part 3: Endpoint Optimization Logic
def optimize_endpoints():
    try:
        print("\\n>>> Optimizing Initial structure...")
        initial_atoms = read(UNOPTIMIZED_INITIAL)
        initial_atoms.calc = copy.deepcopy(ase_calculator)
        optimizer_initial = BFGS(initial_atoms, trajectory=OPTIMIZED_INITIAL)
        optimizer_initial.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Initial structure optimization complete!")

        print("\\n>>> Optimizing Final structure...")
        final_atoms = read(UNOPTIMIZED_FINAL)
        final_atoms.calc = copy.deepcopy(ase_calculator)
        optimizer_final = BFGS(final_atoms, trajectory=OPTIMIZED_FINAL)
        optimizer_final.run(fmax=OPTIMIZE_FMAX)
        print(f">>> Final structure optimization complete!")
    except Exception:
        traceback.print_exc()

if __name__ == "__main__":
    optimize_endpoints()
"""

    # --- 템플릿 2: ML-NEB 스크립트 ---
    mlneb_template = f"""
# {neb_filename} (Auto-generated)
# (필요 시 여기에 ML-NEB 스크립트 내용을 추가)
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

