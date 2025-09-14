# 3_run_mlneb.py (energy 속성 오류를 수정한 최종 버전)

import sys
import copy
import re
import traceback
from ase.io import read
from ase.calculators.espresso import Espresso
from catlearn.optimize.mlneb import MLNEB

# --- 1. 사용자 설정 부분 ---
OPTIMIZED_INITIAL_PWO = 'initial_opt.pwo'
OPTIMIZED_FINAL_PWO = 'final_opt.pwo'
INPUT_QE_FILE = 'espresso.neb.in'
PSEUDO_DIR = './'
N_CORES = 16
N_IMAGES = 5
NEB_FMAX = 0.1
NEB_TRAJECTORY_FILE = 'mlneb_final.traj'

# --- 2. 안정적인 QE 파라미터 파서 ---
def parse_qe_parameters(filename):
    settings = {'control': {}, 'system': {}, 'electrons': {}, 'pseudos': {}, 'kpts': (1, 1, 1)}
    with open(filename, 'r', encoding='utf-8') as f: content = f.read()
    for section in ['CONTROL', 'SYSTEM', 'ELECTRONS']:
        match = re.search(f'&{section}(.*?)/', content, re.IGNORECASE | re.DOTALL)
        if match:
            for line in match.group(1).splitlines():
                line = line.split('!')[0].strip()
                if '=' in line:
                    parts = [p.strip().strip(',') for p in line.split('=')]
                    key, value = parts[0].lower(), parts[1]
                    try: float(value)
                    except ValueError:
                        if not (value.startswith("'") and value.endswith("'")): value = f"'{value}'"
                    settings[section.lower()][key] = value
    match = re.search(r'ATOMIC_SPECIES.*?((?:\n\s*\w+\s+[\d.]+\s+\S+)+)', content, re.IGNORECASE | re.DOTALL)
    if match:
        for line in match.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) >= 3: settings['pseudos'][parts[0]] = parts[2]
    match = re.search(r'K_POINTS\s+.*?[\r\n]+([\s\d\.]+)', content, re.IGNORECASE)
    if match:
        k_values = [int(v) for v in match.group(1).strip().split()[:3]]
        settings['kpts'] = tuple(k_values)
    return settings['pseudos'], settings['control'], settings['system'], settings['electrons'], settings['kpts']

# --- 3. 메인 워크플로우 ---
def run_main_mlneb():
    try:
        print("--- ML-NEB 계산 시작 ---")
        pseudos, control, system, electrons, kpts = parse_qe_parameters(INPUT_QE_FILE)
        
        initial_atoms = read(OPTIMIZED_INITIAL_PWO, format='espresso-out')
        final_atoms = read(OPTIMIZED_FINAL_PWO, format='espresso-out')

        input_data = {'control': control, 'system': system, 'electrons': electrons}
        input_data['control']['calculation'] = "'neb'"
        if 'pseudo_dir' in input_data['control']: del input_data['control']['pseudo_dir']
        
        ase_calculator = Espresso(
            pseudopotentials=pseudos, input_data=input_data, kpts=kpts,
            pseudo_dir=PSEUDO_DIR, nprocs=N_CORES, executable='pw.x')
        ase_calculator.set_label('qe_calc_neb')

        # --- 핵심 수정 부분 ---
        # .pwo에서 읽은 구조에 계산기를 다시 명시적으로 연결
        initial_atoms.calc = copy.deepcopy(ase_calculator)
        final_atoms.calc = copy.deepcopy(ase_calculator)

        print("\n>>> Initializing MLNEB object...")
        mlneb = MLNEB(start=initial_atoms, end=final_atoms,
                      ase_calc=copy.deepcopy(ase_calculator), n_images=N_IMAGES, k=0.1)
        
        print("\n>>> Starting ML-NEB run...")
        mlneb.run(fmax=NEB_FMAX, trajectory=NEB_TRAJECTORY_FILE)
        print(f"\n>>> ML-NEB calculation completed!")

    except Exception:
        print(f"\n!!! ML-NEB 계산 중 오류 발생:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    run_main_mlneb()
