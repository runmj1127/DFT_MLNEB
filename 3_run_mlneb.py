# 3_run_mlneb.py (energy/forces 재계산 로직이 포함된 최종 완성본)

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
        
        # 1. .pwo 파일에서 최적화된 구조 '만' 읽어옵니다.
        print(f">>> Reading optimized structures from .pwo files...")
        initial_atoms = read(OPTIMIZED_INITIAL_PWO, format='espresso-out')
        final_atoms = read(OPTIMIZED_FINAL_PWO, format='espresso-out')

        # 2. 에너지와 힘을 '재계산'하기 위한 계산기 설정 (calculation = 'scf')
        scf_input_data = {'control': control, 'system': system, 'electrons': electrons}
        scf_input_data['control']['calculation'] = "'scf'" # 단일 에너지/힘 계산
        if 'pseudo_dir' in scf_input_data['control']: del scf_input_data['control']['pseudo_dir']
        
        scf_calculator = Espresso(
            pseudopotentials=pseudos, input_data=scf_input_data, kpts=kpts,
            pseudo_dir=PSEUDO_DIR, nprocs=N_CORES, executable='pw.x')
        
        # 3. 각 구조에 계산기를 붙이고 에너지/힘을 명시적으로 재계산하여 저장
        print(">>> Attaching calculator and re-calculating energy/forces for endpoints...")
        initial_atoms.calc = copy.deepcopy(scf_calculator)
        initial_atoms.calc.set_label('qe_calc_initial_scf')
        initial_energy = initial_atoms.get_potential_energy()
        initial_forces = initial_atoms.get_forces()
        print(f">>> Initial structure E={initial_energy:.4f} eV")

        final_atoms.calc = copy.deepcopy(scf_calculator)
        final_atoms.calc.set_label('qe_calc_final_scf')
        final_energy = final_atoms.get_potential_energy()
        final_forces = final_atoms.get_forces()
        print(f">>> Final structure E={final_energy:.4f} eV")
        
        # 4. NEB 계산을 위한 계산기 설정 (calculation = 'neb')
        neb_input_data = {'control': control, 'system': system, 'electrons': electrons}
        neb_input_data['control']['calculation'] = "'neb'"
        neb_calculator = Espresso(
            pseudopotentials=pseudos, input_data=neb_input_data, kpts=kpts,
            pseudo_dir=PSEUDO_DIR, nprocs=N_CORES, executable='pw.x')
        neb_calculator.set_label('qe_calc_neb')
        
        # 5. ML-NEB 실행
        print("\n>>> Initializing MLNEB object...")
        mlneb = MLNEB(start=initial_atoms, end=final_atoms,
                      ase_calc=copy.deepcopy(neb_calculator), n_images=N_IMAGES, k=0.1)
        
        print("\n>>> Starting ML-NEB run...")
        mlneb.run(fmax=NEB_FMAX, trajectory=NEB_TRAJECTORY_FILE)
        print(f"\n>>> ML-NEB calculation completed!")

    except Exception:
        print(f"\n!!! ML-NEB 계산 중 오류 발생:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    run_main_mlneb()
