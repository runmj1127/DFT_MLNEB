# 3_run_mlneb.py (espresso.neb.in 파싱 기능이 포함된 최종 완성본)

import sys
import copy
import re
from ase.io import read
from ase.calculators.espresso import Espresso
from catlearn.optimize.mlneb import MLNEB

# --- 1. 사용자 설정 부분 ---

# 입력 파일: 2_run_optimization.sh 실행으로 생성된 .pwo 파일들
OPTIMIZED_INITIAL_PWO = 'initial_opt.pwo'
OPTIMIZED_FINAL_PWO = 'final_opt.pwo'

# 원본 파라미터 파일
INPUT_QE_FILE = 'espresso.neb.in'

PSEUDO_DIR = './'
N_CORES = 16

# ML-NEB 설정
N_IMAGES = 5
NEB_FMAX = 0.1
NEB_TRAJECTORY_FILE = 'mlneb_final.traj'

# --- 2. 안정적인 QE 파라미터 파서 ---

def parse_qe_parameters(filename):
    """
    espresso.neb.in 파일을 안정적으로 파싱하여 ASE와 호환되는 설정을 추출합니다.
    """
    settings = {
        'control': {}, 'system': {}, 'electrons': {},
        'pseudos': {}, 'kpts': (1, 1, 1)
    }
    
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()

    # Namelist 파싱
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
    
    # ATOMIC_SPECIES 파싱
    match = re.search(r'ATOMIC_SPECIES.*?((?:\n\s*\w+\s+[\d.]+\s+\S+)+)', content, re.IGNORECASE | re.DOTALL)
    if match:
        for line in match.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) >= 3:
                settings['pseudos'][parts[0]] = parts[2]
    
    return settings['pseudos'], settings['control'], settings['system'], settings['electrons'], settings['kpts']


# --- 3. 메인 워크플로우 ---

def run_main_mlneb():
    print("--- ML-NEB 계산 시작 ---")
    
    try:
        # 1. espresso.neb.in에서 계산 파라미터 파싱
        print(f">>> Parsing parameters from {INPUT_QE_FILE}...")
        pseudos, control, system, electrons, kpts = parse_qe_parameters(INPUT_QE_FILE)

        # 2. .pwo 파일에서 최적화된 초기/최종 구조 읽기
        print(f">>> Reading optimized structures from .pwo files...")
        initial_atoms = read(OPTIMIZED_INITIAL_PWO, format='espresso-out')
        final_atoms = read(OPTIMIZED_FINAL_PWO, format='espresso-out')

        # 3. NEB 계산을 위한 input_data 준비
        input_data = {'control': control, 'system': system, 'electrons': electrons}
        input_data['control']['calculation'] = "'neb'" # 계산 종류 변경
        if 'pseudo_dir' in input_data['control']:
            del input_data['control']['pseudo_dir']
        
        # 4. NEB 계산기 설정
        ase_calculator = Espresso(
            pseudopotentials=pseudos,
            input_data=input_data,
            kpts=kpts,
            pseudo_dir=PSEUDO_DIR,
            nprocs=N_CORES,
            executable='pw.x')
        ase_calculator.set_label('qe_calc_neb')

        # 5. ML-NEB 객체 생성 및 실행
        print("\n>>> Initializing MLNEB object...")
        mlneb = MLNEB(start=initial_atoms, end=final_atoms,
                      ase_calc=copy.deepcopy(ase_calculator), n_images=N_IMAGES, k=0.1)
        
        print("\n>>> Starting ML-NEB run...")
        mlneb.run(fmax=NEB_FMAX, trajectory=NEB_TRAJECTORY_FILE)
        print(f"\n>>> ML-NEB calculation completed! Trajectory saved to '{NEB_TRAJECTORY_FILE}'.")

    except FileNotFoundError:
        print(f"\n!!! 오류: 입력 파일(.pwo 또는 .neb.in)을 찾을 수 없습니다.")
        sys.exit(1)
    except Exception as e:
        print(f"\n!!! ML-NEB 계산 중 오류 발생: {e}")
        sys.exit(1)

if __name__ == "__main__":
    run_main_mlneb()
