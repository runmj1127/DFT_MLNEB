# 3_run_mlneb.py (에너지를 직접 읽어와서 문제를 해결하는 최종 버전)

import sys
import copy
import re
import traceback
from ase.io import read
from ase.calculators.espresso import Espresso
from ase.calculators.singlepoint import SinglePointCalculator
from catlearn.optimize.mlneb import MLNEB

# --- (사용자 설정 및 파서 부분은 이전과 동일) ---

# --- 3. 메인 워크플로우 (수정된 부분) ---
def run_main_mlneb():
    try:
        print("--- ML-NEB 계산 시작 ---")
        pseudos, control, system, electrons, kpts = parse_qe_parameters(INPUT_QE_FILE)
        
        # 1. .pwo 파일에서 구조를 읽어옵니다.
        print(f">>> Reading optimized structures from .pwo files...")
        initial_atoms = read(OPTIMIZED_INITIAL_PWO, format='espresso-out')
        final_atoms = read(OPTIMIZED_FINAL_PWO, format='espresso-out')

        # --- 핵심 수정 부분: 에너지와 힘을 직접 추출하여 새로 연결 ---
        print(">>> Re-attaching energy and forces to atoms...")
        initial_energy = initial_atoms.get_potential_energy()
        initial_forces = initial_atoms.get_forces()
        initial_atoms.calc = SinglePointCalculator(initial_atoms, energy=initial_energy, forces=initial_forces)

        final_energy = final_atoms.get_potential_energy()
        final_forces = final_atoms.get_forces()
        final_atoms.calc = SinglePointCalculator(final_atoms, energy=final_energy, forces=final_forces)
        # -----------------------------------------------------------

        # NEB 계산을 위한 input_data 준비
        input_data = {'control': control, 'system': system, 'electrons': electrons}
        input_data['control']['calculation'] = "'neb'"
        if 'pseudo_dir' in input_data['control']: del input_data['control']['pseudo_dir']
        
        # NEB 계산기 설정
        ase_calculator = Espresso(
            pseudopotentials=pseudos, input_data=input_data, kpts=kpts,
            pseudo_dir=PSEUDO_DIR, nprocs=N_CORES, executable='pw.x')
        ase_calculator.set_label('qe_calc_neb')

        print("\n>>> Initializing MLNEB object...")
        mlneb = MLNENEB(start=initial_atoms, end=final_atoms,
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
