# 1_NEB_pwi_converter.py (FIRST/LAST 이미지만 정확히 파싱하는 최종 버전)

import sys
import numpy as np
from ase import Atoms
from ase.io import write

def manual_parse_and_save(input_filename='espresso.neb.in'):
    """
    사용자 정의 NEB 입력 파일을 직접 파싱하여 첫 번째와 마지막 이미지만
    unoptimized .traj 파일로 저장합니다.
    """
    
    images_coords = {'FIRST_IMAGE': [], 'LAST_IMAGE': []}
    cell_params = {}
    
    current_image_key = None
    in_system_section = False
    is_reading_coords = False

    try:
        with open(input_filename, 'r', encoding='utf-8') as f:
            for line in f:
                stripped_line = line.strip()
                line_upper = stripped_line.upper()

                # --- 섹션 탐색 ---
                if '&SYSTEM' in line_upper: in_system_section = True; continue
                if in_system_section and '/' in line_upper: in_system_section = False; continue
                
                if 'FIRST_IMAGE' in line_upper:
                    current_image_key = 'FIRST_IMAGE'
                    is_reading_coords = False
                    continue
                if 'LAST_IMAGE' in line_upper:
                    current_image_key = 'LAST_IMAGE'
                    is_reading_coords = False
                    continue
                if 'INTERMEDIATE_IMAGE' in line_upper:
                    current_image_key = None # 중간 이미지는 무시
                    is_reading_coords = False
                    continue
                
                if 'ATOMIC_POSITIONS' in line_upper and current_image_key:
                    is_reading_coords = True
                    continue
                
                # --- 데이터 파싱 ---
                if in_system_section:
                    parts = stripped_line.replace('=', ' ').replace(',', ' ').split()
                    if len(parts) >= 2 and parts[0].lower() in ['a', 'b', 'c', 'cosab']:
                        try: cell_params[parts[0].lower()] = float(parts[1])
                        except ValueError: pass
                
                if is_reading_coords and current_image_key:
                    parts = stripped_line.split()
                    if len(parts) >= 4 and parts[0].isalpha():
                        coords = {'symbol': parts[0], 'position': [float(p) for p in parts[1:4]]}
                        images_coords[current_image_key].append(coords)

        # --- 셀(cell) 행렬 생성 ---
        a = cell_params.get('a')
        b = cell_params.get('b')
        c = cell_params.get('c')
        cosab = cell_params.get('cosab')

        if not all([a, b, c, cosab is not None]):
            print(f"오류: 셀 파라미터(a, b, c, cosab)를 '{input_filename}'의 &SYSTEM 블록에서 찾을 수 없습니다."); sys.exit(1)
            
        sinab = (1.0 - cosab**2)**0.5
        cell = np.array([
            [a, 0.0, 0.0],
            [b * cosab, b * sinab, 0.0],
            [0.0, 0.0, c]
        ])
        
        # --- ASE Atoms 객체 생성 및 파일 저장 ---
        initial_data = images_coords['FIRST_IMAGE']
        final_data = images_coords['LAST_IMAGE']

        if not initial_data or not final_data:
            print("오류: FIRST_IMAGE 또는 LAST_IMAGE의 좌표를 찾을 수 없습니다."); sys.exit(1)

        initial_atoms = Atoms(symbols=[d['symbol'] for d in initial_data], positions=[d['position'] for d in initial_data], cell=cell, pbc=True)
        write('initial_unoptimized.traj', initial_atoms)
        print("-> 'initial_unoptimized.traj' 파일이 생성되었습니다.")
        
        final_atoms = Atoms(symbols=[d['symbol'] for d in final_data], positions=[d['position'] for d in final_data], cell=cell, pbc=True)
        write('final_unoptimized.traj', final_atoms)
        print("-> 'final_unoptimized.traj' 파일이 생성되었습니다.")

    except FileNotFoundError:
        print(f"오류: 입력 파일 '{input_filename}'을 찾을 수 없습니다.")
        sys.exit(1)
    except Exception as e:
        print(f"파일 처리 중 예기치 않은 오류가 발생했습니다: {e}")
        sys.exit(1)

if __name__ == "__main__":
    manual_parse_and_save()
