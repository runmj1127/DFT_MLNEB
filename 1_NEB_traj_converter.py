# 1_NEB_traj_converter.py (수동 파싱으로 안정성을 확보한 최종 버전)

from ase import Atoms
from ase.io import write
import numpy as np
import sys

def manual_parse_and_save(input_filename='espresso.neb.in'):
    """
    사용자 정의 NEB 입력 파일을 직접 파싱하여 첫 번째와 마지막 이미지를
    unoptimized .traj 파일로 저장합니다. ASE의 자동 파서에 의존하지 않습니다.
    """
    images_coords = {'FIRST_IMAGE': [], 'LAST_IMAGE': []}
    cell_params = {}
    current_image_key = None
    in_system_section = False
    is_reading_coords = False

    with open(input_filename, 'r', encoding='utf-8') as f:
        for line in f:
            stripped_line = line.strip()
            line_upper = stripped_line.upper()
            if '&SYSTEM' in line_upper: in_system_section = True; continue
            if in_system_section and '/' in line_upper: in_system_section = False; continue
            if 'FIRST_IMAGE' in line_upper: current_image_key = 'FIRST_IMAGE'; is_reading_coords = False; continue
            if 'LAST_IMAGE' in line_upper: current_image_key = 'LAST_IMAGE'; is_reading_coords = False; continue
            if 'INTERMEDIATE_IMAGE' in line_upper: current_image_key = None; is_reading_coords = False; continue
            if 'ATOMIC_POSITIONS' in line_upper and current_image_key: is_reading_coords = True; continue
            if in_system_section:
                parts = stripped_line.replace('=', ' ').replace(',', ' ').split()
                if len(parts) >= 2 and parts[0].lower() in ['a', 'b', 'c', 'cosab']:
                    try: cell_params[parts[0].lower()] = float(parts[1])
                    except ValueError: pass
            if is_reading_coords and current_image_key:
                parts = stripped_line.split()
                if len(parts) >= 4 and parts[0].isalpha():
                    images_coords[current_image_key].append({'symbol': parts[0], 'position': [float(p) for p in parts[1:4]]})

    a, b, c, cosab = cell_params.get('a'), cell_params.get('b'), cell_params.get('c'), cell_params.get('cosab')
    if not all([a, b, c, cosab is not None]):
        print(f"오류: 셀 파라미터(a, b, c, cosab)를 '{input_filename}'의 &SYSTEM 블록에서 찾을 수 없습니다."); sys.exit(1)
    sinab = (1.0 - cosab**2)**0.5
    cell = np.array([[a, 0.0, 0.0], [b * cosab, b * sinab, 0.0], [0.0, 0.0, c]])
    
    initial_data, final_data = images_coords['FIRST_IMAGE'], images_coords['LAST_IMAGE']
    if not initial_data or not final_data: print("오류: FIRST_IMAGE 또는 LAST_IMAGE의 좌표를 찾을 수 없습니다."); sys.exit(1)
    
    write('initial_unoptimized.traj', Atoms(symbols=[d['symbol'] for d in initial_data], positions=[d['position'] for d in initial_data], cell=cell, pbc=True))
    print("-> 'initial_unoptimized.traj' 파일이 생성되었습니다.")
    write('final_unoptimized.traj', Atoms(symbols=[d['symbol'] for d in final_data], positions=[d['position'] for d in final_data], cell=cell, pbc=True))
    print("-> 'final_unoptimized.traj' 파일이 생성되었습니다.")

if __name__ == "__main__":
    manual_parse_and_save()
