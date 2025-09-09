# 1_extract_structures.py (완전한 최종 버전)

from ase import Atoms
from ase.io import write
import numpy as np
import sys

def parse_and_save_neb(input_filename='espresso.neb.in'):
    """
    사용자 정의 NEB 입력 파일을 파싱하여 첫 번째와 마지막 이미지를
    unoptimized .traj 파일로 저장합니다.
    """
    
    all_images_data = []
    current_image_atoms = []
    cell_params = {}

    in_positions_section = False
    is_reading_coords = False
    in_system_section = False

    try:
        with open(input_filename, 'r', encoding='utf-8') as f:
            for line in f:
                stripped_line = line.strip()

                if not stripped_line: continue

                # 계산 설정(&SYSTEM 블록) 읽기
                if '&SYSTEM' in stripped_line.upper(): in_system_section = True; continue
                if in_system_section and '/' in stripped_line: in_system_section = False; continue
                if in_system_section:
                    parts = stripped_line.replace('=', ' ').replace(',', ' ').split()
                    if len(parts) >= 2:
                        param_name = parts[0].lower()
                        if param_name in ['a', 'b', 'c', 'cosab']:
                            try:
                                cell_params[param_name] = float(parts[1])
                            except ValueError:
                                print(f"경고: {param_name}의 값을 숫자로 변환할 수 없습니다: {parts[1]}")
                    continue

                # 구조 정보 읽기
                if 'BEGIN_POSITIONS' in stripped_line: in_positions_section = True; continue
                if 'END_POSITIONS' in stripped_line:
                    if current_image_atoms: all_images_data.append(current_image_atoms)
                    break
                if not in_positions_section: continue

                if 'FIRST_IMAGE' in stripped_line or 'INTERMEDIATE_IMAGE' in stripped_line or 'LAST_IMAGE' in stripped_line:
                    if current_image_atoms: all_images_data.append(current_image_atoms)
                    current_image_atoms = []
                    is_reading_coords = False
                    continue

                if 'ATOMIC_POSITIONS' in stripped_line:
                    is_reading_coords = True
                    continue

                if is_reading_coords:
                    parts = stripped_line.split()
                    if len(parts) >= 4 and parts[0].isalpha():
                        symbol = parts[0]
                        position = [float(p) for p in parts[1:4]]
                        current_image_atoms.append({'symbol': symbol, 'position': position})

    except FileNotFoundError:
        print(f"오류: '{input_filename}' 파일을 찾을 수 없습니다.")
        sys.exit(1)
    except Exception as e:
        print(f"파일 처리 중 예기치 않은 오류가 발생했습니다: {e}")
        sys.exit(1)

    # ibrav=12에 해당하는 셀(cell) 행렬 생성
    a = cell_params.get('a', 0)
    b = cell_params.get('b', 0)
    c = cell_params.get('c', 0)
    cosab = cell_params.get('cosab', 0)
    
    if a == 0 or b == 0 or c == 0:
        print("오류: 셀(cell) 파라미터(a, b, c)를 파일에서 찾을 수 없습니다.")
        sys.exit(1)
        
    sinab = (1.0 - cosab**2)**0.5
    cell = np.array([
        [a, 0.0, 0.0],
        [b * cosab, b * sinab, 0.0],
        [0.0, 0.0, c]
    ])
    
    ase_images = []
    for image_data in all_images_data:
        symbols = [atom['symbol'] for atom in image_data]
        positions = [atom['position'] for atom in image_data]
        atoms_obj = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
        ase_images.append(atoms_obj)
        
    if len(ase_images) < 2:
        print("추출된 이미지가 2개 미만입니다. 시작점과 끝점을 모두 추출했는지 확인하세요.")
        sys.exit(1)

    # ASE Atoms 객체를 .traj 파일로 저장
    try:
        write('initial_unoptimized.traj', ase_images[0])
        print("-> 'initial_unoptimized.traj' 파일이 생성되었습니다.")
        
        write('final_unoptimized.traj', ase_images[-1])
        print("-> 'final_unoptimized.traj' 파일이 생성되었습니다.")
    except Exception as e:
        print(f".traj 파일 저장 중 오류가 발생했습니다: {e}")


if __name__ == "__main__":
    parse_and_save_neb()