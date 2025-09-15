# 1_pwi_generator.py (첫 번째/마지막 이미지의 .pwi 파일을 직접 생성)

import sys
import re

def create_pwi_files(input_filename='espresso.neb.in'):
    """
    espresso.neb.in 파일을 파싱하여, 첫 번째와 마지막 이미지를 위한
    별도의 Quantum ESPRESSO 입력 파일(.pwi) 두 개를 생성합니다.
    """
    
    # 파일 전체를 읽어 섹션별로 분리
    with open(input_filename, 'r', encoding='utf-8') as f:
        content = f.read()

    engine_input_match = re.search(r'BEGIN_ENGINE_INPUT(.*?)END_ENGINE_INPUT', content, re.DOTALL)
    if not engine_input_match:
        print("오류: BEGIN_ENGINE_INPUT 또는 END_ENGINE_INPUT을 찾을 수 없습니다.")
        sys.exit(1)
    engine_input = engine_input_match.group(1)
    
    # Namelist, Species, K-points 텍스트 추출
    namelists = ""
    for section in ['CONTROL', 'SYSTEM', 'ELECTRONS']:
        match = re.search(f'&{section}(.*?)/', engine_input, re.IGNORECASE | re.DOTALL)
        if match:
            namelists += match.group(0) + '\n\n'

    species_match = re.search(r'ATOMIC_SPECIES.*?(?=\n\s*BEGIN_POSITIONS|\Z)', engine_input, re.IGNORECASE | re.DOTALL)
    species = species_match.group(0) if species_match else ""

    k_points_match = re.search(r'K_POINTS.*?(?:\n.*)', engine_input, re.IGNORECASE)
    k_points = k_points_match.group(0) if k_points_match else "K_POINTS {gamma}\n 1 1 1 0 0 0"

    # 이미지(좌표) 텍스트 추출
    positions_block_match = re.search(r'BEGIN_POSITIONS(.*?)END_POSITIONS', content, re.DOTALL)
    if not positions_block_match:
        print("오류: BEGIN_POSITIONS 또는 END_POSITIONS를 찾을 수 없습니다.")
        sys.exit(1)
    positions_block = positions_block_match.group(1)
    
    images = re.findall(r'(FIRST_IMAGE|LAST_IMAGE)\s*\n(.*?)(?=FIRST_IMAGE|INTERMEDIATE_IMAGE|LAST_IMAGE|END_POSITIONS)', positions_block, re.DOTALL)
    
    image_data = dict(images)
    first_image_coords = image_data.get('FIRST_IMAGE')
    last_image_coords = image_data.get('LAST_IMAGE')

    if not first_image_coords or not last_image_coords:
        print("오류: FIRST_IMAGE 또는 LAST_IMAGE를 찾을 수 없습니다.")
        sys.exit(1)

    # --- initial_opt.pwi 파일 생성 ---
    with open('initial_opt.pwi', 'w') as f:
        temp_namelists = re.sub(r"calculation\s*=\s*['\"]\s*neb\s*['\"]", "calculation = 'relax'", namelists, flags=re.IGNORECASE)
        f.write(temp_namelists)
        f.write(species + '\n\n')
        f.write(first_image_coords.strip() + '\n\n')
        f.write(k_points + '\n')
    print("-> 'initial_opt.pwi' 파일이 생성되었습니다.")

    # --- final_opt.pwi 파일 생성 ---
    with open('final_opt.pwi', 'w') as f:
        f.write(temp_namelists)
        f.write(species + '\n\n')
        f.write(last_image_coords.strip() + '\n\n')
        f.write(k_points + '\n')
    print("-> 'final_opt.pwi' 파일이 생성되었습니다.")

if __name__ == "__main__":
    create_pwi_files()
