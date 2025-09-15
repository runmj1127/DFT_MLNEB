# 1_pwi_generator.py (파싱 버그를 모두 수정한 최종 버전)

import sys
import re

def create_pwi_files(input_filename='espresso.neb.in'):
    """
    espresso.neb.in 파일을 안정적으로 파싱하여, 첫 번째와 마지막 이미지를 위한
    별도의 Quantum ESPRESSO 입력 파일(.pwi) 두 개를 생성합니다.
    """
    
    # 파일 전체를 읽어 섹션별로 분리
    with open(input_filename, 'r', encoding='utf-8') as f:
        content = f.read()

    engine_input = re.search(r'BEGIN_ENGINE_INPUT(.*?)END_ENGINE_INPUT', content, re.DOTALL).group(1)
    
    # Namelist, Species, K-points 텍스트 추출
    namelists = ""
    for section in ['CONTROL', 'SYSTEM', 'ELECTRONS']:
        match = re.search(f'&{section}(.*?)/', engine_input, re.IGNORECASE | re.DOTALL)
        if match:
            namelists += match.group(0) + '\n\n'

    species = re.search(r'ATOMIC_SPECIES.*?(?:\n\s*\w+\s+[\d.]+\s+\S+)+', engine_input, re.IGNORECASE | re.DOTALL).group(0)
    k_points = re.search(r'K_POINTS.*?(?:\n.*)', engine_input, re.IGNORECASE).group(0)

    # 이미지(좌표) 텍스트 추출
    positions_block = re.search(r'BEGIN_POSITIONS(.*?)END_POSITIONS', content, re.DOTALL).group(1)
    
    # 정규식을 사용하여 각 이미지 블록을 정확히 분리
    images = re.findall(r'(FIRST_IMAGE|INTERMEDIATE_IMAGE|LAST_IMAGE)\s*\n(.*?)(?=FIRST_IMAGE|INTERMEDIATE_IMAGE|LAST_IMAGE|END_POSITIONS)', positions_block, re.DOTALL)
    
    first_image_coords = ""
    last_image_coords = ""

    for image_type, coords in images:
        if image_type == 'FIRST_IMAGE':
            first_image_coords = coords
        elif image_type == 'LAST_IMAGE':
            last_image_coords = coords

    if not first_image_coords or not last_image_coords:
        print("오류: FIRST_IMAGE 또는 LAST_IMAGE를 찾을 수 없습니다.")
        sys.exit(1)

    # --- initial_opt.pwi 파일 생성 ---
    with open('initial_opt.pwi', 'w') as f:
        f.write(namelists.replace("calculation = 'neb'", "calculation = 'relax'"))
        f.write(species + '\n\n')
        f.write(first_image_coords.strip() + '\n\n')
        f.write(k_points + '\n')
    print("-> 'initial_opt.pwi' 파일이 생성되었습니다.")

    # --- final_opt.pwi 파일 생성 ---
    with open('final_opt.pwi', 'w') as f:
        f.write(namelists.replace("calculation = 'neb'", "calculation = 'relax'"))
        f.write(species + '\n\n')
        f.write(last_image_coords.strip() + '\n\n')
        f.write(k_points + '\n')
    print("-> 'final_opt.pwi' 파일이 생성되었습니다.")

if __name__ == "__main__":
    create_pwi_files()
