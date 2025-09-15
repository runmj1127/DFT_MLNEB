# 1_pwi_generator.py (모든 구조의 파싱 오류를 해결한 최종 버전)

import sys
import re

def create_pwi_files(input_filename='espresso.neb.in'):
    """
    어떠한 순서의 espresso.neb.in 파일이라도 안정적으로 파싱하여, 
    첫 번째와 마지막 이미지를 위한 .pwi 파일을 생성합니다.
    """
    try:
        with open(input_filename, 'r', encoding='utf-8') as f:
            content = f.read()

        # --- 모든 섹션을 정규식으로 한번에 추출 ---
        engine_input_match = re.search(r'BEGIN_ENGINE_INPUT(.*?)END_ENGINE_INPUT', content, re.DOTALL)
        if not engine_input_match:
            raise ValueError("BEGIN_ENGINE_INPUT 또는 END_ENGINE_INPUT을 찾을 수 없습니다.")
        engine_input = engine_input_match.group(1)
        
        namelists_text = ""
        for section in ['CONTROL', 'SYSTEM', 'ELECTRONS']:
            match = re.search(f'&{section}(.*?)/', engine_input, re.IGNORECASE | re.DOTALL)
            if match:
                namelists_text += match.group(0) + '\n\n'
        
        species_text = re.search(r'ATOMIC_SPECIES.*?(?=\n\s*(?:&|K_POINTS|BEGIN_POSITIONS|CELL_PARAMETERS|$))', engine_input, re.IGNORECASE | re.DOTALL).group(0)
        k_points_text = re.search(r'K_POINTS.*?(?:\n.*)', engine_input, re.IGNORECASE).group(0)

        positions_block_match = re.search(r'BEGIN_POSITIONS(.*?)END_POSITIONS', content, re.DOTALL)
        if not positions_block_match:
            raise ValueError("BEGIN_POSITIONS 또는 END_POSITIONS를 찾을 수 없습니다.")
        positions_block = positions_block_match.group(1)
        
        # 정규식을 사용하여 각 이미지 블록을 정확히 분리
        images = re.findall(r'(FIRST_IMAGE|LAST_IMAGE)\s*\n\s*ATOMIC_POSITIONS.*?\n(.*?)(?=FIRST_IMAGE|INTERMEDIATE_IMAGE|LAST_IMAGE|END_POSITIONS)', positions_block, re.DOTALL)
        
        image_data = dict(images)
        first_image_coords_text = image_data.get('FIRST_IMAGE')
        last_image_coords_text = image_data.get('LAST_IMAGE')

        if not first_image_coords_text or not last_image_coords_text:
            raise ValueError("FIRST_IMAGE 또는 LAST_IMAGE의 좌표 블록을 찾을 수 없습니다.")

        # --- initial_opt.pwi 파일 생성 ---
        with open('initial_opt.pwi', 'w') as f:
            # calculation 종류를 'relax'로, pseudo_dir를 './'로 강제 설정
            temp_namelists = re.sub(r"calculation\s*=\s*['\"]\s*neb\s*['\"]", "calculation = 'relax'", namelists, flags=re.IGNORECASE)
            if 'pseudo_dir' not in temp_namelists.lower():
                 temp_namelists = temp_namelists.replace('&CONTROL', "&CONTROL\n    pseudo_dir = './'")
            else:
                 temp_namelists = re.sub(r"pseudo_dir\s*=\s*.*", "    pseudo_dir = './'", temp_namelists, flags=re.IGNORECASE)

            f.write(temp_namelists)
            f.write(species_text + '\n\n')
            f.write('ATOMIC_POSITIONS {angstrom}\n')
            f.write(first_image_coords_text.strip() + '\n\n')
            f.write(k_points_text + '\n')
        print("-> 'initial_opt.pwi' 파일이 생성되었습니다.")

        # --- final_opt.pwi 파일 생성 ---
        with open('final_opt.pwi', 'w') as f:
            f.write(temp_namelists)
            f.write(species_text + '\n\n')
            f.write('ATOMIC_POSITIONS {angstrom}\n')
            f.write(last_image_coords_text.strip() + '\n\n')
            f.write(k_points_text + '\n')
        print("-> 'final_opt.pwi' 파일이 생성되었습니다.")

    except Exception as e:
        print(f"파일 처리 중 오류가 발생했습니다: {e}")
        sys.exit(1)

if __name__ == "__main__":
    create_pwi_files()
