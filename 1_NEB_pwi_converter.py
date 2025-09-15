# 1_pwi_generator.py (모든 파싱 문제를 해결한 최종 버전)

import sys
import re

def create_pwi_files(input_filename='espresso.neb.in'):
    """
    어떠한 순서의 espresso.neb.in 파일이라도 안정적으로 파싱하여, 
    첫 번째와 마지막 이미지를 위한 .pwi 파일을 생성합니다.
    """
    try:
        with open(input_filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        # --- 파일 내용을 섹션별로 저장할 변수 ---
        namelist_text = ""
        species_text = ""
        k_points_text = ""
        first_image_coords = ""
        last_image_coords = ""
        
        # --- 상태를 추적하는 변수 ---
        reading_section = None # 현재 읽고 있는 섹션 (eg. 'NAMELIST', 'SPECIES', 'FIRST', 'LAST')

        for line in lines:
            stripped_line = line.strip()
            u_line = stripped_line.upper()

            # 새로운 섹션 시작 감지
            if 'BEGIN_ENGINE_INPUT' in u_line: continue
            if 'END_ENGINE_INPUT' in u_line: break
            if u_line.startswith('&'): reading_section = 'NAMELIST'
            elif 'ATOMIC_SPECIES' in u_line: reading_section = 'SPECIES'
            elif 'K_POINTS' in u_line: reading_section = 'K_POINTS'
            elif 'FIRST_IMAGE' in u_line: reading_section = 'FIRST'
            elif 'LAST_IMAGE' in u_line: reading_section = 'LAST'
            elif 'INTERMEDIATE_IMAGE' in u_line: reading_section = 'INTERMEDIATE'
            
            # 현재 섹션에 따라 내용 저장
            if reading_section == 'NAMELIST':
                namelist_text += line
                if u_line.startswith('/'): reading_section = None
            elif reading_section == 'SPECIES':
                species_text += line
            elif reading_section == 'K_POINTS':
                k_points_text += line
            elif reading_section == 'FIRST' and len(stripped_line.split()) >= 4:
                first_image_coords += line
            elif reading_section == 'LAST' and len(stripped_line.split()) >= 4:
                last_image_coords += line
        
        if not first_image_coords or not last_image_coords:
            raise ValueError("FIRST_IMAGE 또는 LAST_IMAGE의 좌표를 찾을 수 없습니다.")

        # --- initial_opt.pwi 파일 생성 ---
        with open('initial_opt.pwi', 'w') as f:
            temp_namelists = namelist_text.replace("calculation = 'neb'", "calculation = 'relax'")
            if 'pseudo_dir' not in temp_namelists.lower():
                 temp_namelists = temp_namelists.replace('&CONTROL', "&CONTROL\n    pseudo_dir = './'")
            else:
                 temp_namelists = re.sub(r"pseudo_dir\s*=\s*.*", "    pseudo_dir = './'", temp_namelists, flags=re.IGNORECASE)

            f.write(temp_namelists)
            f.write(species_text + '\n')
            f.write('ATOMIC_POSITIONS {angstrom}\n')
            f.write(first_image_coords)
            f.write(k_points_text)
        print("-> 'initial_opt.pwi' 파일이 생성되었습니다.")

        # --- final_opt.pwi 파일 생성 ---
        with open('final_opt.pwi', 'w') as f:
            f.write(temp_namelists)
            f.write(species_text + '\n')
            f.write('ATOMIC_POSITIONS {angstrom}\n')
            f.write(last_image_coords)
            f.write(k_points_text)
        print("-> 'final_opt.pwi' 파일이 생성되었습니다.")

    except Exception as e:
        print(f"파일 처리 중 오류가 발생했습니다: {e}")
        sys.exit(1)

if __name__ == "__main__":
    create_pwi_files()

