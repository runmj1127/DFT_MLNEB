# 1_pwi_generator.py (pseudo_dir를 './'로 고정하는 최종 버전)

import sys
import re

def create_pwi_files(input_filename='espresso.neb.in'):
    """
    espresso.neb.in 파일을 파싱하여, 첫 번째와 마지막 이미지를 위한
    별도의 Quantum ESPRESSO 입력 파일(.pwi) 두 개를 생성합니다.
    pseudo_dir는 './'로 고정됩니다.
    """
    
    namelist_text = ""
    species_text = ""
    initial_coords_text = ""
    final_coords_text = ""
    k_points_text = ""

    in_namelist = False
    in_species = False
    current_image = None

    with open(input_filename, 'r', encoding='utf-8') as f:
        for line in f:
            s_line, u_line = line.strip(), line.strip().upper()

            if u_line.startswith('&'): in_namelist = True
            if in_namelist:
                # 기존 pseudo_dir 라인은 무시
                if 'pseudo_dir' in s_line.lower():
                    continue
                # neb 계산 종류를 relax로 변경
                if 'calculation' in s_line.lower() and 'neb' in s_line.lower():
                    line = "    calculation = 'relax'\n"
                namelist_text += line
                if u_line.startswith('/'):
                    # &CONTROL 섹션이 끝나는 지점에 pseudo_dir 추가
                    if '&CONTROL' in namelist_text:
                         namelist_text = namelist_text.replace('&CONTROL', "&CONTROL\n    pseudo_dir = './'")
                    in_namelist = False
            
            elif 'ATOMIC_SPECIES' in u_line:
                in_species = True
                species_text += line
                continue
            elif in_species and 'BEGIN_POSITIONS' not in u_line:
                species_text += line
            elif 'BEGIN_POSITIONS' in u_line:
                in_species = False
            
            elif 'K_POINTS' in u_line:
                k_points_text += line

            elif 'FIRST_IMAGE' in u_line: current_image = 'FIRST'; continue
            elif 'LAST_IMAGE' in u_line: current_image = 'LAST'; continue
            elif 'ATOMIC_POSITIONS' in u_line: continue
            
            if current_image == 'FIRST' and len(s_line.split()) >= 4:
                initial_coords_text += line
            elif current_image == 'LAST' and len(s_line.split()) >= 4:
                final_coords_text += line

    # --- initial_opt.pwi 파일 생성 ---
    with open('initial_opt.pwi', 'w') as f:
        f.write(namelist_text)
        f.write(species_text)
        f.write('ATOMIC_POSITIONS {angstrom}\n')
        f.write(initial_coords_text)
        f.write(k_points_text)
    print("-> 'initial_opt.pwi' 파일이 생성되었습니다.")

    # --- final_opt.pwi 파일 생성 ---
    with open('final_opt.pwi', 'w') as f:
        f.write(namelist_text)
        f.write(species_text)
        f.write('ATOMIC_POSITIONS {angstrom}\n')
        f.write(final_coords_text)
        f.write(k_points_text)
    print("-> 'final_opt.pwi' 파일이 생성되었습니다.")


if __name__ == "__main__":
    create_pwi_files()
