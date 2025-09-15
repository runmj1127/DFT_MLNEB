# 1_pwi_generator.py (ASE 공식 리더를 사용하는 최종 버전)

import sys
from ase.io import read, write

def convert_qe_to_pwi(input_filename='espresso.neb.in'):
    """
    ASE의 공식 리더를 사용해 espresso.neb.in에서 첫 번째와 마지막 이미지를
    별도의 .pwi 파일로 저장합니다.
    """
    try:
        # --- 초기 구조(.pwi) 생성 ---
        # index=0으로 첫 번째 이미지만 명시적으로 읽기
        initial_image = read(input_filename, index=0, format='espresso-in')
        
        # 계산 종류를 'relax'로 변경
        initial_image.calc.parameters['control']['calculation'] = 'relax'
        # pseudo_dir 경로를 './'로 강제 설정
        initial_image.calc.parameters['control']['pseudo_dir'] = './'
        
        write('initial_opt.pwi', initial_image, format='espresso-in')
        print("-> 'initial_opt.pwi' 파일이 생성되었습니다.")
        
        # --- 최종 구조(.pwi) 생성 ---
        # index=-1로 마지막 이미지만 명시적으로 읽기
        final_image = read(input_filename, index=-1, format='espresso-in')

        # 계산 종류를 'relax'로 변경
        final_image.calc.parameters['control']['calculation'] = 'relax'
        # pseudo_dir 경로를 './'로 강제 설정
        final_image.calc.parameters['control']['pseudo_dir'] = './'
        
        write('final_opt.pwi', final_image, format='espresso-in')
        print("-> 'final_opt.pwi' 파일이 생성되었습니다.")

    except Exception as e:
        print(f"파일 처리 중 오류가 발생했습니다: {e}")
        print("입력 파일의 형식을 다시 확인해 주세요.")
        sys.exit(1)

if __name__ == "__main__":
    convert_qe_to_pwi()
