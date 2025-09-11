# 1_NEB_traj_converter.py (ASE를 사용하여 안정적으로 수정된 최종 버전)

from ase.io import read, write
import sys

def convert_qe_to_traj(input_filename='espresso.neb.in'):
    """
    Quantum ESPRESSO 입력 파일을 읽어 첫 번째와 마지막 이미지를
    unoptimized .traj 파일로 저장합니다.
    """
    try:
        # ASE의 강력한 파서를 사용하여 QE 입력 파일에서 모든 이미지(구조)를 읽어옵니다.
        # index=':'는 모든 이미지를 읽으라는 의미입니다.
        images = read(input_filename, index=':', format='espresso-in')
        
        if len(images) < 2:
            print(f"오류: '{input_filename}'에서 2개 미만의 이미지를 찾았습니다.")
            print("파일에 FIRST_IMAGE와 LAST_IMAGE가 모두 포함되어 있는지 확인하세요.")
            sys.exit(1)

        # 첫 번째 이미지를 저장합니다.
        write('initial_unoptimized.traj', images[0])
        print("-> 'initial_unoptimized.traj' 파일이 생성되었습니다.")
        
        # 마지막 이미지를 저장합니다.
        write('final_unoptimized.traj', images[-1])
        print("-> 'final_unoptimized.traj' 파일이 생성되었습니다.")

    except FileNotFoundError:
        print(f"오류: 입력 파일 '{input_filename}'을 찾을 수 없습니다.")
        sys.exit(1)
    except Exception as e:
        print(f"파일 처리 중 예기치 않은 오류가 발생했습니다: {e}")
        print("입력 파일의 형식이 올바른지 다시 확인해 주세요.")
        sys.exit(1)

if __name__ == "__main__":
    convert_qe_to_traj()
