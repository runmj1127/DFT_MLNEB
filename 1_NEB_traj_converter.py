# 1_NEB_traj_converter.py (Multiple ATOMIC_POSITIONS 에러 해결 버전)

from ase.io import read, write
import sys

def convert_qe_to_traj(input_filename='espresso.neb.in'):
    """
    Quantum ESPRESSO NEB 입력 파일을 읽어 첫 번째와 마지막 이미지를
    .traj 파일로 안정적으로 저장합니다.
    """
    try:
        # ASE를 사용해 첫 번째 이미지(index=0)만 명시적으로 읽어옵니다.
        initial_image = read(input_filename, index=0, format='espresso-in')
        write('initial_unoptimized.traj', initial_image)
        print("-> 'initial_unoptimized.traj' 파일이 생성되었습니다.")
        
        # ASE를 사용해 마지막 이미지(index=-1)만 명시적으로 읽어옵니다.
        final_image = read(input_filename, index=-1, format='espresso-in')
        write('final_unoptimized.traj', final_image)
        print("-> 'final_unoptimized.traj' 파일이 생성되었습니다.")

    except Exception as e:
        print(f"파일 처리 중 오류가 발생했습니다: {e}")
        sys.exit(1)

if __name__ == "__main__":
    convert_qe_to_traj()
