# check_energy.py (에너지 읽기 테스트 스크립트)
from ase.io import read

# 확인할 파일 이름
pwo_file = 'initial_opt.pwo'

try:
    print(f"--- Reading '{pwo_file}' ---")
    
    # .pwo 파일에서 구조와 계산 결과를 읽어옵니다.
    atoms = read(pwo_file, format='espresso-out')
    
    # Atoms 객체에 저장된 에너지를 직접 가져옵니다.
    energy = atoms.get_potential_energy()
    
    print("\n[SUCCESS]")
    print(f"File '{pwo_file}' contains energy information.")
    print(f"Potential Energy: {energy} eV")

except Exception as e:
    print("\n[FAILED]")
    print(f"Could not read energy from '{pwo_file}'.")
    print(f"Error details: {e}")
