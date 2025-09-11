# final_test_with_cell.py (ibrav, a, b, c 등 모든 정보를 포함한 최종 테스트)
import copy
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso
import numpy as np

print(">>> Test script with original cell parameters started.")

# 1. 원본 espresso.neb.in의 셀(cell) 정보를 그대로 사용
a = 13.1337185
b = 13.6537423
c = 28.7501618
cosab = -0.61543451
sinab = (1.0 - cosab**2)**0.5
cell = np.array([
    [a, 0.0, 0.0],
    [b * cosab, b * sinab, 0.0],
    [0.0, 0.0, c]
])

# Atoms 객체를 파일에서 읽는 대신 코드에서 직접 생성
initial_atoms = Atoms(symbols=['O', 'H', 'H'],
                        positions=[[5.0, 5.0, 5.0],
                                   [5.7, 5.0, 5.0],
                                   [4.3, 5.0, 5.0]],
                        cell=cell,
                        pbc=True)

# 2. QE 계산기 설정 (원본 espresso.neb.in과 거의 동일하게)
pseudopotentials = {
    'H': 'H.pbe-rrkjus.UPF',
    'O': 'O.pbe-rrkjus.UPF'
}
input_data = {
    'control': {'calculation': 'relax', 'pseudo_dir': './pseudo/'},
    'system': {
        'ibrav': 12,  # 원본 값
        'a': a,
        'b': b,
        'c': c,
        'cosab': cosab,
        'nat': 3, 
        'ntyp': 2, 
        'ecutwfc': 30.0, 
        'ecutrho': 240.0,
        'occupations': "'smearing'",
        'smearing': "'gaussian'"
    },
    'electrons': {'mixing_beta': 0.7}
}
kpts = (1, 1, 1)

ase_calculator = Espresso(
    label='qe_calc_cell_test',
    nprocs=4,
    executable='pw.x',
    pseudopotentials=pseudopotentials,
    input_data=input_data,
    kpts=kpts)

# 3. 계산기 연결 및 최적화 실행
try:
    print(">>> Attaching calculator and starting optimization...")
    initial_atoms.calc = ase_calculator
    optimizer = BFGS(initial_atoms, trajectory='test_optimization_cell.traj')
    optimizer.run(fmax=0.1)
    print("\\n>>> FINAL TEST SUCCEEDED! Optimization with original cell info finished.")
except Exception as e:
    print(f"\\n>>> FINAL TEST FAILED with error: {e}")









































 







