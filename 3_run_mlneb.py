# run_optimization.py (ASE 계산기를 사용하지 않고 직접 QE를 실행하는 최종 스크립트)

import sys
import numpy as np
import subprocess
from ase import Atoms
from ase.io import read, write

# --- 1. 사용자 설정 ---
INPUT_QE_FILE = 'espresso.neb.in'
PSEUDO_DIR = './'
N_CORES = 16
OPTIMIZE_FMAX = 0.1 # 이 스크립트에서는 직접 사용되지 않음
OPTIMIZED_INITIAL_FILE = 'initial_optimized.traj'
OPTIMIZED_FINAL_FILE = 'final_optimized.traj'

# --- 2. 수동 파서 (구조와 설정 텍스트만 추출) ---
def final_parser(filename):
    cell_params, pseudos, namelist_text = {}, {}, ""
    initial_coords, final_coords = [], []
    current_block, current_image, in_species, in_namelist = None, None, False, False

    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            s_line, u_line = line.strip(), line.strip().upper()
            if not s_line or s_line.startswith('!'): continue
            if 'BEGIN_ENGINE_INPUT' in u_line: continue
            if 'END_ENGINE_INPUT' in u_line: break
            if u_line.startswith('&'): in_namelist = True
            if in_namelist:
                namelist_text += line
                if u_line.startswith('/'): in_namelist = False
                if '=' in s_line:
                    key = s_line.split('=')[0].strip().lower()
                    if key in ['a', 'b', 'c', 'cosab']:
                        cell_params[key] = float(s_line.split('=')[1].split('!')[0].strip())
            elif 'ATOMIC_SPECIES' in u_line: in_species = True; continue
            elif 'BEGIN_POSITIONS' in u_line: in_species = False; continue
            if in_species:
                parts = s_line.split()
                if len(parts) >= 3: pseudos[parts[0]] = parts[2]
            if 'FIRST_IMAGE' in u_line: current_image = initial_coords; continue
            if 'LAST_IMAGE' in u_line: current_image = final_coords; continue
            if 'ATOMIC_POSITIONS' in u_line: continue
            if current_image is not None and len(s_line.split()) >= 4 and s_line.split()[0].isalpha():
                current_image.append(s_line)
    
    a, b, c, cosab = cell_params.get('a'), cell_params.get('b'), cell_params.get('c'), cell_params.get('cosab', 0)
    sinab = (1.0 - cosab**2)**0.5
    cell = np.array([[a, 0.0, 0.0], [b * cosab, b * sinab, 0.0], [0.0, 0.0, c]])
    
    initial_atoms = Atoms(symbols=[line.split()[0] for line in initial_coords], 
                          positions=[[float(p) for p in line.split()[1:4]] for line in initial_coords], 
                          cell=cell, pbc=True)
    final_atoms = Atoms(symbols=[line.split()[0] for line in final_coords], 
                        positions=[[float(p) for p in line.split()[1:4]] for line in final_coords], 
                        cell=cell, pbc=True)
    return initial_atoms, final_atoms, pseudos, namelist_text

# --- 3. QE 실행 함수 ---
def run_qe_optimization(label, atoms, pseudos, namelist_text):
    pwi_filename = f'{label}.pwi'
    pwo_filename = f'{label}.pwo'
    
    with open(pwi_filename, 'w') as f:
        f.write(namelist_text.replace("calculation = 'neb'", "calculation = 'relax'"))
        f.write(f'ATOMIC_SPECIES\n')
        for symbol, pfile in pseudos.items():
            f.write(f' {symbol} 1.0 {PSEUDO_DIR}{pfile}\n')
        f.write('ATOMIC_POSITIONS {angstrom}\n')
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            f.write(f' {symbol} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n')
        f.write('K_POINTS {gamma}\n 1 1 1 0 0 0\n')

    command = f"mpirun -np {N_CORES} pw.x -in {pwi_filename}"
    print(f"\n>>> Running command: {command}")
    
    with open(pwo_filename, 'w') as out_file:
        result = subprocess.run(command, shell=True, stdout=out_file, stderr=subprocess.PIPE, text=True)

    with open(pwo_filename, 'r') as out_file: pwo_content = out_file.read()

    if result.returncode != 0 or "JOB DONE." not in pwo_content:
         print(f"\n!!! QE calculation FAILED for {label} !!!")
         print("--- QE Error Message (from stderr) ---")
         print(result.stderr)
         print("\n--- QE Output File Content (last 50 lines) ---")
         print(''.join(pwo_content.splitlines(True)[-50:]))
         sys.exit(1)
    
    return read(pwo_filename, format='espresso-out')

# --- 4. 메인 워크플로우 ---
def run_workflow():
    print("--- 0단계: espresso.neb.in 파일 직접 파싱 ---")
    initial_atoms, final_atoms, pseudos, namelist_text = final_parser(INPUT_QE_FILE)

    print("\n--- 1단계: 구조 최적화 시작 ---")
    
    optimized_initial = run_qe_optimization('initial_opt', initial_atoms, pseudos, namelist_text)
    write(OPTIMIZED_INITIAL_FILE, optimized_initial)
    print(f">>> Initial structure optimization complete!")

    optimized_final = run_qe_optimization('final_opt', final_atoms, pseudos, namelist_text)
    write(OPTIMIZED_FINAL_FILE, optimized_final)
    print(f">>> Final structure optimization complete!")

if __name__ == "__main__":
    run_workflow()
    print("\n--- 최적화 작업이 성공적으로 완료되었습니다. ---")
