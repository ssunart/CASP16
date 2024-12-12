###PDB 파일들 중에서 리간드가 포함되어 있지 않은 파일들을 찾아내고, 이를 리스트로 출력한 뒤 파일로 저장
import os
from Bio.PDB import PDBParser

# 현재 디렉토리에 있는 모든 pdb 파일을 찾음
pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]

# 리간드가 없는 파일들을 저장할 리스트
no_ligand_files = []

# PDBParser 객체 생성
parser = PDBParser(QUIET=True)

# 물 분자를 나타내는 리지듀 이름 리스트
water_residues = ['HOH', 'WAT', 'H2O']

for pdb_file in pdb_files:
    # PDB 파일 파싱
    structure = parser.get_structure('structure', pdb_file)
    
    has_ligand = False
    
    # 각 모델, 체인, 리지듀를 순회
    for model in structure:
        for chain in model:
            for residue in chain:
                # HETATM (헤테로 원자)이며 물 분자가 아닌 경우 리간드로 간주
                if residue.id[0] != ' ' and residue.resname not in water_residues:
                    has_ligand = True
                    break
            if has_ligand:
                break
        if has_ligand:
            break

    if not has_ligand:
        no_ligand_files.append(pdb_file)

# 리간드가 없는 파일 리스트 출력
print("리간드가 없는 파일 리스트:", no_ligand_files)

# 리간드가 없는 파일 리스트를 텍스트 파일로 저장
with open('no_ligand_files.txt', 'w') as f:
    for file in no_ligand_files:
        f.write(file + '\n')

print("리간드가 없는 파일 리스트가 'no_ligand_files.txt' 파일에 저장되었습니다.")
