###PDB 파일에서 리간드 정보를 추출하고, 이를 SMILES 형식으로 변환하여 CSV 파일로 저장
import os
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem

# 물 분자를 나타내는 리지듀 이름 리스트
water_residues = ['HOH', 'WAT', 'H2O']

# PDB 파일을 읽어 리간드 이름과 SMILES 추출
def extract_ligands(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    ligands = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ' and residue.resname not in water_residues:
                    # 리간드 이름 추출
                    ligand_name = residue.resname
                    
                    # 리간드의 원자 정보 추출
                    ligand_atoms = []
                    for atom in residue:
                        ligand_atoms.append(atom)
                    
                    # 리간드의 PDB 파일 생성
                    io = PDBIO()
                    io.set_structure(residue)
                    
                    class LigandSelect(Select):
                        def accept_residue(self, residue):
                            return True

                    ligand_pdb_file = f"{ligand_name}.pdb"
                    io.save(ligand_pdb_file, LigandSelect())
                    
                    # RDKit으로 SMILES 추출
                    ligand_mol = Chem.MolFromPDBFile(ligand_pdb_file, sanitize=True, removeHs=True)
                    if ligand_mol:
                        smiles = Chem.MolToSmiles(ligand_mol)
                        ligands[ligand_name] = smiles
                    
                    # 생성된 리간드 PDB 파일 삭제
                    os.remove(ligand_pdb_file)

    return ligands

# 현재 디렉토리에 있는 모든 pdb 파일을 찾음
pdb_files = [f for f in os.listdir('.') if f.endswith('.pdb')]

# 리간드 정보 저장을 위한 리스트
ligand_data = []

for pdb_file in pdb_files:
    ligands = extract_ligands(pdb_file)
    for ligand, smiles in ligands.items():
        ligand_data.append({'PDB_File': pdb_file, 'Ligand': ligand, 'SMILES': smiles})

# 데이터프레임으로 변환
df = pd.DataFrame(ligand_data)

# 새로운 CSV 파일로 저장
output_csv = 'ligand_smiles.csv'
df.to_csv(output_csv, index=False)

print(f"리간드 이름과 SMILES 정보가 '{output_csv}' 파일에 저장되었습니다.")
