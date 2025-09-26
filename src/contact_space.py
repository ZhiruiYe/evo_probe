from typing import List, Tuple, Dict
from sklearn.preprocessing import StandardScaler
import numpy as np

def load_mj_matrix(mj_file_path: str,standardize: bool = False) -> Dict[Tuple[str, str], float]:
    """加载MJ矩阵文件"""
    with open(mj_file_path, 'r') as f:
        lines = f.readlines()
    
    header = lines[0].strip().split()
    amino_acids = header
    n_aa = len(amino_acids)
    
    matrix = np.zeros((n_aa, n_aa))
    for i, line in enumerate(lines[1:]):
        values = line.strip().split()
        for j, value in enumerate(values):
            if value != '0.00e+00':
                matrix[i, j] = float(value)
                matrix[j, i] = float(value)
    
    if standardize:
        scaler = StandardScaler()
        matrix = scaler.fit_transform(matrix)
    
    # 创建查询字典
    mj_dict = {}
    for i, aa1 in enumerate(amino_acids):
        for j, aa2 in enumerate(amino_acids):
            mj_dict[(aa1, aa2)] = matrix[i, j]
            mj_dict[(aa2, aa1)] = matrix[i, j]
    
    return mj_dict


def generate_contact_embedding(msa_path: str,
                      target_contacts: List[Tuple[int, int]],
                      mj_dict: Dict,
                      standardize: bool = True) -> np.ndarray:
    from Bio import SeqIO
    from sklearn.preprocessing import StandardScaler    

    # 读取序列
    sequences = []
    with open(msa_path, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences.append(str(record.seq))
    
    m, d = len(sequences), len(target_contacts)
    X = np.zeros((m, d))
    
    # 构建embedding
    for seq_idx, sequence in enumerate(sequences):
        for contact_idx, (i, j) in enumerate(target_contacts):
            if i < len(sequence) and j < len(sequence):
                res_i, res_j = sequence[i], sequence[j]
                if res_i != '-' and res_j != '-':
                        mj_energy = mj_dict.get((res_i.upper(), res_j.upper()), 0)
                        X[seq_idx, contact_idx] = mj_energy  
    
    # 标准化
    if standardize:
        scaler = StandardScaler()
        X = scaler.fit_transform(X)

    return X