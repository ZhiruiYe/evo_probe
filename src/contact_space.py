from typing import List, Tuple, Dict, Optional
from sklearn.preprocessing import StandardScaler
import numpy as np


def load_mj_matrix(mj_file_path: str, standardize: bool = False) -> Dict[Tuple[str, str], float]:
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


class ContactSpace:
    """接触空间类，管理多个node的序列和对应的embedding"""
    
    def __init__(self, target_contacts: List[Tuple[int, int]], mj_dict: Dict[Tuple[str, str], float]):
        self.target_contacts = target_contacts
        self.mj_dict = mj_dict
        self.node_sequences = {}  # node_id -> sequences
        self.node_indices = {}    # node_id -> (start_idx, end_idx)
        self.all_sequences = []
        self.embeddings = None
        self.scaler = None
        
    def add_node_sequences(self, node_id: str, sequences: List[str]):
        """添加节点的序列"""
        start_idx = len(self.all_sequences)
        self.all_sequences.extend(sequences)
        end_idx = len(self.all_sequences)
        
        self.node_sequences[node_id] = sequences
        self.node_indices[node_id] = (start_idx, end_idx)
        
    def add_node_from_fasta(self, node_id: str, fasta_path: str):
        """从FASTA文件添加节点序列"""
        from Bio import SeqIO
        
        sequences = []
        with open(fasta_path, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                sequences.append(str(record.seq))
        
        self.add_node_sequences(node_id, sequences)
        
    def build_embeddings(self):
        """构建所有序列的embedding"""
        if not self.all_sequences:
            raise ValueError("没有添加任何序列")
            
        m, d = len(self.all_sequences), len(self.target_contacts)
        X = np.zeros((m, d))
        
        # 构建embedding
        for seq_idx, sequence in enumerate(self.all_sequences):
            for contact_idx, (i, j) in enumerate(self.target_contacts):
                if i < len(sequence) and j < len(sequence):
                    res_i, res_j = sequence[i], sequence[j]
                    if res_i != '-' and res_j != '-':
                        mj_energy = self.mj_dict.get((res_i.upper(), res_j.upper()), 0)
                        X[seq_idx, contact_idx] = mj_energy
        
        # 标准化
        self.scaler = StandardScaler()
        self.embeddings = self.scaler.fit_transform(X)
        
    def get_node_embeddings(self, node_id: str) -> Optional[np.ndarray]:
        """获取指定节点的embedding"""
        if self.embeddings is None:
            raise ValueError("请先调用build_embeddings()构建embedding")
            
        if node_id not in self.node_indices:
            return None
            
        start_idx, end_idx = self.node_indices[node_id]
        return self.embeddings[start_idx:end_idx]
        
    def get_node_sequences(self, node_id: str) -> Optional[List[str]]:
        """获取指定节点的序列"""
        return self.node_sequences.get(node_id)
        
    def get_all_embeddings(self) -> np.ndarray:
        """获取所有embedding"""
        if self.embeddings is None:
            raise ValueError("请先调用build_embeddings()构建embedding")
        return self.embeddings
        
    def get_node_list(self) -> List[str]:
        """获取所有节点ID列表"""
        return list(self.node_sequences.keys())
        
    def get_sequence_count(self, node_id: Optional[str] = None) -> int:
        """获取序列数量"""
        if node_id is None:
            return len(self.all_sequences)
        return len(self.node_sequences.get(node_id, []))


