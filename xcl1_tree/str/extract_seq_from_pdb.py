#!/usr/bin/env python3
"""
从PDB文件中提取指定链的序列
遵循细菌编码风格：模块化、自包含、可复用
"""

import csv
import os
from typing import Dict, List, Tuple, Optional


def read_pdb_chain_csv(csv_file: str) -> List[Tuple[str, str]]:
    """
    读取CSV文件，获取PDB和链的对应关系
    
    Args:
        csv_file: CSV文件路径，格式为 "pdb, chain"
        
    Returns:
        List[Tuple[str, str]]: (pdb_id, chain_id) 的列表
        
    Example:
        >>> pairs = read_pdb_chain_csv("align_pdb_list.csv")
        >>> print(pairs[0])  # ('1EL0', 'A')
    """
    pdb_chain_pairs = []
    
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader)  # 跳过标题行
        
        for row in reader:
            if len(row) >= 2:
                pdb_id = row[0].strip()
                chain_id = row[1].strip()
                pdb_chain_pairs.append((pdb_id, chain_id))
    
    return pdb_chain_pairs


def extract_sequence_from_pdb(pdb_file: str, chain_id: str = None, return_resids: bool = False) -> Optional[str]:
    """从PDB文件提取指定链的氨基酸序列
    
    Args:
        pdb_file: PDB文件路径
        chain_id: 链标识符，None表示第一条链
        return_resids: 是否返回残基ID信息
    
    Returns:
        如果return_resids=False: 氨基酸序列字符串，失败返回None
        如果return_resids=True: (sequence, resids) 元组，失败返回None
    
    Example:
        seq = extract_sequence_from_pdb('protein.pdb', 'A')
        seq, resids = extract_sequence_from_pdb('protein.pdb', 'A', return_resids=True)
    """
    # 三字母氨基酸代码到单字母的映射
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # 非标准氨基酸
        'MSE': 'M',  # 硒代蛋氨酸
        'SEC': 'U',  # 硒代半胱氨酸
        'PYL': 'O',  # 吡咯赖氨酸
    }
    
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        # 存储每条链的残基
        chain_residues = {}
        current_chain = None
        
        for line in lines:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':  # 只处理Cα原子
                chain = line[21].strip()
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()
                
                if chain not in chain_residues:
                    chain_residues[chain] = {}
                    
                if current_chain is None:
                    current_chain = chain
                
                # 存储残基信息
                if res_name in aa_map:
                    chain_residues[chain][res_num] = aa_map[res_name]
        
        # 选择目标链（自动处理大小写转换）
        if chain_id:
            target_chain = None
            if chain_id.upper() in chain_residues:
                target_chain = chain_id.upper()
            
            if not target_chain:
                available_chains = list(chain_residues.keys())

                if not available_chains:
                    return None
                target_chain = available_chains[0]

        else:
            target_chain = current_chain
        
        # 构建序列
        residues = chain_residues[target_chain]
        if not residues:
            return None
        
        # 按残基编号排序并连接
        sorted_residues = sorted(residues.items())
        sequence = ''.join(aa for _, aa in sorted_residues)
        
        if return_resids:
            resids = [res_num for res_num, _ in sorted_residues]
            return sequence, resids
        else:
            return sequence
        
    except Exception as e:
    
        return None



def create_fasta_header(pdb_id: str, chain_id: str, sequence_length: int) -> str:
    """
    创建FASTA格式的序列头
    
    Args:
        pdb_id: PDB标识符
        chain_id: 链标识符
        sequence_length: 序列长度
        
    Returns:
        str: FASTA头部信息
    """
    return f">{pdb_id}_{chain_id} Chain {chain_id} "


def extract_sequences_from_csv(csv_file: str, pdb_directory: str = None, 
                              output_file: str = None) -> Dict[str, str]:
    """
    根据CSV文件中的PDB-链对应关系，批量提取序列
    
    Args:
        csv_file: 包含PDB和链信息的CSV文件
        pdb_directory: PDB文件所在目录，默认为CSV文件同目录
        output_file: 输出FASTA文件路径，如果不指定则不保存文件
        
    Returns:
        Dict[str, str]: {pdb_chain_id: sequence} 的字典
        
    Example:
        >>> sequences = extract_sequences_from_csv("align_pdb_list.csv")
        >>> print(f"提取了 {len(sequences)} 个序列")
    """
    if pdb_directory is None:
        pdb_directory = os.path.dirname(os.path.abspath(csv_file))
    
    # 读取PDB-链对应关系
    pdb_chain_pairs = read_pdb_chain_csv(csv_file)
    sequences = {}
    
    print(f"开始处理 {len(pdb_chain_pairs)} 个PDB-链对...")
    
    for pdb_id, chain_id in pdb_chain_pairs:
        # 构建PDB文件路径（转换为小写）
        pdb_file = os.path.join(pdb_directory, f"{pdb_id.lower()}.pdb")
        
        # 提取序列
        sequence = extract_sequence_from_pdb(pdb_file, chain_id)
        
        if sequence:
            key = f"{pdb_id}_{chain_id}"
            sequences[key] = sequence
            print(f"✓ {key}: {len(sequence)} 个氨基酸")
        else:
            print(f"✗ {pdb_id}_{chain_id}: 序列提取失败")
    
    # 保存为FASTA文件
    if output_file and sequences:
        write_fasta_file(sequences, output_file)
        print(f"\n序列已保存到: {output_file}")
    
    return sequences


def write_fasta_file(sequences: Dict[str, str], output_file: str) -> None:
    """
    将序列字典写入FASTA文件
    
    Args:
        sequences: {sequence_id: sequence} 的字典
        output_file: 输出文件路径
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        for seq_id, sequence in sequences.items():
            pdb_id, chain_id = seq_id.split('_')
            header = create_fasta_header(pdb_id, chain_id, len(sequence))
            f.write(f"{header}\n")
            
            # 每行60个字符
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")


def main():
    """主函数：处理当前目录下的align_pdb_list.csv"""
    csv_file = "align_pdb_list.csv"
    output_file = "extracted_sequences.fasta"
    
    if not os.path.exists(csv_file):
        print(f"错误: 找不到文件 {csv_file}")
        return
    
    print("🧬 开始从PDB文件提取序列...")
    sequences = extract_sequences_from_csv(csv_file, output_file=output_file)
    
    if sequences:
        print(f"\n✅ 成功提取 {len(sequences)} 个序列")
        print("序列统计:")
        for seq_id, seq in sequences.items():
            print(f"  {seq_id}: {len(seq)} aa")
    else:
        print("❌ 未提取到任何序列")


if __name__ == "__main__":
    main()
