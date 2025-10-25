import numpy as np
import re
from typing import List, Tuple, Dict, Optional

# 氨基酸三字母到单字母的转换字典
AA_THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def parse_pdb_coordinates(pdb_path: str, 
                         chain_id: Optional[str] = None) -> Tuple[Dict[int, np.ndarray], Dict[int, str]]:
    """
    从PDB文件中提取残基的Cβ原子坐标（GLY使用Cα）和氨基酸类型
    
    Args:
        pdb_path: PDB文件路径
        chain_id: 指定链ID，None表示使用第一条链
        
    Returns:
        (coordinates, residue_types): ({残基编号: 坐标数组}, {残基编号: 氨基酸类型}) 元组
    """
    coordinates = {}
    residue_types = {}
    
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"PDB文件不存在: {pdb_path}")
        return {}, {}
    
    target_chain = None
    
    for line in lines:
        if not line.startswith('ATOM'):
            continue
            
        # 解析PDB ATOM行
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        current_chain = line[21].strip()
        res_num = int(line[22:26].strip())
        
        # 确定目标链
        if target_chain is None:
            target_chain = current_chain if chain_id is None else chain_id
        
        # 跳过非目标链
        if current_chain != target_chain:
            continue
            
        # 选择合适的原子：GLY用Cα，其他用Cβ
        if (res_name == 'GLY' and atom_name == 'CA') or \
           (res_name != 'GLY' and atom_name == 'CB'):
            
            # 提取坐标
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            
            coordinates[res_num] = np.array([x, y, z])
            residue_types[res_num] = AA_THREE_TO_ONE.get(res_name, res_name)  # 转换为单字母，未知类型保持原样
    
    print(f"PDB坐标提取完成: {len(coordinates)} 个残基 (链 {target_chain})")
    return coordinates, residue_types

def calculate_distance_matrix(coordinates: Dict[int, np.ndarray]) -> Tuple[np.ndarray, List[int]]:
    """
    计算残基间的欧氏距离矩阵
    
    Args:
        coordinates: {残基编号: 坐标} 字典
        
    Returns:
        (distance_matrix, residue_indices): 距离矩阵和残基编号列表
    """
    if not coordinates:
        return np.array([]), []
    
    # 获取排序的残基编号
    residue_indices = sorted(coordinates.keys())
    n_residues = len(residue_indices)
    
    # 初始化距离矩阵
    distance_matrix = np.zeros((n_residues, n_residues))
    
    # 计算所有残基对的距离
    for i, res_i in enumerate(residue_indices):
        for j, res_j in enumerate(residue_indices):
            if i <= j:  # 只计算上三角，利用对称性
                coord_i = coordinates[res_i]
                coord_j = coordinates[res_j]
                
                # 欧氏距离
                distance = np.linalg.norm(coord_i - coord_j)
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance  # 对称填充
    
    print(f"距离矩阵计算完成: {n_residues} x {n_residues}")
    return distance_matrix, residue_indices


def create_binary_contact_matrix(distance_matrix: np.ndarray, threshold: float = 8.0, remove_diag: int = 3) -> np.ndarray:
    """
    将距离矩阵转换为二元接触矩阵
    
    Args:
        distance_matrix: 残基间距离矩阵
        threshold: 接触距离阈值（Å）
        remove_diag: 忽略主对角线附近的残基对数量
        
    Returns:
        binary_contact_matrix: 二元接触矩阵（1=接触，0=无接触）
    """
    n = distance_matrix.shape[0]
    contact_matrix = np.zeros((n, n), dtype=int)
    
    # 创建距离掩码
    distance_mask = distance_matrix < threshold
    
    # 创建对角线掩码（排除主对角线附近的残基对）
    i_indices, j_indices = np.meshgrid(range(n), range(n), indexing='ij')
    diag_mask = np.abs(i_indices - j_indices) > remove_diag
    
    # 组合掩码并应用到接触矩阵
    contact_matrix = (distance_mask & diag_mask).astype(int)
    
    return contact_matrix





def create_pdb_to_msa_mapping(pdb_residue_dict: Dict[int, str], msa_sequence: str, pdb_start_offset: int = 0) -> Dict[int, int]:
    """
    创建PDB残基ID到MSA位置的映射 
    带残基类型验证和起始偏移量，确保映射正确性
    
    Args:
        pdb_residue_dict: PDB残基字典 {残基ID: 残基类型}，如 {1: 'M', 2: 'K', 5: 'L', ...}
        msa_sequence: MSA中对应的序列（含gap，用-表示）
        pdb_start_offset: PDB序列的起始偏移量（默认0，表示从PDB第1个残基开始）
    
    Returns:
        {pdb_residue_id: msa_position} 映射字典
    """
    # 清理MSA序列中的空格和换行符
    clean_msa = msa_sequence.replace(' ', '').replace('\n', '').strip()

    # 按残基ID排序，获取有序的PDB序列
    sorted_pdb_ids = sorted(pdb_residue_dict.keys())
    pdb_sequence = ''.join(pdb_residue_dict[res_id] for res_id in sorted_pdb_ids)

    # 如果没有指定偏移量，尝试自动找到对齐起始点
    if pdb_start_offset == 0:
        pdb_start_offset = find_alignment_start(pdb_sequence, clean_msa)
        print(f"自动检测到PDB起始偏移量: {pdb_start_offset}")
    
    pdb_to_msa = {}
    pdb_seq_idx = pdb_start_offset  # PDB序列中的位置（考虑偏移量）
    matched_count = 0
    mismatches = []
    
    for msa_pos, msa_char in enumerate(clean_msa):
        if msa_char != '-':  # 非gap位置
            if pdb_seq_idx < len(sorted_pdb_ids):
                # 获取当前PDB残基ID和类型
                pdb_residue_id = sorted_pdb_ids[pdb_seq_idx]
                pdb_char = pdb_residue_dict[pdb_residue_id]
                
                # 验证残基类型是否匹配
                if pdb_char.upper() == msa_char.upper():
                    pdb_to_msa[pdb_residue_id] = msa_pos
                    matched_count += 1
                else:
                    # 记录不匹配的残基
                    mismatches.append((pdb_residue_id, pdb_char, msa_pos, msa_char))
                    if len(mismatches) <= 5:  # 只显示前5个不匹配
                        print(f"警告: 残基不匹配 - PDB残基{pdb_residue_id}({pdb_char}) vs MSA位置{msa_pos}({msa_char})")
                
                pdb_seq_idx += 1
            else:
                print(f"警告: MSA序列比PDB序列长，MSA位置{msa_pos}({msa_char})无对应PDB残基")
                break
                
    # 检查是否还有未映射的PDB残基
    if pdb_seq_idx < len(sorted_pdb_ids):
        unmapped_count = len(sorted_pdb_ids) - pdb_seq_idx
        print(f"警告: PDB序列比MSA序列长，有{unmapped_count}个PDB残基未映射")
        print(f"未映射的PDB残基ID: {sorted_pdb_ids[pdb_seq_idx:pdb_seq_idx+5]}{'...' if unmapped_count > 5 else ''}")
      
    print(f"PDB到MSA映射完成: {len(pdb_to_msa)} 个残基 (匹配: {matched_count}, 不匹配: {len(mismatches)})")
    print(f"PDB序列长度: {len(pdb_sequence)}, MSA序列长度: {len(clean_msa)}")

    
    return pdb_to_msa


def find_alignment_start(pdb_sequence: str, msa_sequence: str, min_match_length: int = 5) -> int:
    """
    自动寻找PDB序列在MSA序列中的对齐起始位置 - 细菌式基因
    
    Args:
        pdb_sequence: PDB序列
        msa_sequence: 清理后的MSA序列
        min_match_length: 最小连续匹配长度
    
    Returns:
        PDB序列的起始偏移量
    """
    # 移除MSA中的gap，得到纯氨基酸序列
    msa_no_gaps = msa_sequence.replace('-', '')
    
    best_offset = 0
    best_match_count = 0
    
    # 尝试不同的PDB起始位置
    for pdb_offset in range(len(pdb_sequence)):
        match_count = 0
        consecutive_matches = 0
        max_consecutive = 0
        
        # 比较从该偏移量开始的序列
        for i in range(min(len(pdb_sequence) - pdb_offset, len(msa_no_gaps))):
            if pdb_sequence[pdb_offset + i].upper() == msa_no_gaps[i].upper():
                match_count += 1
                consecutive_matches += 1
                max_consecutive = max(max_consecutive, consecutive_matches)
            else:
                consecutive_matches = 0
        
        # 如果找到足够长的连续匹配，且总匹配数更多
        if max_consecutive >= min_match_length and match_count > best_match_count:
            best_match_count = match_count
            best_offset = pdb_offset
    
    return best_offset


def convert_contacts_to_msa_indices(contact_pairs: List[Tuple[int, int]], 
                                   pdb_to_msa_map: Dict[int, int]) -> List[Tuple[int, int]]:
    """
    将PDB残基ID的接触对转换为MSA位置索引 - 细菌式基因
    
    Args:
        contact_pairs: PDB残基ID的接触对列表
        pdb_to_msa_map: PDB到MSA的映射字典
        
    Returns:
        MSA位置索引的接触对列表
    """
    msa_contacts = []
    skipped_count = 0
    
    for pdb_res1, pdb_res2 in contact_pairs:
        # 检查两个残基是否都能映射到MSA
        if pdb_res1 in pdb_to_msa_map and pdb_res2 in pdb_to_msa_map:
            msa_pos1 = pdb_to_msa_map[pdb_res1]
            msa_pos2 = pdb_to_msa_map[pdb_res2]
            msa_contacts.append((msa_pos1, msa_pos2))
        else:
            skipped_count += 1
    
    print(f"接触对转换完成: {len(msa_contacts)} 个有效接触对, {skipped_count} 个跳过")
    return msa_contacts


def calculate_contact_difference(chemokine_id: str, alternate_id: str, threshold: float = 8.0,remove_diag: int = 3) -> Dict:
    """
    计算两种蛋白质构象的接触差异分析 - 使用PDB残基ID
    
    Args:
        chemokine_id: 趋化因子折叠的PDB ID
        alternate_id: 替代折叠的PDB ID  
        threshold: 接触距离阈值（Å）
        
    Returns:
        结果字典包含：
        - diff_matrix: 接触差异矩阵 (Cdiff = Calternate - Cchemokine)
        - new_contacts: 新形成的接触 (+1)
        - lost_contacts: 断开的接触 (-1)
        - critical_contacts: 关键接触 (+1 or -1)
        - common_contacts: 保守接触 (0)
        - chemokine_matrix: 趋化因子接触矩阵
        - alternate_matrix: 替代折叠接触矩阵
        - chemokine_residue_types: 趋化因子残基类型信息
        - alternate_residue_types: 替代折叠残基类型信息
    """
    # 获取坐标、残基类型和距离矩阵
    chemokine_coords, chemokine_types = parse_pdb_coordinates(f"/yezhirui/evo_probe/data/{chemokine_id}.pdb", "A")
    alternate_coords, alternate_types = parse_pdb_coordinates(f"/yezhirui/evo_probe/data/{alternate_id}.pdb", "A")
    
    chemokine_dist_matrix, chemokine_residues = calculate_distance_matrix(chemokine_coords)
    alternate_dist_matrix, alternate_residues = calculate_distance_matrix(alternate_coords)
    
    # 创建二元接触矩阵
    chemokine_contacts = create_binary_contact_matrix(chemokine_dist_matrix, threshold, remove_diag)
    alternate_contacts = create_binary_contact_matrix(alternate_dist_matrix, threshold, remove_diag)
    
    # 对齐矩阵（根据残基ID）
    aligned_chemokine, aligned_alternate, common_residues = align_contact_matrices(
        chemokine_contacts, chemokine_residues,
        alternate_contacts, alternate_residues
    )
    
    # 计算接触差异矩阵 Cdiff = Calternate - Cchemokine
    diff_matrix = aligned_alternate - aligned_chemokine
    
    # 清零对角线附近的元素（排除序列相近的残基对）
    diff_matrix = clear_diagonal_band(diff_matrix, remove_diag)
    
    # 分析关键接触
    results = analyze_critical_contacts(diff_matrix, common_residues)
    
    # 添加完整矩阵到结果
    results.update({
        'diff_matrix': diff_matrix,
        'chemokine_matrix': aligned_chemokine, 
        'alternate_matrix': aligned_alternate,
        'residue_indices': common_residues,
        'chemokine_residue_types': chemokine_types,
        'alternate_residue_types': alternate_types,
        
    })
    
    return results


def calculate_contact_difference_msa_id(chemokine_id: str, alternate_id: str, 
                                         msa_seq: str, 
                                         threshold: float = 8.0, remove_diag: int = 3) -> Dict:
    """
    计算两种蛋白质构象的接触差异分析 - 使用MSA统一索引 - 细菌式基因
    
    Args:
        chemokine_id: 趋化因子折叠的PDB ID
        alternate_id: 替代折叠的PDB ID
        chemokine_msa_seq: 趋化因子在MSA中的序列
        alternate_msa_seq: 替代折叠在MSA中的序列
        threshold: 接触距离阈值（Å）
        remove_diag: 忽略对角线附近的残基对数量
        
    Returns:
        结果字典包含MSA索引的接触信息：
        - diff_matrix: 接触差异矩阵 (MSA索引)
        - new_contacts: 新形成的接触 (MSA索引)
        - lost_contacts: 断开的接触 (MSA索引)
        - critical_contacts: 关键接触 (MSA索引)
        - common_contacts: 保守接触 (MSA索引)
        - chemokine_matrix: 趋化因子接触矩阵 (MSA索引)
        - alternate_matrix: 替代折叠接触矩阵 (MSA索引)
        - msa_length: MSA序列长度
    """
    # 首先获取PDB的接触差异分析
    pdb_results = calculate_contact_difference(chemokine_id, alternate_id, threshold, remove_diag)
    
    # 提取PDB序列
    pdb_dict = pdb_results['alternate_residue_types']
    
    # 创建PDB到MSA的映射
    pdb_to_msa = create_pdb_to_msa_mapping(pdb_dict, msa_seq,pdb_start_offset=3)

    
    # 转换接触对到MSA索引
    new_contacts_msa = convert_contacts_to_msa_indices(pdb_results['new_contacts'], pdb_to_msa)
    lost_contacts_msa = convert_contacts_to_msa_indices(pdb_results['lost_contacts'], pdb_to_msa)
    
    
    # # 创建MSA索引的接触矩阵
    # msa_length = max(len(msa_seq.replace(' ', '').replace('\n', '').strip()),
    #                  len(msa_seq.replace(' ', '').replace('\n', '').strip()))
    
    # chemokine_matrix_msa = create_contact_matrix_from_pairs(common_contacts_msa + lost_contacts_msa, msa_length)
    # alternate_matrix_msa = create_contact_matrix_from_pairs(common_contacts_msa + new_contacts_msa, msa_length)
    # diff_matrix_msa = alternate_matrix_msa - chemokine_matrix_msa
    
    # # 清零对角线附近的元素
    # diff_matrix_msa = clear_diagonal_band(diff_matrix_msa, remove_diag)
    
    return {
        # 'diff_matrix': diff_matrix_msa,
        'new_contacts': new_contacts_msa,
        'lost_contacts': lost_contacts_msa,
        'critical_contacts': lost_contacts_msa + new_contacts_msa,
        "pdb_to_msa_map": pdb_to_msa,
        "alternate_residue_types": pdb_results['alternate_residue_types'],
        # 'chemokine_matrix': chemokine_matrix_msa,
        # 'alternate_matrix': alternate_matrix_msa,
        'summary': {
            'new_count': len(new_contacts_msa),
            'lost_count': len(lost_contacts_msa),
            'common_count': len(pdb_results['common_contacts']),
            'critical_count': len(lost_contacts_msa + new_contacts_msa)
        }
    }


def create_contact_matrix_from_pairs(contact_pairs: List[Tuple[int, int]], matrix_size: int) -> np.ndarray:
    """
    从接触对列表创建接触矩阵 - 细菌式基因
    
    Args:
        contact_pairs: 接触对列表
        matrix_size: 矩阵大小
        
    Returns:
        接触矩阵
    """
    matrix = np.zeros((matrix_size, matrix_size), dtype=int)
    
    for pos1, pos2 in contact_pairs:
        if 0 <= pos1 < matrix_size and 0 <= pos2 < matrix_size:
            matrix[pos1, pos2] = 1
            matrix[pos2, pos1] = 1  # 对称矩阵
    
    return matrix


def clear_diagonal_band(matrix: np.ndarray, remove_diag: int = 3) -> np.ndarray:
    """
    清零对角线附近的元素 - 细菌式基因，可重用
    
    Args:
        matrix: 输入矩阵
        remove_diag: 清零对角线附近的范围（i,i+remove_diag以内）
        
    Returns:
        cleared_matrix: 清零后的矩阵
    """
    cleared_matrix = matrix.copy()
    n = matrix.shape[0]
    
    # 清零对角线附近 remove_diag 范围内的所有元素
    for i in range(n):
        for j in range(max(0, i-remove_diag), min(n, i+remove_diag+1)):
            cleared_matrix[i, j] = 0
            
    return cleared_matrix


def align_contact_matrices(matrix1: np.ndarray, residues1: List[int], 
                          matrix2: np.ndarray, residues2: List[int]) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    """
    根据残基ID对齐两个接触矩阵 - 细菌式基因，可重用
    
    Args:
        matrix1, matrix2: 要对齐的接触矩阵
        residues1, residues2: 对应的残基ID列表
        
    Returns:
        (aligned_matrix1, aligned_matrix2, common_residues): 对齐后的矩阵和共同残基列表
    """
    # 找到共同残基
    common_residues = sorted(set(residues1) & set(residues2))
    n_common = len(common_residues)
    
    if n_common == 0:
        print("警告: 没有找到共同残基!")
        return np.array([]), np.array([]), []
    
    print(f"共同残基数量: {n_common} 个 (原始: {len(residues1)} vs {len(residues2)})")
    
    # 创建残基到索引的映射
    res_to_idx1 = {res: idx for idx, res in enumerate(residues1)}
    res_to_idx2 = {res: idx for idx, res in enumerate(residues2)}
    
    # 创建对齐后的矩阵
    aligned_matrix1 = np.zeros((n_common, n_common), dtype=matrix1.dtype)
    aligned_matrix2 = np.zeros((n_common, n_common), dtype=matrix2.dtype)
    
    # 填充对齐后的矩阵
    for i, res_i in enumerate(common_residues):
        for j, res_j in enumerate(common_residues):
            if res_i in res_to_idx1 and res_j in res_to_idx1:
                idx1_i, idx1_j = res_to_idx1[res_i], res_to_idx1[res_j]
                aligned_matrix1[i, j] = matrix1[idx1_i, idx1_j]
            
            if res_i in res_to_idx2 and res_j in res_to_idx2:
                idx2_i, idx2_j = res_to_idx2[res_i], res_to_idx2[res_j]
                aligned_matrix2[i, j] = matrix2[idx2_i, idx2_j]
    
    return aligned_matrix1, aligned_matrix2, common_residues


def analyze_critical_contacts(diff_matrix: np.ndarray, residue_indices: List[int]) -> Dict:
    """
    分析关键接触的细菌式函数 - 单一职责，可重用
    
    Args:
        diff_matrix: 接触差异矩阵
        residue_indices: 残基编号列表
        
    Returns:
        关键接触分析结果
    """
    # 获取上三角矩阵（避免重复计算对称部分）
    i_indices, j_indices = np.triu_indices_from(diff_matrix, k=1)
    
    # 分类接触变化
    new_contacts = []      # +1: 新形成的接触
    lost_contacts = []     # -1: 断开的接触  
    common_contacts = []   # 0: 保守接触
    
    for i, j in zip(i_indices, j_indices):
        value = diff_matrix[i, j]
        residue_pair = (residue_indices[i], residue_indices[j])
        
        if value == 1:
            new_contacts.append(residue_pair)
        elif value == -1:
            lost_contacts.append(residue_pair)
        else:  # value == 0
            common_contacts.append(residue_pair)
    
    return {
        'new_contacts': new_contacts,
        'lost_contacts': lost_contacts, 
        'common_contacts': common_contacts,
        'critical_contacts': lost_contacts + new_contacts,
        'summary': {
            'new_count': len(new_contacts),
            'lost_count': len(lost_contacts),
            'common_count': len(common_contacts),
            'critical_count': len(lost_contacts + new_contacts)
        }
    }



def parse_ancestral_probs(anc_prob_file: str, node_id: str = None) -> Dict[int, Dict[str, float]]:
    """
    解析祖先概率文件，获取每个位点的完整氨基酸概率分布 - 细菌式基因
    
    Args:
        anc_prob_file: 祖先概率文件路径
        node_id: 指定要解析的节点ID（如"499", "507"等），None表示使用第一个节点
        
    Returns:
        {site: {aa: prob, ...}} 字典 - site使用0-based索引，aa为20种氨基酸的概率分布
        例如: {0: {'A': 0.474, 'R': 0.001, ..., 'V': 0.023}, 1: {...}, ...}
    """
    site_distributions = {}
    current_node = None
    target_node_found = False
    parsing_target_node = False
    
    try:
        with open(anc_prob_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # 检查是否是节点标题行
                if line.startswith('Prob distribution at node'):
                    # 提取节点ID
                    node_match = re.search(r'node (\d+)', line)
                    if node_match:
                        current_node = node_match.group(1)  
                        
                        # 如果没指定node_id，使用第一个节点
                        if node_id is None:
                            parsing_target_node = True
                            target_node_found = True
                        # 如果找到了目标节点
                        elif current_node == str(node_id):
                            parsing_target_node = True
                            target_node_found = True
                        else:
                            parsing_target_node = False
                    continue
                
                # 如果不在目标节点中，跳过
                if not parsing_target_node:
                    continue
                    
                # 跳过头部和空行
                if not line or 'site' in line or 'Prob distribution' in line:
                    continue
                
                # 查找包含概率分布的行（含有冒号和氨基酸概率）
                if ':' not in line or 'A(' not in line:
                    continue
                    
                # 解析行：提取site编号和概率分布
                parts = line.split(':')
                if len(parts) != 2:
                    continue
                
                # 提取site编号（在行首）
                site_match = re.match(r'\s*(\d+)', line)
                if not site_match:
                    continue
                site_num = int(site_match.group(1))
                
                # 解析概率分布 "A(0.474) R(0.001) ..."
                prob_part = parts[1].strip()
                aa_probs = {}
                
                # 使用正则表达式提取氨基酸概率
                matches = re.findall(r'([ARNDCQEGHILKMFPSTWYV])\(([0-9.]+)\)', prob_part)
                for aa, prob_str in matches:
                    aa_probs[aa] = float(prob_str)
                
                # 存储完整的概率分布，转换为0-based索引
                if aa_probs:
                    site_distributions[site_num - 1] = aa_probs
        
        if node_id is not None and not target_node_found:
            print(f"警告: 未找到节点 {node_id}，可用节点: 499, 500, 501, 502, 507")
                    
    except FileNotFoundError:
        print(f"祖先概率文件不存在: {anc_prob_file}")
    except Exception as e:
        print(f"解析祖先概率文件时出错: {e}")
        
    return site_distributions


def get_ancestral_probs_max_aa(anc_prob_file: str, node_id: str = None) -> Dict[int, Tuple[str, float]]:
    """
    获取每个位点概率最大的氨基酸类型及其概率 - 复用parse_ancestral_probs
    
    Args:
        anc_prob_file: 祖先概率文件路径
        node_id: 指定要解析的节点ID（如"499", "507"等），None表示使用第一个节点
        
    Returns:
        {site: (max_aa, max_prob)} 字典 - 使用0-based索引
        例如: {0: ('A', 0.474), 1: ('K', 0.823), ...}
    """
    # 调用基础函数获取完整概率分布
    site_distributions = parse_ancestral_probs(anc_prob_file, node_id)
    
    # 提取每个位点的最大概率氨基酸
    site_max_probs = {}
    for site, aa_probs in site_distributions.items():
        if aa_probs:
            max_aa = max(aa_probs, key=aa_probs.get)
            max_prob = aa_probs[max_aa]
            site_max_probs[site] = (max_aa, max_prob)
    
    return site_max_probs


def get_contact_info(contact_list: List[Tuple[int, int]], 
                    contact_index: int, 
                    residue_types: Dict[int, str] = None,
                    node_id: str = None
                    ) -> str:
    """
    根据contact index获取残基信息字符串，支持祖先概率信息
    
    Args:
        contact_list: 接触对列表 [(res1, res2), ...]
        contact_index: 要查询的接触对索引
        residue_types: 残基类型字典 {残基编号: 单字母氨基酸代码}
        node_id: 祖先节点ID，如果提供则从anc_prob.txt读取概率信息
        
    Returns:
        描述字符串，如 "A15-V42" 或 "15A(0.3)-42V(0.4)" (带概率信息)
    """
    if contact_index >= len(contact_list):
        return f"索引 {contact_index} 超出范围 (最大: {len(contact_list)-1})"
    
    res1, res2 = contact_list[contact_index]
    
    # 如果提供了node_id，尝试读取祖先概率信息
    if node_id is not None:
        anc_prob_file = f"/yezhirui/evo_probe/data/anc_prob.txt"  # 假设文件名格式
        site_probs = get_ancestral_probs_max_aa(anc_prob_file,node_id)
        
        if res1 in site_probs and res2 in site_probs:
            aa1, prob1 = site_probs[res1]
            aa2, prob2 = site_probs[res2]
            return f"{res1}{aa1}({prob1:.1f})-{res2}{aa2}({prob2:.1f})"
    
    # 使用传统的残基类型信息
    if residue_types is not None:
        type1 = residue_types.get(res1, 'X')  # 未知类型用X表示
        type2 = residue_types.get(res2, 'X')
        return f"{type1}{res1}-{type2}{res2}"
    else:
        return f"{res1}-{res2}"

def get_contact_info_with_pdb(contact_list, contact_index, msa_to_pdb_map=None, 
                             residue_types=None, node_id=None):
    """
    获取残基信息字符串，同时显示MSA位置和对应的PDB ID
    
    Args:
        contact_list: 接触对列表 [(res1, res2), ...]，这里是MSA位置
        contact_index: 要查询的接触对索引
        msa_to_pdb_map: MSA到PDB的映射字典 {msa_id: pdb_id}
        residue_types: 残基类型字典
        node_id: 祖先节点ID
        
    Returns:
        描述字符串，如 "0(4E)-40(46K)" 表示 "MSA位置(PDB位置+残基类型)"
    """
    if contact_index >= len(contact_list):
        return f"索引 {contact_index} 超出范围 (最大: {len(contact_list)-1})"
    
    msa_res1, msa_res2 = contact_list[contact_index]
    
    # 获取对应的PDB ID
    if msa_to_pdb_map is not None:
        pdb_res1 = msa_to_pdb_map.get(msa_res1, None)
        pdb_res2 = msa_to_pdb_map.get(msa_res2, None)
        
        # 如果提供了node_id，尝试读取祖先概率信息
        if node_id is not None:
            anc_prob_file = f"/yezhirui/evo_probe/data/anc_prob.txt"
            site_probs = get_ancestral_probs_max_aa(anc_prob_file,node_id)
            
            if msa_res1 in site_probs and msa_res2 in site_probs:
                aa1, prob1 = site_probs[msa_res1]
                aa2, prob2 = site_probs[msa_res2]
                
                # 格式: MSA位置(PDB位置+残基类型)
                pdb_part1 = f"{pdb_res1}{aa1}" if pdb_res1 is not None else f"?{aa1}"
                pdb_part2 = f"{pdb_res2}{aa2}" if pdb_res2 is not None else f"?{aa2}"
                return f"{msa_res1}({pdb_part1})-{msa_res2}({pdb_part2})"
        
        # 使用残基类型信息
        if residue_types is not None and pdb_res1 is not None and pdb_res2 is not None:
            type1 = residue_types.get(pdb_res1, 'X')
            type2 = residue_types.get(pdb_res2, 'X')
            return f"{msa_res1}({pdb_res1}{type1})-{msa_res2}({pdb_res2}{type2})"
        
        # 只有PDB ID，没有残基类型
        if pdb_res1 is not None and pdb_res2 is not None:
            return f"{msa_res1}({pdb_res1})-{msa_res2}({pdb_res2})"
    
    # 只显示MSA位置
    return f"{msa_res1}-{msa_res2}"
    

def plot_contacts(seq_len,contact_list, title, contact_list2=None):
    import matplotlib.pyplot as plt

    print(f"num of contacts: {len(contact_list)}")
    # 提取坐标
    res_i = [contact[0] for contact in contact_list]
    res_j = [contact[1] for contact in contact_list]

    plt.figure(figsize=(10, 10))
    plt.scatter(res_j, res_i, s=15, c='blue', marker='o', label=f'Contact 1: {len(contact_list)}') # 注意x,y的对应关系
    
    if contact_list2 is not None:
        print(f"num of contacts 2: {len(contact_list2)}")
        # 提取第二个接触列表的坐标
        res_i2 = [contact[0] for contact in contact_list2]
        res_j2 = [contact[1] for contact in contact_list2]
        plt.scatter(res_j2, res_i2, s=15, c='red', marker='o', label=f'Contact 2: {len(contact_list2)}')

    # 设置坐标轴范围
    plt.xlim(0, seq_len)
    plt.ylim(0, seq_len)
    
    # 翻转y轴方向，使0点从左上角开始
    plt.gca().invert_yaxis()

    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.title(title)
    plt.legend()
    plt.savefig(f"/yezhirui/evo_probe/figs/{title}.png")





if __name__ == "__main__":
    # 示例：分析趋化因子和替代折叠的接触差异
    chemokine_pdb = "1j8i"  # 趋化因子折叠
    alternate_pdb = "2jp1"  # 替代折叠
    
    
    
    print("\n=== 完整分析版本 ===")
    msa_seq = 'EVSDKRT-CVSLTTQRLPVSRIKTYTIT---EGSLRAVIFITKRGLKVCADPQATWVRDVVRSMDRKSNT'
    results = calculate_contact_difference_msa_id(chemokine_pdb, alternate_pdb, msa_seq, threshold=10.0,remove_diag=5)
    
    print(f"\n接触变化统计:")
    print(f"新形成的接触: {results['summary']['new_count']} 个")
    print(f"断开的接触: {results['summary']['lost_count']} 个") 
    print(f"保守接触: {results['summary']['common_count']} 个")
    
    print(f"\n前5个新形成的接触 (+1):")
    for res1, res2 in results['new_contacts'][:5]:
        print(f"  残基 {res1} - {res2}")
        
    print(f"\n前5个断开的接触 (-1):")
    for res1, res2 in results['lost_contacts'][:5]:
        print(f"  残基 {res1} - {res2}")
        
    # print(f"\n差异矩阵形状: {results['diff_matrix'].shape}")
    # unique, counts = np.unique(results['diff_matrix'], return_counts=True)
    # print("矩阵元素统计 (包含对称部分):")
    # for val, count in zip(unique, counts):
    #     if val == 1:
    #         print(f"  +1 (新接触): {count} 个矩阵元素 = {count//2} 个接触对")
    #     elif val == -1:
    #         print(f"  -1 (断开接触): {count} 个矩阵元素 = {count//2} 个接触对") 
    #     else:
    #         print(f"   0 (无变化): {count} 个矩阵元素")



