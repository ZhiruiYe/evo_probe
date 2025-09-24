import numpy as np
from typing import List, Tuple, Dict, Optional

def parse_pdb_coordinates(pdb_path: str, 
                         chain_id: Optional[str] = None) -> Dict[int, np.ndarray]:
    """
    从PDB文件中提取残基的Cβ原子坐标（GLY使用Cα）
    
    Args:
        pdb_path: PDB文件路径
        chain_id: 指定链ID，None表示使用第一条链
        
    Returns:
        coordinates: {残基编号: 坐标数组} 字典
    """
    coordinates = {}
    
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"PDB文件不存在: {pdb_path}")
        return {}
    
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
    
    print(f"PDB坐标提取完成: {len(coordinates)} 个残基 (链 {target_chain})")
    return coordinates

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


def filter_dist_matrix(matrix, threshold=10, remove_diag=3):
    """保持向后兼容性的旧函数"""
    contact_matrix = create_binary_contact_matrix(matrix, threshold, remove_diag)
    i_indices, j_indices = np.where(np.triu(contact_matrix, k=1))
    return list(zip(i_indices, j_indices))

def calculate_contact_difference(chemokine_id: str, alternate_id: str, threshold: float = 8.0,remove_diag: int = 3) -> Dict:
    """
    计算两种蛋白质构象的接触差异分析
    
    Args:
        chemokine_id: 趋化因子折叠的PDB ID
        alternate_id: 替代折叠的PDB ID  
        threshold: 接触距离阈值（Å）
        
    Returns:
        结果字典包含：
        - diff_matrix: 接触差异矩阵 (Cdiff = Calternate - Cchemokine)
        - new_contacts: 新形成的接触 (+1)
        - lost_contacts: 断开的接触 (-1)
        - common_contacts: 保守接触 (0)
        - chemokine_matrix: 趋化因子接触矩阵
        - alternate_matrix: 替代折叠接触矩阵
    """
    # 获取坐标和距离矩阵
    chemokine_coords = parse_pdb_coordinates(f"/yezhirui/evo_probe/data/{chemokine_id}.pdb", "A")
    alternate_coords = parse_pdb_coordinates(f"/yezhirui/evo_probe/data/{alternate_id}.pdb", "A")
    
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
        'residue_indices': common_residues
    })
    
    return results


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
        'summary': {
            'new_count': len(new_contacts),
            'lost_count': len(lost_contacts),
            'common_count': len(common_contacts)
        }
    }


def get_critical_contacts(chemokine_id: str, alternate_id: str, threshold: float = 8.0,remove_diag: int = 3) -> Tuple[List, List]:
    """
    简洁版本：仅返回关键接触 - 完美的细菌式基因
    可以轻松"yoink"到任何项目中
    
    Args:
        chemokine_id: 趋化因子PDB ID
        alternate_id: 替代折叠PDB ID  
        threshold: 接触阈值（Å）
        
    Returns:
        (new_contacts, lost_contacts): 新形成和断开的接触对列表
    """
    results = calculate_contact_difference(chemokine_id, alternate_id, threshold,remove_diag)
    return results['new_contacts'], results['lost_contacts']


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
    
    print("=== 简洁版本：仅获取关键接触 ===")
    new_contacts, lost_contacts = get_critical_contacts(chemokine_pdb, alternate_pdb,threshold=10.0,remove_diag=5)
    print(f"新形成接触: {len(new_contacts)} 个")
    print(f"断开接触: {len(lost_contacts)} 个")
    
    print("\n=== 完整分析版本 ===")
    results = calculate_contact_difference(chemokine_pdb, alternate_pdb, threshold=10.0,remove_diag=5)
    
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
        
    print(f"\n差异矩阵形状: {results['diff_matrix'].shape}")
    unique, counts = np.unique(results['diff_matrix'], return_counts=True)
    print("矩阵元素统计 (包含对称部分):")
    for val, count in zip(unique, counts):
        if val == 1:
            print(f"  +1 (新接触): {count} 个矩阵元素 = {count//2} 个接触对")
        elif val == -1:
            print(f"  -1 (断开接触): {count} 个矩阵元素 = {count//2} 个接触对") 
        else:
            print(f"   0 (无变化): {count} 个矩阵元素")



    seq_len = 65
    
    plot_contacts(seq_len, new_contacts, f"{chemokine_pdb}_{alternate_pdb}_new_contacts")
    plot_contacts(seq_len, lost_contacts, f"{chemokine_pdb}_{alternate_pdb}_lost_contacts")
    plot_contacts(seq_len, new_contacts, f"{chemokine_pdb}_{alternate_pdb}_new_lost_contacts", lost_contacts)

    
