import numpy as np
from typing import List, Tuple, Dict
import os
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
from src.find_contact import get_critical_contacts

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


def auto_dbscan_clustering(X, min_eps=5, max_eps=20, eps_step=1, min_samples=3, 
                          test_fraction=0.25, verbose=True):
    """
    自动选择最优eps参数的DBSCAN聚类
    
    Args:
        X: contact embedding矩阵 (n_samples, n_features)
        min_eps, max_eps, eps_step: eps扫描范围和步长
        min_samples: DBSCAN最小样本数
        test_fraction: 用于参数选择的数据比例
        verbose: 是否打印信息
        
    Returns:
        (clustering, eps_to_select, scan_results): 聚类结果、选择的eps值、扫描结果
    """
    
    eps_test_vals = np.arange(min_eps, max_eps + eps_step, eps_step)
    n_clusters_list = []
    scan_results = []
    
    if verbose:
        print("eps\tn_clusters\tn_noise\tn_samples")
        print("-" * 40)
    
    # eps参数扫描
    for eps in eps_test_vals:
        # 用部分数据测试，提高效率
        n_test = int(len(X) * test_fraction)
        test_indices = np.random.choice(len(X), n_test, replace=False)
        testset = X[test_indices]
        
        clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(testset)
        
        n_clust = len(set(clustering.labels_)) - (1 if -1 in clustering.labels_ else 0)
        n_noise = np.sum(clustering.labels_ == -1)
        
        n_clusters_list.append(n_clust)
        scan_results.append({
            'eps': eps,
            'n_clusters': n_clust,
            'n_noise': n_noise,
            'n_samples': len(testset) - n_noise
        })
        
        if verbose:
            print(f"{eps:.2f}\t{n_clust}\t\t{n_noise}\t\t{len(testset) - n_noise}")
        
        # 提前停止：如果eps很大但只有1个聚类，说明eps过大
        if eps > (max_eps * 0.5) and n_clust <= 1:
            break
    
    # 选择产生最多聚类数的eps值
    eps_to_select = eps_test_vals[np.argmax(n_clusters_list)]
    
    # 用选定的eps在全部数据上进行聚类
    clustering = DBSCAN(eps=eps_to_select, min_samples=min_samples).fit(X)
    
    if verbose:
        print(f"\n选择的eps: {eps_to_select:.2f}")
        n_final_clusters = len(set(clustering.labels_)) - (1 if -1 in clustering.labels_ else 0)
        n_final_noise = np.sum(clustering.labels_ == -1)
        print(f"最终结果: {n_final_clusters}个聚类, {n_final_noise}个噪声点")
    
    return clustering, eps_to_select, scan_results


def write_clusters_to_a3m(msa_path: str, 
                         cluster_labels: np.ndarray,
                         output_dir: str,
                         prefix: str = "cluster") -> List[str]:
    """
    将聚类结果写入多个a3m文件，每个聚类一个文件
    
    Args:
        msa_path: 原始MSA文件路径
        cluster_labels: 聚类标签数组 (包括噪声点-1)
        output_dir: 输出目录
        prefix: 输出文件前缀
        
    Returns:
        生成的a3m文件路径列表
    """
    from Bio import SeqIO
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取所有序列
    sequences = []
    seq_ids = []
    with open(msa_path, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences.append(str(record.seq))
            seq_ids.append(record.id)
    
    # 获取参考序列（第一条序列）
    reference_seq = sequences[0]
    reference_id = seq_ids[0]
    
    # 获取唯一的聚类标签（排除噪声点-1）
    unique_clusters = np.unique(cluster_labels)
    unique_clusters = unique_clusters[unique_clusters != -1]  # 排除噪声点
    
    output_files = []
    
    print(f"发现 {len(unique_clusters)} 个聚类")
    print(f"噪声点数量: {np.sum(cluster_labels == -1)}")
    
    for cluster_id in unique_clusters:
        # 获取该聚类的序列索引
        cluster_indices = np.where(cluster_labels == cluster_id)[0]
        cluster_size = len(cluster_indices)
        
        # 构建输出文件名
        output_file = os.path.join(output_dir, f"{prefix}_{cluster_id}.a3m")
        output_files.append(output_file)
        
        print(f"写入聚类 {cluster_id}: {cluster_size} 条序列 -> {output_file}")
        
        # 写入a3m文件
        with open(output_file, 'w') as f:
            # 首先写入参考序列
            f.write(f">{reference_id}_reference\n")
            f.write(f"{reference_seq}\n")
            
            # 写入该聚类的序列
            for idx in cluster_indices:
                # 跳过参考序列（避免重复）
                if idx == 0:
                    continue
                    
                seq_id = seq_ids[idx] if idx < len(seq_ids) else f"seq_{idx}"
                sequence = sequences[idx]
                
                f.write(f">{seq_id}_cluster_{cluster_id}\n")
                f.write(f"{sequence}\n")
    

    
    return output_files




if __name__ == "__main__":

    chemokine_pdb = "1j8i"  # 趋化因子折叠
    alternate_pdb = "2jp1"  # 替代折叠
    new_contacts, lost_contacts = get_critical_contacts(chemokine_pdb, alternate_pdb,threshold=10.0,remove_diag=5)
    critical_contacts = new_contacts + lost_contacts
    print(f"新形成接触: {len(new_contacts)} 个")
    print(f"断开接触: {len(lost_contacts)} 个")


    seqs_path = f"/yezhirui/evo_probe/data/sample/node499_anc0_samples.fasta"
    mj_dict = load_mj_matrix("/yezhirui/evo_multiconf/data/mj_matrix.txt")
    contact_embedding = generate_contact_embedding(seqs_path, critical_contacts , mj_dict)
    print(f"contact_embedding shape: {contact_embedding.shape}")


    clustering, best_eps, scan_results = auto_dbscan_clustering(contact_embedding, min_eps=5, max_eps=20, eps_step=1, min_samples=3, test_fraction=1, verbose=True)
    
    # 写入聚类结果到a3m文件
    output_dir = f"/yezhirui/evo_multiconf/results/pdb_dbscan_clusters/{uniprot_id}"
    os.makedirs(output_dir, exist_ok=True)
    output_files = write_clusters_to_a3m(msa_path, clustering.labels_, output_dir)
    
