import numpy as np
import sys
sys.path.append("/yezhirui/evo_probe")
from src.find_contact import get_critical_contacts
from src.gmm import *
from src.contact_space import *
from Bio import SeqIO
from sklearn.preprocessing import StandardScaler
from typing import List, Tuple, Dict
from scipy import stats

def generate_orthogonal_vectors(v_target, num_vectors=3, random_state=42):
    """
    生成与目标向量正交的向量
    
    Args:
        v_target: 目标向量 (244,)
        num_vectors: 需要生成的正交向量数量
        random_state: 随机种子
    
    Returns:
        orthogonal_vectors: 正交向量列表，每个向量都与v_target正交
    """
    np.random.seed(random_state)
    
    # 标准化目标向量
    v_target_norm = v_target / np.linalg.norm(v_target)
    
    orthogonal_vectors = []
    
    for i in range(num_vectors):
        # 生成随机向量
        random_vec = np.random.randn(len(v_target))
        
        # 从随机向量中减去在v_target上的投影，使其正交
        projection = np.dot(random_vec, v_target_norm) * v_target_norm
        orthogonal_vec = random_vec - projection
        
        # 对已有的正交向量进行正交化（Gram-Schmidt过程）
        for existing_vec in orthogonal_vectors:
            projection = np.dot(orthogonal_vec, existing_vec) * existing_vec
            orthogonal_vec -= projection
        
        # 标准化
        orthogonal_vec = orthogonal_vec / np.linalg.norm(orthogonal_vec)
        orthogonal_vectors.append(orthogonal_vec)
    
    return orthogonal_vectors


def generate_single_mutants(seq):
    """生成单点突变体"""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    mutants = []
    
    for i, aa in enumerate(seq):
        for new_aa in amino_acids:
            if new_aa != aa:
                mutant = seq[:i] + new_aa + seq[i+1:]
                mutants.append(mutant)
    
    return mutants


def embed_and_project_mutants(mutants: List[str], 
                                     target_contacts: List[Tuple[int, int]],
                                     mj_dict: Dict,
                                     scaler: StandardScaler, # 传入训练好的scaler
                                     v_target: np.ndarray,
                                     v_control: np.ndarray,
                                     origin_point: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    
    d = len(target_contacts)
    m = len(mutants)
    x_raw = np.zeros((m, d)) # (m, d) 的二维数组，m是序列数量

    # 1. 构建原始embedding矩阵 (处理多个序列)
    for seq_idx, sequence in enumerate(mutants):
        for contact_idx, (i, j) in enumerate(target_contacts):
            if i < len(sequence) and j < len(sequence):
                res_i, res_j = sequence[i], sequence[j]
                if res_i != '-' and res_j != '-':
                    mj_energy = mj_dict.get((res_i.upper(), res_j.upper()), 0)
                    x_raw[seq_idx, contact_idx] = mj_energy

    # 2. 使用"已训练"的scaler进行标准化 (最关键的一步)
    # 这里用 .transform() 而不是 .fit_transform()!
    # 这保证了新序列是在旧坐标系下进行标准化
    x_scaled = scaler.transform(x_raw)

    # 3. 计算相对位移的投影
    # origin_point 是我们关心的起点，比如Anc.2最优序列的坐标
    displacement_vectors = x_scaled - origin_point  # 广播，每行都减去origin_point
    
    projection_target = np.dot(displacement_vectors, v_target)  # (m,) 向量
    projection_control = np.dot(displacement_vectors, v_control)  # (m,) 向量
    
    return projection_target, projection_control


def analyze_mutant_projections(anc2_label1_sequences: List[str],
                              critical_contacts: List[Tuple[int, int]],
                              mj_dict: Dict,
                              anc2_scaler: StandardScaler,
                              anc2_contact_embedding: np.ndarray,
                              v_target: np.ndarray,
                              v_controls: List[np.ndarray],
                              anc2_results: Dict) -> Dict:
    """
    分析anc2_label1_sequences的突变体在目标方向和控制方向上的投影
    """
    # 获取label=1的序列在anc2_contact_embedding中的坐标
    anc2_labels = anc2_results['labels']
    label_1_mask = anc2_labels == 1
    label_1_embeddings = anc2_contact_embedding[label_1_mask]
    
    all_results = []
    
    for seq_idx, original_seq in enumerate(anc2_label1_sequences):
        print(f"处理序列 {seq_idx + 1}/{len(anc2_label1_sequences)}")
        
        # 生成该序列的突变体
        mutants = generate_single_mutants(original_seq)
        print(f"  生成 {len(mutants)} 个突变体")
        
        # 获取原始序列的坐标作为origin_point
        origin_point = label_1_embeddings[seq_idx]  # 使用对应的embedding
        
        # 对每个控制向量分析
        seq_results = {
            'original_seq': original_seq,
            'origin_point': origin_point,
            'projections_target': None,
            'projections_controls': [],
            'statistics': {}
        }
        
        # 计算在v_target上的投影
        proj_target, _ = embed_and_project_mutants(
            mutants, critical_contacts, mj_dict, anc2_scaler,
            v_target, v_target, origin_point  # 这里v_control参数不会被使用
        )
        seq_results['projections_target'] = proj_target
        
        # 计算在每个v_control上的投影
        for i, v_control in enumerate(v_controls):
            _, proj_control = embed_and_project_mutants(
                mutants, critical_contacts, mj_dict, anc2_scaler,
                v_target, v_control, origin_point  # 这里v_target参数不会被使用
            )
            seq_results['projections_controls'].append(proj_control)
        
        all_results.append(seq_results)
    
    return all_results


def statistical_analysis(all_results: List[Dict]) -> Dict:
    """
    对所有序列的突变体投影进行统计分析
    """
    # 收集所有投影数据
    all_target_projections = []
    all_control_projections = [[] for _ in range(3)]  # 3个控制方向
    
    for seq_result in all_results:
        all_target_projections.extend(seq_result['projections_target'])
        for i in range(3):
            all_control_projections[i].extend(seq_result['projections_controls'][i])
    
    all_target_projections = np.array(all_target_projections)
    all_control_projections = [np.array(proj) for proj in all_control_projections]
    
    print(f"\n=== 统计分析结果 ===")
    print(f"总突变体数量: {len(all_target_projections)}")
    
    # 计算基本统计量
    target_mean = np.mean(all_target_projections)
    target_std = np.std(all_target_projections)
    target_abs_mean = np.mean(np.abs(all_target_projections))
    
    print(f"\nv_target方向:")
    print(f"  平均位移: {target_mean:.4f}")
    print(f"  标准差: {target_std:.4f}")
    print(f"  绝对位移平均: {target_abs_mean:.4f}")
    
    control_stats = []
    for i, control_proj in enumerate(all_control_projections):
        control_mean = np.mean(control_proj)
        control_std = np.std(control_proj)
        control_abs_mean = np.mean(np.abs(control_proj))
        
        print(f"\nv_control_{i}方向:")
        print(f"  平均位移: {control_mean:.4f}")
        print(f"  标准差: {control_std:.4f}")
        print(f"  绝对位移平均: {control_abs_mean:.4f}")
        
        control_stats.append({
            'mean': control_mean,
            'std': control_std,
            'abs_mean': control_abs_mean
        })
    
    # 统计显著性检验
    print(f"\n=== 显著性检验 ===")
    
    # 比较绝对位移
    target_abs = np.abs(all_target_projections)
    p_values = []
    
    for i, control_proj in enumerate(all_control_projections):
        control_abs = np.abs(control_proj)
        
        # Wilcoxon秩和检验（Mann-Whitney U检验）
        statistic, p_value = stats.mannwhitneyu(
            target_abs, control_abs, alternative='greater'
        )
        p_values.append(p_value)
        
        print(f"v_target vs v_control_{i}:")
        print(f"  统计量: {statistic}")
        print(f"  p值: {p_value:.2e}")
        print(f"  显著性: {'是' if p_value < 0.05 else '否'} (α=0.05)")
    
    # 计算效应量 (Cohen's d)
    print(f"\n=== 效应量分析 ===")
    effect_sizes = []
    for i, control_proj in enumerate(all_control_projections):
        control_abs = np.abs(all_control_projections[i])
        target_abs = np.abs(all_target_projections)
        
        pooled_std = np.sqrt(((len(target_abs) - 1) * np.var(target_abs) + 
                             (len(control_abs) - 1) * np.var(control_abs)) / 
                            (len(target_abs) + len(control_abs) - 2))
        
        cohens_d = (np.mean(target_abs) - np.mean(control_abs)) / pooled_std
        effect_sizes.append(cohens_d)
        
        print(f"v_target vs v_control_{i}: Cohen's d = {cohens_d:.4f}")
    
    return {
        'target_projections': all_target_projections,
        'control_projections': all_control_projections,
        'target_stats': {'mean': target_mean, 'std': target_std, 'abs_mean': target_abs_mean},
        'control_stats': control_stats,
        'p_values': p_values,
        'effect_sizes': effect_sizes
    }


def get_sequences_by_label(gmm_labels, label_id, sequences):
    """根据GMM标签筛选序列"""
    mask = gmm_labels == label_id
    filtered_sequences = [sequences[i] for i in range(len(sequences)) if mask[i]]
    return filtered_sequences


if __name__ == "__main__":
    print("=== 初始化数据和模型 ===")
    
    # 1. 获取关键接触点
    chemokine_pdb = "1j8i"  # 趋化因子折叠
    alternate_pdb = "2jp1"  # 替代折叠
    new_contacts, lost_contacts = get_critical_contacts(chemokine_pdb, alternate_pdb,threshold=10.0,remove_diag=5)
    critical_contacts = new_contacts + lost_contacts
    print(f"关键接触点数量: {len(critical_contacts)}")

    # 2. 加载MJ矩阵
    mj_dict = load_mj_matrix("/yezhirui/evo_probe/data/mj_matrix.txt")
    
    # 3. 生成contact embedding
    anc2_seqs_path = f"/yezhirui/evo_probe/data/sample/node501_anc2_samples.fasta"
    anc2_contact_embedding, anc2_scaler = generate_contact_embedding(anc2_seqs_path, critical_contacts, mj_dict)

    anc3_seqs_path = f"/yezhirui/evo_probe/data/sample/node502_anc3_samples.fasta"
    anc3_contact_embedding, anc3_scaler = generate_contact_embedding(anc3_seqs_path, critical_contacts, mj_dict)

    anc4_seqs_path = f"/yezhirui/evo_probe/data/sample/node507_anc4_samples.fasta"
    anc4_contact_embedding, anc4_scaler = generate_contact_embedding(anc4_seqs_path, critical_contacts, mj_dict)

    print(f"ANC2 embedding shape: {anc2_contact_embedding.shape}")
    print(f"ANC3 embedding shape: {anc3_contact_embedding.shape}")
    print(f"ANC4 embedding shape: {anc4_contact_embedding.shape}")

    # 4. 进行GMM聚类分析
    datasets = {
        'ANC2': anc2_contact_embedding,
        'ANC3': anc3_contact_embedding,
        'ANC4': anc4_contact_embedding
    }

    all_results = compare_gmm_guided_projections(
        datasets,
        k=2,
    )

    anc2_results = all_results['ANC2']
    anc3_results = all_results['ANC3']
    anc4_results = all_results['ANC4']

    # 5. 定义投影方向
    v_target = anc3_results['projection_vector']
    v_controls = generate_orthogonal_vectors(v_target, num_vectors=3)
    
    print(f"v_target shape: {v_target.shape}")
    print(f"v_controls数量: {len(v_controls)}")

    # 6. 获取anc2_label1_sequences
    anc2_sequences = []
    with open(anc2_seqs_path, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            anc2_sequences.append(str(record.seq))

    anc2_label1_sequences = get_sequences_by_label(anc2_results['labels'], 1, anc2_sequences)
    print(f"ANC2标签1的序列数量: {len(anc2_label1_sequences)}")

    print("\n=== 开始突变体分析 ===")
    
    # 7. 分析突变体投影
    mutant_results = analyze_mutant_projections(
        anc2_label1_sequences,
        critical_contacts,
        mj_dict,
        anc2_scaler,
        anc2_contact_embedding,
        v_target,
        v_controls,
        anc2_results
    )
    
    # 8. 统计分析
    stats_results = statistical_analysis(mutant_results)
    
    print("\n=== 分析完成 ===")
    print("结果总结:")
    print(f"- 处理了 {len(anc2_label1_sequences)} 个原始序列")
    print(f"- 生成了 {len(stats_results['target_projections'])} 个突变体")
    print(f"- v_target方向平均绝对位移: {stats_results['target_stats']['abs_mean']:.4f}")
    
    significant_controls = sum(1 for p in stats_results['p_values'] if p < 0.05)
    print(f"- 显著高于控制方向的数量: {significant_controls}/3")
    
    if significant_controls == 3:
        print("✓ 结论：突变体在v_target方向的位移显著高于所有控制方向")
    elif significant_controls > 0:
        print(f"✓ 结论：突变体在v_target方向的位移显著高于{significant_controls}个控制方向")
    else:
        print("✗ 结论：突变体在v_target方向的位移未显著高于控制方向")



