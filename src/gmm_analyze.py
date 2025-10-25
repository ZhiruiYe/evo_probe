
from typing import List, Tuple, Dict
import numpy as np


def analyze_centroid_evolution(gmm_results: Dict[str, Dict]) -> Dict[str, any]:
    """
    分析GMM质心的演化轨迹
    
    Args:
        gmm_results: 包含各ANC的GMM结果字典
        
    Returns:
        Dict: 质心分析结果
    """
    # 提取质心
    centroids = {}
    for anc_name, result in gmm_results.items():
        gmm_model = result[0]  # (gmm_model, labels, probabilities)
        centroids[anc_name] = gmm_model.means_[0]  # k=1时只有一个质心
    
    # 计算质心间距离
    anc_names = sorted(centroids.keys())  # 确保顺序: ANC0, ANC1, ANC2, ANC3, ANC4
    distances = {}
    
    for i in range(len(anc_names) - 1):
        curr_anc = anc_names[i]
        next_anc = anc_names[i + 1]
        
        distance = np.linalg.norm(centroids[curr_anc] - centroids[next_anc])
        distances[f"{curr_anc}->{next_anc}"] = distance
    
    # 计算累积距离
    cumulative_distance = 0
    cumulative_distances = {}
    
    for i, anc_name in enumerate(anc_names):
        if i == 0:
            cumulative_distances[anc_name] = 0
        else:
            prev_anc = anc_names[i-1]
            step_distance = distances[f"{prev_anc}->{anc_name}"]
            cumulative_distance += step_distance
            cumulative_distances[anc_name] = cumulative_distance
    
    return {
        'centroids': centroids,
        'step_distances': distances,
        'cumulative_distances': cumulative_distances,
        'anc_order': anc_names,
        'total_distance': cumulative_distance
    }


def plot_centroid_evolution(evolution_results: Dict[str, any], figsize: Tuple[int, int] = (15, 5)):
    """
    可视化质心演化轨迹
    """
    import matplotlib.pyplot as plt
    
    anc_names = evolution_results['anc_order']
    step_distances = list(evolution_results['step_distances'].values())
    cumulative_distances = [evolution_results['cumulative_distances'][anc] for anc in anc_names]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # 子图1: 步进距离
    transitions = list(evolution_results['step_distances'].keys())
    ax1.bar(range(len(step_distances)), step_distances, alpha=0.7, color='skyblue')
    ax1.set_xlabel('evolutionary transition')
    ax1.set_ylabel('centroid distance')
    ax1.set_title('step-wise centroid changes')
    ax1.set_xticks(range(len(transitions)))
    ax1.set_xticklabels(transitions, rotation=45)
    ax1.grid(True, alpha=0.3)
    
    # 子图2: 累积距离
    ax2.plot(range(len(anc_names)), cumulative_distances, 'o-', linewidth=2, markersize=8, color='orange')
    ax2.set_xlabel('ancestral node')
    ax2.set_ylabel('cumulative distance from ANC0')
    ax2.set_title('cumulative centroid drift')
    ax2.set_xticks(range(len(anc_names)))
    ax2.set_xticklabels(anc_names)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # 打印数值结果
    print("=== 质心演化分析 ===")
    print(f"总演化距离: {evolution_results['total_distance']:.4f}")
    print("\n步进距离:")
    for transition, distance in evolution_results['step_distances'].items():
        print(f"  {transition}: {distance:.4f}")
    print("\n累积距离:")
    for anc, distance in evolution_results['cumulative_distances'].items():
        print(f"  {anc}: {distance:.4f}")


def compare_centroid_dimensions(gmm_results: Dict[str, Dict], top_n: int = 10) -> Dict[str, any]:
    """
    比较质心在各维度上的变化
    
    Args:
        gmm_results: GMM结果字典
        top_n: 显示变化最大的前N个维度
        
    Returns:
        Dict: 维度变化分析结果
    """
    # 提取质心
    centroids = {}
    for anc_name, result in gmm_results.items():
        gmm_model = result[0]
        centroids[anc_name] = gmm_model.means_[0]
    
    anc_names = sorted(centroids.keys())
    n_features = len(centroids[anc_names[0]])
    
    # 计算每个维度的总变化量
    dimension_changes = np.zeros(n_features)
    
    for i in range(len(anc_names) - 1):
        curr_centroid = centroids[anc_names[i]]
        next_centroid = centroids[anc_names[i + 1]]
        
        # 累积每个维度的绝对变化
        dimension_changes += np.abs(next_centroid - curr_centroid)
    
    # 找出变化最大的维度
    top_indices = np.argsort(dimension_changes)[-top_n:][::-1]
    
    return {
        'dimension_changes': dimension_changes,
        'top_changing_dimensions': top_indices,
        'top_changes': dimension_changes[top_indices],
        'centroids': centroids,
        'anc_order': anc_names
    }