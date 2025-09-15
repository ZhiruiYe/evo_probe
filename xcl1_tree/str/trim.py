#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
序列裁剪工具 - 细菌式编程风格

遵循细菌编程原则：小巧、模块化、可复制粘贴
专门用于从对齐序列中提取指定位置的片段
"""

import os
from pathlib import Path


def parse_aligned_fasta(fasta_file_path):
    """
    解析对齐后的FASTA文件
    
    Args:
        fasta_file_path (str): FASTA文件路径
    
    Returns:
        list: [(header, sequence), ...] 序列列表
    
    Example:
        sequences = parse_aligned_fasta("aligned.fasta")
        for header, seq in sequences:
            print(f"{header}: {len(seq)} bp")
    """
    sequences = []
    current_header = None
    current_sequence = []
    
    with open(fasta_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # 保存前一个序列
                if current_header and current_sequence:
                    sequences.append((current_header, ''.join(current_sequence)))
                
                # 开始新序列
                current_header = line
                current_sequence = []
            else:
                # 添加序列行
                current_sequence.append(line)
        
        # 保存最后一个序列
        if current_header and current_sequence:
            sequences.append((current_header, ''.join(current_sequence)))
    
    return sequences


def trim_sequences(sequences, start_pos, end_pos):
    """
    从序列中提取指定位置的片段
    
    Args:
        sequences (list): [(header, sequence), ...] 序列列表
        start_pos (int): 起始位置 (1-based)
        end_pos (int): 结束位置 (1-based, inclusive)
    
    Returns:
        list: [(header, trimmed_sequence), ...] 裁剪后的序列列表
    
    Example:
        trimmed = trim_sequences(sequences, 136, 273)
        print(f"提取了 {len(trimmed)} 个序列的 {start_pos}-{end_pos} 位片段")
    """
    trimmed_sequences = []
    
    # 转换为0-based索引
    start_idx = start_pos - 1
    end_idx = end_pos  # end_pos是inclusive的，所以不需要减1
    
    for header, sequence in sequences:
        # 检查序列长度
        if len(sequence) >= end_pos:
            trimmed_seq = sequence[start_idx:end_idx]
            trimmed_sequences.append((header, trimmed_seq))
        else:
            print(f"⚠️  序列 {header} 长度不足 ({len(sequence)} < {end_pos})，跳过")
    
    return trimmed_sequences


def remove_gap_only_positions(sequences):
    """
    移除在所有序列中都是gap(-)的位置
    
    Args:
        sequences (list): [(header, sequence), ...] 序列列表
    
    Returns:
        list: [(header, cleaned_sequence), ...] 清理后的序列列表
    
    Example:
        cleaned = remove_gap_only_positions(sequences)
        print(f"移除了gap-only位置后的序列")
    """
    if not sequences:
        return sequences
    
    # 获取序列长度（假设所有序列长度相同）
    seq_length = len(sequences[0][1])
    
    # 找出非gap-only的位置
    valid_positions = []
    for pos in range(seq_length):
        has_non_gap = False
        for header, sequence in sequences:
            if pos < len(sequence) and sequence[pos] != '-':
                has_non_gap = True
                break
        if has_non_gap:
            valid_positions.append(pos)
    
    # 构建清理后的序列
    cleaned_sequences = []
    for header, sequence in sequences:
        cleaned_seq = ''.join([sequence[pos] for pos in valid_positions if pos < len(sequence)])
        cleaned_sequences.append((header, cleaned_seq))
    
    return cleaned_sequences


def write_fasta_file(sequences, output_path, line_width=80):
    """
    将序列写入FASTA文件
    
    Args:
        sequences (list): [(header, sequence), ...] 序列列表
        output_path (str): 输出文件路径
        line_width (int): 每行字符数，默认80
    
    Returns:
        int: 写入的序列数量
    
    Example:
        count = write_fasta_file(sequences, "trimmed.fasta")
        print(f"已写入 {count} 个序列")
    """
    with open(output_path, 'w', encoding='utf-8') as file:
        for header, sequence in sequences:
            file.write(f"{header}\n")
            # 按指定宽度换行
            for i in range(0, len(sequence), line_width):
                file.write(f"{sequence[i:i+line_width]}\n")
    
    return len(sequences)


def calculate_statistics(original_sequences, trimmed_sequences):
    """
    计算裁剪统计信息
    
    Args:
        original_sequences (list): 原始序列列表
        trimmed_sequences (list): 裁剪后序列列表
    
    Returns:
        dict: 统计信息字典
    
    Example:
        stats = calculate_statistics(orig_seqs, trimmed_seqs)
        print(f"保留率: {stats['retention_rate']:.1%}")
    """
    if not original_sequences or not trimmed_sequences:
        return {}
    
    original_length = len(original_sequences[0][1]) if original_sequences else 0
    trimmed_length = len(trimmed_sequences[0][1]) if trimmed_sequences else 0
    
    stats = {
        'original_count': len(original_sequences),
        'trimmed_count': len(trimmed_sequences),
        'original_length': original_length,
        'trimmed_length': trimmed_length,
        'retention_rate': len(trimmed_sequences) / len(original_sequences) if original_sequences else 0,
        'length_ratio': trimmed_length / original_length if original_length > 0 else 0
    }
    
    return stats


def trim_aligned_fasta(input_file, output_file, start_pos, end_pos, remove_gaps=True):
    """
    完整的序列裁剪流程
    
    Args:
        input_file (str): 输入FASTA文件路径
        output_file (str): 输出FASTA文件路径
        start_pos (int): 起始位置 (1-based)
        end_pos (int): 结束位置 (1-based, inclusive)
        remove_gaps (bool): 是否移除gap-only位置
    
    Returns:
        dict: 处理结果统计
    
    Example:
        result = trim_aligned_fasta("aligned.fasta", "trimmed.fasta", 136, 273)
        print(f"成功处理 {result['trimmed_count']} 个序列")
    """
    print(f"📖 正在读取对齐文件: {input_file}")
    original_sequences = parse_aligned_fasta(input_file)
    print(f"   找到 {len(original_sequences)} 个序列")
    
    if original_sequences:
        print(f"   原始序列长度: {len(original_sequences[0][1])} bp")
    
    print(f"✂️  正在提取位置 {start_pos}-{end_pos}")
    trimmed_sequences = trim_sequences(original_sequences, start_pos, end_pos)
    print(f"   成功提取 {len(trimmed_sequences)} 个序列")
    
    if remove_gaps and trimmed_sequences:
        print(f"🧹 正在移除gap-only位置")
        original_trimmed_length = len(trimmed_sequences[0][1]) if trimmed_sequences else 0
        trimmed_sequences = remove_gap_only_positions(trimmed_sequences)
        final_length = len(trimmed_sequences[0][1]) if trimmed_sequences else 0
        print(f"   序列长度: {original_trimmed_length} → {final_length} bp")
    
    print(f"💾 正在写入文件: {output_file}")
    written_count = write_fasta_file(trimmed_sequences, output_file)
    print(f"   已写入 {written_count} 个序列")
    
    # 计算统计信息
    stats = calculate_statistics(original_sequences, trimmed_sequences)
    
    print(f"\n📊 处理统计:")
    print(f"   原始序列数: {stats.get('original_count', 0)}")
    print(f"   输出序列数: {stats.get('trimmed_count', 0)}")
    print(f"   序列保留率: {stats.get('retention_rate', 0):.1%}")
    print(f"   长度比例: {stats.get('length_ratio', 0):.1%}")
    
    return stats


def main():
    """主函数 - 执行序列裁剪"""

    # 裁剪参数
    start_position = 137
    end_position = 273

    # 文件路径配置
    input_file = "/yezhirui/evo_probe/xcl1_tree/str/Promals3d_aligned.fasta"
    output_file = f"/yezhirui/evo_probe/xcl1_tree/str/Promals3d_aligned_{start_position}_{end_position}.fasta"
    
    
    print(f"🧬 序列裁剪工具")
    print(f"输入文件: {input_file}")
    print(f"输出文件: {output_file}")
    print(f"提取位置: {start_position}-{end_position}")
    print(f"{'='*60}")
    
    # 执行裁剪
    result = trim_aligned_fasta(
        input_file, 
        output_file, 
        start_position, 
        end_position,
        remove_gaps=True
    )
    
    print(f"{'='*60}")
    print(f"✅ 裁剪完成！")
    
    return result


if __name__ == "__main__":
    main()
