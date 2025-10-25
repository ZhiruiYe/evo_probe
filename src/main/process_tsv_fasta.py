#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生物体序列过滤器 - 细菌式编程风格

遵循细菌编程原则：小巧、模块化、可复制粘贴
每个函数都是独立的"基因"，可以单独使用
"""

import csv
from collections import Counter
from pathlib import Path


def count_organisms(tsv_file_path):
    """
    统计TSV文件中每个生物体(Organism)的序列数量
    
    Args:
        tsv_file_path (str): TSV文件路径
    
    Returns:
        dict: {organism_name: count} 的字典
    
    Example:
        counts = count_organisms("data.tsv")
        print(counts["Homo sapiens"])  # 输出该物种的序列数量
    """
    organism_counts = Counter()
    
    with open(tsv_file_path, 'r', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            organism = row.get('Organism', '').strip()
            if organism:
                organism_counts[organism] += 1
    
    return dict(organism_counts)


def filter_organisms_by_count(tsv_file_path, min_count=10):
    """
    过滤掉序列数量少于指定阈值的生物体
    
    Args:
        tsv_file_path (str): 输入TSV文件路径
        min_count (int): 最小序列数量阈值，默认10
    
    Returns:
        tuple: (filtered_rows, valid_organisms)
            - filtered_rows: 过滤后的数据行列表
            - valid_organisms: 保留的生物体集合
    
    Example:
        rows, orgs = filter_organisms_by_count("data.tsv", 5)
        print(f"保留了 {len(orgs)} 个物种")
    """
    # 先统计每个生物体的数量
    organism_counts = count_organisms(tsv_file_path)
    
    # 找出符合条件的生物体
    valid_organisms = {
        organism for organism, count in organism_counts.items() 
        if count >= min_count
    }
    
    # 过滤数据
    filtered_rows = []
    with open(tsv_file_path, 'r', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        header = reader.fieldnames
        
        for row in reader:
            organism = row.get('Organism', '').strip()
            if organism in valid_organisms:
                filtered_rows.append(row)
    
    return filtered_rows, valid_organisms, header


def write_filtered_tsv(filtered_rows, header, output_path):
    """
    将过滤后的数据写入新的TSV文件
    
    Args:
        filtered_rows (list): 过滤后的数据行
        header (list): 表头字段名
        output_path (str): 输出文件路径
    
    Returns:
        int: 写入的行数
    
    Example:
        count = write_filtered_tsv(rows, header, "filtered_data.tsv")
        print(f"已写入 {count} 行数据")
    """
    with open(output_path, 'w', encoding='utf-8', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=header, delimiter='\t')
        writer.writeheader()
        writer.writerows(filtered_rows)
    
    return len(filtered_rows)


def print_organism_summary(valid_organisms, organism_counts):
    """
    打印生物体统计摘要
    
    Args:
        valid_organisms (set): 保留的生物体集合
        organism_counts (dict): 生物体计数字典
    
    Example:
        print_organism_summary(valid_orgs, counts)
    """
    print(f"\n=== 过滤结果摘要 ===")
    print(f"保留的生物体数量: {len(valid_organisms)}")
    print(f"\n保留的生物体列表:")
    
    # 按序列数量排序显示
    sorted_organisms = sorted(
        [(org, organism_counts[org]) for org in valid_organisms],
        key=lambda x: x[1], reverse=True
    )
    
    for i, (organism, count) in enumerate(sorted_organisms, 1):
        print(f"{i:3d}. {organism:<50} ({count:4d} 序列)")


def parse_fasta_file(fasta_file_path):
    """
    解析FASTA文件，建立Entry ID到序列的映射
    
    Args:
        fasta_file_path (str): FASTA文件路径
    
    Returns:
        dict: {entry_id: (header, sequence)} 的字典
    
    Example:
        fasta_dict = parse_fasta_file("sequences.fasta")
        header, seq = fasta_dict["A0A0B4J1N3"]
    """
    fasta_dict = {}
    current_header = None
    current_sequence = []
    
    with open(fasta_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # 保存前一个序列
                if current_header and current_sequence:
                    # 从header中提取Entry ID (格式: >sp|ENTRY|NAME 或 >tr|ENTRY|NAME)
                    entry_id = current_header.split('|')[1] if '|' in current_header else current_header[1:]
                    fasta_dict[entry_id] = (current_header, ''.join(current_sequence))
                
                # 开始新序列
                current_header = line
                current_sequence = []
            else:
                # 添加序列行
                current_sequence.append(line)
        
        # 保存最后一个序列
        if current_header and current_sequence:
            entry_id = current_header.split('|')[1] if '|' in current_header else current_header[1:]
            fasta_dict[entry_id] = (current_header, ''.join(current_sequence))
    
    return fasta_dict


def extract_sequences_from_tsv(tsv_file_path, fasta_dict):
    """
    从TSV文件中提取Entry列表，并从FASTA字典中获取对应序列
    
    Args:
        tsv_file_path (str): TSV文件路径
        fasta_dict (dict): FASTA序列字典
    
    Returns:
        tuple: (found_sequences, missing_entries)
            - found_sequences: [(entry_id, header, sequence), ...] 找到的序列列表
            - missing_entries: [entry_id, ...] 未找到的Entry列表
    
    Example:
        sequences, missing = extract_sequences_from_tsv("data.tsv", fasta_dict)
        print(f"找到 {len(sequences)} 个序列，缺失 {len(missing)} 个")
    """
    found_sequences = []
    missing_entries = []
    
    with open(tsv_file_path, 'r', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            entry_id = row.get('Entry', '').strip()
            if entry_id:
                if entry_id in fasta_dict:
                    header, sequence = fasta_dict[entry_id]
                    found_sequences.append((entry_id, header, sequence))
                else:
                    missing_entries.append(entry_id)
    
    return found_sequences, missing_entries


def write_fasta_file(sequences, output_path):
    """
    将序列列表写入FASTA文件
    
    Args:
        sequences (list): [(entry_id, header, sequence), ...] 序列列表
        output_path (str): 输出FASTA文件路径
    
    Returns:
        int: 写入的序列数量
    
    Example:
        count = write_fasta_file(sequences, "filtered.fasta")
        print(f"已写入 {count} 个序列")
    """
    with open(output_path, 'w', encoding='utf-8') as file:
        for entry_id, header, sequence in sequences:
            file.write(f"{header}\n")
            # 将序列按80个字符一行进行格式化
            for i in range(0, len(sequence), 80):
                file.write(f"{sequence[i:i+80]}\n")
    
    return len(sequences)


def extract_and_filter_fasta(tsv_file_path, fasta_file_path, output_fasta_path):
    """
    根据TSV文件中的Entry从FASTA文件中提取序列并写入新文件
    
    Args:
        tsv_file_path (str): TSV文件路径
        fasta_file_path (str): 原始FASTA文件路径
        output_fasta_path (str): 输出FASTA文件路径
    
    Returns:
        tuple: (found_count, missing_count, missing_entries)
    
    Example:
        found, missing, missing_list = extract_and_filter_fasta(
            "filtered.tsv", "all.fasta", "filtered.fasta"
        )
    """
    print(f"正在解析FASTA文件: {fasta_file_path}")
    fasta_dict = parse_fasta_file(fasta_file_path)
    print(f"解析完成，共找到 {len(fasta_dict)} 个序列")
    
    print(f"正在从TSV文件提取序列: {tsv_file_path}")
    found_sequences, missing_entries = extract_sequences_from_tsv(tsv_file_path, fasta_dict)
    
    print(f"写入FASTA文件: {output_fasta_path}")
    written_count = write_fasta_file(found_sequences, output_fasta_path)
    
    return len(found_sequences), len(missing_entries), missing_entries


def main():

    tsv_path ="/yezhirui/evo_probe/xcl1_tree/seq/fig_28_or_info.tsv"
    fasta_input = "/yezhirui/evo_probe/xcl1_tree/seq/uniprotkb_taxonomy_id_7742.fasta"
    fasta_output = "/yezhirui/evo_probe/xcl1_tree/seq/fig_28_or_seqs.fasta"
    
    # 提取对应的FASTA序列
    print(f"\n=== 开始提取FASTA序列 ===")
    found_count, missing_count, missing_entries = extract_and_filter_fasta(
        tsv_path, fasta_input, fasta_output
    )
    
    print(f"\n=== FASTA提取结果 ===")
    print(f"输出FASTA文件: {fasta_output}")
    print(f"成功找到序列: {found_count}")
    print(f"缺失序列数量: {missing_count}")
    
    


if __name__ == "__main__":
    main()
