#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA到PHY格式转换工具 - 细菌式编程风格

遵循细菌编程原则：小巧、模块化、可复制粘贴
专门用于将FASTA格式转换为简单的PHY格式
"""


def parse_fasta(fasta_file_path):
    """
    解析FASTA文件
    
    Args:
        fasta_file_path (str): FASTA文件路径
    
    Returns:
        list: [(header, sequence), ...] 序列列表
    
    Example:
        sequences = parse_fasta("input.fasta")
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
                    sequences.append((current_header[1:], ''.join(current_sequence)))  # 去掉'>'
                
                # 开始新序列
                current_header = line
                current_sequence = []
            else:
                # 添加序列行
                current_sequence.append(line)
        
        # 保存最后一个序列
        if current_header and current_sequence:
            sequences.append((current_header[1:], ''.join(current_sequence)))  # 去掉'>'
    
    return sequences


def write_phy_format(sequences, output_file_path, header_width=50):
    """
    将序列写成PHY格式
    
    Args:
        sequences (list): [(header, sequence), ...] 序列列表
        output_file_path (str): 输出文件路径
        header_width (int): header字段宽度，默认50
    
    Returns:
        int: 写入的序列数量
    
    Example:
        count = write_phy_format(sequences, "output.phy")
        print(f"已写入 {count} 个序列")
    """
    with open(output_file_path, 'w', encoding='utf-8') as file:
        for header, sequence in sequences:
            # 确保header不超过指定宽度，用空格填充到指定宽度
            formatted_header = header[:header_width].ljust(header_width)
            file.write(f"{formatted_header}\t{sequence}\n")
    
    return len(sequences)


def fasta_to_phy(input_fasta_file, output_phy_file, header_width=50):
    """
    完整的FASTA到PHY转换流程
    
    Args:
        input_fasta_file (str): 输入FASTA文件路径
        output_phy_file (str): 输出PHY文件路径
        header_width (int): header字段宽度
    
    Returns:
        dict: 转换结果统计
    
    Example:
        result = fasta_to_phy("input.fasta", "output.phy")
        print(f"转换了 {result['count']} 个序列")
    """
    print(f"📖 正在读取FASTA文件: {input_fasta_file}")
    sequences = parse_fasta(input_fasta_file)
    print(f"   找到 {len(sequences)} 个序列")
    
    if sequences:
        seq_length = len(sequences[0][1]) if sequences else 0
        print(f"   序列长度: {seq_length} bp")
    
    print(f"💾 正在写入PHY文件: {output_phy_file}")
    written_count = write_phy_format(sequences, output_phy_file, header_width)
    print(f"   已写入 {written_count} 个序列")
    
    result = {
        'count': written_count,
        'seq_length': seq_length if sequences else 0
    }
    
    print(f"\n📊 转换统计:")
    print(f"   序列数量: {result['count']}")
    print(f"   序列长度: {result['seq_length']} bp")
    print(f"   Header宽度: {header_width} 字符")
    
    return result


def main():
    """主函数 - 执行FASTA到PHY转换"""
    # 文件路径配置
    input_file = "/yezhirui/evo_probe/xcl1_tree/str/final_trimmed.fasta"
    output_file = "/yezhirui/evo_probe/xcl1_tree/str/final_trimmed.phy"
    
    # 格式参数
    header_field_width = 50  # header字段宽度
    
    print(f"🔄 FASTA到PHY格式转换工具")
    print(f"输入文件: {input_file}")
    print(f"输出文件: {output_file}")
    print(f"Header宽度: {header_field_width}")
    print(f"{'='*60}")
    
    # 执行转换
    result = fasta_to_phy(input_file, output_file, header_field_width)
    
    print(f"{'='*60}")
    print(f"✅ 转换完成！")
    
    return result


if __name__ == "__main__":
    main()
