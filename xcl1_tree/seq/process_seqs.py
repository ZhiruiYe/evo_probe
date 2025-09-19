def extract_species_name(header):
    """
    从FASTA头部提取物种名称（[]中的内容）
    
    Args:
        header: FASTA头部行
    
    Returns:
        物种名称，如果没找到则返回'Unknown'
    """
    import re
    match = re.search(r'\[([^\]]+)\]', header)
    return match.group(1) if match else 'Unknown'

def extract_isoform_sequences(input_fasta, output_fasta):
    """
    从FASTA文件中提取所有头部包含'isoform'的序列并写入新文件
    
    Args:
        input_fasta: 输入FASTA文件路径
        output_fasta: 输出FASTA文件路径
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        current_header = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            
            if line.startswith('>'):
                # 如果之前的序列包含isoform，写入输出文件
                if current_header and 'isoform' in current_header.lower():
                    outfile.write(current_header + '\n')
                    outfile.write(''.join(current_sequence) + '\n')
                
                # 重置当前序列信息
                current_header = line
                current_sequence = []
            else:
                # 累积序列数据
                current_sequence.append(line)
        
        # 处理最后一个序列
        if current_header and 'isoform' in current_header.lower():
            outfile.write(current_header + '\n')
            outfile.write(''.join(current_sequence) + '\n')

def extract_non_isoform_sequences(input_fasta, output_fasta):
    """
    从FASTA文件中提取所有头部不包含'isoform'的序列并写入新文件
    
    Args:
        input_fasta: 输入FASTA文件路径
        output_fasta: 输出FASTA文件路径
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        current_header = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            
            if line.startswith('>'):
                # 如果之前的序列不包含isoform，写入输出文件
                if current_header and 'isoform' not in current_header.lower():
                    outfile.write(current_header + '\n')
                    outfile.write(''.join(current_sequence) + '\n')
                
                # 重置当前序列信息
                current_header = line
                current_sequence = []
            else:
                # 累积序列数据
                current_sequence.append(line)
        
        # 处理最后一个序列
        if current_header and 'isoform' not in current_header.lower():
            outfile.write(current_header + '\n')
            outfile.write(''.join(current_sequence) + '\n')

def filter_sequences_by_length(input_fasta, output_fasta, min_length=50, max_length=150):
    """
    移除序列长度小于指定长度的序列并写入新文件
    
    Args:
        input_fasta: 输入FASTA文件路径
        output_fasta: 输出FASTA文件路径
        min_length: 最小序列长度阈值，默认50
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        current_header = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            
            if line.startswith('>'):
                # 如果之前的序列长度满足要求，写入输出文件
                if current_header and len(''.join(current_sequence)) >= min_length and len(''.join(current_sequence)) <= max_length:
                    outfile.write(current_header + '\n')
                    outfile.write(''.join(current_sequence) + '\n')
                else :
                    print(current_header)
                
                # 重置当前序列信息
                current_header = line
                current_sequence = []
            else:
                # 累积序列数据
                current_sequence.append(line)
        
        # 处理最后一个序列
        if current_header and len(''.join(current_sequence)) >= min_length:
            outfile.write(current_header + '\n')
            outfile.write(''.join(current_sequence) + '\n')

def group_isoforms_by_species(input_fasta, output_fasta):
    """
    按物种聚合isoform序列并写入新文件
    
    Args:
        input_fasta: 输入FASTA文件路径
        output_fasta: 输出FASTA文件路径
    """
    species_groups = {}
    
    # 读取所有isoform序列并按物种分组
    with open(input_fasta, 'r') as infile:
        current_header = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            
            if line.startswith('>'):
                # 处理前一个序列
                if current_header and 'isoform' in current_header.lower():
                    species = extract_species_name(current_header)
                    if species not in species_groups:
                        species_groups[species] = []
                    species_groups[species].append((current_header, ''.join(current_sequence)))
                
                # 重置当前序列信息
                current_header = line
                current_sequence = []
            else:
                # 累积序列数据
                current_sequence.append(line)
        
        # 处理最后一个序列
        if current_header and 'isoform' in current_header.lower():
            species = extract_species_name(current_header)
            if species not in species_groups:
                species_groups[species] = []
            species_groups[species].append((current_header, ''.join(current_sequence)))
    
    # 按物种排序并写入输出文件
    with open(output_fasta, 'w') as outfile:
        for species in sorted(species_groups.keys()):
            for header, sequence in species_groups[species]:
                outfile.write(header + '\n')
                outfile.write(sequence + '\n')
            outfile.write('\n')

# 使用示例
if __name__ == "__main__":
    input_file = "blast_nr_598.fasta"
    # isoform_output = "isoform_seqs.fasta"
    # grouped_output = "isoform_seqs_by_species.fasta"
    non_isoform_output = "non_isoform_seqs.fasta"
    
    # # 第一步：提取所有isoform序列
    # extract_isoform_sequences(input_file, isoform_output)
    # print(f"已将包含isoform的序列提取到 {isoform_output}")
    
    # # 第二步：按物种聚合isoform序列
    # group_isoforms_by_species(input_file, grouped_output)
    # print(f"已按物种聚合isoform序列到 {grouped_output}")
    
    # 第三步：提取所有非isoform序列
    extract_non_isoform_sequences(input_file, non_isoform_output)
    print(f"已将不包含isoform的序列提取到 {non_isoform_output}")
    
    # 第四步：移除长度小于50的序列
    filtered_output = "filtered_seqs_min50_max120.fasta"
    filter_sequences_by_length(non_isoform_output, filtered_output, 50, 120)
    print(f"已移除长度<50的序列，结果保存到 {filtered_output}")