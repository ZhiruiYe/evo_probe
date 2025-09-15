"""PDB下载工具模块 

自包含的PDB文件下载器，可独立使用
"""

import os
import requests
import time
from pathlib import Path
from typing import List, Optional, Tuple
import logging


def download_pdb_file(pdb_id: str, output_dir: str, max_retries: int = 3) -> Tuple[bool, str]:
    """下载单个PDB文件
    
    Args:
        pdb_id: 4字符PDB ID (如 '1abc')
        output_dir: 输出目录路径
        max_retries: 最大重试次数
    
    Returns:
        (success: bool, filepath: str) - 成功状态和文件路径
    
    Example:
        success, path = download_pdb_file('1abc', './pdbs')
    """
    pdb_id = pdb_id.lower().strip()
    if len(pdb_id) != 4:
        return False, f"Invalid PDB ID: {pdb_id}"
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    filepath = output_dir / f"{pdb_id}.pdb"
    
    # 如果文件已存在且大小>0，跳过下载
    if filepath.exists() and filepath.stat().st_size > 0:
        return True, str(filepath)
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                with open(filepath, 'w') as f:
                    f.write(response.text)
                return True, str(filepath)
            elif response.status_code == 404:
                return False, f"PDB {pdb_id} not found"
            else:
                logging.warning(f"HTTP {response.status_code} for {pdb_id}, attempt {attempt + 1}")
                
        except Exception as e:
            logging.warning(f"Download error for {pdb_id}: {e}, attempt {attempt + 1}")
            
        if attempt < max_retries - 1:
            time.sleep(2 ** attempt)  # 指数退避
    
    return False, f"Failed to download {pdb_id} after {max_retries} attempts"


def batch_download_pdbs(pdb_ids: List[str], output_dir: str, 
                       delay: float = 0.1) -> Tuple[List[str], List[str]]:
    """批量下载PDB文件
    
    Args:
        pdb_ids: PDB ID列表
        output_dir: 输出目录
        delay: 下载间隔(秒)，避免服务器过载
    
    Returns:
        (successful_ids, failed_ids) - 成功和失败的PDB ID列表
    
    Example:
        success, failed = batch_download_pdbs(['1abc', '2def'], './pdbs')
    """
    successful = []
    failed = []
    
    for i, pdb_id in enumerate(pdb_ids):
        success, path = download_pdb_file(pdb_id, output_dir)
        if success:
            successful.append(pdb_id)
            logging.info(f"Downloaded {pdb_id} -> {path}")
        else:
            failed.append(pdb_id)
            logging.error(f"Failed to download {pdb_id}: {path}")
        
        # 进度日志
        if (i + 1) % 10 == 0:
            logging.info(f"Progress: {i + 1}/{len(pdb_ids)} PDBs processed")
        
        time.sleep(delay)
    
    return successful, failed


def extract_unique_pdb_ids(csv_file: str, pdb_columns: List[str] = None) -> List[str]:
    """从CSV文件提取唯一的PDB ID
    
    Args:
        csv_file: CSV文件路径
        pdb_columns: PDB ID列名列表，默认['pdb_id_1', 'pdb_id_2']
    
    Returns:
        去重的PDB ID列表
    
    Example:
        pdb_ids = extract_unique_pdb_ids('pairs.csv')
    """
    import pandas as pd
    
    if pdb_columns is None:
        pdb_columns = ['pdb_id_1', 'pdb_id_2']
    
    df = pd.read_csv(csv_file)
    
    all_pdb_ids = []
    for col in pdb_columns:
        if col in df.columns:
            all_pdb_ids.extend(df[col].dropna().astype(str).tolist())
    
    # 去重并过滤有效的PDB ID
    unique_ids = list(set(pdb_id.strip().lower() 
                         for pdb_id in all_pdb_ids 
                         if len(pdb_id.strip()) == 4))
    
    return sorted(unique_ids)


if __name__ == "__main__":
    # 独立运行示例
    logging.basicConfig(level=logging.INFO)
    
    # 从CSV提取PDB ID并下载
    pdb_ids = extract_unique_pdb_ids("/yezhirui/evo_probe/xcl1_tree/str/align_pdb_list.csv",['pdb'])
    print(f"Found {len(pdb_ids)} unique PDB IDs")
    
    successful, failed = batch_download_pdbs(
        pdb_ids, 
        "/yezhirui/evo_probe/xcl1_tree/str",
        delay=0.1
    )
    
    print(f"Successfully downloaded: {len(successful)}")
    print(f"Failed downloads: {len(failed)}")
    if failed:
        print("Failed IDs:", failed[:10])  # 显示前10个失败的
