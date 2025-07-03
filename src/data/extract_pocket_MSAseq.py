import pickle
import numpy as np
import os
import argparse

def extract_pocket_sequences(pkl_file):
    """
    从pkl文件中提取口袋序列
    
    参数:
    - pkl_file: pkl文件路径
    
    返回:
    - 包含每条链口袋序列的字典
    """
    # 读取pkl文件
    with open(pkl_file, 'rb') as f:
        data = pickle.load(f)
    
    # 提取需要的信息
    seq_info = data['seq_info']
    pocket_mask = data['pocket_mask']
    protein_mask = data['protein_mask']
    
    # 创建结果字典
    result = {}
    
    # 获取所有链的信息
    chains = seq_info['chains']
    
    # 跟踪当前序列的索引位置
    current_idx = 0
    
    for chain in chains:
        chain_id = chain['chain_id']
        sequence = chain['sequence']
        chain_type = chain['chain_type']
        
        # 只处理蛋白质链
        if protein_mask[current_idx]:
            # 提取该链对应的pocket_mask部分
            chain_length = len(sequence)
            chain_pocket_mask = pocket_mask[current_idx:current_idx+chain_length]
            
            # 根据pocket_mask构建口袋序列
            pocket_sequence = ''
            for i, (res, mask) in enumerate(zip(sequence, chain_pocket_mask)):
                if mask == 1:
                    pocket_sequence += res
                else:
                    pocket_sequence += '-'
            
            # 将结果添加到字典
            result[chain_id] = {
                'sequence': sequence,
                'pocket_sequence': pocket_sequence,
                'chain_type': chain_type
            }
        
        # 更新索引位置
        current_idx += len(sequence)
    
    return result

def main():
    parser = argparse.ArgumentParser(description='提取蛋白质口袋序列')
    parser.add_argument('--input', type=str, required=True, help='输入的pkl文件路径')
    parser.add_argument('--output', type=str, required=True, help='输出的序列文件路径')
    args = parser.parse_args()
    
    # 提取口袋序列
    pocket_sequences = extract_pocket_sequences(args.input)
    
    # 将结果写入文件
    with open(args.output, 'w') as f:
        for chain_id, info in pocket_sequences.items():
            f.write(f"Chain {chain_id} ({info['chain_type']}):\n")
            f.write(f"Original sequence: {info['sequence']}\n")
            f.write(f"Pocket sequence: {info['pocket_sequence']}\n\n")
    
    print(f"已成功提取口袋序列并保存到 {args.output}")

if __name__ == "__main__":
    main()
