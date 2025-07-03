#!/usr/bin/env python3
import os
import sys
import argparse
import gemmi
import glob
import warnings

# 抑制BioPython的PDB构建警告
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def get_asym_id_to_chain_id_mapping(cif_file, use_label_asym_id=False):
    """
    从CIF文件中获取asym_id到chain_id的映射
    
    参数:
        cif_file (str): CIF文件的路径
        use_label_asym_id (bool): 是否强制使用label_asym_id作为chain_id
    
    返回:
        dict: asym_id到chain_id的映射字典
    """
    doc = gemmi.cif.read_file(cif_file)
    block = doc.sole_block()
    
    mapping = {}
    try:
        # 尝试读取label_asym_id和auth_asym_id
        atom_sites = block.find(['_atom_site.label_asym_id', '_atom_site.auth_asym_id'])
        
        # 收集所有的asym_id和对应的auth_asym_id
        asym_to_auth = {}
        asym_to_label = {}
        for label_asym_id, auth_asym_id in atom_sites:
            asym_to_auth[label_asym_id] = auth_asym_id
            asym_to_label[label_asym_id] = label_asym_id
        
        # 检查auth_asym_id的唯一性
        auth_chain_ids = set(asym_to_auth.values())
        label_chain_ids = set(asym_to_label.values())
        
        print(f"发现的label_asym_ids: {sorted(label_chain_ids)}")
        print(f"发现的auth_asym_ids: {sorted(auth_chain_ids)}")
        
        # 如果强制使用label_asym_id，或者auth_asym_id不够唯一，则使用label_asym_id
        if use_label_asym_id or len(auth_chain_ids) < len(asym_to_auth):
            print("使用label_asym_id作为chain_id")
            mapping = asym_to_label
        else:
            print("使用auth_asym_id作为chain_id")
            mapping = asym_to_auth
            
    except Exception as e:
        print(f"无法获取asym_id到chain_id的映射: {e}")
    
    # 如果找不到映射，假设asym_id就是chain_id
    if not mapping:
        print("警告: 无法获取asym_id到chain_id的映射，假设它们相同")
    
    return mapping

def convert_cif_to_pdb(cif_file, output_dir=None, use_label_asym_id=False):
    """
    将CIF文件转换为PDB文件，只提取assembly ID为1的结构
    
    参数:
        cif_file (str): CIF文件的路径
        output_dir (str, optional): 输出目录，如果不指定则使用与输入文件相同的目录
        use_label_asym_id (bool): 是否强制使用label_asym_id作为chain_id
    
    返回:
        bool: 转换是否成功
    """
    # 获取输入文件的目录和文件名
    input_dir = os.path.dirname(cif_file)
    base_name = os.path.splitext(os.path.basename(cif_file))[0]
    
    # 设置输出目录
    if output_dir is None:
        output_dir = input_dir
    
    # 创建输出目录（如果不存在）
    os.makedirs(output_dir, exist_ok=True)
    
    # 构建输出文件路径
    output_file = os.path.join(output_dir, f"{base_name}.pdb")
    
    try:
        print(f"处理文件: {cif_file}")
        
        # 获取asym_id到chain_id的映射
        asym_to_chain = get_asym_id_to_chain_id_mapping(cif_file, use_label_asym_id)
        
        # 读取CIF文件获取assembly信息
        doc = gemmi.cif.read_file(cif_file)
        block = doc.sole_block()
        
        # 找到assembly 1的asym_ids
        assembly_1_asym_ids = None
        try:
            # 尝试查找assembly信息
            assembly_ids = list(block.find_values('_pdbx_struct_assembly.id'))
            print(f"找到的assembly IDs: {assembly_ids}")
            
            # 如果有多个assembly并且含有assembly 1
            if len(assembly_ids) > 1 and '1' in assembly_ids:
                try:
                    # 查找assembly 1包含的链的asym_id
                    assembly_gens = block.find(['_pdbx_struct_assembly_gen.assembly_id', 
                                             '_pdbx_struct_assembly_gen.asym_id_list'])
                    
                    assembly_1_asym_ids = []
                    for gen in assembly_gens:
                        if gen[0] == '1':  # 找到assembly_id=1的记录
                            assembly_1_asym_ids.extend(gen[1].split(','))
                    
                    if assembly_1_asym_ids:
                        print(f"Assembly 1包含的asym_ids: {','.join(assembly_1_asym_ids)}")
                        
                        # 转换为最终的chain_ids
                        assembly_1_final_chain_ids = set()
                        for asym_id in assembly_1_asym_ids:
                            chain_id = asym_to_chain.get(asym_id, asym_id)
                            assembly_1_final_chain_ids.add(chain_id)
                        
                        print(f"Assembly 1映射后的chain_ids: {','.join(sorted(assembly_1_final_chain_ids))}")
                except Exception as e:
                    print(f"提取assembly 1的链信息时出错: {e}")
        except Exception as e:
            print(f"检查assembly信息时出错: {e}")
            print(f"将使用完整结构")
        
        # 直接从CIF数据构建PDB内容
        pdb_lines = []
        pdb_lines.append("HEADER    CONVERTED FROM CIF")
        
        # 读取原子数据
        try:
            atom_columns = [
                '_atom_site.group_PDB',
                '_atom_site.id', 
                '_atom_site.type_symbol',
                '_atom_site.label_atom_id',
                '_atom_site.label_alt_id',
                '_atom_site.label_comp_id',
                '_atom_site.label_asym_id',
                '_atom_site.label_entity_id',
                '_atom_site.label_seq_id',
                '_atom_site.pdbx_PDB_ins_code',
                '_atom_site.Cartn_x',
                '_atom_site.Cartn_y', 
                '_atom_site.Cartn_z',
                '_atom_site.occupancy',
                '_atom_site.B_iso_or_equiv',
                '_atom_site.auth_asym_id'
            ]
            
            atom_data = block.find(atom_columns)
            
            for row in atom_data:
                (group_PDB, atom_id, element, atom_name, alt_loc, residue_name, 
                 label_asym_id, entity_id, seq_id, ins_code, x, y, z, occupancy, 
                 b_factor, auth_asym_id) = row
                
                # 过滤替代位点：只保留A或空的替代位点
                # 如果alt_loc不是'.'、'?'或'A'，则跳过这个原子
                if alt_loc not in ['.', '?', 'A']:
                    continue
                
                # 决定使用哪个asym_id作为chain_id
                if use_label_asym_id or len(set(asym_to_chain.keys())) > len(set(asym_to_chain.values())):
                    final_chain_id = asym_to_chain.get(label_asym_id, label_asym_id)
                else:
                    final_chain_id = asym_to_chain.get(auth_asym_id, auth_asym_id)
                
                # 如果指定了assembly 1，检查是否应该包含这个原子
                if assembly_1_asym_ids is not None:
                    if label_asym_id not in assembly_1_asym_ids:
                        continue
                
                # 处理缺失值和替代位点
                # 对于替代位点，如果我们只保留A替代位点，则在PDB中应设为空格
                if alt_loc in ['.', '?']:
                    alt_loc = ' '
                elif alt_loc == 'A':
                    alt_loc = ' '  # 将A替代位点也设为空格，因为我们只保留A
                    
                ins_code = ins_code if ins_code not in ['.', '?'] else ' '
                seq_id = seq_id if seq_id not in ['.', '?'] else '1'
                occupancy = occupancy if occupancy not in ['.', '?'] else '1.00'
                b_factor = b_factor if b_factor not in ['.', '?'] else '0.00'
                
                # 去除原子名称中的双引号
                atom_name = atom_name.strip('"')
                
                # 构建PDB行
                if group_PDB == 'ATOM':
                    record_type = 'ATOM  '
                elif group_PDB == 'HETATM':
                    record_type = 'HETATM'
                else:
                    record_type = 'ATOM  '
                
                # 格式化PDB行 - 修复格式以确保正确的PDB格式
                # PDB格式：列13-16为原子名称（左对齐），列17为替代位点，列18-20为残基名称（右对齐）
                atom_name_field = f"{atom_name:<4}"  # 原子名称左对齐4位
                residue_name_field = f"{residue_name:>3}"  # 残基名称右对齐3位
                pdb_line = f"{record_type}{int(atom_id):>5} {atom_name_field}{alt_loc}{residue_name_field} {final_chain_id}{int(seq_id):>4}{ins_code}   {float(x):>8.3f}{float(y):>8.3f}{float(z):>8.3f}{float(occupancy):>6.2f}{float(b_factor):>6.2f}          {element:>2}"
                pdb_lines.append(pdb_line)
                
        except Exception as e:
            print(f"处理原子数据时出错: {e}")
            return False
        
        pdb_lines.append("END")
        
        # 写入PDB文件
        with open(output_file, 'w') as f:
            f.write('\n'.join(pdb_lines) + '\n')
        
        # 统计最终的链ID
        final_chain_ids = set()
        for line in pdb_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21]
                final_chain_ids.add(chain_id)
        
        print(f"最终PDB文件包含的链IDs: {','.join(sorted(final_chain_ids))}")
        print(f"成功将 {cif_file} 转换为 {output_file}")
        return True
        
    except Exception as e:
        print(f"转换 {cif_file} 过程中出现错误: {str(e)}")
        return False

def batch_convert(cif_dir, output_dir=None, filter_file=None, use_label_asym_id=False):
    """
    批量转换CIF文件为PDB文件
    
    参数:
        cif_dir (str): CIF文件所在的目录
        output_dir (str, optional): 输出目录
        filter_file (str, optional): 包含要处理文件信息的文本文件，如果不指定则处理目录中所有CIF文件
        use_label_asym_id (bool): 是否强制使用label_asym_id作为chain_id
    """
    # 如果未指定输出目录，则在CIF目录下创建pdb_output子目录
    if output_dir is None:
        output_dir = os.path.join(cif_dir, "pdb_output")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 统计转换结果
    success_count = 0
    failed_count = 0
    not_found_count = 0
    
    try:
        # 定义要处理的文件列表
        pdb_ids = []
        
        # 如果提供了filter_file，则从文件中读取PDB ID
        if filter_file:
            if not os.path.exists(filter_file):
                print(f"过滤文件不存在: {filter_file}")
                sys.exit(1)
                
            with open(filter_file, 'r') as f:
                for line in f:
                    # 分割行内容，取第一部分作为PDB ID
                    parts = line.strip().split()
                    if parts:
                        pdb_ids.append(parts[0])  # 只取第一列作为PDB ID
                        
            print(f"从文件 {filter_file} 中读取了 {len(pdb_ids)} 个PDB ID")
            
            # 处理从filter_file获取的PDB ID
            for pdb_id in pdb_ids:
                # 构建CIF文件路径
                cif_file = os.path.join(cif_dir, f"{pdb_id}.cif")
                
                # 检查文件是否存在
                if os.path.exists(cif_file):
                    # 转换文件
                    if convert_cif_to_pdb(cif_file, output_dir, use_label_asym_id):
                        success_count += 1
                    else:
                        failed_count += 1
                else:
                    print(f"文件不存在: {cif_file}")
                    not_found_count += 1
        
        # 如果没有提供filter_file，则处理目录中所有的CIF文件
        else:
            if not os.path.exists(cif_dir):
                print(f"CIF目录不存在: {cif_dir}")
                sys.exit(1)
                
            print(f"未提供过滤文件，将处理目录 {cif_dir} 中的所有CIF文件")
            cif_files = glob.glob(os.path.join(cif_dir, "*.cif"))
            
            if not cif_files:
                print(f"目录 {cif_dir} 中没有找到CIF文件")
                return
                
            print(f"找到 {len(cif_files)} 个CIF文件")
            
            # 处理所有CIF文件
            for cif_file in cif_files:
                if convert_cif_to_pdb(cif_file, output_dir, use_label_asym_id):
                    success_count += 1
                else:
                    failed_count += 1
        
        print(f"\n批量转换完成！")
        print(f"成功转换: {success_count} 个文件")
        print(f"转换失败: {failed_count} 个文件")
        if filter_file:
            print(f"未找到: {not_found_count} 个文件")
        
    except Exception as e:
        print(f"批量转换过程中出现错误: {str(e)}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='批量将CIF文件转换为PDB文件')
    parser.add_argument('cif_dir', help='CIF文件所在目录')
    parser.add_argument('--filter_file', help='包含要处理的PDB ID的文本文件（可选）', default=None)
    parser.add_argument('--output_dir', help='输出目录（可选）', default=None)
    parser.add_argument('--use_label_asym_id', action='store_true', 
                       help='强制使用label_asym_id作为chain_id（适用于auth_asym_id都相同的情况）')
    
    args = parser.parse_args()
    
    # 检查目录是否存在
    if not os.path.exists(args.cif_dir):
        print(f"错误: CIF目录不存在: {args.cif_dir}")
        sys.exit(1)
        
    if args.filter_file and not os.path.exists(args.filter_file):
        print(f"错误: 过滤文件不存在: {args.filter_file}")
        sys.exit(1)
    
    batch_convert(args.cif_dir, args.output_dir, args.filter_file, args.use_label_asym_id)

if __name__ == '__main__':
    main() 