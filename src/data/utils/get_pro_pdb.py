# 从蛋白核酸复合物中获取蛋白质的PDB文件

import os
import argparse
import pandas as pd
from tqdm import tqdm # 用于显示进度条

# 标准氨基酸三字母代码集合
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", 
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", 
    "TYR", "VAL"
}

def get_pro_pdb(input_pdb_path, output_pdb_path):
    """
    从蛋白核酸复合物PDB文件中提取蛋白质链，并保存为新的PDB文件。
    只保留记录类型为 'ATOM' 且残基名称为标准氨基酸的行。
    """
    try:
        with open(input_pdb_path, 'r') as infile, open(output_pdb_path, 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM"):
                    # PDB格式中，残基名称在第18-20列 (1-indexed) 或 17-19 (0-indexed)
                    residue_name = line[17:20].strip()
                    if residue_name in STANDARD_AMINO_ACIDS:
                        outfile.write(line)
                elif line.startswith("TER") or line.startswith("ENDMDL"): # 保留链终止符和模型结束符
                     # 可选：检查上一个写入的原子是否是蛋白质，以决定是否写入 TER
                     # 为简单起见，暂时保留所有TER和ENDMDL
                    outfile.write(line)
    except FileNotFoundError:
        print(f"错误：输入文件未找到 {input_pdb_path}")
    except Exception as e:
        print(f"处理文件 {input_pdb_path} 时发生错误: {e}")


def main():
    parser = argparse.ArgumentParser(description='从蛋白核酸复合物 PDB 文件中提取蛋白质链')
    parser.add_argument('--pdb_file', type=str, help='单个蛋白核酸复合物的 PDB 文件路径')
    parser.add_argument('--input_dir', type=str, help='包含 PDB 文件的输入目录路径')
    parser.add_argument('--output_dir', type=str, required=True, help='保存蛋白质 PDB 文件的输出目录')

    args = parser.parse_args()

    # 检查输入参数
    if not args.pdb_file and not args.input_dir:
        parser.error("必须提供 --pdb_file 或 --input_dir 参数之一")
    if args.pdb_file and args.input_dir:
        parser.error("不能同时提供 --pdb_file 和 --input_dir 参数")

    # 确保输出目录存在
    os.makedirs(args.output_dir, exist_ok=True)

    if args.pdb_file:
        # 处理单个文件
        if not os.path.isfile(args.pdb_file):
             print(f"错误: 文件不存在 {args.pdb_file}")
             return
        base_name = os.path.basename(args.pdb_file)
        output_filename = os.path.splitext(base_name)[0] + '_pro.pdb'
        output_path = os.path.join(args.output_dir, output_filename)
        print(f"正在处理文件: {args.pdb_file} -> {output_path}")
        get_pro_pdb(args.pdb_file, output_path)
        print("处理完成.")

    elif args.input_dir:
        # 处理目录中的所有 PDB 文件
        if not os.path.isdir(args.input_dir):
             print(f"错误: 目录不存在 {args.input_dir}")
             return
        pdb_files = [f for f in os.listdir(args.input_dir) if f.endswith('.pdb')]
        if not pdb_files:
             print(f"警告: 在目录 {args.input_dir} 中未找到 PDB 文件")
             return

        print(f"在目录 {args.input_dir} 中找到 {len(pdb_files)} 个 PDB 文件，开始处理...")
        for filename in tqdm(pdb_files, desc="处理 PDB 文件"):
            input_path = os.path.join(args.input_dir, filename)
            output_filename = os.path.splitext(filename)[0] + '_pro.pdb'
            output_path = os.path.join(args.output_dir, output_filename)
            get_pro_pdb(input_path, output_path)
        print("所有文件处理完成.")


if __name__ == '__main__':
    main()