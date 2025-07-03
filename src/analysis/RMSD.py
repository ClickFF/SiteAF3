#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# 读入mmcif文件，计算生成的结构中配体与真实结构的配体的RMSD
# 要求：
# 1.可以指定受体和配体的种类：蛋白质，DNA，RNA，其它
# 2. align受体
# 3. 计算配体与真实结构的配体的RMSD

import argparse
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# 忽略BioPython的PDB构造警告，例如不连续的末端等。
warnings.simplefilter('ignore', PDBConstructionWarning)

# 标准氨基酸和核酸残基名称（用于类型检查）
STANDARD_AA_NAMES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
    'TYR', 'VAL'
}
NUCLEIC_ACID_RESIDUES = {
    "A", "C", "G", "T", "U",  # RNA/DNA five main bases
    "DA", "DC", "DG", "DT",   # DNA specific
    "RA", "RC", "RG", "RU"    # RNA specific
}
# 核酸骨架原子（用于默认对齐）
NA_BACKBONE_ATOMS_DEFAULT = {"P", "O5\'", "C5\'", "C4\'", "C3\'", "O3\'"}


def get_corresponding_atoms(atoms_list_ref, atoms_list_model):
    """
    从两组原子列表中找出并返回对应的原子对。
    对应关系基于残基名、残基ID (hetflag, seqnum, icode) 和原子名。
    返回 (paired_ref_atoms, paired_model_atoms)
    """
    model_atoms_map = {}
    for atom_m in atoms_list_model:
        # 使用 (残基名, 残基ID元组, 原子名) 作为唯一键
        key = (atom_m.get_parent().get_resname(), 
               atom_m.get_parent().get_id(), 
               atom_m.get_name())
        model_atoms_map[key] = atom_m

    paired_ref_atoms = []
    paired_model_atoms = []

    for atom_r in atoms_list_ref:
        key = (atom_r.get_parent().get_resname(), 
               atom_r.get_parent().get_id(), 
               atom_r.get_name())
        if key in model_atoms_map:
            paired_ref_atoms.append(atom_r)
            paired_model_atoms.append(model_atoms_map[key])
            
    return paired_ref_atoms, paired_model_atoms

def get_atoms_from_selection(structure, chain_ids_str=None, molecule_type=None, atom_names_str="all"):
    """
    从结构中根据链ID、分子类型和原子名称选择原子。
    链ID和分子类型中至少一个应该被指定（由调用者保证）。
    """
    selected_atoms = []
    
    target_chain_ids = None
    if chain_ids_str:
        target_chain_ids = {c.strip() for c in chain_ids_str.split(',') if c.strip()}

    target_atom_names_set = set()
    if atom_names_str and atom_names_str.lower() != "all":
        target_atom_names_set = {a.strip().upper() for a in atom_names_str.split(',') if a.strip()}

    mol_type_lower = molecule_type.lower() if molecule_type else None

    for model_struct in structure: # Bio.PDB.Structure.Model object
        for chain in model_struct:
            # 1. 链ID筛选 (如果提供了链ID)
            if target_chain_ids and chain.id not in target_chain_ids:
                continue # 如果指定了目标链，但当前链不在其中，则跳过

            for residue in chain:
                resname_upper = residue.get_resname().strip().upper()
                
                keep_residue = False
                
                # 2. 分子类型筛选 (如果提供了分子类型)
                if mol_type_lower:
                    if mol_type_lower == "protein":
                        if resname_upper in STANDARD_AA_NAMES:
                            keep_residue = True
                    elif mol_type_lower in ["dna", "rna"]:
                        if resname_upper in NUCLEIC_ACID_RESIDUES:
                            keep_residue = True
                    elif mol_type_lower == "smallmolecule":
                        if residue.id[0].startswith('H_') and \
                           resname_upper not in STANDARD_AA_NAMES and \
                           resname_upper not in NUCLEIC_ACID_RESIDUES:
                            keep_residue = True
                    elif mol_type_lower == "other":
                        keep_residue = True # 'other' 类型表示接受该链/残基，不进行特定化学类型检查
                else:
                    # 如果没有指定分子类型 (mol_type_lower is None), 
                    # 那么只要链ID匹配 (或没有指定链ID，即选择所有链), 残基就应该被考虑
                    keep_residue = True 
                
                if keep_residue:
                    for atom in residue:
                        # 3. 原子名称筛选
                        if not target_atom_names_set or atom.name.strip().upper() in target_atom_names_set:
                            selected_atoms.append(atom)
    return selected_atoms


def calculate_rmsd_manual(atoms1, atoms2):
    """
    手动计算两组原子之间的RMSD。
    原子列表必须长度相同且原子一一对应。
    """
    if len(atoms1) != len(atoms2):
        raise ValueError(
            f"原子列表长度必须相同才能计算RMSD。列表1: {len(atoms1)}个原子, 列表2: {len(atoms2)}个原子。"
        )
    if not atoms1: # 两个列表都为空
        return 0.0
    
    diff = np.array([atom1.coord - atom2.coord for atom1, atom2 in zip(atoms1, atoms2)])
    return np.sqrt(np.sum(diff * diff) / len(atoms1))


def main():
    parser = argparse.ArgumentParser(
        description="计算两个分子结构（参考和模型）中配体的RMSD，首先对齐它们的受体。",
        formatter_class=argparse.RawTextHelpFormatter
    )

    req_group = parser.add_argument_group('必要文件参数')
    req_group.add_argument("--ref_file", required=True, help="参考结构mmCIF文件路径。")
    req_group.add_argument("--model_file", required=True, help="模型结构mmCIF文件路径。")

    sel_group = parser.add_argument_group('受体和配体选择参数 (每个实体至少指定链或类型之一)')
    sel_group.add_argument("--receptor-chains", dest="receptor_chains", help="受体的链ID (逗号分隔, 例如: A,B)。应用于参考和模型结构。")
    sel_group.add_argument("--receptor-type", dest="receptor_type", choices=["protein", "dna", "rna", "other"], help="受体的分子类型。应用于参考和模型结构。")
    sel_group.add_argument("--ligand-chains", dest="ligand_chains", help="配体的链ID (逗号分隔)。应用于参考和模型结构。")
    sel_group.add_argument("--ligand-type", dest="ligand_type", choices=["protein", "dna", "rna", "smallmolecule", "other"], help="配体的分子类型。应用于参考和模型结构。")

    opt_group = parser.add_argument_group('原子名称可选参数')
    opt_group.add_argument("--receptor-align-atoms", dest="receptor_align_atoms", default=None,
                           help="用于受体对齐的原子名称 (逗号分隔, 例如: CA,C,N)。\\n"
                                "默认为: 'CA' (protein), 'P,O5\\\'\',C5\\\'\',C4\\\'\',C3\\\'\',O3\\\'\'\' (dna/rna), 'all' (other 或未指定受体类型时)。")
    opt_group.add_argument("--ligand-rmsd-atoms", dest="ligand_rmsd_atoms", default="all",
                           help="用于配体RMSD计算的原子名称 (逗号分隔, 例如: C1,C2,N1 或 'all')。\\n"
                                "默认为: 'all'。")

    args = parser.parse_args()

    # 验证选择参数：对于受体和配体，必须提供链或类型中的至少一个
    if not args.receptor_chains and not args.receptor_type:
        parser.error("对于受体选择，必须通过 --receptor-chains 或 --receptor-type 指定选择标准。")
    if not args.ligand_chains and not args.ligand_type:
        parser.error("对于配体选择，必须通过 --ligand-chains 或 --ligand-type 指定选择标准。")

    # 为receptor_align_atoms设置智能默认值
    effective_receptor_align_atoms = args.receptor_align_atoms
    if effective_receptor_align_atoms is None: # 用户未指定
        if args.receptor_type: # 仅当提供了受体类型时，才使用基于类型的默认值
            receptor_type_lower = args.receptor_type.lower()
            if receptor_type_lower == "protein":
                effective_receptor_align_atoms = "CA"
            elif receptor_type_lower in ["dna", "rna"]:
                effective_receptor_align_atoms = ",".join(NA_BACKBONE_ATOMS_DEFAULT)
            else: # 'other' 类型
                effective_receptor_align_atoms = "all"
        else: # 如果没有提供受体类型，则默认为 "all"
            effective_receptor_align_atoms = "all"
    
    # cif_parser = MMCIFParser(QUIET=True) # Will be determined dynamically

    def get_structure_parser(filename):
        if filename.lower().endswith(('.cif', '.mmcif')):
            return MMCIFParser(QUIET=True)
        elif filename.lower().endswith('.pdb'):
            return PDBParser(PERMISSIVE=1) # PERMISSIVE=1 is good for PDB files
        else:
            raise ValueError(f"不支持的文件格式: {filename}. 请使用 .pdb, .cif, 或 .mmcif 文件扩展名。")

    try:
        print(f"加载参考结构: {args.ref_file}")
        ref_parser = get_structure_parser(args.ref_file)
        ref_struct = ref_parser.get_structure("reference", args.ref_file)
        
        print(f"加载模型结构: {args.model_file}")
        model_parser = get_structure_parser(args.model_file)
        model_struct = model_parser.get_structure("model", args.model_file)
    except Exception as e:
        print(f"错误：无法加载结构文件。 {e}")
        return

    print("\n--- 原子选择 ---")

    def format_selection_criteria(chains, m_type):
        parts = []
        if chains:
            parts.append(f"链 '{chains}'")
        if m_type:
            parts.append(f"类型 '{m_type}'")
        return " 和 ".join(parts) if parts else "未指定 (将尝试选择所有原子，依赖原子名)"

    receptor_criteria_desc = format_selection_criteria(args.receptor_chains, args.receptor_type)
    print(f"受体选择标准: {receptor_criteria_desc}")
    print(f"  用于对齐的原子: '{effective_receptor_align_atoms}'")

    ref_receptor_atoms = get_atoms_from_selection(
        ref_struct, args.receptor_chains, args.receptor_type, effective_receptor_align_atoms
    )
    print(f"  参考结构中选中的受体原子数: {len(ref_receptor_atoms)}")
    
    model_receptor_atoms = get_atoms_from_selection(
        model_struct, args.receptor_chains, args.receptor_type, effective_receptor_align_atoms
    )
    print(f"  模型结构中选中的受体原子数: {len(model_receptor_atoms)}")

    ligand_criteria_desc = format_selection_criteria(args.ligand_chains, args.ligand_type)
    print(f"配体选择标准: {ligand_criteria_desc}")
    print(f"  用于RMSD的原子: '{args.ligand_rmsd_atoms}'")

    ref_ligand_atoms_initial = get_atoms_from_selection(
        ref_struct, args.ligand_chains, args.ligand_type, args.ligand_rmsd_atoms
    )
    print(f"  参考结构中初始选中的配体原子数: {len(ref_ligand_atoms_initial)}")

    model_ligand_atoms_initial = get_atoms_from_selection(
        model_struct, args.ligand_chains, args.ligand_type, args.ligand_rmsd_atoms
    )
    print(f"  模型结构中初始选中的配体原子数: {len(model_ligand_atoms_initial)}")

    # 检查原子选择结果 (受体部分)
    if not ref_receptor_atoms or not model_receptor_atoms:
        print("错误：受体原子选择导致参考或模型中的原子列表为空。请检查链ID、类型和原子名称。")
        return
    if len(ref_receptor_atoms) != len(model_receptor_atoms):
        print(f"错误：受体对齐的原子数量不匹配。参考: {len(ref_receptor_atoms)}, 模型: {len(model_receptor_atoms)}。对齐需要数量相等的对应原子。")
        return

    # 处理配体原子并找到共同原子
    common_ref_ligand_atoms = []
    common_model_ligand_atoms = []
    num_common_ligand_atoms = 0
    can_calculate_ligand_rmsd = True

    if not ref_ligand_atoms_initial and not model_ligand_atoms_initial:
        print("  警告：参考和模型中的配体均未选中任何原子。")
        can_calculate_ligand_rmsd = False
    elif not ref_ligand_atoms_initial or not model_ligand_atoms_initial:
        print(f"  错误：配体原子初始选择导致参考 ({len(ref_ligand_atoms_initial)} 个原子) 或模型 ({len(model_ligand_atoms_initial)} 个原子) 中的列表有一个为空。")
        can_calculate_ligand_rmsd = False
    else:
        common_ref_ligand_atoms, common_model_ligand_atoms = get_corresponding_atoms(
            ref_ligand_atoms_initial, model_ligand_atoms_initial
        )
        num_common_ligand_atoms = len(common_ref_ligand_atoms)
        print(f"  在初始选择的配体原子之间，找到 {num_common_ligand_atoms} 对共同原子用于RMSD计算。")
        if num_common_ligand_atoms == 0:
            print("  警告：未找到共同的配体原子。")
            can_calculate_ligand_rmsd = False

    print("\n--- 受体对齐 ---")
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_receptor_atoms, model_receptor_atoms) # fixed (ref), moving (model)
    
    # 将变换应用于整个模型结构中的所有原子
    # 这会就地更新 model_struct 中所有原子的坐标
    super_imposer.apply(model_struct.get_atoms()) 
    
    print(f"受体对齐 RMSD (基于 {len(ref_receptor_atoms)} 对原子): {super_imposer.rms:.4f} Å")

    # model_ligand_atoms 列表中的原子对象的坐标已被 super_imposer.apply 更新
    # 无需重新从 model_struct 获取
    # common_model_ligand_atoms 包含的是 model_ligand_atoms_initial 中的原子引用，
    # 而 model_ligand_atoms_initial 的原子来自 model_struct，所以其坐标也已更新。

    print("\n--- 配体 RMSD 计算 ---")
    ligand_rmsd_value = 0.0  # 使用不同的变量名以避免与args中的同名参数混淆

    if can_calculate_ligand_rmsd and num_common_ligand_atoms > 0:
        try:
            ligand_rmsd_value = calculate_rmsd_manual(common_ref_ligand_atoms, common_model_ligand_atoms)
            print(f"配体 RMSD (基于 {num_common_ligand_atoms} 对共同原子): {ligand_rmsd_value:.4f} Å")
        except ValueError as e:
            print(f"错误：计算配体RMSD失败。{e}")
            # ligand_rmsd_value 保持 0.0
    elif not can_calculate_ligand_rmsd:
        print(f"配体 RMSD: 由于原子选择问题或未找到共同原子，RMSD计算跳过，设为 0.0 Å。")
        # ligand_rmsd_value 保持 0.0
    else: # num_common_ligand_atoms is 0 (and can_calculate_ligand_rmsd might have been true initially but then set to false)
        print(f"配体 RMSD (基于 0 对共同原子): 0.0000 Å")
        # ligand_rmsd_value 保持 0.0

if __name__ == "__main__":
    main()
