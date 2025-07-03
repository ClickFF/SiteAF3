"""
定义与残基相关的常量。
"""
import numpy as np
# 从 alphafold3.constants.residue_names 导入必要的常量
from alphafold3.constants import residue_names

# 标准氨基酸残基的单字母代码
restype_1to3 = residue_names.PROTEIN_COMMON_ONE_TO_THREE

# 核酸残基的单字母代码
# residue_names.DNA_COMMON_ONE_TO_TWO maps 'A':'DA', 'T':'DT', 'G':'DG', 'C':'DC'
dna_1to3 = residue_names.DNA_COMMON_ONE_TO_TWO
# residue_constants.py 原来有 'U': 'DU'。
# AlphaFold 3 将 U 视为 RNA (residue_names.RNA_TYPES)，DT/T 视为 DNA。
# CCD_NAME_TO_ONE_LETTER 包含 '0AU': 'U', '2AU': 'U' 等，将 'DU' 映射回 'U' 的情况较少。
# 为了与 AlphaFold 3 的定义更一致，这里不显式添加 'U':'DU'。

# 标准氨基酸残基的三字母代码到单字母代码的映射
# 以及广泛的CCD名称到单字母代码的映射
restype_3to1 = dict(residue_names.CCD_NAME_TO_ONE_LETTER)
# 确保residue_names中定义的DNA和RNA类型（作为键）正确映射到其单字母形式
for dna_two_letter_val in residue_names.DNA_TYPES:  # e.g., DA, DG, DC, DT
    one_letter = residue_names.letters_three_to_one(dna_two_letter_val, default='X')
    if one_letter != 'X':
        restype_3to1[dna_two_letter_val] = one_letter
for rna_one_letter_val in residue_names.RNA_TYPES:  # e.g., A, G, C, U
    restype_3to1[rna_one_letter_val] = rna_one_letter_val


# 标准氨基酸残基的列表（三字母代码）
amino_acid_3_list = list(residue_names.PROTEIN_TYPES)

# 核酸残基的列表（residue_names.py中RNA是单字母, DNA是双字母）
nucleic_acid_3_list = list(residue_names.DNA_TYPES) + list(residue_names.RNA_TYPES)

# 标准氨基酸残基的集合
amino_acid_3_set = set(amino_acid_3_list)

# 核酸残基的集合
nucleic_acid_3_set = set(nucleic_acid_3_list)

# 核酸残基名称到单字母代码的映射
_NUCLEIC_ONE_LETTERS_FOR_MAPPING = set(residue_names.DNA_TYPES_ONE_LETTER_WITH_UNKNOWN) | {'U'} | set(residue_names.RNA_TYPES)
resname_to_nucleic_acid_1letter = {
    name: one_letter
    for name, one_letter in residue_names.CCD_NAME_TO_ONE_LETTER.items()
    if one_letter in _NUCLEIC_ONE_LETTERS_FOR_MAPPING
}
# 确保标准DNA双字母和RNA单字母的映射存在
for dna_dl in residue_names.DNA_TYPES:
    resname_to_nucleic_acid_1letter[dna_dl] = residue_names.letters_three_to_one(dna_dl, default=residue_names.UNK_NUCLEIC_ONE_LETTER)
for rna_sl in residue_names.RNA_TYPES:
    resname_to_nucleic_acid_1letter[rna_sl] = rna_sl
# 'DU' 和 'DI' 等特殊情况：
# residue_names.CCD_NAME_TO_ONE_LETTER 包含 '0AU': 'U' (代表Uridine), '2AU': 'U'
# 不直接有 'DU': 'U'。 如果需要，可以手动添加。
if 'DU' in residue_names.CCD_NAME_TO_ONE_LETTER:
    resname_to_nucleic_acid_1letter['DU'] = residue_names.CCD_NAME_TO_ONE_LETTER['DU']
elif 'U' in _NUCLEIC_ONE_LETTERS_FOR_MAPPING: # fallback if DU not in CCD but U is a target
    resname_to_nucleic_acid_1letter['DU'] = 'U'

if 'DI' in residue_names.CCD_NAME_TO_ONE_LETTER: # Deoxyriboinosine
     resname_to_nucleic_acid_1letter['DI'] = residue_names.CCD_NAME_TO_ONE_LETTER['DI']


def is_nucleic_acid(resname):
    """
    检查残基名称是否为核酸残基
    
    Args:
        resname: 残基名称（通常是三字母或双字母代码）
        
    Returns:
        bool: 如果是核酸残基则为True，否则为False
    """
    # 检查是否在预定义的核酸符号集合中（例如 'DA', 'A'）
    if resname in nucleic_acid_3_set:
        return True
    # 检查是否在更广泛的名称到单字母映射中，并且映射为已知的核酸单字母
    if resname in resname_to_nucleic_acid_1letter:
        # letters_three_to_one(resname, default='?') in _NUCLEIC_ONE_LETTERS_FOR_MAPPING
        # This check is implicitly covered if resname_to_nucleic_acid_1letter is built correctly
        return True
    return False

def is_aa(resname):
    """
    检查残基名称是否为氨基酸残基
    
    Args:
        resname: 残基名称（三字母代码）
        
    Returns:
        bool: 如果是氨基酸残基则为True，否则为False
    """
    # 标准氨基酸残基
    if resname in amino_acid_3_set: # amino_acid_3_set contains standard 3-letter codes like ALA, ARG
        return True
    
    # 尝试通过 residue_names.letters_three_to_one 将 resname 转换为单字母代码
    # 并检查该单字母是否为已知的蛋白质单字母代码
    one_letter = residue_names.letters_three_to_one(resname, default='_') # Use a default not in protein letters
    if one_letter in residue_names.PROTEIN_TYPES_ONE_LETTER:
        return True
        
    # 非标准残基，但仍被处理为氨基酸 (根据原有的判断逻辑，可以保留作为补充)
    # residue_names.CCD_NAME_TO_ONE_LETTER 应该已经覆盖了这些
    # 例如 'MSE' -> 'M', 'HYP' -> 'P'
    # if resname in ['MSE', 'HYP', 'TPO', 'SEP', 'PTR', 'CSO', 'CSD', 'CME']:
    #     return True
        
    # 原有的通用三字母判断规则（作为最后的补充）
    if len(resname) == 3 and resname.isalpha() and resname.upper() == resname:
        # 排除已明确为核酸的
        if not is_nucleic_acid(resname):
            # 排除一些常见的辅因子和小分子 (可以从residue_names.py中的非聚合物类型获取更权威的列表，但暂时保留)
            if resname not in ['HOH', 'WAT', 'H2O', 'NAG', 'MAN', 'BMA', 'FUC', 'GDP', 'GTP', 'ADP', 'ATP', 'NAD', 'FAD', residue_names.UNL]:
                # 如果CCD能把它映射到一个蛋白质单字母，那它就是蛋白质
                if residue_names.letters_three_to_one(resname, default='_') in residue_names.PROTEIN_TYPES_ONE_LETTER:
                    return True
    
    return False

# 残基类型的顺序（索引到one-hot编码）
# 使用 AlphaFold 3 的 POLYMER_TYPES_ORDER_WITH_ALL_UNKS_AND_GAP 作为权威顺序
# 它的键是聚合物的标准表示 (ALA, A (RNA), DA (DNA), UNK, GAP 等)
restype_order = residue_names.POLYMER_TYPES_ORDER_WITH_ALL_UNKS_AND_GAP

# 残基类型总数
restype_num = residue_names.POLYMER_TYPES_NUM_ORDER_WITH_ALL_UNKS_AND_GAP

# 核酸残基集合，用于判断残基类型 (供第二个 is_nucleic_acid 函数使用)
# 从 CCD_NAME_TO_ONE_LETTER 中提取所有值为核酸单字母 ('A', 'C', 'G', 'T', 'U', 'N') 的键
nucleic_acid_residues = {
    name for name, one_letter in residue_names.CCD_NAME_TO_ONE_LETTER.items()
    if one_letter in _NUCLEIC_ONE_LETTERS_FOR_MAPPING
}
# 确保标准DNA和RNA类型包含在内
nucleic_acid_residues.update(residue_names.DNA_TYPES)
nucleic_acid_residues.update(residue_names.RNA_TYPES)
# 添加一些原始文件中存在的特殊核酸（如果CCD中没有）
# e.g. 'ADE', 'THY', 'GUA', 'CYT', 'URA'
# 这些通常会被CCD映射，例如 'ADE' -> 'A'
for common_na_name in ['ADE', 'THY', 'GUA', 'CYT', 'URA']:
    if common_na_name in residue_names.CCD_NAME_TO_ONE_LETTER and \
       residue_names.CCD_NAME_TO_ONE_LETTER[common_na_name] in _NUCLEIC_ONE_LETTERS_FOR_MAPPING:
        nucleic_acid_residues.add(common_na_name)


def is_nucleic_acid(resname): # Note: This redefines the previous one.
    """
    判断残基是否为核酸残基
    
    Args:
        resname: 残基的三字母代码 (或双字母/单字母)
        
    Returns:
        True 如果是核酸残基，False 如果不是
    """
    # 检查残基名称是否在广泛的核酸残基集合中
    if resname in nucleic_acid_residues:
        return True
    
    # 检查残基名称是否以D开头（DNA残基常见命名方式, 作为补充规则）
    # 但要小心如 ASP, DNS 等非核酸。CCD覆盖应更准确。
    # if resname.startswith('D') and len(resname) > 1 and len(resname) <=3:
    #    one_letter_type = residue_names.letters_three_to_one(resname, default='_')
    #    if one_letter_type in ('A', 'G', 'C', 'T', 'N'): # Check if it's a DNA type
    #        return True
            
    return False

# 非标准残基到标准残基的映射
# 目标: 非标准名称 -> 标准聚合物名称 (ALA, MET, DA, A, UNK, N, DN, UNL)
substitute_non_standard_restype = {}
_all_standard_forms_for_subst = set(residue_names.PROTEIN_TYPES) | \
                               set(residue_names.RNA_TYPES) | \
                               set(residue_names.DNA_TYPES) | \
                               {residue_names.UNK, residue_names.UNK_RNA, residue_names.UNK_DNA, 
                                residue_names.UNL, residue_names.GAP}

for ccd_name, one_letter in residue_names.CCD_NAME_TO_ONE_LETTER.items():
    if ccd_name in _all_standard_forms_for_subst:
        continue

    standard_form = None
    if one_letter in residue_names.PROTEIN_COMMON_ONE_TO_THREE:
        standard_form = residue_names.PROTEIN_COMMON_ONE_TO_THREE[one_letter]
    elif one_letter in residue_names.DNA_COMMON_ONE_TO_TWO: # one_letter is A,T,G,C for DNA
        standard_form = residue_names.DNA_COMMON_ONE_TO_TWO[one_letter]
    elif one_letter in residue_names.RNA_TYPES: # one_letter is A,U,G,C for RNA
        standard_form = one_letter # Standard form is the one-letter itself
    elif one_letter == 'X': # Protein unknown
        standard_form = residue_names.UNK
    elif one_letter == residue_names.UNK_NUCLEIC_ONE_LETTER: # 'N'
        # Heuristic: if ccd_name starts with 'D' or contains 'deoxy', maybe DNA?
        # This is hard. Default to UNK_RNA ('N') as it's a single letter 'N'.
        # AlphaFold 3 uses UNK_RNA = 'N' and UNK_DNA = 'DN'.
        # If ccd_name is 'DN' or looks like DNA, map to UNK_DNA.
        if ccd_name == 'DN' or 'D' == ccd_name[0] and ccd_name[1:] in residue_names.RNA_TYPES : # e.g. DA, DG, DC, DU (if U is RNA type)
             standard_form = residue_names.UNK_DNA # 'DN'
        else:
             standard_form = residue_names.UNK_RNA # 'N'
    # else: could be a ligand not mapping to UNL via CCD, or other unhandled cases.

    if standard_form and ccd_name != standard_form:
        substitute_non_standard_restype[ccd_name] = standard_form

# Ensure 'UNK' itself (if considered non-standard by some tool expecting 'X') maps to residue_names.UNK
# However, the loop already skips if ccd_name is in _all_standard_forms_for_subst.
# residue_names.UNK ('UNK') is standard.
# The original dict had 'UNK': 'X'. Our target standard_form for 'X' is residue_names.UNK ('UNK').
# So, if some `some_unk_code` maps to 'X' via CCD, then `substitute_non_standard_restype[some_unk_code] = residue_names.UNK`.
# This is more consistent with using full polymer names as targets.

# 从 alphafold3.constants 导入 atom_types
from alphafold3.constants import atom_types

# 合法原子列表基于 alphafold3.constants.atom_types.ATOM37 和 atom_types.ATOM29 构建
_protein_atoms = list(atom_types.ATOM37)
_nucleic_atoms = list(atom_types.ATOM29)
_combined_atom_list = _protein_atoms + [atom for atom in _nucleic_atoms if atom not in _protein_atoms]
atom_order = {name: i for i, name in enumerate(_combined_atom_list)}
atom_type_num = len(_combined_atom_list)