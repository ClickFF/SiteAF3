"""
PDB结构解析器。
"""
import dataclasses
import numpy as np
from typing import List, Tuple, Dict, Any, Optional

from Bio.PDB import Chain
from src.data import residue_constants

@dataclasses.dataclass
class ChainFeatures:
    """PDB链的特征。"""
    aatype: np.ndarray  # 残基类型，形状为(残基数)
    atom_positions: np.ndarray  # 原子坐标，形状为(残基数, 原子数, 3)
    atom_mask: np.ndarray  # 原子掩码，形状为(残基数, 原子数)
    residue_index: np.ndarray  # 残基索引，形状为(残基数)
    chain_id: int  # 链ID的整数表示
    b_factors: np.ndarray  # B因子，形状为(残基数, 原子数)

def process_chain(chain: Chain, chain_id: int) -> ChainFeatures:
    """
    处理Bio.PDB.Chain并提取特征。
    
    Args:
        chain: Bio.PDB.Chain对象
        chain_id: 链ID的整数表示
        
    Returns:
        ChainFeatures对象
    """
    # 获取残基数量
    residues = list(chain.get_residues())
    n_res = len(residues)
    
    # 初始化特征数组
    aatype = np.zeros(n_res, dtype=np.int32)
    atom_positions = np.zeros((n_res, residue_constants.atom_type_num, 3), dtype=np.float32)
    atom_mask = np.zeros((n_res, residue_constants.atom_type_num), dtype=np.float32)
    residue_index = np.zeros(n_res, dtype=np.int32)
    b_factors = np.zeros((n_res, residue_constants.atom_type_num), dtype=np.float32)
    
    # 遍历所有残基
    for i, residue in enumerate(residues):
        # 记录残基索引
        residue_index[i] = residue.id[1]
        
        # 获取残基类型
        res_name = residue_constants.substitute_non_standard_restype.get(residue.resname, residue.resname)
        res_shortname = residue_constants.restype_3to1.get(res_name, 'X')
        aatype[i] = residue_constants.restype_order.get(res_shortname, residue_constants.restype_num)
        
        # 处理残基中的原子
        for atom in residue.get_atoms():
            atom_name = atom.name
            # 查找原子类型索引
            if atom_name in residue_constants.atom_order:
                atom_idx = residue_constants.atom_order[atom_name]
                atom_positions[i, atom_idx] = atom.coord
                atom_mask[i, atom_idx] = 1.0
                if hasattr(atom, 'bfactor'):
                    b_factors[i, atom_idx] = atom.bfactor
    
    return ChainFeatures(
        aatype=aatype,
        atom_positions=atom_positions,
        atom_mask=atom_mask,
        residue_index=residue_index,
        chain_id=chain_id,
        b_factors=b_factors
    ) 