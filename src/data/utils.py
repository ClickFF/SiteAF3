"""
数据处理工具函数。
"""
import pickle
import os
import numpy as np
from typing import Dict, Any, List, Optional


def write_pkl(file_path: str, data_dict: Dict[str, Any]) -> None:
    """
    将数据字典写入pickle文件。
    
    Args:
        file_path: 文件路径
        data_dict: 要保存的数据字典
    """
    with open(file_path, 'wb') as f:
        pickle.dump(data_dict, f, protocol=4)


def read_pkl(file_path: str) -> Dict[str, Any]:
    """
    从pickle文件读取数据字典。
    
    Args:
        file_path: 文件路径
        
    Returns:
        数据字典
    """
    with open(file_path, 'rb') as f:
        return pickle.load(f)


def chain_str_to_int(chain_id: str) -> int:
    """
    将链ID字符串转换为整数表示。
    
    Args:
        chain_id: 链ID字符串
        
    Returns:
        链ID的整数表示
    """
    # 将链ID转换为ASCII值
    if len(chain_id) == 1:
        return ord(chain_id.upper())
    else:
        # 处理多字符链ID
        value = 0
        for i, c in enumerate(chain_id.upper()):
            value += ord(c) * (256 ** i)
        return value


def parse_chain_feats(
        chain_dict: Dict[str, Any],
        center_pos: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    解析链特征并居中到指定位置。
    
    Args:
        chain_dict: 链特征字典
        center_pos: 中心位置，默认为None
        
    Returns:
        处理后的链特征字典
    """
    # 如果给定了中心位置，则移动原子坐标
    if center_pos is not None:
        chain_dict['atom_positions'] = chain_dict['atom_positions'] - center_pos
    
    return chain_dict


def concat_np_features(
        feat_list: List[Dict[str, Any]],
        add_batch_dim: bool = True,
) -> Dict[str, Any]:
    """
    合并多个特征字典为一个。
    
    Args:
        feat_list: 特征字典列表
        add_batch_dim: 是否添加批次维度
        
    Returns:
        合并后的特征字典
    """
    concat_feats = {}
    # 合并所有特征
    for feat_name in feat_list[0].keys():
        if isinstance(feat_list[0][feat_name], np.ndarray):
            concat_feats[feat_name] = np.concatenate(
                [feats[feat_name] for feats in feat_list], axis=0)
        else:
            # 非数组类型的特征取第一个
            concat_feats[feat_name] = feat_list[0][feat_name]
    
    # 如果需要，添加批次维度
    if add_batch_dim:
        for feat_name, feat in concat_feats.items():
            if isinstance(feat, np.ndarray):
                concat_feats[feat_name] = feat[np.newaxis]
    
    return concat_feats 