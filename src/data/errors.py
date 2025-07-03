"""
定义处理PDB数据时可能遇到的错误类型。
"""

class DataError(Exception):
    """处理数据时发生的一般错误。"""
    pass

class LengthError(DataError):
    """序列长度相关的错误。"""
    pass

class ChainError(DataError):
    """链ID或链结构相关的错误。"""
    pass

class ResidueError(DataError):
    """残基类型或结构相关的错误。"""
    pass

class FormatError(DataError):
    """文件格式相关的错误。"""
    pass 