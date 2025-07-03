# AlphaFold3嵌入生成工具

## 概述
本工具用于使用AlphaFold3模型生成蛋白质和核酸分子的**序列嵌入**或**结构嵌入**表示。这些嵌入可用于后续的机器学习任务，例如核酸适配体设计、分子对接预测和结构分析。
注意：如果输入PDB文件是经过预处理的（例如，仅包含口袋和配体），则生成的嵌入将反映该预处理后的结构。

## 功能特点
- 直接从PDB文件生成**序列嵌入**（包含单残基和配对嵌入）
- 生成**结构嵌入**，额外包含原子坐标和原子掩码
- 支持蛋白质、RNA和DNA混合结构
- 自动检测链类型（蛋白质、RNA或DNA）
- 可选使用AlphaFold3官方MSA搜索工具构建高质量MSA
- 可选生成**固定掩码 (fixed mask)**，用于在生成任务（如diffusion模型）中固定分子结构的特定部分（例如一条链）
- 提供详细的嵌入统计信息
- （注意：批量处理PKL文件的嵌入生成逻辑当前未完全实现）

## 安装要求
- Python 3.11+
- AlphaFold3环境 (包含相关依赖)
- JAX（GPU版本推荐）
- BioPython
- NumPy
- Pandas
- MDTraj
- Haiku
- （可选）MSA搜索工具（jackhmmer, nhmmer, hmmalign, hmmsearch, hmmbuild）和相应数据库

## 使用方法

### 单个PDB文件生成序列嵌入
```bash
# 默认生成序列嵌入
python -m embeddings.embed_af3 --pdb_file /path/to/structure.pdb --output_dir /path/to/output --weight_dir /path/to/af3/weights --verbose
```

### 生成结构嵌入（包含坐标和可选的固定掩码）
```bash
# 生成结构嵌入，但不固定任何部分
python -m embeddings.embed_af3 --pdb_file /path/to/structure.pdb --output_dir /path/to/output --weight_dir /path/to/af3/weights --generate_struct_emb --verbose

# 生成结构嵌入，并固定链A
python -m embeddings.embed_af3 --pdb_file /path/to/structure.pdb --output_dir /path/to/output --weight_dir /path/to/af3/weights --generate_struct_emb --fixed_chain A --verbose
```

### 使用MSA生成嵌入
```bash
# 生成序列嵌入，并使用MSA（需要配置数据库）
python -m embeddings.embed_af3 --pdb_file /path/to/structure.pdb --output_dir /path/to/output --weight_dir /path/to/af3/weights --use_msa --db_dir /path/to/databases --verbose

# 生成结构嵌入，固定链A，并使用MSA
python -m embeddings.embed_af3 --pdb_file /path/to/structure.pdb --output_dir /path/to/output --weight_dir /path/to/af3/weights --generate_struct_emb --fixed_chain A --use_msa --db_dir /path/to/databases --verbose
```

### （实验性）批量处理PKL文件目录
```bash
# 注意：此功能当前主要用于文件复制/加载，嵌入逻辑需自行添加
python -m embeddings.embed_af3 --pkl_dir /path/to/pkl/files --output_dir /path/to/output --weight_dir /path/to/af3/weights --verbose
```

## 参数说明
- `--pdb_file`: （必需，除非提供`--pkl_dir`）待处理的PDB文件路径。
- `--pkl_dir`: （可选）待处理的PKL文件目录（嵌入逻辑未实现）。
- `--output_dir`: 输出目录路径（默认为"/work/hat170/aptamer/test_output"）。
- `--weight_dir`: AlphaFold3模型权重目录路径（默认为"/work/hat170/aptamer/alphafold3/weight"）。
- `--generate_struct_emb`: （可选）生成结构嵌入（包含坐标和原子掩码），而不是仅生成序列嵌入。如果设置，输出文件名将包含 `_struct_emb`。
- `--fixed_chain`: （可选，仅与`--generate_struct_emb`一起使用）要固定的链的ID。如果提供，该链的所有标准残基将在`fixed_mask`中标记为1。
- `--db_dir`: （可选，与`--use_msa`一起使用）MSA数据库目录路径。
- `--msa_dir`: （可选，与`--use_msa`一起使用）预计算的MSA文件目录路径。如果提供，将跳过MSA计算。
- `--use_msa`: （可选）是否使用MSA特征生成嵌入（默认为False）。如果设置为True但未提供`--msa_dir`，则需要`--db_dir`和已安装的MSA工具。
- `--verbose`: （可选）显示详细输出。

## 输出文件结构

输出为 Pickle (`.pkl`) 文件。

### 序列嵌入输出 (`*_seq_emb.pkl`)
```python
{
  'embeddings': {
    'single': numpy.ndarray(shape=(num_tokens, ...)),  # 单残基/token嵌入
    'pair': numpy.ndarray(shape=(num_tokens, num_tokens, ...)),  # 配对嵌入
    'target_feat': numpy.ndarray(...) # 目标特征
  },
  'seq_info': {
    'name': 'structure_name',
    'chains': [
      {
        'chain_id': 'A',
        'sequence': 'SEQUENCE...',
        'chain_type': 'protein/rna/dna/...'
      },
      # ... 其他链信息
    ]
  },
  'input_features': ...  # 用于模型输入的原始特征字典
}
```

### 结构嵌入输出 (`*_struct_emb.pkl`)
结构嵌入输出包含序列嵌入的所有内容，并额外添加了结构坐标、原子掩码和固定信息：
```python
{
  'embeddings': {
    'single': numpy.ndarray(shape=(num_tokens, ...)),
    'pair': numpy.ndarray(shape=(num_tokens, num_tokens, ...)),
    'target_feat': numpy.ndarray(...),
    'structure_atom_coords': numpy.ndarray(shape=(num_tokens, MAX_ATOMS_PER_TOKEN, 3)), # 原子坐标 (根据AF3原子顺序)
    'structure_atom_mask': numpy.ndarray(shape=(num_tokens, MAX_ATOMS_PER_TOKEN)) # 原子存在掩码 (1表示存在, 0表示不存在)
  },
  'fixed_info': {
    'fixed_chain_id': 'A' or None,  # 请求固定的链ID
    'fixed_residues': [(chain_id, res_id), ...],  # 最终被映射并用于固定的PDB残基列表
    'fixed_mask': numpy.ndarray(shape=(num_tokens,)),  # 固定部分的二进制掩码 (1表示固定, 0表示可变)
    'fixed_embed_indices': [idx1, idx2, ...] # fixed_mask中为1的嵌入索引列表
  },
  'seq_info': { ... }, # 同序列嵌入输出
  'input_features': { ... } # 同序列嵌入输出
}
```

## 结构嵌入与固定掩码功能说明
当使用 `--generate_struct_emb` 标志时，工具会计算并包含额外的结构信息：
1.  **`structure_atom_coords`**: 每个 token（通常对应一个残基）的原子坐标。坐标按照 AlphaFold3 定义的原子顺序排列（最多 `MAX_ATOMS_PER_TOKEN` 个原子）。
2.  **`structure_atom_mask`**: 一个掩码，指示 `structure_atom_coords` 中的哪些原子坐标是有效的（即原子存在于输入结构中）。

此外，可以通过 `--fixed_chain` 参数指定一个链ID。如果指定：
1.  该链的所有标准残基将被识别。
2.  这些残基在最终输出的 `fixed_mask` 中对应的位置会被标记为 `1`。
3.  `fixed_info` 字典会记录请求的链ID、实际固定的残基列表（PDB ID）、固定掩码本身以及掩码中为1的嵌入索引。

**`fixed_mask` 的用途**：
此掩码主要用于下游的生成任务，特别是基于diffusion的模型。通过将此掩码应用于模型的特定层或损失函数，可以强制模型在生成新结构（例如，设计新的适配体链）时保持指定部分（例如蛋白质受体链）的结构不变。

**适用场景**：
- **条件化分子生成**：在生成新的分子（如RNA、DNA、肽）时，将已知结构（如蛋白质受体）固定。
- **结构修复或编辑**：保持分子的大部分结构固定，只允许特定区域变化。
- **研究分子相互作用**：分析固定一部分结构对另一部分结构或相互作用的影响。

## 预处理文件说明（可选）
如果输入 PDB 文件是经过预处理的（例如，使用 `NAcutoff.py` 提取了口袋区域），本工具仍然可以处理。生成的嵌入将反映这个预处理后的结构。如果使用了 `--fixed_chain`，它将固定预处理后 PDB 文件中存在的指定链的残基。

## 错误处理
此工具会在以下情况下抛出错误：
- **权重加载错误**：无法找到或加载模型权重 (`FileNotFoundError`, `RuntimeError`)。
- **文件不存在错误**：指定的 PDB 文件不存在 (`FileNotFoundError`)。
- **PDB解析错误**：无法解析输入 PDB 文件 (`ValueError`)。
- **特征生成错误**：无法从 PDB 文件或序列创建模型输入特征 (`ValueError`)。
- **MSA工具/数据库错误**：如果`use_msa`为True，但相关工具或数据库缺失/配置错误 (`RuntimeError`)。
- **GPU内存错误**：GPU内存不足 (`jaxlib.xla_extension.XlaRuntimeError: RESOURCE_EXHAUSTED`)。
- **API兼容性错误**：如果 AlphaFold3 或其依赖库 API 发生变化，可能导致不兼容 (`AttributeError`, `TypeError` 等)。

错误信息会尽量提供上下文，便于调试。

## 注意事项
1.  **模型权重**: 确保 AlphaFold3 模型权重已正确下载并放置在 `--weight_dir` 指定的目录。
2.  **内存需求**: 生成嵌入（特别是使用MSA时）可能需要大量 CPU 和 GPU 内存。
3.  **API兼容性**: 代码基于特定版本的 AlphaFold3 API。未来 API 更新可能需要代码调整。
4.  **固定掩码准确性**: `fixed_mask` 的生成依赖于 PDB 文件结构与 AlphaFold3 内部处理后序列/token 的正确映射。如果 PDB 文件复杂或存在不标准残基，映射可能不完美，建议检查 `fixed_info` 中的 `fixed_embed_indices` 和 `fixed_residues`。
5.  **PKL 处理**: `--pkl_dir` 功能目前不执行嵌入生成，主要用于文件管理。

## 工作流程
1.  解析命令行参数。
2.  初始化 `Embedder` 对象，加载 AlphaFold3 模型权重。
3.  **如果处理 PDB 文件**:
    a.  调用 `seq_emb_af3` 生成基础序列嵌入（包括 MSA 处理，如果请求）。
    b.  **如果请求了结构嵌入 (`--generate_struct_emb`)**:
        i.  调用 `struct_emb_af3`。
        ii. 解析 PDB 获取原子坐标。
        iii. 根据 `--fixed_chain` 或（未来可能实现的）`--fixed_residues` 确定要固定的残基。
        iv. 创建 PDB 残基 ID 到嵌入索引的映射。
        v.  生成 `structure_atom_coords`, `structure_atom_mask` 和 `fixed_mask`。
        vi. 将结构和固定信息添加到结果字典中。
    c.  保存结果到 `.pkl` 文件（文件名反映是序列嵌入还是结构嵌入）。
4.  **如果处理 PKL 目录** (当前逻辑):
    a.  遍历目录中的 `.pkl` 文件。
    b.  加载每个文件并（可选地）执行某些操作（当前仅为加载/保存）。
    c.  保存文件到输出目录。

## 内存优化建议
- 对于小型GPU（≤16GB内存）：优先处理单链或较小复合物；考虑禁用MSA（`--use_msa False`）；调整JAX内存分配（`XLA_PYTHON_CLIENT_MEM_FRACTION`）。
- 对于中型GPU（24-32GB内存）：可以处理中等大小的复合物和MSA。

## 常见问题
- **CUDA/JAX错误**: 确认 JAX, CUDA, cuDNN 版本兼容，GPU驱动最新，内存足够。
- **模型权重错误**: 检查 `--weight_dir` 路径是否正确且包含所有权重文件。
- **API兼容性错误**: 参考错误信息和最新的 AlphaFold3 文档更新代码。
- **固定掩码不准确**: 检查 PDB 文件质量和 `--fixed_chain` ID 是否正确。查看 verbose 输出中的警告信息。

## 示例运行
```bash
# 生成序列嵌入 (输出: ..._seq_emb.pkl)
python -m embeddings.embed_af3 --pdb_file /work/hat170/aptamer/test_input/1A1V_processed.pdb --output_dir /work/hat170/aptamer/test_output --verbose

# 生成结构嵌入，不固定 (输出: ..._struct_emb.pkl)
python -m embeddings.embed_af3 --pdb_file /work/hat170/aptamer/test_input/1A1V_processed.pdb --output_dir /work/hat170/aptamer/test_output --generate_struct_emb --verbose

# 生成结构嵌入，固定链 A (输出: ..._struct_emb.pkl, fixed_mask 将包含链A的残基)
python -m embeddings.embed_af3 --pdb_file /work/hat170/aptamer/test_input/1A1V_processed.pdb --output_dir /work/hat170/aptamer/test_output --generate_struct_emb --fixed_chain A --verbose

# 生成结构嵌入，固定链 A，并使用MSA
python -m embeddings.embed_af3 --pdb_file /work/hat170/aptamer/test_input/1A1V_processed.pdb --output_dir /work/hat170/aptamer/test_output --generate_struct_emb --fixed_chain A --use_msa --db_dir /path/to/databases --verbose
```

## MSA处理说明
本工具**可选地**利用AlphaFold3官方的数据管道功能为蛋白质和RNA链生成高质量多序列比对(MSA)：

1. **蛋白质链MSA与模板搜索**：使用AlphaFold3数据管道自动运行：
   - 使用Jackhmmer搜索UniRef90、Mgnify和Small BFD数据库生成unpaired MSA
   - 使用Jackhmmer搜索UniProt数据库生成paired MSA
   - 使用Hmmsearch进行模板搜索，找到结构相似的同源模板

2. **RNA链MSA**：使用AlphaFold3数据管道自动运行：
   - 使用Nhmmer搜索NT-RNA、Rfam和RNACentral数据库
   - 合并搜索结果生成RNA的MSA

注意：所有MSA工具和数据库必须正确安装和配置，并且通过`--db_dir`参数指定数据库目录，**同时**设置`--use_msa`标志，才会执行MSA搜索。 