# 数据与计算入口

## 原始输入

数据来自 [LMFDB 虚二次域类群表](https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups)。四个 Git LFS 文件位于 `data/raw/`：

| 文件 | 同余条件 |
|---|---|
| cl3mod8.0.gz | $\lvert D\rvert\equiv3\pmod8$ |
| cl7mod8.0.gz | $\lvert D\rvert\equiv7\pmod8$ |
| cl4mod16.0.gz | $\lvert D\rvert\equiv4\pmod{16}$ |
| cl8mod16.0.gz | $\lvert D\rvert\equiv8\pmod{16}$ |

文件名一般为 `cl{r}mod{m}.{k}.gz`；块 k 覆盖 $k2^{28}\le|D|<(k+1)2^{28}$。当前四个文件均为 k=0。

每行记录为 `a h c1 ... ct`。从 $D_0=-k2^{28}-r$ 开始，更新 $D_i=D_{i-1}-ma_i$；该域的类群为 $\prod_j C_{c_j}$，类数必须等于 $\prod_j c_j$。

例如 `cl3mod8.0.gz` 开头四行：

```text
0 1 1
1 1 1
1 1 1
2 2 2
```

依次对应 D = −3、−11、−19、−35，类数为 1、1、1、2。

`src/dataset.py` 的每条 `FieldRecord` 包含 `d`、`D_K`、`h_K`、`Cl_K_invariants`、`sylow3_invariants`，并提供 `rank3`、`torsion3_size`。其中 $K=\mathbb Q(\sqrt d)$，d 为无平方因子的负数；当 $D\equiv1\pmod4$ 时 d=D，否则 d=D/4。

## 本轮汇总文件

| 文件 | 内容 |
|---|---|
| p3_statistics.csv | 三个核心统计、所选四种 Sylow 类型与 other，以及两种误差 |
| sylow3_types.csv | 每个检查点的全部实际出现的 Sylow 类型及比例 |
| rank3_counts.csv | 每个检查点的 3-rank 分布 |
| p3_statistics.manifest.json | 输入 SHA-256、文件大小、数据范围、完整检查点列表 |

以上文件均位于 `data/processed/`。`N_divisible = N_fields - N_trivial`；`S_torsion3 = sum_r 3^r N_rank`；`Avg_torsion3 = S_torsion3 / N_fields`。类型比例的分母始终为全部域，`other` 为四种所选类型的补集。

类型键按因子降序，例如 `C9xC3`；平凡群为 `trivial`。完整计算过程、当前数字、所有拟合系数见 [自动生成的计算说明](results/calculations.md)。

## 重建与检查

```bash
python -m scripts.list_fields --max-D 100
python -m scripts.build_stats
python -m scripts.reproduce --rebuild --check
```

默认使用 C++ gzip 扫描器；`--engine python` 使用参考实现。测试比较两者在类型分类、检查点边界和块边界的结果。扫描器的内存占用由检查点和群类型数量决定，不随记录总数线性增长。

每个所需数据块的四个同余类必须齐全；程序拒绝缺失块、LFS 指针、错误类数或越界记录。LMFDB 原表有 k=0,…,4095。扩展命令见 README；额外原始文件被 Git 忽略，不自动上传。
