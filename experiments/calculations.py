"""Generate a mathematical walkthrough from the current data and actual fits."""
import json
import os
from pathlib import Path

import pandas as pd

from src.statistics import DENSITY, ETA3, P_DIVISIBLE, TARGETS, aut_size, cl_probability

LABELS = {
    ("divisible", None): "3 整除类数",
    ("torsion", None): "3-torsion 总和",
    ("sylow", "trivial"): "Sylow 平凡",
    ("sylow", "C3"): "Sylow C3",
    ("sylow", "C9"): "Sylow C9",
    ("sylow", "C3xC3"): "Sylow C3 × C3",
    ("sylow", "other"): "Sylow 其他类型",
}


def write_calculations(data_dir, output_dir, results, sensitivity, requested_xmin, requested_xmax):
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    manifest = json.loads((data_dir/"processed/p3_statistics.manifest.json").read_text())
    summaries = pd.read_csv(data_dir/"processed/p3_statistics.csv").set_index("X")
    first = results[0]["fit"]
    X = first["xmax"]
    row = summaries.loc[X]
    n, ndiv, torsion = (int(row[c]) for c in ("N_fields", "N_divisible", "S_torsion3"))
    ranks = pd.read_csv(data_dir/"processed/rank3_counts.csv")
    ranks = ranks[ranks.X == X]
    types = pd.read_csv(data_dir/"processed/sylow3_types.csv")
    ntypes = int(((types.X == X) & (types.N_group > 0)).sum())
    fence = chr(96)*3
    lines = ["# 本轮计算说明", "",
             "本说明由复现程序根据当前汇总表和实际拟合结果自动生成。所有数字来自原始 LMFDB 记录。",
             "重新运行复现命令会同时更新本说明、系数表和图像。", "",
             "## 1. 数据范围与输入", "",
             f"- 输入数据覆盖上限：{manifest['xmax']:,}；累计检查点：{len(summaries)} 个。",
             f"- 本报告统计端点：X = {X:,}，包含 {n:,} 个域、{ntypes} 种实际出现的 Sylow 3 类型。",
             f"- 请求拟合下界：{requested_xmin:g}；上界："
             f"{requested_xmax if requested_xmax is not None else '全部可用检查点'}。",
             f"- 实际拟合区间：{first['xmin']:,} ≤ X ≤ {first['xmax']:,}，共 {first['points']} 个点。",
             "", "每个数据块的四个同余类必须齐全。程序检查 gzip 数据、判别式次序和块范围，并验证类数等于循环因子阶数的乘积。",
             "", "| 输入文件 | 字节数 | SHA-256 |", "|---|---:|---|"]
    for item in manifest["inputs"]:
        lines.append(f"| {item['filename']} | {item['bytes']:,} | {item['sha256']} |")
    lines.extend(["", "## 2. 从原始记录到三个统计对象", "",
                  r"令 \(D<0\) 为基本判别式，\(A_D=\mathrm{Cl}(\mathbb Q(\sqrt D))\)，\(h_D=|A_D|\)。",
                  "文件名给出 r、m、k；每行是 a、h、c₁、…、cₜ。", "",
                  r"$$D_0=-k\,2^{28}-r,\qquad D_i=D_{i-1}-ma_i,\qquad "
                  r"A_{D_i}\cong\prod_j C_{c_j},\qquad h_{D_i}=\prod_j c_j.$$",
                  "", "提取每个 cⱼ 中的最大 3 次幂，删除等于 1 的因子，按降序排列：", "",
                  r"$$q_j=3^{v_3(c_j)},\qquad G_D=A_D[3^\infty]\cong\prod_{q_j>1}C_{q_j}.$$",
                  "", "令 r₃(D) 是保留下来的因子个数，则", "",
                  r"$$A_D[3]\cong C_3^{\,r_3(D)},\qquad |A_D[3]|=3^{r_3(D)},\qquad"
                  r"\mathbf1_{\{3\mid h_D\}}=\mathbf1_{\{r_3(D)>0\}}.$$",
                  "", "| 示例类群 | Sylow 3 类型 | 3-rank | 3-torsion 大小 | 3 整除类数 |",
                  "|---|---|---:|---:|---:|",
                  "| C18 | C9 | 1 | 3 | 1 |",
                  "| C3 × C3 | C3 × C3 | 2 | 9 | 1 |",
                  "| C18 × C12 | C9 × C3 | 2 | 9 | 1 |",
                  "| C2 | 平凡 | 0 | 1 | 0 |",
                  "", "C9 和 C3 × C3 的阶都为 9，但对 torsion 总和的贡献分别为 3 与 9。单位元也计入 torsion。",
                  "", "## 3. 累计计数与实际代入", "",
                  r"$$N(X)=\sum_{\substack{-X\le D<0\\D\ {\rm fundamental}}}1,\quad "
                  r"N_{\rm div}(X)=\sum_D\mathbf1_{\{3\mid h_D\}},\quad "
                  r"S_3(X)=\sum_D3^{r_3(D)},\quad N_G(X)=\sum_D\mathbf1_{\{G_D\cong G\}}.$$",
                  "",
                  r"$$P_{\rm div}(X)=N_{\rm div}(X)/N(X),\qquad "
                  r"M_3(X)=S_3(X)/N(X),\qquad P_G(X)=N_G(X)/N(X).$$",
                  "", "所有比例的分母均为全部域的数量。其他类型是所选四种类型的补集，不对四类比例重新归一化。", "",
                  f"在 X = {X:,}：", "",
                  f"- N(X) = {n:,}。",
                  f"- N_div(X) = {n:,} − {int(row['N_trivial']):,} = {ndiv:,}。",
                  f"- P_div(X) = {ndiv:,} / {n:,} = {ndiv/n:.12f} = {100*ndiv/n:.8f}%。",
                  f"- S₃(X) = {torsion:,}；M₃(X) = {torsion:,} / {n:,} = {torsion/n:.12f}。",
                  "", "| Sylow 类型 | 实测数量 | 实测比例 |", "|---|---:|---:|"])
    for group in [*TARGETS, "other"]:
        count = int(row[f"N_{group}"])
        lines.append(f"| {group} | {count:,} | {count/n:.12f} |")
    lines.extend(["", "完整秩分布提供加权校验：", "",
                  "| r | 秩为 r 的域数量 M_r(X) | 对 S₃ 的贡献 3ʳ M_r(X) |",
                  "|---:|---:|---:|"])
    contributions = []
    for rank in ranks.itertuples():
        contribution = 3**int(rank.rank3)*int(rank.N_rank)
        contributions.append(contribution)
        lines.append(f"| {int(rank.rank3)} | {int(rank.N_rank):,} | {contribution:,} |")
    if sum(contributions) != torsion:
        raise ValueError("Rank counts do not reconstruct the torsion sum")
    lines.extend(["", f"加权总和：{' + '.join(f'{v:,}' for v in contributions)} = {torsion:,}。",
                  "", "## 4. 主项与 CL 常数", "",
                  r"$$\rho=\frac3{\pi^2},\qquad N(X)\sim\rho X,\qquad"
                  r"\eta_3=\prod_{j=1}^{\infty}(1-3^{-j}),\qquad p_3=1-\eta_3.$$",
                  "", f"- ρ = {DENSITY:.15f}。",
                  f"- η₃ = {ETA3:.15f}。",
                  f"- p₃ = {P_DIVISIBLE:.15f}。",
                  "- 代码使用 j = 1,…,99 的乘积，以双精度计算；不把 p₃ 取整为 0.44。",
                  "", r"$$P_{\rm CL}(G)=\eta_3/|\mathrm{Aut}(G)|.$$",
                  "",
                  r"若 \(G=\prod_i C_{3^{\lambda_i}}\)，令 \(\lambda'_j=\#\{i:\lambda_i\ge j\}\)，"
                  r"\(m_a=\#\{i:\lambda_i=a\}\)，则",
                  "",
                  r"$$|\mathrm{Aut}(G)|=3^{\sum_j(\lambda'_j)^2}"
                  r"\prod_{a\ge1}\prod_{j=1}^{m_a}(1-3^{-j}).$$",
                  "", "| 类型 | 自同构群大小 | CL 主项比例 |", "|---|---:|---:|"])
    for group, parts in TARGETS.items():
        lines.append(f"| {group} | {aut_size(parts)} | {cl_probability(parts):.12f} |")
    other_mass = 1-sum(cl_probability(parts) for parts in TARGETS.values())
    lines.extend([f"| other | 不对应单个群 | {other_mass:.12f} |", "",
                  "CL 整除概率和固定 Sylow 类型分布是猜想输入，不是从样本拟合出的极限。"
                  "3-torsion 平均值趋于 2，以及总和中的 X^(5/6) 次主项，是已知定理。",
                  "定理依据：[Bhargava–Shankar–Tsimerman，定理 2 与 4](https://arxiv.org/abs/1005.0672)。",
                  "",
                  r"$$N_{\rm div}(X)\text{ 的主项为 }p_3\rho X,\quad "
                  r"S_3(X)\text{ 的主项为 }2\rho X=\frac6{\pi^2}X,\quad "
                  r"N_G(X)\text{ 的主项为 }P_{\rm CL}(G)\rho X.$$",
                  "",
                  "减掉单位元后，非平凡 3-torsion 元素的平均数量趋于 1。这个加权平均不能直接换成整除概率。",
                  "", "## 5. 两种误差的定义和实际数值", "",
                  "为每个对象写 T(X)，每域主项常数写 L：整除计数用 p₃，torsion 总和用 2，Sylow 计数用 P_CL(G)。", "",
                  r"$$E_{\rm linear}(X)=T(X)-L\rho X,\qquad E_{\rm empirical}(X)=T(X)-LN(X).$$",
                  "", "CLI 中 xuzhu 对应 linear，conditional 对应 empirical；conditional 不表示以 3 整除类数为条件的概率。",
                  "",
                  r"$$E_{\rm linear}-E_{\rm empirical}=L(N(X)-\rho X).$$",
                  "", f"本端点 N(X) − ρX = {n-DENSITY*X:.9f}。", "",
                  "| 对象 | L | T(X) | E_linear(X) | E_empirical(X) |",
                  "|---|---:|---:|---:|---:|"])
    primary = [r for r in results if r["error_convention"] == "xuzhu" and r["fit"]["model"] == "power"]
    for entry in primary:
        observable, group = entry["observable"], entry["group"]
        count = ndiv if observable == "divisible" else torsion if observable == "torsion" else int(row[f"N_{group}"])
        L = entry["ratio_leading_term"]
        lines.append(f"| {LABELS[(observable, group)]} | {L:.12f} | {count:,} | "
                     f"{count-L*DENSITY*X:.6f} | {count-L*n:.6f} |")
    lines.extend(["", "## 6. 最小二乘、数值方法与残差", "",
                  r"实际选中的检查点为 \((X_i,E_i)\)，最小化",
                  "",
                  r"$$\mathrm{RSS}=\sum_i(E_i-f(X_i))^2,\qquad"
                  r"\mathrm{RMSE}=\sqrt{\mathrm{RSS}/n_{\rm points}}.$$",
                  "",
                  "这是未加权的计数误差最小二乘，不是在 log(E) 上拟合，也不是直接对比例或平均值做最小二乘。"
                  "较大的 X 通常有较大的计数误差，因此更能影响此损失。", "",
                  "| 模型 | f(X) | 待拟合参数 |", "|---|---|---|",
                  "| power | A X^θ | A、θ |",
                  "| five_sixths | A X^(5/6) | A |",
                  "| five_sixths_plus_two_thirds | A X^(5/6) + B X^(2/3) | A、B |",
                  "", "自由幂次模型先消去振幅：令 $X_*$ 为拟合区间最大值，$t_i=X_i/X_*$，则",
                  "",
                  r"$$b(\theta)=\frac{\sum_i E_it_i^\theta}{\sum_i t_i^{2\theta}},\qquad "
                  r"A(\theta)=b(\theta)X_*^{-\theta}.$$",
                  "",
                  "代码在 θ∈[−1,2] 上扫描 601 个点，再在最佳网格点两侧做有界最小化。"
                  "若最优值位于搜索边界则报错。固定幂次模型对缩放后的基函数做线性最小二乘。",
                  "", r"残差为 \(r_i=E_i-\widehat f(X_i)\)；图中的归一化误差为 \(E_i/X_i^{5/6}\)。",
                  "", "## 7. 本轮拟合系数", "",
                  f"所有拟合均使用 {first['points']} 个点，区间为 [{first['xmin']:,}, {first['xmax']:,}]。"])
    for convention, title in [("xuzhu", "减去线性主项"), ("conditional", "减去实际域总数主项")]:
        entries = [r for r in results if r["error_convention"] == convention]
        index = {(r["observable"], r["group"], r["fit"]["model"]): r["fit"] for r in entries}
        lines.extend(["", f"### {title}", "",
                      "| 对象 | 自由 A | 自由 θ | 自由 RMSE | 固定 5/6 的 A | 固定 RMSE |",
                      "|---|---:|---:|---:|---:|---:|"])
        for entry in primary:
            key = entry["observable"], entry["group"]
            free, fixed = index[(*key, "power")], index[(*key, "five_sixths")]
            lines.append(f"| {LABELS[key]} | {free['parameters']['A']:.10g} | "
                         f"{free['parameters']['theta']:.9f} | {free['RMSE']:.6f} | "
                         f"{fixed['parameters']['A']:.10g} | {fixed['RMSE']:.6f} |")
        lines.extend(["", "| 对象 | 双幂次 A | 双幂次 B | 双幂次 RMSE |", "|---|---:|---:|---:|"])
        for entry in primary:
            key = entry["observable"], entry["group"]
            fit = index[(*key, "five_sixths_plus_two_thirds")]
            lines.append(f"| {LABELS[key]} | {fit['parameters']['A']:.10g} | "
                         f"{fit['parameters']['B']:.10g} | {fit['RMSE']:.6f} |")
    lines.extend(["", "## 8. 从 E 还原计数、比例与平均值", "",
                  r"线性主项下，拟合计数为 \(\widehat T(X)=L\rho X+\widehat E(X)\)。"
                  r"图中每域值严格使用 \(\widehat T(X)/N(X)\)。精确恒等式为",
                  "",
                  r"$$\frac{T(X)}{N(X)}-L=\frac{E_{\rm linear}(X)}{N(X)}"
                  r"+L\frac{\rho X-N(X)}{N(X)}.$$",
                  "",
                  r"用 \(N(X)\sim\rho X\) 作渐近近似，自由幂次模型对应",
                  "",
                  r"$$\frac{T(X)}{N(X)}\approx L+\frac A\rho X^{\theta-1}.$$",
                  "",
                  "固定 θ=5/6 时，每域修正指数是 −1/6，系数是 A/ρ，与计数修正系数 A 不同。"
                  "如果使用 empirical 误差，则 T/N − L = E_empirical/N 是精确等式。"
                  "torsion 平均值由总和拟合再除以 N(X) 得到，没有另外优化平均值误差。",
                  "", "四个核心对象的自由幂次模型，使用本轮 linear 拟合系数："])
    symbols = {("divisible", None): r"N_{\rm div}", ("torsion", None): r"S_3",
               ("sylow", "C9"): r"N_{C_9}", ("sylow", "C3xC3"): r"N_{C_3\times C_3}"}
    for entry in primary:
        key = entry["observable"], entry["group"]
        if key not in symbols:
            continue
        A, theta = entry["fit"]["parameters"]["A"], entry["fit"]["parameters"]["theta"]
        L, symbol = entry["ratio_leading_term"], symbols[key]
        lines.extend(["", f"**{LABELS[key]}**", "",
                      r"$$"+symbol+rf"(X)\approx {L*DENSITY:.12g}X {A:+.10g}X^{{{theta:.9f}}}."+r"$$",
                      "",
                      rf"$$\frac{{{symbol}(X)}}{{N(X)}}\approx {L:.12g}"
                      rf" {A/DENSITY:+.10g}X^{{{theta-1:.9f}}}.$$"])
    lines.extend(["", "这些是有限范围内的数值模型，不能据此写出已证明的 little-o 余项。", "",
                  "## 9. 拟合区间检查", "",
                  "下表是 linear 误差下的自由指数；全部截点结果见 fit_sensitivity.csv。", "",
                  "| 对象 | 请求下界 | 实际首个 X | 点数 | θ | RMSE |",
                  "|---|---:|---:|---:|---:|---:|"])
    for item in sensitivity:
        key = item["observable"], item["group"] or None
        if item["error"] != "xuzhu" or key not in symbols or item["requested_xmin"] not in (10**6, 10**7, 10**8):
            continue
        lines.append(f"| {LABELS[key]} | {item['requested_xmin']:,} | {item['actual_xmin']:,} | "
                     f"{item['points']} | {item['theta']:.9f} | {item['RMSE']:.6f} |")
    lines.extend(["", "累计检查点相关，因此不把普通最小二乘标准误当作独立样本置信区间。"
                  "AIC/BIC 仅作描述性诊断。比较模型时还应检查残差形状和指数随截点的漂移，"
                  "不能仅凭一次拟合接近 5/6 就断言存在该次主项。",
                  "", "## 10. 复现与文件对应", "", fence+"bash",
                  "source .venv/bin/activate",
                  "python -m scripts.reproduce --check",
                  "python -m scripts.reproduce --rebuild --check",
                  "",
                  "python -m experiments.fit --observable divisible --error xuzhu --model power",
                  "python -m experiments.fit --observable torsion --error xuzhu --model five_sixths",
                  "python -m experiments.fit --group C9 --error conditional --model power",
                  "python -m experiments.fit --group C3xC3 --error conditional --model power",
                  fence, "",
                  "第一条复现命令从汇总表生成图像、拟合结果和本说明；第二条先从原始 gzip 记录重新累计。"
                  "独立拟合还输出逐点 fitted_error、residual、fitted_ratio 的 CSV。", "",
                  "| 文件 | 内容 |", "|---|---|"])
    for filename, purpose in [
        ("processed/p3_statistics.csv", "三个统计的总表及两种误差"),
        ("processed/sylow3_types.csv", "每个检查点的完整 Sylow 同构类型分布"),
        ("processed/rank3_counts.csv", "按 3-rank 统计，重建 torsion 总和"),
        ("processed/p3_statistics.manifest.json", "输入哈希和检查点列表")]:
        link = os.path.relpath(data_dir/filename, output_dir)
        lines.append(f"| [{Path(filename).name}]({link}) | {purpose} |")
    lines.extend(["| [fit_summary.csv](fit_summary.csv) | 42 组模型比较及 RMSE |",
                  "| [fits.json](fits.json) | 未舍入的参数、损失诊断和输入哈希 |",
                  "| [fit_sensitivity.csv](fit_sensitivity.csv) | 拟合截点结果 |",
                  "| [final_statistics.csv](final_statistics.csv) | 本报告端点的实测值 |", ""])
    output = output_dir/"calculations.md"
    text = "\n".join(lines).replace(r"\(", "$").replace(r"\)", "$")
    output.write_text(text, encoding="utf-8")
    return output
