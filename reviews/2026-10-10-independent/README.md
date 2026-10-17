# 等离子体教材独立评审与 TODO

> 说明：本轮的机器索引（`*.json`）、`checks/*.txt` 与评审方复算脚本（`recheck_*.py`）
> 已移出仓库，只保留这份人类可读的记录；需要时用 `check/` 下的检查脚本重跑复现。

日期：2026-10-10。基线：`f6317daf8b4203b96674304e08471caccb9f895f`（v1.6.0）。

本轮确认 **13 项待处理问题：P1 × 1、P2 × 10、P3 × 2**。最优先的是第 9 章对称双流模型：正文、例题和速查卡使用了与已声明总密度不一致的解析根，导致稳定阈值、最大增长率和数值答案成组错误。第 12 章诊断定义、图形遮挡和数值验证流程也有可复现问题。

本目录只增加评审、TODO 与证据；未修正文、检查工具或发布 PDF。既有 `processed.md` 的结案记录不代替本轮独立核验。

| ID | 级别 | 问题 | 报告 |
| --- | --- | --- | --- |
| PP-C01 | P1 | 对称双流模型的频率归一化及参考系混用 | [内容](01-content.md) |
| PP-C02 | P2 | Buneman 共振增长率的质量比幂次错误 | [内容](01-content.md) |
| PP-C03 | P2 | 冷分布检验把 δ′ 积分算反，负平方仍写成正频率 | [内容](01-content.md) |
| PP-C04 | P2 | EEPF 被定义成只除密度的 EEDF | [内容](01-content.md) |
| PP-C05 | P2 | “漏 √ε 使温度低估一倍”不成立 | [内容](01-content.md) |
| PP-C06 | P2 | 集体汤姆逊散射与强耦合／直接离子散射混淆 | [内容](01-content.md) |
| PP-C07 | P2 | 斜入射反射条件把 sec θ 写成 cos θ | [内容](01-content.md) |
| PP-C08 | P3 | 两种漂移相等对应的俯仰角写成余角 | [内容](01-content.md) |
| PP-C09 | P3 | 十倍振幅与十个 e-folding 的时间比写错 | [内容](01-content.md) |
| PP-U01 | P2 | 图 9.1 注释遮挡连线，图 9.2 分类节点重叠 | [阅读与 UI](02-ui-reading.md) |
| PP-U02 | P2 | 附录 F 的章号导航沿用旧章节映射 | [阅读与 UI](02-ui-reading.md) |
| PP-W01 | P2 | 编译检查把任意单个缺字及单次严重溢出判 PASS | [验证流程](03-workflow.md) |
| PP-W02 | P2 | 附录 G.2 求解未收敛，却把 status=1 解释为收敛 | [验证流程](03-workflow.md) |

执行入口：[TODO](TODO.md)。机器索引：`findings.json`。每项均有源码行号、影响、反例及验收要求。

## 实测与范围

- 干净克隆后运行 `python check/full_check.py`：编译成功、218 页，FAIL 0 / WARN 1；WARN 是首次编译前尚无 `.toc`，散引用检查被跳过。编译后再运行 `--no-compile`：FAIL 0 / WARN 0。
- `ui_readability.py`：FAIL 0 / WARN 10；`consistency_audit.py`：FAIL 0 / WARN 0。日志归档于 `checks/*.txt`。这三套脚本覆盖结构和部分统计，不验证下列所有科学结论或实际图形遮挡。
- 已跟踪 PDF 为 **218 页、187 书签**，SHA-256 `d19ab39b11babce9d963b14fa07ede917a6fcec9f40317617de31a7e41daac9a`。重编译 PDF 哈希不同，但同环境下全部 218 页的提取文字、72 dpi 渲染逐页一致，书签也一致；详见 `制品比较`。这不等于全书科学内容通过。
- 接触表抽看 47 页，另放大检查关键页；本轮内容深查集中在第 6、7、9、12 章、附录 F/G 及相关速查卡。其余章节以结构、导航和版面抽样为主，**未声称 12 章逐式审完**。页清单见 `audit.json`。
- 独立脚本复现了色散根、分布函数、几何条件和两类门禁负例；附录 G 的 3 段原样代码也实际执行，G.2 的进程退出 0 与求解器 `success=False` 同时出现。证据见 `verification.json` 与 `example-runs.json`。
- 所有 PDF 页码为该 218 页文件的一基物理页。在线依据均于 2026-10-10 查阅。没有对外部线上交互资源、历史 release 全集或全部参考引文作认证。

## 复现入口

在仓库根目录运行：

```powershell
python check/full_check.py
python check/full_check.py --no-compile
python check/ui_readability.py
python check/consistency_audit.py
python reviews/2026-10-10-independent/recheck_review.py
python reviews/2026-10-10-independent/recheck_examples.py
```

评审脚本依赖 NumPy、SciPy、SymPy、PyMuPDF；示例执行另需 Matplotlib。现有完整检查需 LaTeX、latexmk、pdfinfo、pdftotext。环境版本记录见 `environment.json`。两个评审脚本只写 `tmp/review-20261010/`；它们的退出 0 表示取证／进程完成，**不是教材合格证**。修复后的数值应按 TODO 的正确关系验收，不能要求保留错误来让本轮历史反例继续通过。
