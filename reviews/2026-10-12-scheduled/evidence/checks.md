# 检查记录

基线 `69605fdd6e77e1cf702ccefd1a9992bc2dadee77`，Windows Python 3.14；检查器须用 UTF-8 标准输出，否则 GBK 无法编码数学字符并在打印汇总时异常退出。设置 `PYTHONIOENCODING=utf-8` 后复跑。

| 检查 | 结果 |
|---|---|
| `python check/physics_regression.py` | 退出 0；PASS 1351，INFO 1（两处 α 判据表述需要按上下文阅读，不是失败项）。 |
| `python check/figure_layout.py` | 退出 0；FAIL 0 / WARN 0 / 15 项。三张插图节点布局、遮挡与 6.5 pt 最小字号检查通过；统一色板静态检查通过。 |
| `python check/pdf_outline.py` | 退出 0；PASS 6 / WARN 2 / 11 项。页数 391、书签 392、目录匹配和链接数量自洽；警告分别为 O3b 17/392 的 −1 页偏移及 O4 附录 F 31/35 目标页与 aux 页差异。 |
| 上一轮 `reviews/2026-10-11-scheduled/CLOSURE.md` | 六项修复记录均有当前版本为 1.8.1 的验收过程；本轮回归未重现相同数值错误。 |

本轮检查覆盖数值回归、插图几何及 PDF 书签/目录数据。没有声称已人工逐页阅读 391 页；O3b/O4 仍需按待办完成页面级目标核验。
