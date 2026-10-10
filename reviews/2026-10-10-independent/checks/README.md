# 现有检查实测日志

全部命令从仓库根目录运行。以下为原始 stdout/stderr 的 UTF-8 归档，仅把本机工作区路径替换为 `$ROOT`。

| 命令 | 退出码 | 结果 | 日志 |
| --- | --- | --- | --- |
| `python check/full_check.py` | 0 | FAIL 0 / WARN 1，含完整构建 | [首次](full-first.txt) |
| `python check/full_check.py --no-compile` | 0 | FAIL 0 / WARN 0，已有 toc 后补验 | [复跑](full-recheck.txt) |
| `python check/ui_readability.py` | 0 | FAIL 0 / WARN 10 | [UI](ui.txt) |
| `python check/consistency_audit.py` | 0 | FAIL 0 / WARN 0 | [一致性](consistency.txt) |

报告条目数取自实际日志：首次 31、复跑 30、UI 29、一致性 14；不将 README 宣称的检查数当作实测数。日志门禁负例与数值反例另见 [verification.json](../verification.json)。
