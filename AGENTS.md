# AGENTS.md

## 长任务进度上报

预计运行超过约 30 秒的命令，一律用进度 CLI 包起来执行：

    node "$env:DSH_PROGRESS_CLI" run --task <任务id> -- <你的命令>

- CLI 会自动宣布任务开始、从命令输出里读取进度（行内出现 `12%` 或 `12/88` 都能识别）、
  退出时按退出码写终态（0 = done，非 0 = failed），不需要手写上报脚本。
- `--task` 的任务 id 必须作为完整单词出现在命令行里（面板靠它把这行进度和 job 对上）；
  id 只用字母、数字、`.`、`_`、`-`，1–40 个字符。
- 包装器表达不了的情况（脚本自己算进度），把进度事件追加到
  `$DSH_PROGRESS_DIR/<任务id>.jsonl`（一行一条 JSON，例如
  `{"v":1,"task":"<id>","state":"running","pct":42,"msg":"linking"}`），
  或用 `node "$env:DSH_PROGRESS_CLI" emit --task <id> --pct N --msg "..."` 上报、
  `done` / `fail` 收尾。
- 一个任务 id 只留一个写入者；结束必须写终态，不要把进度文件读回来——它是给人看的。
