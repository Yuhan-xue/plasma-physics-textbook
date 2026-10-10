# 构建与数值验证流程

## PP-W01／P2：用数量白名单放过真实缺字和严重溢出

**定位**：`check/full_check.py:490–506`，其中 `504` 把 missing≤1 判 PASS，`506` 把 overfull≤5 判 PASS；`check/ui_readability.py:474` 附近也只按 overfull 个数分级。

当前构建确实只有一个 U+000A 缺字记录，但检查并未验证字符身份。注释“1 = 仅 U+000A 换行符属正常”被实现成“任意一个缺字都正常”；溢出则不论宽度有多大，只要次数少就通过。

本轮不改生产日志，在临时目录给 `check_C_compile(False)` 注入以下独立夹具：

| 输入 | C 类检查结果 | 应表达的结果 |
| --- | --- | --- |
| 无错误记录 | PASS | 正常对照 |
| 一个已知 U+000A 记录 | PASS | 可按精确白名单豁免 |
| 一个真实 Ω（U+03A9）缺字 | PASS，无 WARN/FAIL | 真实数学字符缺失，应阻断或显式升级 |
| 一次 200 pt overfull | PASS，仅 INFO 列出宽度 | 严重溢出，不能靠次数少豁免 |
| 一条 `! Undefined control sequence.` | FAIL | 错误检测正对照有效 |

详细结果见 [verification.json](verification.json) 的 `compile_gate_controls`。这是日志分类函数的隔离实验，不声称本次正式 PDF 实际缺 Ω 或实际溢出 200 pt。

**建议与验收**：按码点、字体及明确原因豁免特定缺字；按溢出幅度、位置和内容类型设阈值。未知缺字、单次严重溢出均必须被识别；已知换行记录仍允许。对编译日志不存在、检查项未执行等状态保留独立的“未验证”，不要合并成通过。

## PP-W02／P2：附录 G.2 未收敛，状态码却标成收敛

**定位**：主文件 `5570–5580`；PDF p.208。原样示例与结果见 [example-runs.json](example-runs.json)。

本例调用 `solve_bvp(..., tol=1e-8, max_nodes=50000)` 后打印：`solver status (1=最近一次网格上收敛)`。SciPy 的实际约定是 status=0 才收敛、status=1 表示超出允许网格节点数，`success` 仅在达到精度时为真。[SciPy 官方返回值说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_bvp.html)

当前示例原样运行的结果：

```text
process exit_code = 0
status = 1
success = false
message = The maximum number of mesh nodes is exceeded.
mesh_nodes = 45151
max_rms_residual = 6.110961375200931e-06
requested tol = 1e-08
max relative deviation from Debye = 4.6635732387168005e-05
```

45151 是最后一次允许的网格数；下一轮需要更多节点会超过上限，并不要求最终已保存网格恰好达到 50000。残差约为设定容差的 611 倍。与线性解的某个误差指标较小，不能覆盖求解器明确失败的事实；这也不意味着当前曲线每个点都完全错误。

**建议**：修正状态解释，以 `sol.success`、`sol.message` 和残差判定是否可接受结果；未收敛时停止合格结论。随后根据无量纲化、尺度、网格和容差检查原因，不通过无限增加节点掩盖问题。教材应说明可复现的接受标准，并留一个故意限制节点数的失败样例。

**验收**：正常例必须 `success=True` 且达到声明容差，再对边界条件及线性极限作独立验证；受限节点负例必须被标为失败，不能进而打印“收敛”或返回代表验收成功的状态。

## 本轮已验证的正面结果与边界

- 完整 LaTeX 构建退出 0；现有检查结果见 [checks](checks/README.md)。首轮干净构建的散引用检查先于 `.toc` 生成，复跑后才实际覆盖，报告分别保留两次结果。
- 重编译 PDF 与仓库 PDF 在本机同渲染器、72 dpi 下逐页相同，文字／书签一致；PDF 字节哈希不同，不据此报告附件不一致。
- G.1 输出 E×B 漂移 1000 m/s、磁镜转向约 −0.016522／0.016505 m、相对动能漂移 1.026×10⁻⁷；G.3 输出三组热修正 0.03／0.27／3.0。三个示例进程都完成，但只有 G.2 在本轮被发现有求解状态错误。此处没有宣称 G.1 的磁矩或 G.3 的全区间物理有效性已经完整验收。
- 示例在 Matplotlib Agg 模式执行，`plt.show()` 的无交互窗口警告属于运行方式差异；原始代码未改写。
