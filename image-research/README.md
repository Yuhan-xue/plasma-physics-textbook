# 配图从哪里来（image-research/）

这本书的插图有两条来源：**手绘示意图**（TikZ，随正文编译）与**真实图像**（装置照片、天文观测、
机构示意图，放在 `figures/`）。这份文档只讲后者是怎么选出来的。

## 选图标准

| 项目 | 要求 |
|---|---|
| 许可 | 只要 CC0 / 公有领域 / CC BY / CC BY-SA；NC、ND、「保留所有权利」一律不用 |
| 分辨率 | 长边 ≥ 1200 px（矢量图不限） |
| 内容 | 能说明一个具体的物理事实，而不是装饰性配图 |
| 来源 | Wikimedia Commons（API 检索，许可与作者取 `extmetadata`，不靠人工推断） |

已采用的图像共 113 张，逐张的作者、许可与来源页见 [`figures/CREDITS.md`](../figures/CREDITS.md)。
正文里每张图都在图注下方署名（`\figsource{作者｜许可｜来源页}`），并由
`check/figure_layout.py` 强制检查：文件存在、有署名、已登记、单张 ≤ 1.5 MB。

## 排除的许可

| 许可 | 处理 |
|---|---|
| CC BY-NC-*、CC BY-ND-* | 排除（禁止商用 / 禁止演绎） |
| FAL（Free Art License）、GFDL 1.2（无双重许可） | 排除 |
| Copyrighted free use | 排除（条款不明确） |
| CC BY-SA 3.0 igo（ESA/IGO 变体） | 保留，按 CC BY-SA 署名 |

因此未采用：`Oxygen discharge tube.jpg`、`Nitrogen discharge tube.jpg`、
`Hydrogen deuterium glow comparison.png`（FAL）；`Argon discharge tube.jpg`、`Neon discharge tube.jpg`、
`Lightning strike jan 2007.jpg`（GFDL 1.2）；`TorchCuttingCloseup.jpg`（Copyrighted free use）。
MIT OpenCourseWare、DoITPoMS（CC BY-NC-SA）、PhET（商用付费）、HyperPhysics 与
IAEA《Fusion Physics》（保留所有权利）同样不用。

## 建议自绘、不必外找的主题

- **磁重联示意图**：Commons 上该主题几乎只有视频，清晰标注版不如 TikZ 重绘；
- **磁镜 / 磁瓶位形**、**∇B、曲率与 E×B 漂移**：没有 ≥1200 px 的合适位图，继续用 TikZ。

## 采集时踩过的坑

- `upload.wikimedia.org` 对不符合 Wikimedia User-Agent 政策的客户端一律返回 **HTTP 429**，
  与请求频率无关；必须用 `客户端名/版本 (联系方式)` 形式的 UA。
- 用 `urllib` 取任意尺寸缩略图会被拒（HTTP 400），需选 Commons 列出的标准档位；
  稳定做法是用 `Special:FilePath/<文件名>?width=1500`。
- 检索与下载脚本在 `tools/`（`commons_pick.py`、`fetch_figure.py`），
  候选池表格可用 `tools/gen_table.py` 从 `tools/verified.json` 重新生成。
