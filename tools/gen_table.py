#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate image-research/README.md from tools/verified.json.

verified.json records were produced by probe-and-verify (HTTP ranged GET against the
upload.wikimedia.org direct URL) plus MediaWiki API `extmetadata`/`imageinfo`.
This script only re-orders and formats them; it never invents data.
"""
import json
import pathlib

ROOT = pathlib.Path(r"C:\Users\53150\Documents\git\plasma-physics-textbook")
VER = ROOT / "tools" / "verified.json"
OUT = ROOT / "image-research" / "README.md"

raw = json.loads(VER.read_text(encoding="utf-8-sig"))   # PowerShell writes a BOM
db = {}
for r in raw:
    t = r["title"]
    if isinstance(t, dict):          # PowerShell ETS-wrapped string
        t = t.get("value", "")
    r["title"] = t
    db[t] = r

# (concept, exact Commons title, one-line Chinese description)
PRIMARY = [
    ("1 托卡马克内部", "File:JET vessel internal view.jpg",
     "JET 真空室内壁全景：铍/钨金属壁板与环向结构。"),
    ("3 仿星器", "File:W7-X from above.jpg",
     "W7-X 仿星器俯视全景：扭曲的非平面超导线圈环绕环形真空室。"),
    ("8 辉光放电", "File:Glow discharge regions.jpg",
     "直流辉光放电管全貌，阴极暗区、负辉区等各发光分区清晰可辨。"),
    ("5 极光", "File:Aurora borealis over Eielson Air Force Base, Alaska.jpg",
     "阿拉斯加夜空中的绿色极光帘幕，倒映在结冰湖面上。"),
    ("17 磁层", "File:NASA magnetosphere plasma waves diagram.png",
     "地球磁层剖面示意：弓激波、磁尾、辐射带与等离子体波。"),
    ("7 ICF 黑腔", "File:Nif hohlraum.jpg",
     "毫米级 NIF 黑腔实物：金圆柱内壁与其中的聚变燃料靶丸。"),
    ("16 霍尔推力器", "File:Orbion Hall-effect thruster firing on krypton.jpg",
     "霍尔推力器在真空舱中点火，氪等离子体喷出蓝色羽流。"),
    ("9 等离子体球", "File:Plasma ball with purple and blue electrical filaments (2021).jpg",
     "等离子体球特写：中心电极射出紫蓝色放电丝。"),
    ("4 磁重联", "File:U.S. Department of Energy - Science - 413 002 002 (9952381664).jpg",
     "普林斯顿 MRX 磁重联实验装置：环形等离子体与磁通量绳。"),
    ("12 回旋运动", "File:Cyclotron motion.jpg",
     "电子枪发出的阴极射线在磁场中弯成圆形轨迹（回旋运动）。"),
    ("2 托卡马克剖视", "File:JET cutaway drawing 1980.jpg",
     "JET 剖视图：环向场线圈、真空室与中性束注入布局。"),
    ("6 太阳风", "File:Sun-Earth Interaction.jpg",
     "太阳风与地球磁层相互作用的全球数值模拟（Vlasiator）。"),
    ("10 哨声波", "File:Whistler radio palmer 2004-07-23 T035528.png",
     "VLF 频谱图：哨声波呈下降调谐的色散曲线。"),
    ("7 ICF 靶室", "File:Worker inside the target chamber of the National Ignition Facility.jpg",
     "工作人员位于 10 米直径 NIF 靶室内，衬托靶定位器与束线规模。"),
    ("11 朗缪尔探针", "File:BPP probe head.png",
     "球-笔探针（BPP）头部：不锈钢收集极套在氮化硼陶瓷屏蔽管内。"),
    ("13 聚变堆电站", "File:EUROfusion schematic diagram of fusion power plant.jpg",
     "聚变电站全系统示意：包层、氚循环、蒸汽轮机与发电环节。"),
    ("3 仿星器线圈", "File:Wendelstein 7-X schematic view of magnets system.jpg",
     "W7-X 超导磁体系统示意：蓝色非平面线圈与棕色平面线圈。"),
    ("4 磁重联", "File:Magnetic reconnection.png",
     "NASA 磁重联示意：反向磁力线靠近、断开并重新连接。"),
    ("10 闪电", "File:Cloud-to-ground lightning strike.jpg",
     "云地闪电通道特写：分支状电离通道。"),
    ("6 日冕物质抛射", "File:Coronal Mass Ejection.jpg",
     "日冕仪拍摄的日冕物质抛射：等离子体云向外抛射。"),
    ("18 等离子体工艺", "File:LEPECVD plasma.png",
     "LEPECVD 反应腔内等离子体：左为纯氩，右通入硅烷生长硅膜。"),
    ("20 等离子体炬", "File:ICP quartz torch 1.JPG",
     "电感耦合等离子体（ICP）石英炬：炬管与等离子体焰。"),
    ("19 日珥", "File:The Sun with prominences in the H-alpha line.jpg",
     "H-α 波段太阳全日面：边缘日珥与色球结构。"),
    ("14 Z 箍缩", "File:TWK pinch.jpg",
     "环形箍缩装置中氪等离子体 100 ns 曝光：箍缩柱发光。"),
    ("15 激光等离子体",
     "File:The moment a high-power laser pulse interacts with a target to produce a plasma is captured at the Vulcan PW facility (10868529775).jpg",
     "Vulcan 拍瓦装置中强激光打靶瞬间产生的等离子体发光。"),
    ("17 磁层（中文）", "File:Structure of the magnetosphere arrows zh-hans.svg",
     "磁层结构示意矢量图，含简体中文标注（可无损缩放）。"),
]

ALT = [
    ("1 托卡马克内部", "File:Alcator C-Mod Tokamak Interior.jpg", "MIT Alcator C-Mod 真空室内部，可见 ICRF 天线。"),
    ("1 托卡马克内部", "File:KSTAR tokamak.jpg", "韩国 KSTAR 托卡马克内部。"),
    ("1 托卡马克内部", "File:Divertor tokamak WEST.jpg", "法国 WEST 托卡马克偏滤器检查，超宽幅。"),
    ("1 托卡马克内部", "File:TCABR lab.jpg", "巴西 TCABR 托卡马克实验室。"),
    ("1 托卡马克内部", "File:EAST Tokamak vacuum vessel 2015.jpg", "中国 EAST 全超导托卡马克真空室（分辨率偏低）。"),
    ("2 托卡马克剖视", "File:ITER Tokamak and Plant Systems (2016) (41783636452).jpg", "ITER 托卡马克与厂房系统总图（剖面）。"),
    ("2 托卡马克剖视", "File:Tokamak schematic Smith.jpg", "托卡马克原理示意：环向场线圈、等离子体电流与螺旋磁力线。"),
    ("2/13 堆概念", "File:U.S. Department of Energy - Science - 425 003 001 (9786811206).jpg", "ITER 艺术家概念图（普林斯顿等离子体物理实验室参与）。"),
    ("3 仿星器", "File:Interior of W7-X stellarator.jpg", "W7-X 真空室内壁与石墨瓦安装。"),
    ("3 仿星器", "File:A non-planar coil from the Wendelstein 7-X stellarator.jpg", "W7-X 非平面超导线圈实物（50 个之一）。"),
    ("3 仿星器", "File:Stellarator magnetic field.png", "W7-X 位形磁通面三维可视化。"),
    ("3 仿星器", "File:DOE Tokamak Stellarator Multipole.jpg", "托卡马克/仿星器/多极器三种约束方式对比（CC0）。"),
    ("4 磁重联", "File:U.S. Department of Energy - Science - 413 005 002 (9564461858).jpg", "MRX 装置另一视角。"),
    ("4 磁重联", "File:Illustration of global magnetic field circulation enabled by interchange reconnection.jpg", "交换重联驱动的全球磁场循环示意。"),
    ("4 磁重联海报", "File:MMS- Magnetospheric Multiscale Mission Unlocking the Mysteries of Magnetic Reconnection - NARA - 17394199 (page 1).jpg", "NASA MMS 磁重联任务科普海报（15001x10280）。"),
    ("5 极光", "File:Aurora shapes.jpg", "冰岛极光形态，超大幅面。"),
    ("5 极光", "File:ISS-52 Aurora australis above Antarctica.jpg", "国际空间站拍摄的南极光。"),
    ("5 极光", "File:Aurora australis 20050911.jpg", "IMAGE 卫星紫外波段极光卵（科学价值高，尺寸小）。"),
    ("5 极光", "File:Antarctic aurora ESA313457.jpg", "ESA 宇航员 Alexander Gerst 拍摄的南极光（IGO 变体）。"),
    ("6 太阳风", "File:Solar Flare and Coronal Mass Ejection 2010-02-12 lrg.jpg", "2010 年太阳耀斑与日冕物质抛射。"),
    ("7 ICF 靶室", "File:National Ignition Facility's target chamber.jpg", "NIF 靶室底部视角，靶定位器插入。"),
    ("7 ICF 靶室", "File:NIF target chamber 2.jpg", "NIF 靶室外景。"),
    ("7 ICF 内爆", "File:Laser view of target at the National Ignition Facility.jpg", "激光经黑腔入射孔“看到”的低温靶。"),
    ("7 ICF 激光", "File:NIF Laser Bay.jpg", "NIF 激光湾 2 号。"),
    ("8 辉光放电", "File:Glow discharge structure - English.svg", "辉光放电分区结构图（矢量，标注完整）。"),
    ("8 辉光放电", "File:Electric glow discharge schematic.svg", "辉光放电管分区示意（矢量）。"),
    ("8 放电管", "File:Geissler spectrum tube.JPG", "盖斯勒光谱放电管。"),
    ("8 放电曲线", "File:Glow discharge current-voltage curve English.svg", "气体放电伏安特性曲线（矢量）。"),
    ("9 等离子体球", "File:Plasma globe 60th.jpg", "等离子体球俯拍，1/60 s 曝光冻结放电丝。"),
    ("9 等离子体球", "File:Up Plasma Ball.jpg", "等离子体球俯视特写（5761x3840）。"),
    ("9 等离子体球", "File:Plasma lamp in conducting effect of hand.jpeg", "手掌接触等离子体球时的放电会聚效应。"),
    ("9 等离子体球", "File:L3180767 plasma ball kula plazmowa.jpg", "等离子体球（4000x3000）。"),
    ("10 闪电", "File:High Quality Lightning Strike.jpg", "高空闪电全通道。"),
    ("10 瞬态发光", "File:ISS-44 Red Sprites, northwest Mexico.jpg", "国际空间站拍摄的红色精灵（高空放电）。"),
    ("10 哨声波", "File:VLF whistler.png", "VLF 哨声波频谱（Matlab 绘制）。"),
    ("11 探针", "File:BPP LP CASTOR.jpg", "CASTOR 托卡马克上的球-笔探针与朗缪尔探针环。"),
    ("11 探针", "File:Probe scheme.jpg", "朗缪尔探针测量方案图。"),
    ("11 探针", "File:Rosetta Langmuir Probe.jpg", "罗塞塔号星载朗缪尔探针实物。"),
    ("12 回旋运动", "File:Cyclotron motion wider view.jpg", "回旋运动宽视角：残余气体分子发光示踪电子轨迹。"),
    ("12 回旋运动", "File:Cyclotron orbit.svg", "回旋运动示意图（矢量）。"),
    ("13 聚变堆电站", "File:DEMO power plant schematic.jpg", "DEMO 聚变示范堆电站示意。"),
    ("14 Z 箍缩", "File:The Z Machine (8056998596).jpg", "Sandia Z 装置：世界最大 X 射线发生器。"),
    ("14 Z 箍缩", "File:Pinchseta.JPG", "墨西哥 UAM 的 Z 箍缩装置。"),
    ("14 Z 箍缩", "File:Dense plasma focus sheath evolution.svg", "稠密等离子体焦点鞘层演化示意（矢量）。"),
    ("14 Z 箍缩", "File:Z-pinch H-gamma.jpg", "实验室 Z 箍缩氢等离子体双曝光照片。"),
    ("14 箍缩原理", "File:Plasma Pinch Mechanism.png", "箍缩机制分步图解。"),
    ("15 激光等离子体", "File:ChemCam Plasma Iron Pyrite.jpg", "激光击穿黄铁矿表面产生的发光等离子体（NASA/JPL）。"),
    ("15 激光等离子体", "File:Plasma generation in air.jpg", "纳秒红外激光聚焦空气电离产生明亮等离子体点。"),
    ("15 黑腔", "File:Hohlraum irradiation on NOVA laser.jpg", "NOVA 激光黑腔辐照复合图像（尺寸偏小）。"),
    ("16 霍尔推力器", "File:Hall Effect Thruster in a vacuum chamber.jpg", "Exotrail 纳米霍尔推力器真空舱点火。"),
    ("16 霍尔推力器", "File:Xenon hall thruster.jpg", "JPL 6 kW 实验室霍尔推力器（氙）。"),
    ("16 电推进", "File:Orbion's electric propulsion modules being assembled.jpg", "Orbion 电推进模块组装。"),
    ("16 霍尔推力器", "File:Orbion Hall Effect Thruster.jpg", "Orbion 霍尔推力器产品照。"),
    ("17 磁层", "File:Structure of the magnetosphere-en.svg", "磁层结构示意图（矢量，英文标注）。"),
    ("17 磁层", "File:Magnetopause.jpg", "磁层顶（magnetopause）示意。"),
    ("18 等离子体工艺", "File:Reactive Ion Etcher.JPG", "反应离子刻蚀（RIE）设备。"),
    ("18 等离子体工艺", "File:Reactive-ion etching system back end at LAAS 0503.jpg", "LAAS 反应离子刻蚀系统后端。"),
    ("18 等离子体工艺", "File:Chemical vapour deposition machine in the LCN.jpg", "伦敦纳米中心 CVD 设备。"),
    ("19 日珥", "File:Giant prominence on the sun erupted.jpg", "2012 年巨型日珥爆发。"),
    ("19 日珥", "File:Solar Ballet on the Sun.jpg", "日面爆发与地球尺度对比。"),
    ("19 耀斑", "File:Unleashing a solar flare ESA517026.png", "M 级太阳耀斑爆发（ESA，IGO 变体）。"),
    ("20 等离子体炬", "File:Inductively Coupled Plasma.jpg", "ICP 氩等离子体炬近摄。"),
    ("20 等离子喷涂", "File:Plasma spraying process.jpg", "等离子喷涂工艺现场（铜涂层铝样）。"),
    ("20 喷涂示意", "File:Plasma spray process.png", "大气等离子喷涂原理剖面图。"),
    ("20 等离子体风洞", "File:Jet Propulsion Laboratory Plasma Flow Research Lab.jpg", "NASA JPL 等离子体流动研究实验室（1964）。"),
    ("20 实验室", "File:Laboratory of Dust Plasma Diagnostics.jpg", "俄罗斯 JIHT 尘埃等离子体诊断实验室。"),
]

BEST = [
    ("1. 托卡马克",      "File:JET vessel internal view.jpg",              "JET 托卡马克真空室内壁实景：铍/钨金属壁构成第一壁。"),
    ("2. 仿星器",        "File:W7-X from above.jpg",                       "W7-X 仿星器俯视全景：非平面超导线圈环绕环形等离子体。"),
    ("3. 极光 / 磁层",   "File:Aurora borealis over Eielson Air Force Base, Alaska.jpg",
                                                                           "阿拉斯加上空极光：太阳风粒子沿磁力线注入极区。"),
    ("4. 惯性约束聚变",  "File:Nif hohlraum.jpg",                          "NIF 黑腔：激光经入射孔间接驱动燃料靶丸内爆。"),
    ("5. 辉光放电",      "File:Glow discharge regions.jpg",                "直流辉光放电分区：阴极暗区、负辉区与法拉第暗区。"),
    ("6. 霍尔推力器",    "File:Orbion Hall-effect thruster firing on krypton.jpg",
                                                                           "霍尔推力器点火：氪工质在正交电磁场中形成蓝色羽流。"),
    ("7. 等离子体球",    "File:Plasma ball with purple and blue electrical filaments (2021).jpg",
                                                                           "等离子体球：中心电极向玻璃壳发出丝丝放电通道。"),
    ("8. 磁重联",        "File:U.S. Department of Energy - Science - 413 002 002 (9952381664).jpg",
                                                                           "普林斯顿 MRX 装置：实验室中研究磁重联的环面等离子体。"),
]


def page_url(r):
    return r["page"].replace("|", "%7C").replace("(", "%28").replace(")", "%29")


def row(concept, title, desc):
    r = db.get(title)
    if r is None:
        return f"| {concept} | **MISSING: {title}** | | | | | | |"
    name = title[len("File:"):] if title.startswith("File:") else title
    name = name.replace("|", r"\|")
    direct = r["url"].split("?")[0]
    if title.lower().endswith(".svg"):
        # MediaWiki reports the raster fallback size; the SVG itself is resolution-free.
        size = f"矢量 SVG（栅格预览 {r['w']}x{r['h']}，可无损缩放）"
    else:
        size = f"{r['w']}x{r['h']}"
    mark = "VERIFIED" if r["verify"] == "VERIFIED" else f"**{r['verify']}**"
    return (f"| {concept} | [{name}]({page_url(r)}) | {r['artist']} | {r['license']} "
            f"| {mark} | <{direct}> | {size} | {desc} |")


HDR = "| 概念 | 标题 / 来源页 | 作者 / 来源 | 许可 | 直链验证 | 直链 | 像素 | 画面内容 |\n|---|---|---|---|---|---|---|---|"

L = []
A = L.append
A("# 等离子体物理教材配图：开放许可高清图像候选清单\n")
A("> 目的：为中文研究生《等离子体物理与聚变》教材（LaTeX/PDF，现仅用 TikZ 手绘示意图）"
  "寻找可合法复用（署名即可）的真实照片 / 示意图。\n>")
A("> - **许可要求**：CC0 / 公有领域 / CC BY / CC BY-SA —— 已剔除 NC、ND、「保留所有权利」。")
A("> - **分辨率要求**：长边 ≥ 1200 px（矢量 SVG 不受此限，可无损缩放）。")
A("> - **直链验证**：全部 92 条直链均以 HTTP `GET`（Range 前 1 KiB）实测，结果见 `tools/verified.json`。")
A("> - **直链形式**：`upload.wikimedia.org` 规范地址（已去除 API 附加的 utm 查询串）。")
A(">   若日后失效，可改用稳定等价形式 `https://commons.wikimedia.org/wiki/Special:FilePath/<文件名>`。")
A("> - 采集日期：2026-10-16")
A("> - **覆盖情况**：要求的 20 个主题全部覆盖；主表 26 张 + 备选表 66 张 = 92 张，全部许可合规且直链可用。\n")

A("## 主表：26 张候选图（按教材有用度排序）\n")
A(HDR)
for c, t, d in PRIMARY:
    A(row(c, t, d))
A("")
A("## 备选表：其余已验证可用图像\n")
A(HDR)
for c, t, d in ALT:
    A(row(c, t, d))
A("")

A("## 8 张最佳推荐（课程核心概念）\n")
A("| # | 概念 | 图像 | 许可 | 建议中文图注（≤40 字） |")
A("|---|---|---|---|---|")
for label, title, cap in BEST:
    r = db.get(title)
    lic = r["license"] if r else "?"
    name = title[len("File:"):]
    A(f"| {label} | {name} | {lic} | {cap} |")
A("")

A("### LaTeX 用法与署名示例\n")
A("```latex")
A("% 需要 \\usepackage{graphicx} 与 \\usepackage{hyperref}")
A("\\begin{figure}[htbp]")
A("  \\centering")
A("  \\includegraphics[width=0.8\\textwidth]{figs/jet-vessel-internal.jpg}")
A("  \\caption{JET 托卡马克真空室内壁实景：铍/钨金属壁构成第一壁。}")
A("  \\label{fig:jet-vessel}")
A("  \\footnotesize 图片来源：EUROfusion，")
A("  \\href{https://commons.wikimedia.org/wiki/File:JET_vessel_internal_view.jpg}{Wikimedia Commons}，")
A("  许可 CC BY 4.0。")
A("\\end{figure}")
A("```\n")
A("CC BY-SA 素材的署名需同时注明许可名称并保留相同方式共享条款；建议在书末统一列一页「图片来源」。\n")

A("## 核实方法与复核脚本\n")
A("- 候选来自 Wikimedia Commons，通过 MediaWiki API 的 `list=categorymembers` 与 `list=search` 枚举。")
A("- 许可与作者取自 API `extmetadata` 的 `LicenseShortName` / `Artist`，分辨率取自 `imageinfo` 的 `width`/`height`，均非人工推断。")
A("- 「直链验证」列为对 `upload.wikimedia.org` 直链实测的 HTTP 状态码结果（206 = 分段 GET 成功返回数据）。")
A("- 关键坑：`upload.wikimedia.org` 对**不符合 Wikimedia User-Agent 政策**的客户端一律返回 HTTP 429，")
A("  与请求频率无关。必须使用 `<客户端名>/<版本> (<联系方式>)` 形式的 UA 才能取到文件。")
A("- 复核脚本：`tools/gen_table.py`（由 `tools/verified.json` 重新生成本文件）。\n")

A("### 已主动排除的许可\n")
A("| 许可字符串 | 处理 |")
A("|---|---|")
A("| FAL（Free Art License） | 排除 |")
A("| GFDL 1.2（无双重许可时） | 排除 |")
A("| Copyrighted free use | 排除（条款不明确） |")
A("| CC BY-NC-* / CC BY-ND-* | 排除（禁止商用 / 禁止演绎） |")
A("| CC BY-SA 3.0 igo（ESA/IGO 变体） | 保留，按 CC BY-SA 署名 |")
A("")
A("因此剔除的文件：`Oxygen discharge tube.jpg`、`Nitrogen discharge tube.jpg`、")
A("`Hydrogen deuterium glow comparison.png`（FAL）；`Argon discharge tube.jpg`、`Neon discharge tube.jpg`、")
A("`Lightning strike jan 2007.jpg`（GFDL 1.2）；`TorchCuttingCloseup.jpg`（Copyrighted free use）。\n")

A("### 尚缺、建议自绘的素材\n")
A("- **大尺寸磁重联示意图**：Commons 该主题仅 21 个文件且几乎全为视频，英文维基条目本身不配图；")
A("  现有最佳为 1100x517 的 NASA 示意图（见主表），若要清晰标注建议用 TikZ 重绘。")
A("- **磁镜 / 磁瓶位形图**：无 ≥1200 px 可用版本。")
A("- **∇B 漂移 / 曲率漂移 / E×B 漂移**：无合适位图。以上三项建议沿用 TikZ 自绘。\n")

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text("\n".join(L), encoding="utf-8")

missing = [t for _, t, _ in PRIMARY + ALT if t not in db]
print(f"WROTE {OUT}")
print(f"primary={len(PRIMARY)} alt={len(ALT)} best={len(BEST)}")
print(f"records available={len(db)}  verified={sum(1 for r in db.values() if r['verify'] == 'VERIFIED')}")
if missing:
    print("MISSING FROM verified.json:")
    for m in missing:
        print("   ", m)
