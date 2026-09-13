# 评审收尾与可重复验证记录

日期：2026-09-13。教材基线：3cfa48e6052a562d558b8567aea914d2b0d07a53；首轮评审提交：eb0f72a4fab965dc4b60eeda8d48cf50269c0c1d。远端该评审提交之后暂无教材修改。本文补齐首轮评审留在临时工作文件中的验证方法，入口见[总评](00_overall.md)。

## 本轮收尾范围

完成问题分级和全部章节入口；纠正总评把 4.1–4.3 写成不存在的误判，改为可核对的傅里叶和玻尔兹曼小节内容错配；将第 2 章的“太阳风”改回题干实际的“星际介质”；校正第 2、5 章源文件锚点；补记 Bennett 例题中间指数错误。

更新内容属于评审文档。下面程序用于复核评审中的具体断言；运行成功不表示教材中的全部错误已修正，也不表示通过了全书逐页出版校对。

## 数值复算

以下程序仅依赖 Python 标准库。保存为临时文件后运行即可，常数采用足以核对本轮数量级的 SI 数值。程序中的合理范围用于容忍常数取整；它核验表中所用输入下的结果，不为未知实验条件作标定。

~~~python
from math import pi, sqrt

e, eps0, me = 1.602176634e-19, 8.8541878188e-12, 9.1093837139e-31
kb, amu, c = 1.380649e-23, 1.66053906892e-27, 299792458.0
mu0 = 4*pi*1e-7  # 本轮数值精度下的近似
def debye(n, temperature_eV):
    return sqrt(eps0*temperature_eV/(n*e))
def omega_pe(n):
    return sqrt(n*e*e/(eps0*me))

ld = debye(1e3, kb*1e6/e)
diffusion = 3*e/(70*amu*1e5)
re = e*e/(4*pi*eps0*me*c*c)
checks = [
    ("rare_gas_lambda_m", ld, 2180, 2185),
    ("rare_gas_ND_sphere", 4*pi/3*1e3*ld**3, 4.34e13, 4.37e13),
    ("interstellar_lambda_m", debye(1e6, .86), 6.89, 6.90),
    ("interstellar_omega_s-1", omega_pe(1e6), 5.63e4, 5.65e4),
    ("fusion_omega_s-1", omega_pe(1e20), 5.63e11, 5.65e11),
    ("thermal_magnetic_drift_m_s", 2*1e4/(5*2), 1999, 2001),
    ("one_loss_cone_fraction", (1-sqrt(1-1/5))/2, .0527, .0529),
    ("two_loss_cones_fraction", 1-sqrt(1-1/5), .1054, .1057),
    ("collision_main_coefficient",
     e**4/(12*pi**1.5*eps0**2*sqrt(me)*e**1.5), 2.05e-12, 2.06e-12),
    ("collision_appendix_coefficient",
     e**4/(3*(2*pi)**1.5*eps0**2*sqrt(me)*e**1.5), 2.90e-12, 2.91e-12),
    ("two_stream_ku2_gamma_ratio", sqrt(sqrt(5)-2), .485, .487),
    ("Lawson_no_radiation_m-3_s", 12*1e4/(3.5e6*1.1e-22), 3.11e20, 3.12e20),
    ("ICF_internal_energy_J", 3*(.3e-6/4.15e-27)*1e4*e, 3.47e5, 3.48e5),
    ("diffusion_m2_s", diffusion, 41.3, 41.4),
    ("radial_loss_time_s", .06**2/diffusion, 8.70e-5, 8.72e-5),
    ("RF_Debye_length_m", debye(1e16, 3), 1.28e-4, 1.30e-4),
    ("RF_sheath_estimate_m",
     sqrt(2)/3*debye(1e16, 3)*(400/3)**.75, .00237, .00239),
    ("probe_Debye_length_m", debye(1.1e17, .29), 1.20e-5, 1.22e-5),
    ("phase_coefficient_m", re, 2.81e-15, 2.83e-15),
    ("cutoff_50GHz_m-3", eps0*me*(2*pi*50e9)**2/e**2, 3.10e19, 3.11e19),
    ("phase_density_2mm_10pi_m-3", 10*pi/(re*.002), 5.57e18, 5.58e18),
    ("coil_inductance_H", mu0*50**2*2e-4/(2*pi), .999e-7, 1.001e-7),
    ("Bennett_thermal_energy_J", mu0*1e5**2/(8*pi*1e19), 4.99e-17, 5.01e-17),
]
for name, value, low, high in checks:
    assert low < value < high, (name, value)
    print(f"{name}: {value:.9g}")
print(f"{len(checks)} numerical checks passed")
~~~

复核结果：23 项通过。数值对应第 2、3、5、8–12 章的问题，推导与适用假设仍须同时阅读逐章评审。例如 RF 鞘层只是沿用教材的简化估算，不能把其复算结果当作真实 RF 鞘层的精确解。

## 附录 G 原始代码与探针拟合复现

在仓库根目录运行以下程序，需要 NumPy、SciPy、Matplotlib。程序读取源文件中的三个原始 verbatim 代码块，在系统临时目录内执行，图像和临时目录随退出清理，不改教材。复核环境：Python 3.14.7、NumPy 2.5.3、SciPy 1.18.1；以无窗口绘图后端执行。

~~~python
from pathlib import Path
from tempfile import TemporaryDirectory
import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

source = Path("src/plasma_physics_textbook_v2.tex").read_text(encoding="utf-8")
appendix = source[source.index(r"\chapter*{附录G"):]
blocks = re.findall(r"\\begin\{verbatim\}\n(.*?)\\end\{verbatim\}", appendix, re.S)
assert len(blocks) == 3
original_cwd = Path.cwd()
with TemporaryDirectory(prefix="plasma-review-") as temp_dir:
    try:
        os.chdir(temp_dir)
        states = []
        for i, code in enumerate(blocks, 1):
            state = {}
            try:
                exec(compile(code, f"G{i}_original", "exec"), state)
                states.append("OK")
            except (SyntaxError, ValueError) as error:
                states.append(type(error).__name__)
                print(f"G{i}: {type(error).__name__}: {error}")
            finally:
                plt.close("all")
            if i == 1:
                t = state["t"]
                w = state["q"]/state["m"]
                exact = np.array([
                    1000*t + 99000/w*np.sin(w*t),
                    99000/w*(np.cos(w*t)-1),
                    1e4*t,
                ]).T
                error = np.max(np.abs(state["sol"][:, :3]-exact))
                assert error < 1e-8
                print(f"G1: max position error = {error:.6g} m")
        assert states == ["OK", "SyntaxError", "SyntaxError"], states
        # 只去掉重复续行符，进一步定位 G.2 的边界条件问题。
        adapted = blocks[1].replace("\\\\" + "\n", "\\" + "\n")
        try:
            exec(compile(adapted, "G2_continuation_only", "exec"), {})
        except ValueError as error:
            assert "shape (2,)" in str(error) and "(3,)" in str(error)
            print("G2 continuation-only:", error)
        else:
            raise AssertionError("Expected the baseline boundary-shape error")
    finally:
        plt.close("all")
        os.chdir(original_cwd)

v = np.array([-5., -2., 0., 2.])
electron_current = np.array([.10, .25, .55, 1.05])  # mA；扣除 -0.50 mA 离子电流
a, b = np.polyfit(v, np.log(electron_current), 1)
te = 1/a
vpl = (np.log(2.10)-b)/a  # 总平台 1.60 mA + 离子份额 0.50 mA
assert 2.94 < te < 2.95 and 4.05 < vpl < 4.06
print(f"probe fit: Te={te:.6g} eV, Vpl={vpl:.6g} V")
~~~

复核结果：

- G.1 最大位置差约 3.52×10⁻¹⁰ m；
- G.2、G.3 原始块触发续行语法错误；
- G.2 仅改续行后的边界残差维度为 3，预期维度为 2；
- 按第 12 章评审明确的四点、恒定离子电流和电子平台假设，拟合得 2.94419 eV、4.05372 V。

这些断言特意描述当前基线的已知失败。当教材代码被修正时，应更新为求解成功、边界残差和物理极限检查，而不是继续期待这些错误。

## 文档与范围检查

核对所有逐章评审、附录和本文件均可从总评到达；本地相对文件链接和源文件行号有效。关键数值复算和 G.1–G.3 的复现方法已随评审提交。源教材、PDF 和外部设备配置不在本次评审收尾的修改范围内。

首轮 PDF 仅抽查第 16、25、32、100、148、189 页，本次没有新增逐页视觉检查。原子数据、经验定标与书目仍按逐章评审所列条件要求作者核对，不把这些后续教材修订工作记为已完成。
