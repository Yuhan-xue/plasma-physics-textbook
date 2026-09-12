# 附录 A–G 专项评审

评审日期：2026-09-12；基线及分级见[总体评审](00_overall.md)。附录是自学者求助与校对的入口，应该至少与正文同样可靠。

## P1｜附录 A 的拉普拉斯算子重复使用错误坐标

位置：[直角坐标](../src/plasma_physics_textbook_v2.tex#L4311)、[柱坐标](../src/plasma_physics_textbook_v2.tex#L4348)、[球坐标](../src/plasma_physics_textbook_v2.tex#L4390)。

直角坐标三项都写成对 \(x\) 的二阶导数；应分别为 \(x,y,z\)。柱坐标第二、三项应为 \(r^{-2}\partial_\theta^2\phi+\partial_z^2\phi\)，球坐标最后一项应为 \((r^2\sin^2\theta)^{-1}\partial_\varphi^2\phi\)。例如 \(\phi=y^2\) 的拉普拉斯为 2，错误式给 0。建议用同一标量场在三种坐标下交叉验证。

## P1｜附录 B 与正文的单位、章节引用和碰撞公式不一致

位置：[参数表](../src/plasma_physics_textbook_v2.tex#L4476)、[温度单位说明](../src/plasma_physics_textbook_v2.tex#L4534)、[碰撞频率](../src/plasma_physics_textbook_v2.tex#L4519)。

若其余量使用 SI 数值，\(T_e\) 的 eV 数值必须先换成 \(\Theta_e=1.602176634\times10^{-19}T_{e,\mathrm{eV}}\) J；“\(k_BT\to T\)”不能作为无条件数值替换。附录的碰撞频率系数与正文相差 \(\sqrt2\)，耦合参数、德拜球计数也未统一。表中“章节”仍指向旧章序，建议改用实际标签而不是手写章节号。

## P1｜附录 C 的欧姆定律、阿尔芬波和 Grad–Shafranov 推导仍有基本错误

位置：[C.1](../src/plasma_physics_textbook_v2.tex#L4567)、[C.3](../src/plasma_physics_textbook_v2.tex#L4618)、[C.4](../src/plasma_physics_textbook_v2.tex#L4669)。

由 \(j=en(v_i-v_e)\) 和附录自己的电子动量方程，忽略电子惯性应得 \(E+v_i\times B=\eta j+j\times B/(en)-\nabla p_e/(en)\)，而附录两个非电阻项都反号。C.3 声明 \(\exp[i(kz-\omega t)]\)，却在时间导数处使用 \(+i\omega\)；第 4629 行还把依赖 \(z\) 的量提出旋度，不能得到下一行的非零结果。C.4 算子最后一项写成 \(\partial_x^2\phi\)，应为 \(\partial_Z^2\psi\)。

## P1｜附录 D 的 SI–CGS 电荷密度和电流密度换算方向、数量级错误

位置：[单位对照表](../src/plasma_physics_textbook_v2.tex#L4731)。由 \(1\,\mathrm{statC}=3.335641\times10^{-10}\,\mathrm C\)，逐步处理面积和体积后，\(\rho[\mathrm{statC/cm^3}]\approx2.9979\times10^3\rho[\mathrm{C/m^3}]\)，\(j[\mathrm{statA/cm^2}]\approx2.9979\times10^5j[\mathrm{A/m^2}]\)。

教材写成分别除以 \(3\times10^9\) 和除以 \(3\times10^5\)，前者方向反了，后者也与常见 statA 定义不符。[NIST 单位换算表](https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b9)可作为基础电荷换算核对入口。建议给出 1 C/m³ 和 1 A/m² 的完整例子。

## P1｜附录 F 的“易错点”自己引入符号和条件错误

位置：[漂移误区](../src/plasma_physics_textbook_v2.tex#L4886)、[傅里叶约定](../src/plasma_physics_textbook_v2.tex#L4930)。

附录把 \(\exp[-i(k\cdot r-\omega t)]\) 与 \(\exp[+i(\omega t-k\cdot r)]\) 当作不同约定，其实两式完全相同；它们又与正文使用的 \(\exp[i(k\cdot r-\omega t)]\) 相反。另称真空场中梯度漂移与曲率漂移“总是大小相等”，但单粒子贡献分别含 \(v_\perp^2/2\) 和 \(v_\parallel^2\)，只有特定热平均才可相等。应从统一指数约定和完整速度项重写整组易错点。

## P1｜附录 G.2、G.3 代码不能按原样运行

位置：[G.2](../src/plasma_physics_textbook_v2.tex#L5079)、[G.3](../src/plasma_physics_textbook_v2.tex#L5153)。

在 Python 3.14.7、NumPy 2.5.3、SciPy 1.18.1 下，直接执行 LaTeX 中的 verbatim 代码：G.1 运行成功，与解析轨道最大位置差约 \(3.52\times10^{-10}\) m；G.2 和 G.3 在双反斜杠续行处触发 SyntaxError。临时把 G.2 续行改成单反斜杠后，求解器又报边界残差应为 2 项却返回 3 项。SciPy 官方接口要求边界残差数等于未知函数数加待求参数数（[接口说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_bvp.html)）。

建议先修复代码再提供期望图和收敛检查。G.2 只能选两个相容边界条件；G.3 应补回溯测试，报告 \(k\lambda_D\) 的近似范围。

## P1｜附录 G.1 的磁镜扩展不满足 \(\nabla\cdot B=0\)

位置：[G.1 扩展说明](../src/plasma_physics_textbook_v2.tex#L5065)。

只把 \(B_z\) 换成 \(B_0(1+(z/L)^2)\) 会使 \(\nabla\cdot B=\partial_zB_z\ne0\)，且当前洛伦兹力函数接收常量数组，不能仅改参数就产生磁镜反射。建议给满足散度为零的近轴场（例如配套 \(B_r\simeq-rB_z'(z)/2\)），把场改为位置函数，并检查磁矩、总能量和平行动能。

## P2｜附录 E 与 G 的学习入口不可直接复现

位置：[书单](../src/plasma_physics_textbook_v2.tex#L4776)、[G.4](../src/plasma_physics_textbook_v2.tex#L5194)、[电子资源说明](../src/plasma_physics_textbook_v2.tex#L869)。

仓库没有说明中承诺的独立脚本、PIC 教程和论文导读入口；书单缺少版次、ISBN/DOI 或出版方链接。建议逐条补齐可检索元数据，为每个外部资源指定“解决本书哪个疑问、读哪一节、读完能算什么”，并让 G.1–G.3 形成可下载、可运行、带期望输出的最小学习包。
