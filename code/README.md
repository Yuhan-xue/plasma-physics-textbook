# code/ —— 书中会用到的脚本

正文里不贴整段程序：需要动手的地方只写"用哪个脚本、跑出来应该看到什么"，
脚本本体放在这里，可以直接运行。

## 三个脚本

| 脚本 | 对应位置 | 做什么 | 你会看到 |
| --- | --- | --- | --- |
| [`appendix_g/orbit_tracking.py`](appendix_g/orbit_tracking.py) | 附录 G「单粒子轨道追踪」 | 均匀磁场里的回旋 + $E\times B$ 漂移 + 近轴磁镜反射 | 三联图：轨道、能量守恒曲线、磁镜反射点 |
| [`appendix_g/debye_bvp.py`](appendix_g/debye_bvp.py) | 附录 G「德拜屏蔽电势的数值求解」 | 球坐标泊松方程的边值问题，与线性德拜解逐点对比 | 电势曲线对比图 + 求解器自报的收敛信息（`success`、`status`、残差范数） |
| [`appendix_g/langmuir_dispersion.py`](appendix_g/langmuir_dispersion.py) | 附录 G「朗缪尔波色散曲线」 | Bohm–Gross 色散关系与冷等离子体对照 | 色散曲线 + 适用范围（$k\lambda_D$）回溯检查 |

## 怎么跑

装好依赖，在仓库根目录执行即可：

```
pip install numpy scipy matplotlib
python code/appendix_g/debye_bvp.py
```

三个脚本都会把图存到当前目录（`single_particle_orbit.png`、`debye_shielding.png`、
`langmuir_dispersion.png`）。想固定输出位置，改一下脚本里保存图片的文件名即可。

## 与正文的关系

- 正文负责讲清物理目标、方程，以及"应当验证什么"；脚本负责给出可直接运行的实现。
- 附录 G 里印出的运行信息（是否收敛、迭代情况、残差大小）就是这些脚本的真实输出，
  版本为 NumPy 2.5 / SciPy 1.18；换环境后数值可能在小数位上有差别。
- 想在正文里新增脚本引用，只需写清路径；全量检查会核对"正文提到的每个脚本都存在且能通过语法编译"。
