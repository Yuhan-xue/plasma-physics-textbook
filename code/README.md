# 书中脚本（可直接运行）

教材里出现的完整程序放在这里，正文只留「骨架 + 路径」——这样正文更短，
而想动手的读者能拿到真正能跑的版本。

## 目录

| 文件 | 对应章节 | 做什么 |
| --- | --- | --- |
| `appendix_g/orbit_tracking.py` | 附录 G.1 单粒子轨道追踪 | 均匀磁场回旋 + $E\times B$ 漂移 + 近轴磁镜反射；含能量守恒检查与三联图 |
| `appendix_g/debye_bvp.py` | 附录 G.2 德拜屏蔽电势的数值求解 | 无量纲球坐标泊松方程边值问题；与线性德拜解逐点对比、幅值扫描、`max_nodes` 负例 |
| `appendix_g/langmuir_dispersion.py` | 附录 G.3 朗缪尔波色散曲线 | Bohm--Gross 色散关系、冷等离子体对照、$k\lambda_D$ 适用范围回溯检查 |

## 怎么跑

```
pip install numpy scipy matplotlib
python code/appendix_g/debye_bvp.py
```

脚本用相对路径保存图片（`single_particle_orbit.png`、`debye_shielding.png`、
`langmuir_dispersion.png`），在仓库根目录运行时图片会落在根目录；想固定输出位置，
把 `savefig` 里的文件名换成带目录的路径即可。

## 与正文的关系

- 正文给出物理目标、方程与「应当验证什么」，脚本给出完整实现；
- 正文里印出的运行输出（如 `success/status/message` 与残差）就是这些脚本在
  NumPy 2.5.3 / SciPy 1.18.1 下的实测输出；
- 路径守卫在 `check/full_check.py` 的「代码清单」一类里：正文提到的每个 `code/**.py`
  必须存在且能通过语法编译，同时正文不再出现超过 20 行的代码清单。
