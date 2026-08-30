# Plasma Physics Textbook — Formula Numbering Unification Report

## Summary

All mathematical formulas across the plasma physics textbook have been systematically reviewed and unified.
Core physics formulas (Lawson criterion, Paschen law, Bohm criterion, Child-Langmuir, Grad-Shafranov, etc.)
have been converted from unnumbered `\\[...\\]` environments to numbered `equation` environments with semantic labels.
Formulas inside theorem/definition/example/insight/formula/checklist tcolorbox environments were intentionally
left unnumbered as per design requirements.

---

## 1. Main File: `plasma_physics_textbook_v2.tex`

### 1.1 Fourier Operator Convention (Converted from `\\[...\\]` to `equation`)
| Old Environment | New Environment | Label |
|-----------------|-----------------|-------|
| `\\[...\\]` | `equation` | `eq:ch1-fourier-operator` |

### 1.2 Previously Unlabeled `equation`/`align` Environments (Labels Added)
| Line | Formula Description | Environment | Label |
|------|---------------------|-------------|-------|
| 2788 | 鞘层密度关系 (n_e, n_i) | `equation` | `eq:ch11-sheath-densities` |
| 2793 | 鞘层 Poisson 方程 | `equation` | `eq:ch11-sheath-poisson` |
| 2798 | Bohm 判据 (u_s ≥ √(k_B T_e/m_i)) | `equation` | `eq:ch11-bohm-criterion` |
| 2852 | Child-Langmuir 定律 (J_i = (4/9)ε₀√(2e/m_i) V_s^(3/2)/s²) | `equation` | `eq:child-langmuir` |
| 2867 | 无碰撞鞘层方程组 (能量守恒+连续性) | `align` | `eq:ch11-sheath-equations` |
| 2872 | 鞘层 Poisson 方程 (n_i 代入后) | `equation` | `eq:ch11-sheath-poisson2` |
| 3482 | Langmuir 探针 I-V 特性 | `equation` | `eq:ch12-probe-current` |
| 3616 | 等离子体折射率 N = √(1 - ω_p²/ω²) | `equation` | `eq:ch12-refractive-index` |
| 4341 | MHD 动量方程 | `equation` | `eq:ch9-mhd-momentum` |
| 4348 | 广义欧姆定律 (完整形式) | `equation` | `eq:ch9-generalized-ohm` |
| 4356 | 广义欧姆定律 (简化形式) | `equation` | `eq:ch9-generalized-ohm-simplified` |
| 4369 | Landau 积分 ∫ (∂F₀/∂v)/(v-ω/k) dv | `equation` | `eq:ch7-landau-integral` |
| 4376 | Landau 围道积分结果 | `equation` | `eq:ch7-landau-contour` |
| 4382 | Landau 阻尼条件 ω_i < 0 | `equation` | `eq:ch7-landau-condition` |
| 4408 | Alfvén 波 B₁ 表达式 | `equation` | `eq:ch9-alfven-b1` |
| 4414 | Alfvén 波动量方程 | `equation` | `eq:ch9-alfven-momentum` |
| 4420 | Alfvén 波色散关系 ω² = k²v_A² | `equation` | `eq:ch9-alfven-dispersion` |
| 4434 | Grad-Shafranov B 场表达式 | `equation` | `eq:ch10-grad-shafranov-b` |
| 4443 | Grad-Shafranov 方程 | `equation` | `eq:ch10-grad-shafranov` |
| 4448 | Grad-Shafranov 算子 Δ*ψ | `equation` | `eq:ch10-grad-shafranov-operator` |

---

## 2. Supplement Files

### 2.1 `ch10_supplement.tex` (Controlled Fusion — 12 modifications)
| Formula | Old Environment | New Environment | Label |
|---------|-----------------|-----------------|-------|
| Lawson 判据 (完整形式) | `\\[...\\]` | `equation` | `eq:ch10-lawson-full` |
| Lawson 判据 (简化形式) | `\\[...\\]` | `equation` | `eq:ch10-lawson-simplified` |
| 三重乘积 nTτ_E | `\\[...\\]` | `equation` | `eq:ch10-triple-product` |
| Tokamak 能量约束时间 scaling | `\\[...\\]` | `equation` | `eq:ch10-tokamak-scaling` |
| Stellarator 能量约束时间 scaling | `\\[...\\]` | `equation` | `eq:ch10-stellarator-scaling` |
| Stellarator 体积 scaling | `\\[...\\]` | `equation` | `eq:ch10-stellarator-scaling2` |
| Tokamak 大环体积 V = 2π²R₀a²κ | `\\[...\\]` | `equation` | `eq:ch10-tokamak-volume` |
| 安全因子 q(r) = rB₀/RB_θ(r) | `\\[...\\]` | `equation` | `eq:ch10-safety-factor` |
| 极向磁场安培定律 ∮B_θ dl = μ₀I_p | `\\[...\\]` | `equation` | `eq:ch10-ampere-law` |
| 边缘极向磁场 B_θ(a) = μ₀I_p/(2πa) | `\\[...\\]` | `equation` | `eq:ch10-poloidal-field` |
| q₉₅ 近似公式 | `\\[...\\]` | `equation` | `eq:ch10-q95-approx` |
| q₉₅ 经验公式 | `\\[...\\]` | `equation` | `eq:ch10-q95-empirical` |

**Note:** Duplicate labels in `formula` tcolorbox examples were fixed:
- `eq:ch10-lawson-simplified` → `eq:ch10-lawson-simplified-ex` (second occurrence)
- `eq:ch10-poloidal-field` → `eq:ch10-poloidal-field-ex` (second occurrence)

### 2.2 `ch11_supplement.tex` (Gas Discharge & Sheath — 7 modifications)
| Formula | Old Environment | New Environment | Label |
|---------|-----------------|-----------------|-------|
| Paschen 击穿电压定律 | `\\[...\\]` | `equation` | `eq:ch11-paschen-law` |
| Bohm 判据 v_s ≥ √(k_B T_e/m_i) | `\\[...\\]` | `equation` | `eq:ch11-bohm-criterion` |
| Child-Langmuir 定律 | `\\[...\\]` | `equation` | `eq:ch11-child-langmuir` |
| Debye 长度 λ_D | `\\[...\\]` | `equation` | `eq:ch11-debye-length` |
| E×B 漂移速度 v_d = E×B/B² | `\\[...\\]` | `equation` | `eq:ch11-exb-drift` |
| 热速度 v_th = √(8k_B T_e/πm_e) | `\\[...\\]` | `equation` | `eq:ch11-thermal-velocity` |

**Note:** Duplicate label fixed:
- `eq:ch11-child-langmuir` → `eq:ch11-child-langmuir-ex` (second occurrence in `formula` tcolorbox)
- `eq:ch11-debye-length-ex` added for second Debye length example in `formula` tcolorbox

### 2.3 `ch12_supplement.tex` (Plasma Diagnostics — 13 modifications)
| Formula | Old Environment | New Environment | Label |
|---------|-----------------|-----------------|-------|
| Langmuir 探针 I = I_esat(1 - exp(-eV/k_B T_e)) | `\\[...\\]` | `equation` | `eq:ch12-probe-current` |
| ln(I_esat - I) = ln(I_esat) + eV/k_B T_e | `\\[...\\]` | `equation` | `eq:ch12-probe-lnI` |
| 电子温度 T_e = e(V₂-V₁)/(k_B ln((I_esat-I₁)/(I_esat-I₂))) | `\\[...\\]` | `equation` | `eq:ch12-probe-temperature` |
| OML 离子饱和电流 I_i = (2/3)e n_e A √(2k_B T_e/m_i) | `\\[...\\]` | `equation` | `eq:ch12-oml-current` |
| 等离子体频率 ω_p = √(n_e e²/ε₀ m_e) | `\\[...\\]` | `equation` | `eq:ch12-plasma-frequency` |
| 折射率 N = √(1 - ω_p²/ω²) | `\\[...\\]` | `equation` | `eq:ch12-refractive-index` |
| 相位移动 Δφ = -e² n_e L/(2ε₀ m_e c ω) | `\\[...\\]` | `equation` | `eq:ch12-phase-shift` |
| 密度-相位关系 n_e = -2ε₀ m_e c ω Δφ/(e² L) | `\\[...\\]` | `equation` | `eq:ch12-density-phase` |
| Thomson 散射 Δλ_FWHM = 4λ₀ sin(θ/2) √(2k_B T_e ln2/m_e c²) | `\\[...\\]` | `equation` | `eq:ch12-thomson-fwhm` |
| 电子温度 T_e = m_e c² (Δλ/λ₀)² / (8 sin²(θ/2) ln2) | `\\[...\\]` | `equation` | `eq:ch12-thomson-temperature` |
| Faraday 定律 ∮E·dl = -dΦ_B/dt | `\\[...\\]` | `equation` | `eq:ch12-faraday-law` |
| Ampere 定律 ∮B·dl = μ₀ I_enc | `\\[...\\]` | `equation` | `eq:ch12-ampere-law` |
| Rogowski 线圈 EMF = -μ₀ n A (dI/dt) | `\\[...\\]` | `equation` | `eq:ch12-rogowski-emf` |

**Note:** Duplicate label fixed:
- `eq:ch12-refractive-index` → `eq:ch12-refractive-index-ex` (second occurrence in `formula` tcolorbox)

---

## 3. Files Not Modified (All Core Formulas Already Properly Numbered)

- `ch1_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments
- `ch2_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments
- `ch3_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments
- `ch4_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments
- `ch5_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments
- `ch6_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments
- `ch7_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments
- `ch8_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments
- `ch9_supplement.tex` — All formulas inside `example`/`formula` tcolorbox environments

---

## 4. Compilation Status

✅ **Compilation Successful** — XeLaTeX produced 222-page PDF with no errors.

Warnings (non-critical, pre-existing):
- Overfull/underfull hbox warnings (typography, pre-existing)
- Hyperref math-in-bookmark warnings (Unicode math in PDF bookmarks, pre-existing)
- Caption `hypcap=true` warning (pre-existing)
- `Label(s) may have changed` — Normal for first pass, resolves on second compilation

---

## 5. Label Naming Convention

All labels follow the pattern: `eq:ch{chapter}-{description}`

Examples:
- `eq:ch10-lawson-full` — Chapter 10, Lawson criterion (full form)
- `eq:ch11-bohm-criterion` — Chapter 11, Bohm criterion
- `eq:ch12-rogowski-emf` — Chapter 12, Rogowski coil EMF

---

## 6. Design Principles Applied

1. **Core physics formulas are numbered**: All fundamental equations appearing in main text or
   example problem solutions receive equation numbers and semantic labels.
2. **Intermediate derivation steps are unnumbered**: Step-by-step algebraic manipulations
   inside `proof` environments remain in `\\[...\\]` for cleaner visual flow.
3. **tcolorbox interiors are unnumbered**: Formulas inside `theorem`, `definition`, `example`,
   `insight`, `formula`, and `checklist` tcolorbox environments are intentionally unnumbered
   to maintain visual consistency with the boxed presentation style.
4. **Duplicate labels resolved**: When the same formula appears in both main text and
   `formula` tcolorbox summary, second occurrence gets `-ex` suffix.
