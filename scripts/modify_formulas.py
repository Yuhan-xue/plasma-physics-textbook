import re
import os

base_dir = r"C:\Users\53150\Documents\kimi\workspace\plasma_physics_textbook"

# ============================================================
# MODIFICATIONS FOR plasma_physics_textbook_v2.tex
# ============================================================

main_file_modifications = []

# Unnumbered displaymath outside tcolorbox that should be numbered
main_file_modifications.append({
    'type': 'replace',
    'old': r"""\
\[
\frac{\partial}{\partial t} \to -\ii\omega, \qquad \nabla \to +\ii\vect{k}.
\]""",
    'new': r"""\begin{equation}
\frac{\partial}{\partial t} \to -\ii\omega, \qquad \nabla \to +\ii\vect{k}.
\label{eq:ch1-fourier-operator}
\end{equation}"""
})

# Numbered equations without labels - add labels
# Line 383 (align): Fourier transform definition
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'align',
    'content_snippet': ' f(\\vect{r},t) &= \\int \\frac{\\diff^3 k\\,\\diff\\omega}{(2\\pi)^4}\\; \\tilde{f}(\\vect{k},\\omega)\\,\\ee^{+\\ii(\\vect{k}\\cdot\\vect{r}-\\omega t)}, \\\\',
    'label': 'eq:ch1-fourier-transform'
})

# Line 2784 (align): Bohm sheath densities
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'align',
    'content_snippet': ' n_i(x) &= n_0 \\left(1 - \\frac{2e\\phi(x)}{m_i u_0^2}\\right)^{-1/2} \\\\\\\\ n_e(x) &= n_0 \\exp\\!\\left(\\frac{e\\phi(x)}{k_B T_e}\\right)',
    'label': 'eq:ch11-sheath-densities'
})

# Line 2789 (equation): Sheath Poisson
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\eps_0 \\frac{\\diff^2\\phi}{\\diff x^2} = e(n_e - n_i)',
    'label': 'eq:ch11-sheath-poisson'
})

# Line 2793 (equation): Bohm criterion
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\frac{1}{k_B T_e} \\ge \\frac{m_i u_0^2}{k_B T_e} \\quad\\Rightarrow\\quad u_0 \\ge \\sqrt{\\frac{k_B T_e}{',
    'label': 'eq:ch11-bohm-criterion'
})

# Line 2860 (align): Sheath equations
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'align',
    'content_snippet': ' \\frac{1}{2}m_i u^2(x) + e\\phi(x) &= \\frac{1}{2}m_i u_B^2 \\quad (\\text{能量守恒}) \\\\\\\\ n_i(x) u(x) &= n_s',
    'label': 'eq:ch11-sheath-equations'
})

# Line 2865 (equation): Sheath Poisson 2
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\frac{\\diff^2\\phi}{\\diff x^2} = -\\frac{e n_i}{\\eps_0}',
    'label': 'eq:ch11-sheath-poisson2'
})

# Line 3474 (equation): Langmuir probe
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\ln I = \\ln I_{\\mathrm{e,sat}} + \\frac{V - V_p}{T_e}.',
    'label': 'eq:ch12-probe-current'
})

# Line 3606 (equation): Refractive index
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' n_r \\approx 1 - \\frac{1}{2}\\frac{\\omega_{pe}^2}{\\omega^2} = 1 - \\frac{e^2}{2\\eps_0 m_e \\omega^2} n_',
    'label': 'eq:ch12-refractive-index'
})

# Line 4331 (equation): MHD momentum
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\rho_m \\left(\\frac{\\partial \\vect{v}}{\\partial t} + (\\vect{v}\\cdot\\grad)\\vect{v}\\right) = \\vect{j}\\',
    'label': 'eq:ch9-mhd-momentum'
})

# Line 4337 (equation): Generalized Ohm's law
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\vect{E} + \\vect{v}_i\\times\\vect{B} = \\frac{1}{n_e e}\\vect{R} + \\frac{m_e}{n_e e^2}\\frac{\\partial \\',
    'label': 'eq:ch9-generalized-ohm'
})

# Line 4342 (equation): Generalized Ohm's law simplified
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\boxed{ \\vect{E} + \\vect{v}\\times\\vect{B} = \\eta\\vect{j} + \\frac{1}{n_e e}\\grad p_e - \\frac{1}{n_e ',
    'label': 'eq:ch9-generalized-ohm-simplified'
})

# Line 4356 (equation): Landau integral
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\int_{-\\infty}^{+\\infty} \\frac{\\partial F_0/\\partial v}{v-\\omega/k}\\,\\diff v',
    'label': 'eq:ch7-landau-integral'
})

# Line 4360 (equation): Landau contour
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\int_{-\\infty}^{+\\infty} \\frac{\\partial F_0/\\partial v}{v-\\omega/k}\\,\\diff v = \\mathcal{P}\\int_{-\\i',
    'label': 'eq:ch7-landau-contour'
})

# Line 4367 (equation): Landau condition
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\left.\\frac{\\partial F_0}{\\partial v}\\right|_{v=\\omega/k} < 0 \\quad (\\omega/k>0)',
    'label': 'eq:ch7-landau-condition'
})

# Line 4392 (equation): Alfven B1
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\vect{B}_1 = -\\frac{1}{\\ii\\omega}\\curl(v_1\\unitvec{x}\\times B_0\\unitvec{z}) = -\\frac{B_0 v_1}{\\ii\\o',
    'label': 'eq:ch9-alfven-b1'
})

# Line 4397 (equation): Alfven momentum
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\ii\\omega\\rho_0 v_1 = \\frac{1}{\\mu_0}(-\\ii k_z B_{1x})(-B_0) - 0 = \\frac{\\ii k_z B_0^2 k_z v_1}{\\mu',
    'label': 'eq:ch9-alfven-momentum'
})

# Line 4402 (equation): Alfven dispersion
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\boxed{\\omega = \\pm k_z v_A}',
    'label': 'eq:ch9-alfven-dispersion'
})

# Line 4415 (equation): Grad-Shafranov B
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\vect{B} = \\frac{1}{R}\\grad\\psi\\times\\unitvec{\\varphi} + \\frac{F(\\psi)}{R}\\unitvec{\\varphi}',
    'label': 'eq:ch10-grad-shafranov-b'
})

# Line 4423 (equation): Grad-Shafranov equation
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\nabla^{2*}\\psi = -\\mu_0 R^2 \\frac{\\diff p}{\\diff\\psi} - F\\frac{\\diff F}{\\diff\\psi}',
    'label': 'eq:ch10-grad-shafranov'
})

# Line 4427 (equation): Grad-Shafranov operator
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\nabla^{2*}\\psi \\equiv R^2\\divg\\left(\\frac{1}{R^2}\\grad\\psi\\right) = R\\frac{\\partial }{\\partial R}\\',
    'label': 'eq:ch10-grad-shafranov-operator'
})

# Line 4745 (equation): Lorentz force
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' m\\frac{\\diff\\vect{v}}{\\diff t} = q(\\vect{E}+\\vect{v}\\times\\vect{B}), \\qquad \\frac{\\diff\\vect{r}}{\\d',
    'label': 'eq:app-lorentz-force'
})

# Line 4830 (equation): Debye radial
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\frac{1}{r^2}\\frac{\\diff}{\\diff r}\\left(r^2\\frac{\\diff\\phi}{\\diff r}\\right) = \\frac{e n_0}{\\eps_0}\\',
    'label': 'eq:app-debye-radial'
})

# Line 4905 (equation): Bohm-Gross
main_file_modifications.append({
    'type': 'add_label_before_end',
    'env_type': 'equation',
    'content_snippet': ' \\omega^2 = \\omega_{pe}^2 + \\frac{3k_\\mathrm{B}T_e}{m_e}k^2 = \\omega_{pe}^2 + 3k^2v_{te}^2',
    'label': 'eq:app-bohm-gross'
})


# ============================================================
# MODIFICATIONS FOR ch10_supplement.tex
# ============================================================

ch10_modifications = []

# Lawson criterion (full)
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    n\tau \geq \frac{3k_B T}{\frac{E_{\alpha}}{4}\langle\sigma v\rangle - \frac{P_{rad}}{n^2} \cdot n}
\]""",
    'new': r"""\begin{equation}
    n\tau \geq \frac{3k_B T}{\frac{E_{\alpha}}{4}\langle\sigma v\rangle - \frac{P_{rad}}{n^2} \cdot n}
\label{eq:ch10-lawson-full}
\end{equation}"""
})

# Lawson criterion (simplified)
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    n\tau \geq \frac{12k_B T}{E_{\alpha}\langle\sigma v\rangle}
\]""",
    'new': r"""\begin{equation}
    n\tau \geq \frac{12k_B T}{E_{\alpha}\langle\sigma v\rangle}
\label{eq:ch10-lawson-simplified}
\end{equation}"""
})

# Triple product definition
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    nT\tau_E = n_e \times T_e \times \tau_E
\]""",
    'new': r"""\begin{equation}
    nT\tau_E = n_e \times T_e \times \tau_E
\label{eq:ch10-triple-product}
\end{equation}"""
})

# Tokamak scaling law
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \tau_E^{Goldston} = 0.04 \cdot H \cdot I_p^{1.25} \cdot B^{0.5} \cdot n^{0.25} \cdot a^{1.5} \cdot R^{1.75}
\]""",
    'new': r"""\begin{equation}
    \tau_E^{Goldston} = 0.04 \cdot H \cdot I_p^{1.25} \cdot B^{0.5} \cdot n^{0.25} \cdot a^{1.5} \cdot R^{1.75}
\label{eq:ch10-tokamak-scaling}
\end{equation}"""
})

# Stellarator scaling law
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \tau_E^{Stellarator} \approx 0.17 \cdot a^2 \cdot R^{0.5} \cdot n^{0.5} \cdot B^2
\]""",
    'new': r"""\begin{equation}
    \tau_E^{Stellarator} \approx 0.17 \cdot a^2 \cdot R^{0.5} \cdot n^{0.5} \cdot B^2
\label{eq:ch10-stellarator-scaling}
\end{equation}"""
})

# Another stellarator scaling law
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \tau_E \approx 0.1 \cdot a^2 \cdot R^{0.5} \cdot B^{0.5} \cdot n^{0.3}
\]""",
    'new': r"""\begin{equation}
    \tau_E \approx 0.1 \cdot a^2 \cdot R^{0.5} \cdot B^{0.5} \cdot n^{0.3}
\label{eq:ch10-stellarator-scaling2}
\end{equation}"""
})

# Tokamak volume
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    V_p \approx 2\pi R \times \pi a^2 \times \kappa
\]""",
    'new': r"""\begin{equation}
    V_p \approx 2\pi R \times \pi a^2 \times \kappa
\label{eq:ch10-tokamak-volume}
\end{equation}"""
})

# Safety factor q
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    q = \frac{r B_{\phi}}{R B_{\theta}}
\]""",
    'new': r"""\begin{equation}
    q = \frac{r B_{\phi}}{R B_{\theta}}
\label{eq:ch10-safety-factor}
\end{equation}"""
})

# Ampere's law for tokamak
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \oint \vect{B}_{\theta}\cdot\diff\vect{l} = \mu_0 I_{enc} \quad\Rightarrow\quad 2\pi a B_{\theta} = \mu_0 I_p
\]""",
    'new': r"""\begin{equation}
    \oint \vect{B}_{\theta}\cdot\diff\vect{l} = \mu_0 I_{enc} \quad\Rightarrow\quad 2\pi a B_{\theta} = \mu_0 I_p
\label{eq:ch10-ampere-law}
\end{equation}"""
})

# Edge poloidal field
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    B_{\theta}(a) = \frac{\mu_0 I_p}{2\pi a}
\]""",
    'new': r"""\begin{equation}
    B_{\theta}(a) = \frac{\mu_0 I_p}{2\pi a}
\label{eq:ch10-poloidal-field}
\end{equation}"""
})

# q_95 approximation
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    q_{95} \approx q_a \times \frac{(1+\kappa^2)}{2} \times \text{三角修正}
\]""",
    'new': r"""\begin{equation}
    q_{95} \approx q_a \times \frac{(1+\kappa^2)}{2} \times \text{三角修正}
\label{eq:ch10-q95-approx}
\end{equation}"""
})

# q_95 empirical formula
ch10_modifications.append({
    'type': 'replace',
    'old': r"""\[
    q_{95} \approx \frac{5 a^2 B_0}{R I_p} \cdot \frac{(1+\kappa^2)}{2}
\]""",
    'new': r"""\begin{equation}
    q_{95} \approx \frac{5 a^2 B_0}{R I_p} \cdot \frac{(1+\kappa^2)}{2}
\label{eq:ch10-q95-empirical}
\end{equation}"""
})

# ============================================================
# MODIFICATIONS FOR ch11_supplement.tex
# ============================================================

ch11_modifications = []

# Paschen's law
ch11_modifications.append({
    'type': 'replace',
    'old': r"""\[
    V_B = \frac{B \cdot p d}{\ln(A p d) - \ln\left(\ln\left(1+\frac{1}{\gamma}\right)\right)}
\]""",
    'new': r"""\begin{equation}
    V_B = \frac{B \cdot p d}{\ln(A p d) - \ln\left(\ln\left(1+\frac{1}{\gamma}\right)\right)}
\label{eq:ch11-paschen-law}
\end{equation}"""
})

# Bohm criterion
ch11_modifications.append({
    'type': 'replace',
    'old': r"""\[
    v_B = \sqrt{\frac{k_B T_e}{m_i}}
\]""",
    'new': r"""\begin{equation}
    v_B = \sqrt{\frac{k_B T_e}{m_i}}
\label{eq:ch11-bohm-criterion}
\end{equation}"""
})

# Child-Langmuir law
ch11_modifications.append({
    'type': 'replace',
    'old': r"""\[
    J_i = \frac{4}{9}\eps_0 \sqrt{\frac{2e}{m_i}} \frac{V_s^{3/2}}{s^2}
\]""",
    'new': r"""\begin{equation}
    J_i = \frac{4}{9}\eps_0 \sqrt{\frac{2e}{m_i}} \frac{V_s^{3/2}}{s^2}
\label{eq:ch11-child-langmuir}
\end{equation}"""
})

# Debye length
ch11_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \lambda_D = \sqrt{\frac{\eps_0 k_B T_e}{n_e e^2}} = 743\sqrt{\frac{T_e\,[\mathrm{eV}]}{n_e\,[\mathrm{cm^{-3}}]}}\,\mathrm{cm}
\]""",
    'new': r"""\begin{equation}
    \lambda_D = \sqrt{\frac{\eps_0 k_B T_e}{n_e e^2}} = 743\sqrt{\frac{T_e\,[\mathrm{eV}]}{n_e\,[\mathrm{cm^{-3}}]}}\,\mathrm{cm}
\label{eq:ch11-debye-length}
\end{equation}"""
})

# E×B drift
ch11_modifications.append({
    'type': 'replace',
    'old': r"""\[
    v_{E\times B} = \frac{E}{B}
\]""",
    'new': r"""\begin{equation}
    v_{E\times B} = \frac{E}{B}
\label{eq:ch11-exb-drift}
\end{equation}"""
})

# Electron thermal velocity
ch11_modifications.append({
    'type': 'replace',
    'old': r"""\[
    v_{th,e} = \sqrt{\frac{2 k_B T_e}{m_e}} = 5.93\times10^5 \sqrt{T_e\,[\mathrm{eV}]}
\]""",
    'new': r"""\begin{equation}
    v_{th,e} = \sqrt{\frac{2 k_B T_e}{m_e}} = 5.93\times10^5 \sqrt{T_e\,[\mathrm{eV}]}
\label{eq:ch11-thermal-velocity}
\end{equation}"""
})


# ============================================================
# MODIFICATIONS FOR ch12_supplement.tex
# ============================================================

ch12_modifications = []

# Langmuir probe current
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    I = I_{esat} \exp\left(\frac{e(V_p - V_{pl})}{k_B T_e}\right) - |I_{isat}|
\]""",
    'new': r"""\begin{equation}
    I = I_{esat} \exp\left(\frac{e(V_p - V_{pl})}{k_B T_e}\right) - |I_{isat}|
\label{eq:ch12-probe-current}
\end{equation}"""
})

# ln I formula
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \ln I = \ln I_{esat} + \frac{e(V_p - V_{pl})}{k_B T_e}
\]""",
    'new': r"""\begin{equation}
    \ln I = \ln I_{esat} + \frac{e(V_p - V_{pl})}{k_B T_e}
\label{eq:ch12-probe-lnI}
\end{equation}"""
})

# T_e from slope
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \frac{d(\ln I)}{dV_p} = \frac{e}{k_B T_e} = \frac{1}{T_e\,[\mathrm{eV}]}
\]""",
    'new': r"""\begin{equation}
    \frac{d(\ln I)}{dV_p} = \frac{e}{k_B T_e} = \frac{1}{T_e\,[\mathrm{eV}]}
\label{eq:ch12-probe-temperature}
\end{equation}"""
})

# OML current
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    I_{esat} = A_p n_e e \sqrt{\frac{k_B T_e}{2\pi m_e}}
\]""",
    'new': r"""\begin{equation}
    I_{esat} = A_p n_e e \sqrt{\frac{k_B T_e}{2\pi m_e}}
\label{eq:ch12-oml-current}
\end{equation}"""
})

# Plasma frequency
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \omega_p = \sqrt{\frac{n_e e^2}{\eps_0 m_e}}
\]""",
    'new': r"""\begin{equation}
    \omega_p = \sqrt{\frac{n_e e^2}{\eps_0 m_e}}
\label{eq:ch12-plasma-frequency}
\end{equation}"""
})

# Refractive index
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    N = \sqrt{1 - \frac{\omega_p^2}{\omega^2}} \approx 1 - \frac{\omega_p^2}{2\omega^2}
\]""",
    'new': r"""\begin{equation}
    N = \sqrt{1 - \frac{\omega_p^2}{\omega^2}} \approx 1 - \frac{\omega_p^2}{2\omega^2}
\label{eq:ch12-refractive-index}
\end{equation}"""
})

# Phase shift
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    |\Delta\phi| = \frac{\omega_p^2 L}{2\omega c}
\]""",
    'new': r"""\begin{equation}
    |\Delta\phi| = \frac{\omega_p^2 L}{2\omega c}
\label{eq:ch12-phase-shift}
\end{equation}"""
})

# Density from phase shift
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    n_e = \frac{4\pi\eps_0 m_e f c \cdot |\Delta\phi|}{e^2 L}
\]""",
    'new': r"""\begin{equation}
    n_e = \frac{4\pi\eps_0 m_e f c \cdot |\Delta\phi|}{e^2 L}
\label{eq:ch12-density-phase}
\end{equation}"""
})

# Thomson scattering FWHM
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \Delta\lambda_{1/2} = 2\lambda_0 \sqrt{\frac{8k_B T_e \ln 2}{m_e c^2}} \sin\frac{\theta}{2}
\]""",
    'new': r"""\begin{equation}
    \Delta\lambda_{1/2} = 2\lambda_0 \sqrt{\frac{8k_B T_e \ln 2}{m_e c^2}} \sin\frac{\theta}{2}
\label{eq:ch12-thomson-fwhm}
\end{equation}"""
})

# T_e from Thomson scattering
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    T_e = \frac{m_e c^2}{8k_B \ln 2} \left(\frac{\Delta\lambda_{1/2}}{2\lambda_0 \sin(\theta/2)}\right)^2
\]""",
    'new': r"""\begin{equation}
    T_e = \frac{m_e c^2}{8k_B \ln 2} \left(\frac{\Delta\lambda_{1/2}}{2\lambda_0 \sin(\theta/2)}\right)^2
\label{eq:ch12-thomson-temperature}
\end{equation}"""
})

# Faraday's law (Rogowski coil)
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \mathcal{E} = -N \frac{d\Phi}{dt} = -N A \frac{dB_{\theta}}{dt}
\]""",
    'new': r"""\begin{equation}
    \mathcal{E} = -N \frac{d\Phi}{dt} = -N A \frac{dB_{\theta}}{dt}
\label{eq:ch12-faraday-law}
\end{equation}"""
})

# Ampere's law
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \oint B_{\theta} dl = \mu_0 I_p \quad\Rightarrow\quad B_{\theta} = \frac{\mu_0 I_p}{2\pi R}
\]""",
    'new': r"""\begin{equation}
    \oint B_{\theta} dl = \mu_0 I_p \quad\Rightarrow\quad B_{\theta} = \frac{\mu_0 I_p}{2\pi R}
\label{eq:ch12-ampere-law}
\end{equation}"""
})

# Rogowski coil EMF
ch12_modifications.append({
    'type': 'replace',
    'old': r"""\[
    \mathcal{E} = -N A \frac{d}{dt}\left(\frac{\mu_0 I_p}{2\pi R}\right) = -\frac{N A \mu_0}{2\pi R} \frac{dI_p}{dt}
\]""",
    'new': r"""\begin{equation}
    \mathcal{E} = -N A \frac{d}{dt}\left(\frac{\mu_0 I_p}{2\pi R}\right) = -\frac{N A \mu_0}{2\pi R} \frac{dI_p}{dt}
\label{eq:ch12-rogowski-emf}
\end{equation}"""
})


def apply_modifications(fpath, modifications):
    with open(fpath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    changes = []
    for mod in modifications:
        if mod['type'] == 'replace':
            if mod['old'] in content:
                content = content.replace(mod['old'], mod['new'], 1)
                changes.append(f"Replaced: {mod['old'][:60]}... -> {mod['new'][:60]}...")
            else:
                changes.append(f"FAILED to find: {mod['old'][:60]}...")
        elif mod['type'] == 'add_label_before_end':
            env_type = mod['env_type']
            label = mod['label']
            content_snippet = mod['content_snippet']
            
            # Find the environment containing the content snippet
            # Search for the snippet within the environment
            pattern = re.compile(r'\\begin\{' + env_type + r'\}(.*?)\\end\{' + env_type + r'\}', re.DOTALL)
            found = False
            for m in pattern.finditer(content):
                inner = m.group(1)
                if content_snippet in inner:
                    # Check if already has label
                    if '\\label' not in inner:
                        # Replace \end{env_type} with \label{...}\n\end{env_type}
                        old_end = r'\end{' + env_type + r'}'
                        new_end = r'\label{' + label + r'}' + '\n' + old_end
                        # Need to replace the specific occurrence
                        start = m.start()
                        end = m.end()
                        env_text = content[start:end]
                        new_env_text = env_text.replace(old_end, new_end, 1)
                        content = content[:start] + new_env_text + content[end:]
                        changes.append(f"Added label {label} to {env_type} at position {start}")
                        found = True
                        break
                    else:
                        changes.append(f"Already has label in {env_type} containing {content_snippet[:40]}...")
                        found = True
                        break
            if not found:
                changes.append(f"FAILED to find {env_type} containing {content_snippet[:60]}...")
    
    with open(fpath, 'w', encoding='utf-8') as f:
        f.write(content)
    
    return changes


# Apply modifications
all_changes = {}

print("Processing plasma_physics_textbook_v2.tex...")
all_changes['plasma_physics_textbook_v2.tex'] = apply_modifications(
    os.path.join(base_dir, 'plasma_physics_textbook_v2.tex'),
    main_file_modifications
)

print("Processing ch10_supplement.tex...")
all_changes['ch10_supplement.tex'] = apply_modifications(
    os.path.join(base_dir, 'ch10_supplement.tex'),
    ch10_modifications
)

print("Processing ch11_supplement.tex...")
all_changes['ch11_supplement.tex'] = apply_modifications(
    os.path.join(base_dir, 'ch11_supplement.tex'),
    ch11_modifications
)

print("Processing ch12_supplement.tex...")
all_changes['ch12_supplement.tex'] = apply_modifications(
    os.path.join(base_dir, 'ch12_supplement.tex'),
    ch12_modifications
)

# Print summary
for fname, changes in all_changes.items():
    print(f"\n{'='*60}")
    print(f"{fname}: {len(changes)} modifications attempted")
    for c in changes:
        print(f"  {c}")

print("\nDone!")
