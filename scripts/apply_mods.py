import json
import re
import os

base_dir = r"C:\Users\53150\Documents\kimi\workspace\plasma_physics_textbook"

with open(os.path.join(base_dir, 'modifications.json'), 'r', encoding='utf-8') as f:
    file_mods = json.load(f)

def add_label_to_env(content, env_type, content_snippet, label):
    """Find an environment containing the snippet and add label before \end{env_type}"""
    pattern = re.compile(r'\\begin\{' + env_type + r'\}(.*?)\\end\{' + env_type + r'\}', re.DOTALL)
    for m in pattern.finditer(content):
        inner = m.group(1)
        if content_snippet in inner:
            if '\\label' not in inner:
                old_end = r'\end{' + env_type + r'}'
                new_end = r'\label{' + label + r'}' + '\n' + old_end
                start = m.start()
                end = m.end()
                env_text = content[start:end]
                new_env_text = env_text.replace(old_end, new_end, 1)
                content = content[:start] + new_env_text + content[end:]
                return content, True
            else:
                return content, False
    return content, False

def main(ctx):
    results = {}
    
    # Process JSON modifications (replace type)
    for item in file_mods:
        fname = item['file']
        fpath = os.path.join(base_dir, fname)
        with open(fpath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        changes = []
        for mod in item['modifications']:
            if mod['type'] == 'replace':
                old = mod['old']
                new = mod['new']
                if old in content:
                    content = content.replace(old, new, 1)
                    changes.append(f"Replaced: {old[:50]}...")
                else:
                    changes.append(f"FAILED: {old[:50]}...")
        
        with open(fpath, 'w', encoding='utf-8') as f:
            f.write(content)
        
        results[fname] = changes
    
    # Now add labels to numbered equations in main file
    main_path = os.path.join(base_dir, 'plasma_physics_textbook_v2.tex')
    with open(main_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    label_additions = [
        ('align', ' f(\\vect{r},t) &= \\int \\frac{\\diff^3 k\\,\\diff\\omega}{(2\\pi)^4}\\; \\tilde{f}(\\vect{k},\\omega)\\,\\ee^{+\\ii(\\vect{k}\\cdot\\vect{r}-\\omega t)}, \\\\\\\\', 'eq:ch1-fourier-transform'),
        ('align', ' n_i(x) &= n_0 \\left(1 - \\frac{2e\\phi(x)}{m_i u_0^2}\\right)^{-1/2} \\\\\\\\ n_e(x) &= n_0 \\exp\\!\\left(\\frac{e\\phi(x)}{k_B T_e}\\right)', 'eq:ch11-sheath-densities'),
        ('equation', ' \\eps_0 \\frac{\\diff^2\\phi}{\\diff x^2} = e(n_e - n_i)', 'eq:ch11-sheath-poisson'),
        ('equation', ' \\frac{1}{k_B T_e} \\ge \\frac{m_i u_0^2}{k_B T_e} \\quad\\Rightarrow\\quad u_0 \\ge \\sqrt{\\frac{k_B T_e}{', 'eq:ch11-bohm-criterion'),
        ('align', ' \\frac{1}{2}m_i u^2(x) + e\\phi(x) &= \\frac{1}{2}m_i u_B^2 \\quad (\\text{能量守恒}) \\\\\\\\ n_i(x) u(x) &= n_s', 'eq:ch11-sheath-equations'),
        ('equation', ' \\frac{\\diff^2\\phi}{\\diff x^2} = -\\frac{e n_i}{\\eps_0}', 'eq:ch11-sheath-poisson2'),
        ('equation', ' \\ln I = \\ln I_{\\mathrm{e,sat}} + \\frac{V - V_p}{T_e}.', 'eq:ch12-probe-current'),
        ('equation', ' n_r \\approx 1 - \\frac{1}{2}\\frac{\\omega_{pe}^2}{\\omega^2} = 1 - \\frac{e^2}{2\\eps_0 m_e \\omega^2} n_', 'eq:ch12-refractive-index'),
        ('equation', ' \\rho_m \\left(\\frac{\\partial \\vect{v}}{\\partial t} + (\\vect{v}\\cdot\\grad)\\vect{v}\\right) = \\vect{j}\\', 'eq:ch9-mhd-momentum'),
        ('equation', ' \\vect{E} + \\vect{v}_i\\times\\vect{B} = \\frac{1}{n_e e}\\vect{R} + \\frac{m_e}{n_e e^2}\\frac{\\partial \\', 'eq:ch9-generalized-ohm'),
        ('equation', ' \\boxed{ \\vect{E} + \\vect{v}\\times\\vect{B} = \\eta\\vect{j} + \\frac{1}{n_e e}\\grad p_e - \\frac{1}{n_e ', 'eq:ch9-generalized-ohm-simplified'),
        ('equation', ' \\int_{-\\infty}^{+\\infty} \\frac{\\partial F_0/\\partial v}{v-\\omega/k}\\,\\diff v', 'eq:ch7-landau-integral'),
        ('equation', ' \\int_{-\\infty}^{+\\infty} \\frac{\\partial F_0/\\partial v}{v-\\omega/k}\\,\\diff v = \\mathcal{P}\\int_{-\\i', 'eq:ch7-landau-contour'),
        ('equation', ' \\left.\\frac{\\partial F_0}{\\partial v}\\right|_{v=\\omega/k} < 0 \\quad (\\omega/k>0)', 'eq:ch7-landau-condition'),
        ('equation', ' \\vect{B}_1 = -\\frac{1}{\\ii\\omega}\\curl(v_1\\unitvec{x}\\times B_0\\unitvec{z}) = -\\frac{B_0 v_1}{\\ii\\o', 'eq:ch9-alfven-b1'),
        ('equation', ' \\ii\\omega\\rho_0 v_1 = \\frac{1}{\\mu_0}(-\\ii k_z B_{1x})(-B_0) - 0 = \\frac{\\ii k_z B_0^2 k_z v_1}{\\mu', 'eq:ch9-alfven-momentum'),
        ('equation', ' \\boxed{\\omega = \\pm k_z v_A}', 'eq:ch9-alfven-dispersion'),
        ('equation', ' \\vect{B} = \\frac{1}{R}\\grad\\psi\\times\\unitvec{\\varphi} + \\frac{F(\\psi)}{R}\\unitvec{\\varphi}', 'eq:ch10-grad-shafranov-b'),
        ('equation', ' \\nabla^{2*}\\psi = -\\mu_0 R^2 \\frac{\\diff p}{\\diff\\psi} - F\\frac{\\diff F}{\\diff\\psi}', 'eq:ch10-grad-shafranov'),
        ('equation', ' \\nabla^{2*}\\psi \\equiv R^2\\divg\\left(\\frac{1}{R^2}\\grad\\psi\\right) = R\\frac{\\partial }{\\partial R}\\', 'eq:ch10-grad-shafranov-operator'),
        ('equation', ' m\\frac{\\diff\\vect{v}}{\\diff t} = q(\\vect{E}+\\vect{v}\\times\\vect{B}), \\qquad \\frac{\\diff\\vect{r}}{\\d', 'eq:app-lorentz-force'),
        ('equation', ' \\frac{1}{r^2}\\frac{\\diff}{\\diff r}\\left(r^2\\frac{\\diff\\phi}{\\diff r}\\right) = \\frac{e n_0}{\\eps_0}\\', 'eq:app-debye-radial'),
        ('equation', ' \\omega^2 = \\omega_{pe}^2 + \\frac{3k_\\mathrm{B}T_e}{m_e}k^2 = \\omega_{pe}^2 + 3k^2v_{te}^2', 'eq:app-bohm-gross'),
    ]
    
    main_labels = []
    for env_type, snippet, label in label_additions:
        content, success = add_label_to_env(content, env_type, snippet, label)
        if success:
            main_labels.append(f"Added label {label}")
        else:
            main_labels.append(f"FAILED to add label {label}")
    
    with open(main_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    results['plasma_physics_textbook_v2.tex'] = results.get('plasma_physics_textbook_v2.tex', []) + main_labels
    
    for fname, changes in results.items():
        print(f"\n{'='*60}")
        print(f"{fname}: {len(changes)} modifications")
        for c in changes:
            print(f"  {c}")
    
    return results
