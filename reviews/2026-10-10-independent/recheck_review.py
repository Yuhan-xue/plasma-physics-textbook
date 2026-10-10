"""Independent review evidence, not a textbook release gate.
Run from any directory: python reviews/2026-10-10-independent/recheck_review.py
Requires numpy, scipy, sympy, pymupdf. Writes only to --out (default tmp/review-20261010).
Exit 0 means the documented baseline counterexamples were reproduced, NOT that the book is correct.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
import sympy as s
from scipy.integrate import quad
import pymupdf

ROOT=Path(__file__).resolve().parents[2]
BASE='f6317daf8b4203b96674304e08471caccb9f895f'

def complex_pair(z): return [float(z.real),float(z.imag)]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=ROOT/'tmp/review-20261010')
    out=parser.parse_args().out.resolve();out.mkdir(parents=True,exist_ok=True)
    data={'baseline':BASE,'meaning':'evidence reproduction; exit 0 is not release approval'}
    w,a=s.symbols('w a',real=True)
    D=1-s.Rational(1,2)/(w-a)**2-s.Rational(1,2)/(w+a)**2
    numerator=s.factor(s.together(D)).as_numer_denom()[0]
    expected=w**4-(2*a*a+1)*w*w+a**4-a*a
    assert s.expand(numerator-expected)==0
    cases=[]
    for av in [.5,math.sqrt(3/8),1.,1.2,1.5]:
        roots=np.roots([1,0,-(2*av**2+1),0,av**4-av**2])
        cases.append({'a':av,'roots':[complex_pair(x) for x in roots],
                      'max_growth':float(max(x.imag for x in roots))})
    av=1.2;claimed_growth=math.sqrt(math.sqrt(1+4*av*av)-av*av-1)
    z=av+1j*claimed_growth
    residual=abs(1-.5/(z-av)**2-.5/(z+av)**2)
    wp=5e10;u=2e7
    data['two_stream']={'cleared_polynomial':str(expected),'cases':cases,
        'book_a_1_2_growth':claimed_growth,'book_root_residual_in_original_D':residual,
        'correct_k_threshold':2*wp/u,'correct_lambda_min':2*math.pi/(2*wp/u),
        'correct_gamma_max':wp/math.sqrt(8),'correct_k_opt':math.sqrt(1.5)*wp/u,
        'correct_growth_time_s':math.sqrt(8)/wp}
    assert cases[3]['max_growth']<1e-10 and claimed_growth>.39 and residual>.1
    bun=[]
    for mass in [100,1836,10000,1000000,100000000]:
        mu=1/mass;roots=np.roots([1,-2,-mu,2*mu,-mu]);root=max(roots,key=lambda z:z.imag)
        asym=math.sqrt(3)/2*(mu/2)**(1/3)
        bun.append({'mi_me':mass,'root':[float(root.real),float(root.imag)],
                    'asymptotic_gamma_wp':asym,'gamma_over_asymptotic':float(root.imag/asym),
                    'residual':float(abs(1-1/(root-1)**2-mu/root**2))})
    data['buneman']={'at':'k u_e / omega_pe = 1','cases':bun,
        'log_mass_slope_last_two':math.log(bun[-1]['root'][1]/bun[-2]['root'][1])/math.log(100)}
    assert max(x['residual'] for x in bun)<1e-10
    assert abs(data['buneman']['log_mass_slope_last_two']+1/3)<.002
    v,c=s.symbols('v c',nonzero=True)
    cold=s.simplify(-s.diff(1/(v-c),v).subs(v,0))
    assert cold==1/c**2
    data['cold_delta']={'integral_delta_prime_over_v_minus_c':str(cold),
                         'book_sign':'negative','correct_sign':'positive'}
    T=4.;norm=2/(math.sqrt(math.pi)*T**1.5)
    p=lambda e:norm*math.sqrt(e)*math.exp(-e/T)
    g=lambda e:norm*math.exp(-e/T)
    V=np.linspace(-12,-1,101);I2=np.exp(V/T)/T**2
    recovered=1/np.polyfit(V,np.log(I2),1)[0]
    data['energy_distributions']={'T_eV':T,'integral_normalized_EEDF':quad(p,0,np.inf)[0],
        'integral_sqrt_energy_times_EEPF':quad(lambda e:math.sqrt(e)*g(e),0,np.inf)[0],
        'integral_EEPF_unweighted':quad(g,0,np.inf)[0],
        'temperature_from_log_I_second_derivative':float(recovered),
        'energy_eV':[2,4,8],'log_EEDF_slopes_per_eV':[1/(2*e)-1/T for e in [2,4,8]]}
    assert abs(recovered-4)<1e-12
    eps0=8.8541878128e-12;ec=1.602176634e-19;n=1e20;temp=10.;wave=1064e-9
    ld=math.sqrt(eps0*temp/(n*ec));nd=4*math.pi/3*n*ld**3
    data['thomson']={'n_m3':n,'T_eV':temp,'wavelength_m':wave,'lambda_D_m':ld,'N_D':nd,
        'angles':[{'angle_deg':ang,'alpha':1/(4*math.pi/wave*math.sin(math.radians(ang)/2)*ld)} for ang in [1,90]]}
    assert nd>1000 and data['thomson']['angles'][0]['alpha']>1 and data['thomson']['angles'][1]['alpha']<1
    theta=math.radians(60)
    data['oblique_reflection']={'incident_angle_from_vertical_deg':60,
        'correct_omega_over_local_wp':1/math.cos(theta),'book_omega_over_local_wp':math.cos(theta),
        'correct_kz_squared_normalized':(1/math.cos(theta))**2*math.cos(theta)**2-1,
        'book_kz_squared_normalized':math.cos(theta)**4-1}
    data['pitch_angle']={'correct_tan_squared':2,'correct_angle_degrees':math.degrees(math.atan(math.sqrt(2))),
        'book_tan_squared':.5,'book_angle_degrees':math.degrees(math.atan(math.sqrt(.5))),
        'gradB_over_curvature_at_book_angle':.5/2}
    data['efolding']={'t_10_efolds_over_t_tenfold':10/math.log(10),'book_ratio':math.log(10)}
    # Isolate the production log reader by pointing its SRC at a disposable fixture directory.
    # No production log/source is overwritten and no compiler is invoked for these controls.
    spec=importlib.util.spec_from_file_location('review_target_full_check',ROOT/'check/full_check.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    fixtures=out/'log-fixtures';fixtures.mkdir(exist_ok=True);mod.SRC=str(fixtures)
    variants={'clean':'This is XeTeX\n',
              'known_newline':'Missing character: There is no ^^J (U+000A) in font Test!\n',
              'real_missing_glyph':'Missing character: There is no Ω (U+03A9) in font Test!\n',
              'huge_overfull':'Overfull \\hbox (200.0pt too wide) in paragraph at lines 1--2\n',
              'real_tex_error':'This is XeTeX\n! Undefined control sequence.\n'}
    data['compile_gate_controls']={}
    for name,body in variants.items():
        (fixtures/mod.MAIN.replace('.tex','.log')).write_text(body,encoding='utf-8')
        mod.RESULTS.clear();mod.check_C_compile(False)
        data['compile_gate_controls'][name]=list(mod.RESULTS)
    assert any(r[0]=='FAIL' for r in data['compile_gate_controls']['real_tex_error'])
    for name in ['real_missing_glyph','huge_overfull']:
        assert not any(r[0] in ('FAIL','WARN') for r in data['compile_gate_controls'][name])
    pdf=ROOT/'pdf/plasma_physics_textbook_v2.pdf';doc=pymupdf.open(pdf)
    data['published_pdf']={'pages':len(doc),'bookmarks':len(doc.get_toc()),
                            'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()}
    terms=['冷分布极限','总电子密度','恰好一半','后者归一到1','当密度升高时','斜入射时有效判据','投掷角满足','第3–4章']
    data['pdf_text_hits']={term:[i+1 for i,pag in enumerate(doc) if term in ''.join(pag.get_text().split())] for term in terms}
    doc.close()
    data['source_hashes']={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [ROOT/'src/plasma_physics_textbook_v2.tex',ROOT/'src/ch6_supplement.tex',ROOT/'src/ch9_supplement.tex',ROOT/'check/full_check.py']}
    (out/'verification.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS: documented mathematical counterexamples and isolated gate controls reproduced')
    print('This is NOT a release approval. Evidence:',out/'verification.json')

if __name__=='__main__':main()
