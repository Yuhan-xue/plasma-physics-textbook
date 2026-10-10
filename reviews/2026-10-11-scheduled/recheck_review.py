"""只读独立复核，门禁变异在内存进行；退出 0 不等于教材已修正。"""
from pathlib import Path
import importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
text=(ROOT/'src/appendix_k_problems.tex').read_text(encoding='utf-8')
e,eps0,mu0=1.602176634e-19,8.8541878128e-12,1.25663706212e-6
# 柱坐标取 e_R=(1,0,0), e_phi=(0,1,0)，曲率指向 -e_R。
b=np.array([0.,1.,0.]); curvature=np.array([-1.,0.,0.])
drift=np.cross(b,curvature)
assert np.dot(drift,curvature)==0 and np.allclose(drift,[0,0,1])
# 用例题自身的 Bohm 通量和 Child-Langmuir 公式闭合计算。
T,n0,ns,M,V=3.,1e16,.61e16,6.68e-26,100.
cs=math.sqrt(e*T/M); ld=math.sqrt(eps0*T/(n0*e)); j=e*ns*cs
s_book=math.sqrt(2)/3*ld*(2*V/T)**.75
s_cl=math.sqrt(4/9*eps0*math.sqrt(2*e/M)*V**1.5/j)
j_from_book=4/9*eps0*math.sqrt(2*e/M)*V**1.5/s_book**2
assert math.isclose(s_cl/s_book,math.sqrt(n0/ns),rel_tol=1e-12)
assert math.isclose(j_from_book/j,n0/ns,rel_tol=1e-12)
# N_D 不含 k,omega，不能替代流体闭合；非相对论性不含 eta,L，也不能替代理想性。
eta=5.2e-5*10/10**1.5
vA=7.69e6
# 静止导电介质中 B_y=B0*cos(x/L)*exp(-t/tau) 的电阻扩散解。
# 流速为零仍非相对论，但磁通不冻结在静止流体上。
tau=mu0/eta
fraction=math.exp(-(tau*math.log(2))/tau)
assert math.isclose(fraction,.5) and vA/299792458<.03
# 现有 P01 丢符号且跨整篇检索：把阻尼答案变为增长，检查仍 PASS。
spec=importlib.util.spec_from_file_location('gate',ROOT/'check/physics_regression.py')
gate=importlib.util.module_from_spec(spec); spec.loader.exec_module(gate)
orig_read=gate.read
def run(mutant):
    rows=[]
    gate.rec=lambda *args: rows.append(args)
    def read(path):
        s=orig_read(path)
        if mutant and Path(path).name=='appendix_k_problems.tex':
            s=s.replace(r'=-2.00\times10^{-2}',r'=+2.00\times10^{-2}')
        return s
    gate.read=read
    gate.check_P01_new_appendix_numbers()
    return rows
baseline,mutated=run(False),run(True)
assert '=-2.00\\times10^{-2}' in text
assert baseline==mutated and not any(r[0]=='FAIL' for r in mutated)
# 上轮关键问题的独立数值正对照。
a=.75; w2=a*a+.5-.5*math.sqrt(1+8*a*a)
roots=np.roots([1,-2,-1/1836,2/1836,-1/1836])
root=max(roots,key=lambda z:z.imag)
assert abs(root.imag-.05431458)<1e-8
assert abs(1-1/(root-1)**2-(1/1836)/root**2)<1e-12
# 上轮 PP-C04 的归一 EEDF 与 EEPF 仍不能用常数互换。
energies=np.array([.25,1.,4.]); temperature=4.
p=2/math.sqrt(math.pi)*temperature**-1.5*np.sqrt(energies)*np.exp(-energies/temperature)
g=p/np.sqrt(energies)
assert np.allclose(p/g,[.5,1.,2.])
main=(ROOT/'src/plasma_physics_textbook_v2.tex').read_text(encoding='utf8')
result={
 'meaning':'PASS = 当前缺陷反例与旧项正对照复现；未修改教材',
 'P-C01':dict(source_marker='曲率漂移沿磁力线的曲率中心方向' in text,curvature=curvature.tolist(),drift=drift.tolist(),dot=float(drift@curvature)),
 'P-C02':dict(s_book_m=s_book,s_consistent_CL_m=s_cl,wall_current_Am2=j,current_implied_by_book_Am2=j_from_book,thickness_ratio=s_cl/s_book),
 'P-C03':dict(source_marker='磁冻结图像在非相对论等离子体中非常好' in text,book_vA_over_c=vA/299792458,resistive_example_eta=eta,diffusion_time_L1m_s=tau,static_fluid_flux_fraction_after_tau_ln2=fraction,note='电阻感应方程的静止流体反例；题目自身缺少温度、电阻率及流动尺度，不能据 vA/c 判定理想性'),
 'P-W01':dict(mutation='K.5 最终阻尼率 -2.00e-2 改成 +2.00e-2',baseline=baseline,mutated=mutated),
 'PP-C04_reopened':dict(source_still_equates_EEPF_to_normalized_EEDF='f_\\eps=F/n_e' in main,epsilon_eV=energies.tolist(),normalized_EEDF=p.tolist(),EEPF=g.tolist(),ratio=(p/g).tolist(),reference='https://cantera.org/3.2/cxx/d6/ded/PlasmaPhase_8h_source.html'),
 'PP-C06_partial':dict(wrong_appendix_marker_present=r'\alpha\ll1$ 为集体散射' in main,main_text_correct_marker_present=r'\alpha\ll1$ 为\textbf{非相干散射}' in main,definition='alpha=1/(q*lambda_D)',alpha_small_means='q*lambda_D much larger than 1, noncollective'),
 'previous':dict(two_stream_w2=w2,two_stream_growth=math.sqrt(-w2),buneman_ku_over_wpe=1,buneman_root=[root.real,root.imag],pitch_angle_deg=math.degrees(math.atan(math.sqrt(2))),efolding_ratio=10/math.log(10))
}
(HERE/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
