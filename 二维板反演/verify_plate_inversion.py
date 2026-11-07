"""数值验收：无噪声恢复、起伏测点、分量/单位、非唯一性及跨语言输出。"""
from pathlib import Path
import json
import numpy as np
from two_dimensional_plate_inversion import TRUE, INITIAL, xasu, maqt

base=Path(__file__).parent
x=np.arange(0.,2001.,20.)
checks=[]
for component in ['ppt','dx','dz','dt']:
    for terrain in [False,True]:
        z=35*np.sin(x/400) if terrain else np.zeros_like(x)
        obs=xasu(TRUE,x,z,component)
        p,pred,h,status,sv=maqt(x,z,obs,np.ones_like(x),component=component)
        err=float(np.max(np.abs(p[:5]-TRUE[:5])))
        assert err<.01,(component,terrain,err)
        assert np.all(np.diff(h[:,1])<0)
        np.testing.assert_allclose(xasu(TRUE,x,z,component,'si'),obs*4*np.pi*100)
        checks.append(dict(case=component,terrain=terrain,max_parameter_error=err,status=status))
# 七参数试验只要求解释数据，不能把小残差当作唯一几何恢复。
z=np.zeros_like(x); obs=xasu(TRUE,x,z)
p,pred,h,status,sv=maqt(x,z,obs,np.ones_like(x),active=range(7))
assert np.max(np.abs(obs-pred))<1e-5
checks.append(dict(case='seven_parameters',max_signal_error=float(np.max(np.abs(obs-pred))),smallest_scaled_singular_value=float(sv[-1]),parameters=p.tolist()))
# 宽板浅埋时，验证 atan2 仍连续（旧 ATAN 会有象限问题）。
shallow=np.array([1000.,400.,1000.,400.,60.,60.,2000.])
fine=np.linspace(-1000,3000,10001)
signal=xasu(shallow,fine,np.zeros_like(fine))
assert np.max(np.abs(np.diff(signal))) < .01*np.ptp(signal)
checks.append(dict(case='shallow_wide_atan2',max_adjacent_change=float(np.max(np.abs(np.diff(signal))))))
matdir=base/'inversion_results_matlab'
needed=[matdir/name for name in ['fit.csv','observations.csv','result.json']]
if all(path.exists() for path in needed):
    # 必须复用 MATLAB 本次实际使用的数据和设置，不能假定它用了示例CSV。
    data=np.loadtxt(matdir/'observations.csv',delimiter=',',skiprows=1,ndmin=2)
    mf=np.loadtxt(matdir/'fit.csv',delimiter=',',skiprows=1,ndmin=2)
    config=json.loads((matdir/'result.json').read_text(encoding='utf-8-sig'))
    np.testing.assert_allclose(mf[:,:3],data[:,:3],rtol=1e-12,atol=1e-10,
        err_msg='MATLAB输出文件来自不同运行，请重新运行MATLAB后再验收')
    p,pred,h,status,sv=maqt(*data[:,:3].T,data[:,3],
        p0=config['initial'],active=np.asarray(config['active']).reshape(-1)-1,
        component=config['component'],unit=config['unit'],
        background=config.get('background',0.),max_iter=config.get('max_iter',200))
    gap=float(np.max(np.abs(pred-mf[:,3])))
    assert gap<1e-5*max(1.,np.max(np.abs(mf[:,3]))),gap
    checks.append(dict(case='matlab_python_same_data',input='MATLAB observations.csv',max_prediction_difference=gap))
else:
    checks.append(dict(case='matlab_python_same_data',status='skipped',
        reason='尚无完整MATLAB结果；Python本身的数值检查已执行'))
out=base/'verification_report.json'
out.write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
print(out.read_text(encoding='utf-8'))
