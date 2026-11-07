
from pathlib import Path
import argparse
import json
import sys
import os
import numpy as np

NAMES = ['x0_m', 'z0_m', 'width_2b_m', 'length_2L_m', 'alpha_deg', 'i_deg', 'M']
INITIAL = np.array([900., 1100., 500., 700., 50., 60., 2000.])
TRUE = np.array([1000., 1000., 400., 800., 60., 60., 2000.])
LOW = np.array([-2000., 10., 10., 10., 5., -89., 1.])
HIGH = np.array([4000., 5000., 2500., 4000., 175., 89., 10000.])


def corners(p):
    u = np.array([np.cos(np.deg2rad(p[4])), np.sin(np.deg2rad(p[4]))])
    t, b = p[:2] - p[3]/2*u, p[:2] + p[3]/2*u
    return np.array([t+[-p[2]/2, 0], t+[p[2]/2, 0],
                     b+[-p[2]/2, 0], b+[p[2]/2, 0]])


def valid(p, z):
    return np.all(p >= LOW) and np.all(p <= HIGH) and corners(p)[:, 1].min() > np.max(z)+1.


def xasu(p, x, z, component='ppt', unit='ppt'):
    """PPT 第15页 XASU，atan2 修复 ATAN 的象限歧义。
    ppt 单位直接沿用教材 M 的幅度约定；SI 模式 M 为 A/m，输出 nT。
    dx/dz/dt 对应原 two_dimensional_plate_forward.m 的几何与公式。
    """
    alpha, inc = np.deg2rad(p[4:6])
    if component == 'ppt':
        a = p[1]-z+p[3]*np.sin(alpha)/2
        b = p[1]-z-p[3]*np.sin(alpha)/2
        c = x-p[0]+p[3]*np.cos(alpha)/2
        d = x-p[0]-p[3]*np.cos(alpha)/2
        e, f, g, h = c+p[2]/2, c-p[2]/2, d+p[2]/2, d-p[2]/2
        logu = np.log(b*b+f*f)+np.log(a*a+g*g)-np.log(b*b+e*e)-np.log(a*a+h*h)
        v = np.arctan2(p[2]*b, b*b+c*c-p[2]**2/4)
        w = np.arctan2(p[2]*a, a*a+d*d-p[2]**2/4)
        out = p[6]/(2*np.pi)*np.sin(alpha)*(np.sin(alpha-inc)*logu/2+np.cos(alpha-inc)*(v-w))
    else:
        pts = corners(p)
        dx, dz = pts[:, 0, None]-x, pts[:, 1, None]-z
        r, phi = np.hypot(dx, dz), np.arctan2(dz, dx)
        q = phi[0]-phi[1]-phi[2]+phi[3]
        lr = np.log(r[1])+np.log(r[2])-np.log(r[0])-np.log(r[3])
        hx = p[6]/(2*np.pi)*np.sin(alpha)*(lr*np.cos(alpha-inc)-np.sin(alpha-inc)*q)
        hz = p[6]/(2*np.pi)*np.sin(alpha)*(lr*np.sin(alpha-inc)+np.cos(alpha-inc)*q)
        out = {'dx': hx, 'dz': hz, 'dt': hx*np.cos(inc)+hz*np.sin(inc)}[component]
    return out * (4*np.pi*1e2 if unit == 'si' else 1.)


def detu(p, x, z, active, component, unit):
    cols = []
    for j in active:
        h = 1e-4*max(abs(p[j]), 1.)
        pp, pm = p.copy(), p.copy()
        pp[j] += h
        pm[j] -= h
        cols.append((xasu(pp, x, z, component, unit)-xasu(pm, x, z, component, unit))/(2*h))
    return np.column_stack(cols)


def qf(obs, pred, sigma, background):
    r = (obs-pred-background)/sigma
    return float(r@r), r


def maqt(x, z, obs, sigma, p0=INITIAL, active=(0,1,2,3,4),
         component='ppt', unit='ppt', background=0., max_iter=200):
    """规格化后 (A*+lambda I)y=b*，dp=y/sqrt(diag(J'J))。"""
    p = np.asarray(p0, float).copy()
    active = np.asarray(active, int)
    if not valid(p, z):
        raise ValueError('初始板体必须位于全部测点下方且在参数边界内')
    if len(obs) <= len(active) or len(active) == 0:
        raise ValueError('测点数必须大于自由参数数，且至少有一个自由参数')
    lam, history, status = .1, [], 'max_iter'
    pred = xasu(p, x, z, component, unit)
    f, r = qf(obs, pred, sigma, background)
    history.append([0, f, lam, *p])
    for it in range(1, max_iter+1):
        j = detu(p, x, z, active, component, unit)/sigma[:, None]
        scale = np.maximum(np.linalg.norm(j, axis=0), 1e-15)
        js = j/scale
        numerical_floor = 100*np.finfo(float).eps*max(1., np.max(np.abs(obs-background)))
        if np.max(np.abs(js.T@r)) < 1e-9 or f < 1e-20 or np.max(np.abs(obs-pred-background)) <= numerical_floor:
            status = 'gradient_or_exact_fit'
            break
        accepted = False
        for _ in range(24):
            a = js.T@js+lam*np.eye(len(active))
            chol = np.linalg.cholesky(a)
            y = np.linalg.solve(chol.T, np.linalg.solve(chol, js.T@r))
            trial = p.copy()
            trial[active] += y/scale
            if valid(trial, z):
                tp = xasu(trial, x, z, component, unit)
                tf, tr = qf(obs, tp, sigma, background)
                if tf < f:
                    relative = (f-tf)/max(f, 1e-30)
                    p, pred, f, r = trial, tp, tf, tr
                    history.append([it, f, lam, *p])
                    lam = max(lam/10, 1e-12)
                    accepted = True
                    break
            lam *= 10
            if lam > 1e16:
                break
        if not accepted:
            status = 'no_acceptable_step'
            break
        if relative < 1e-9:
            status = 'relative_objective'
            break
    j = detu(p, x, z, active, component, unit)/sigma[:, None]
    sv = np.linalg.svd(j/np.maximum(np.linalg.norm(j, axis=0), 1e-15), compute_uv=False)
    return p, pred+background, np.array(history), status, sv


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data', type=Path, help='CSV: x_m,z_m,signal[,std]，首行为表头')
    ap.add_argument('--component', choices=['ppt','dx','dz','dt'], default='ppt')
    ap.add_argument('--unit', choices=['ppt','si'], default='ppt')
    ap.add_argument('--initial', nargs=7, type=float, default=INITIAL)
    ap.add_argument('--active', nargs='+', type=int, default=[1,2,3,4,5], help='自由参数编号1至7')
    ap.add_argument('--noise', type=float, default=.01, help='模拟噪声标准差/信号峰值')
    ap.add_argument('--background', type=float, default=0.)
    ap.add_argument('--out', type=Path, default=Path(__file__).parent/'inversion_results_python')
    ap.add_argument('--no-plot', action='store_true')
    args = ap.parse_args()
    active = np.array(args.active)-1
    if len(set(active)) != len(active) or np.any((active < 0)|(active > 6)):
        ap.error('active 必须为不重复的1至7')
    if args.noise < 0:
        ap.error('noise 必须非负')
    truth = None
    if args.data:
        data = np.loadtxt(args.data, delimiter=',', skiprows=1, ndmin=2)
        if data.shape[1] not in (3,4):
            ap.error('CSV 必须有3或4列')
        x, z, obs = data[:, :3].T
        sigma = data[:, 3] if data.shape[1] == 4 else np.ones(len(x))
    else:
        truth = TRUE.copy()
        x, z = np.arange(0., 2001., 20.), np.zeros(101)
        clean = xasu(truth, x, z, args.component, args.unit)
        sd = args.noise*np.max(np.abs(clean))
        sigma = np.full(len(x), max(sd, 1e-10))
        obs = clean+args.background+np.random.default_rng(2026).normal(0, sd, len(x))
    if not all(np.all(np.isfinite(a)) for a in (x,z,obs,sigma)) or np.any(sigma <= 0):
        ap.error('输入必须有限，std 必须为正')
    p, pred, hist, status, sv = maqt(x,z,obs,sigma,args.initial,active,args.component,args.unit,args.background)
    args.out.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.out/'observations.csv', np.c_[x,z,obs,sigma], delimiter=',', header='x_m,z_m,signal,std', comments='')
    np.savetxt(args.out/'fit.csv', np.c_[x,z,obs,pred,obs-pred], delimiter=',', header='x_m,z_m,observed,predicted,residual', comments='')
    np.savetxt(args.out/'history.csv', hist, delimiter=',', header='iteration,weighted_SSE,lambda,'+','.join(NAMES), comments='')
    result = dict(parameters=dict(zip(NAMES,p.tolist())), initial=list(args.initial), truth=None if truth is None else truth.tolist(), active=args.active, component=args.component, unit=args.unit, status=status, rmse=float(np.sqrt(np.mean((obs-pred)**2))), scaled_jacobian_singular_values=sv.tolist())
    (args.out/'result.json').write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(result,indent=2,ensure_ascii=False))
    print('\n反演完成。结果目录：', args.out.resolve())
    print('课堂只需运行本文件；verify_plate_inversion.py 是可选验收脚本。')
    if not args.no_plot:
        os.environ.setdefault('MPLCONFIGDIR', str(args.out/'.matplotlib'))
        import matplotlib
        # Spyder/IPython 保留其交互后端；命令行使用文件输出后端。
        interactive = 'IPython' in sys.modules
        if not interactive:
            matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2,2,figsize=(12,8))
        ax[0,0].plot(x,obs,'.',label='Observed'); ax[0,0].plot(x,pred,label='Inverted')
        ax[0,0].plot(x,xasu(np.array(args.initial),x,z,args.component,args.unit)+args.background,'--',label='Initial')
        ax[0,0].set_ylabel('nT' if args.unit == 'si' else 'Textbook amplitude'); ax[0,0].legend()
        ax[0,1].plot(x,obs-pred); ax[0,1].set_ylabel('Residual')
        for model,label in [(np.array(args.initial),'Initial'),(p,'Inverted')]+([] if truth is None else [(truth,'True')]):
            c=corners(model)[[0,1,3,2,0]]
            ax[1,0].plot(c[:,0],c[:,1],label=label)
        ax[1,0].plot(x,z,'k.',ms=2); ax[1,0].invert_yaxis(); ax[1,0].axis('equal'); ax[1,0].legend(); ax[1,0].set_ylabel('Depth z (m)')
        ax[1,1].semilogy(hist[:,0],np.maximum(hist[:,1],1e-30),'o-'); ax[1,1].set_xlabel('Iteration'); ax[1,1].set_ylabel('Weighted SSE')
        for a in ax.flat: a.grid(True)
        for a in [ax[0,0],ax[0,1],ax[1,0]]: a.set_xlabel('x (m)')
        fig.tight_layout(); fig.savefig(args.out/'inversion.png',dpi=180)
        print('结果图：', (args.out/'inversion.png').resolve())
        if interactive:
            plt.show()
        else:
            plt.close(fig)


if __name__ == '__main__':
    main()
