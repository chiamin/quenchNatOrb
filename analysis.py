import pylab as pl
from collections import OrderedDict
import sys, glob
from math import pi, acos
sys.path.append ('/home/chiamin/mypy/')
import plotsetting as ps
import numpy as np
import fitfun as ff
import matplotlib.colors as colors
import cmasher as cmr

def exactG (V, tp):
    if V != 0:
        return 0.5 * pi / acos(-0.5*V)
    else:
        return 4*tp**2/(1+tp**2)**2

def get_para (fname, key, typ, last=False, n=1):
    with open(fname) as f:
        for line in f:
            if key in line:
                val = list(map(typ,line.split()[-n:]))
                if n == 1: val = val[0]
                if not last:
                    return val
        return val

'''def get_hop_t (fname):
    L_lead = get_para (fname, 'L_lead', int)
    L_device = get_para (fname, 'L_device', int)
    L = 2*L_lead + L_device
    ts = np.full (L, np.nan)
    with open(fname) as f:
        for line in f:
            if 'H left lead' in line:
                part = 'L'
            elif 'H right lead' in line:
                part = 'R'
            elif 'H dev' in line:
                part = 'S'
            elif 'Hk, t' in line:
                if part == 'L': offset = 0
                elif part == 'S': offset = L_lead
                elif part == 'R': offset = L_lead + L_device
                tmp = line.split()
                i = int(tmp[-3])
                t = float(tmp[-1])
                ts[i+1+offset] = t
            elif 't_contactL' in line:
                t = float(line.split()[-1])
                ts[L_lead] = t
            elif 't_contactR' in line:
                t = float(line.split()[-1])
                ts[L_lead+L_device] = t
    return ts'''

def get_hop_t (fname):
    L_lead = get_para (fname, 'L_lead', int)
    L_device = get_para (fname, 'L_device', int)
    L = 2*L_lead + L_device
    ts = np.full (L, np.nan)
    t_lead = get_para (fname, 't_lead', float)
    t_device = get_para (fname, 't_device', float)
    t_contactL = get_para (fname, 't_contactL', float)
    t_contactR = get_para (fname, 't_contactR', float)
    idevL, idevR = get_para (fname, 'device site', int, n=2)
    for ilink in range(L):
        if ilink == idevL-1:
            ts[ilink] = t_contactL
        elif ilink == idevR:
            ts[ilink] = t_contactR
        elif ilink >= idevL and ilink < idevR:
            ts[ilink] = t_device
        else:
            ts[ilink] = t_lead
    return ts

def get_data (fname):
    L_lead = get_para (fname, 'L_lead', int)
    L_device = get_para (fname, 'L_device', int)
    L = 2*L_lead + L_device
    Nstep = get_para (fname, 'step =', int, last=True)
    idevL, idevR = get_para (fname, 'device site', int, n=2)
    hopt = get_hop_t (fname)

    Is = np.full ((Nstep,L), np.nan)
    ns = np.full ((Nstep,L+1), np.nan)
    Ss = np.full ((Nstep,L+1), np.nan)
    dims = np.full ((Nstep,L), np.nan)
    with open(fname) as f:
        for line in f:
            line = line.lstrip()
            if line.startswith ('step ='):
                tmp = line.split()
                step = int(tmp[-1])
            elif line.startswith('current'):
                tmp = line.split()
                ilink = int(tmp[-4])
                I = float(tmp[-1])
                Is[step-1,ilink-1] = I * 2*pi * hopt[ilink]
            elif line.startswith('n ') and 'step' in globals():
                tmp = line.split()
                i = int(tmp[1])
                n = float(tmp[-1])
                ns[step-1,i-1] = n
            elif line.startswith('entang entropy'):
                tmp = line.split()
                i = int(tmp[2])
                S = float(tmp[-1])
                Ss[step-1,i-1] = S
            elif line.startswith('bond dim '):
                tmp = line.split()
                i = int(tmp[-3])
                m = float(tmp[-1])
                dims[step-1,i-1] = m
    return Is, ns, Ss, dims

def plot_prof (ax, data, dt, label=''):
    Nstep, L = np.shape(data)
    sc = ax.imshow (data, origin='lower', extent=[1, L, dt, Nstep*dt], aspect='auto')
    cb = pl.colorbar (sc)
    ax.set_xlabel ('site')
    ax.set_ylabel ('time')
    cb.ax.set_ylabel (label)

def plot_time_slice (ax, data, n, xs=[], **args):
    Nstep, L = np.shape(data)
    itv = Nstep // n
    if xs == []:
        xs = range(1,L+1)
    for d in data[::itv,:]:
        ax.plot (xs, d, marker='.', **args)
    ax.plot (xs, data[-1,:], marker='.', **args)

def get_basis_en (fname):
    ens = []
    with open(fname) as f:
        for line in f:
            if 'orbitals, segment, ki, energy' in line:
                for line in f:
                    tmp = line.split()
                    if not tmp[0].isdigit():
                        return ens
                    ens.append (float(tmp[-1]))

if __name__ == '__main__':
    fname = glob.glob ('out')[0]

    files = [i for i in sys.argv[1:] if i[0] != '-']
    fi,axi = pl.subplots()
    for fname in files:
        print (fname)

        en_basis = get_basis_en (fname)

        # Get data
        Is, ns, Ss, dims = get_data (fname)
        dt = get_para (fname, 'dt', float)
        m = get_para (fname, 'Largest link dim', int)
        Nstep, L = np.shape(Is)

        # I profile
        f,ax = pl.subplots()
        plot_prof (ax, Is, dt, 'current')
        ax.set_title ('$m='+str(m)+'$')
        ps.set(ax)

        '''f8,ax8 = pl.subplots()
        plot_time_slice (ax8, Is, n=3)
        ps.set(ax8)'''

        # n profile
        '''f2,ax2 = pl.subplots()
        plot_prof (ax2, ns, dt, 'density')
        ax2.set_title ('$m='+str(m)+'$')
        ps.set(ax2)

        f,ax = pl.subplots()
        plot_time_slice (ax, ns, n=5, xs=en_basis, ls='None')
        #ax.plot (range(1,len(en_basis)+1), en_basis)
        ax.set_xlabel ('energy')
        ax.set_ylabel ('occupasion')
        ps.set(ax)'''

        # S profile
        f5,ax5 = pl.subplots()
        plot_prof (ax5, Ss, dt, 'entropy')
        ax5.set_title ('$m='+str(m)+'$')
        ps.set(ax5)

        f6,ax6 = pl.subplots()
        plot_time_slice (ax6, Ss, n=3)
        ax6.set_xlabel ('site')
        ax6.set_ylabel ('entropy')
        ps.set(ax6)

        f,ax = pl.subplots()
        plot_time_slice (ax, dims, n=5)
        ax.set_xlabel ('site')
        ax.set_ylabel ('bond dimension')
        ps.set(ax)

        ts = dt * np.arange(1,Nstep+1)
        idevL, idevR = get_para (fname, 'device site', int, n=2)
        print (idevL, idevR)
        #
        '''N_dev = np.sum (ns[:, idevL-1:idevR], axis=1)
        f3,ax3 = pl.subplots()
        ax3.plot (ts, N_dev, marker='.')
        ps.set(ax3)'''

        # current vs time
        muL = get_para (fname, 'mu_biasL', float)
        muR = get_para (fname, 'mu_biasR', float)
        Vb = muL - muR
        Il = Is[:, idevL-2] / Vb
        Ir = Is[:, idevR-1] / Vb
        axi.plot (ts, Il, marker='.', label='left')
        axi.plot (ts, Ir, marker='.', label='right')
        #axi.plot (ts, Is[:, idevL-3] / Vb, marker='.', label='l2')
        #axi.plot (ts, Is[:, idevL] / Vb, marker='.', label='r2')
        axi.set_xlabel ('time')
        axi.set_ylabel ('conductance')
        axi.legend()

        V_lead = get_para (fname, 'V_lead', float)
        tp = 0.8#get_para (fname, 't_device', float)
        axi.axhline (exactG(V_lead, tp), ls='--', c='gray')

        ps.set(axi)

    pl.show()
