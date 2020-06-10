import numpy as np
import matplotlib
matplotlib.use('Agg')
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import csv
import scipy
from strfuncs import *
from collections import defaultdict
from scipy.interpolate import make_interp_spline, BSpline

def load_DF(filename):
    df = pd.read_csv(filename, compression = 'zip')
    return df

def write_to_DF(inputtype, df, d, En, index, model = 'Nilsson'):
    if inputtype in df.columns:
        print(inputtype)
    else:
        df.insert(-1, inputtype, index)
        print(df)
    sys.exit()
    return

def parse_WoodsSaxon(resfile):
    df = pd.DataFrame()
    resultfile = open(resfile).read()
    results = resultfile.replace('*','').split('\n1\n')[1:]
    defs = []
    ZN = ''
    states = []
    E = defaultdict(list)
    parity = {}
    omega = {}
    for r in results:
        ZN = r[find_nth(r, "NUCLEUS", 0):find_nth(r, "A=", 0)]
        Z = ZN[find_nth(ZN, "Z=",0)+2:find_nth(ZN, "N=", 0)]
        N = ZN[find_nth(ZN, "N= ", 0) + 2:]
        params = r[find_nth(r, "CENTRAL", 0):find_nth(r, ", D E", 0)]
        defo = r[find_nth(r, "BET ", 0):find_nth(r, "BETA3", 0)].replace("BET",'').replace(" ", '').replace(",",'')
        defs.append(float(defo))
        aux = r[find_nth(r, "AUXILIARY PARAMETERS", 0):find_nth(r, "SINGLE PARTICLE NUCLEAR", 0)]
        spe = [glfr(x.strip()) for x in r[find_nth(r, "1)", 0):find_nth(r, "MAXIMA", 0)].split("      ") if len(x) > 2][:-1]
        for i, e , j, p in spe:
            ind = int(i[:find_nth(i, ")",0)])
            E[ind].append(float(e))
            if "E" in p:
                parity[ind] = p[:find_nth(p, "E", 0)]
            else:
                parity[ind] = p
            omega[ind] = j
    df.insert(0, "Omega", omega.values())
    df.insert(1, "Parity", parity.values())
    df.index = omega
    Ed = [defaultdict(list) for d in range(len(defs))]
    for n in df.index:
        for d in range(len(E[n])):
            Ed[d][n] = E[n][d]
    for d in range(len(defs)):
        df.insert(len(df.columns), defs[d], Ed[d].values())
    df.insert(len(df.columns), 'ZN', str(Z) + str(N))
    df = df.sort_values(by = 0.35)
    return df

def spaghetti(df):
    oddplot = True
    defcols = np.array([x for x in df.columns[2:-1].values])
    defs = np.array([float(x) for x in df.columns[2:-1].values])
    dnew = np.linspace(defs.min(), defs.max(), 300)
    fig, ax = plt.subplots(figsize = (10,10))
    ax.set_title('Z, N = ' + df['ZN'][0], fontsize = 24)
    ax.set_xlabel(r'$\beta_2$', fontsize = 24) 
    ax.set_ylabel('MeV', fontsize = 24) 
    for i in df.index:
        if '+' in df.loc[[i]]['Parity'].values:
            linestyle = 'solid'
        else:
            linestyle = 'dotted'
            #linestyle = (0, (3, 1, 1, 1, 1, 1))
        spl = make_interp_spline(defs, df.loc[[i]][defcols].values[0], k = 3)
        Esmooth = spl(dnew)
        ax.plot(dnew, Esmooth, c = 'k', linestyle = linestyle, lw = 1)
        if oddplot:
            plt.text(dnew[-1]*1.005, Esmooth[-1],"[" + df.loc[[i]]['Omega'].values[0] + df.loc[[i]]['Parity'].values[0]+ "]", fontsize = 2)
            oddplot = False
        else:
            plt.text(dnew[-1]*1.025, Esmooth[-1], "[" + df.loc[[i]]['Omega'].values[0] + df.loc[[i]]['Parity'].values[0] + "]", fontsize = 2)
            oddplot = True
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.savefig('nilssonplot.png', format = 'png')
    plt.savefig('nilssonplot.eps', format = 'eps')
    sys.exit()

df = load_DF('WoodsSaxonDF.zip')
spaghetti(df)

sys.exit()

WSdf = parse_WoodsSaxon('swbeta.out')
print(WSdf)
WSdf.to_csv('WoodsSaxonDF.zip', index = False, compression = 'zip')

