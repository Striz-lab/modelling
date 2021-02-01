import csv
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)



def Diffusion():
    with open('Diffuziya.csv'.format(1), 'r') as f:
        reader = csv.reader(f)
        rdr = [list(read) for read in reader][3::]

    fig, ax = plt.subplots()
    ax.plot([float(r[0]) for r in rdr], [float(r[1]) for r in rdr], 'r.-')
    ax.grid(linestyle='--', linewidth=0.5)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    ax.set_ylabel(r'$<\bigtriangleup r>$')
    ax.set_xlabel(r'$t$')
    ax.set_title('Диффузия')
        
    plt.show()
    fig.savefig('Diffusion.png')
    
def Maxwell():
    mu = 100  # mean of distribution
    sigma = 15  # standard deviation of distribution
    with open('Maxwell.csv'.format(1), 'r') as f:
        reader = csv.reader(f)
        rdr = [list(read) for read in reader][3::]

    fig, ax = plt.subplots(1, 1, figsize=(8, 4), sharey=True)
    ax.bar([float(r[0]) for r in rdr], [float(r[1]) for r in rdr])
    #ax.plot( [(float(r[0])**2) for r in rdr],[np.log1p((float(r[1])/(float(r[0])**2))) for r in rdr], '.')
    
    #y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    #ax.plot(bins, y, '--')
    #ax[1].grid(linestyle='--', linewidth=0.5)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylabel(r'$N$')
    ax.set_xlabel(r'$V$')
    ax.set_title('Распределение Максвелла по модулям скоростей')
        
    plt.show()
    fig.savefig('Maxwell.png')
    
    
def Energy():
    with open('Energy.csv'.format(1), 'r') as f:
        reader = csv.reader(f)
        rdr = [list(read) for read in reader][3::]

    fig, ax = plt.subplots()
    ax.plot([float(r[0]) for r in rdr], [float(r[1]) for r in rdr], '.', color='tab:blue')
    ax.plot([float(r[0]) for r in rdr], [float(r[2]) for r in rdr], '.', color='tab:red')
    ax.plot([float(r[0]) for r in rdr], [float(r[3]) for r in rdr], '.', color='tab:purple')
    ax.grid(linestyle='--', linewidth=0.5)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylabel(r'$N$')
    ax.set_xlabel(r'$V$')
    ax.set_title('Энергии')
        
    plt.show()
    fig.savefig('Energy.png')

Diffusion()
Maxwell()
Energy()

