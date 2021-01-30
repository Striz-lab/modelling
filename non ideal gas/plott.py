import csv
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)



def Diffusion():
    with open('Diffuziya.csv'.format(1), 'r') as f:
        reader = csv.reader(f)
        rdr = [list(read) for read in reader][3::]

    fig, ax = plt.subplots()
    ax.plot([float(r[0]) for r in rdr], [float(r[1])**0.5 for r in rdr], 'r.-')
    ax.grid(linestyle='--', linewidth=0.5)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    ax.set_ylabel(r'$<\bigtriangleup r>$')
    ax.set_xlabel(r'$t$')
    ax.set_title('Диффузия')
        
    plt.show()
    fig.savefig('Diffusion.png')
    
def Maxwell():
    with open('Maxwell.csv'.format(1), 'r') as f:
        reader = csv.reader(f)
        rdr = [list(read) for read in reader][3::]

    fig, ax = plt.subplots()
    ax.plot([float(r[0]) for r in rdr], [float(r[1]) for r in rdr], '.')
    ax.grid(linestyle='--', linewidth=0.5)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

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
    ax.plot([float(r[0]) for r in rdr], [0.5*float(r[2]) for r in rdr], '.', color='tab:red')
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

