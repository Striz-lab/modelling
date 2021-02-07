import csv
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
'''
Рисует картинку
file -- название файла без расширения
n -- количество кривых на рисунке
Ox -- подпись к оси Ох
Oy -- подпись к оси Оу
Caption -- подпись к графику
name -- название сохраненной картинки с расширением
a -- костыль, позволяющий построить график в нелинейных осях
'''

def picture(file, n, Ox, Oy, Caption, name, a = 0):
    with open(file + '.csv'.format(1), 'r') as f:
        reader = csv.reader(f)
        rdr = [list(read) for read in reader][3::]
        
    fig, ax = plt.subplots()
    if a == 1:
        ax.plot( [np.log1p((float(r[1]))/(((float(r[0])**2)))) for r in rdr],[(float(r[0])**2) for r in rdr], '.')
    else:
        for i in range(n):
            ax.plot([float(r[0]) for r in rdr], [float(r[i + 1]) for r in rdr], 'r.-', color='#'+ str(int((i - 1)*(i - 2)*4.5))+ str(0) + str(i*(2 - i)*8) + str(0) + str(i*(i - 1)*4) + str(0))
            
    ax.grid(linestyle='--', linewidth=0.5)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel(Oy)
    ax.set_xlabel(Ox)
    ax.set_title(Caption)
            
    plt.show()
    fig.savefig(name)


picture('Diffusion', 1, r'$t$', r'$<\Delta r^2>$', 'Диффузия', 'Diffusion.png')
picture('Velocities', 1, r'$V$', r'$N$', 'Распределение Максвелла по модулям скоростей', 'Maxwell.png')
picture('Velocities', 1, r'$V^2$', r'$\ln{\frac{N}{V^2}}$', 'Линеаризованное распределение Максвелла по модулям скоростей', 'Maxwell_lin.png', 1)
picture('Energy', 3, r'$t$', r'$E$', 'Зависимость энергии от времени', 'Energy.png')
picture('Temperature', 1, r'$t$', r'$T$', 'Зависимость температуры от времени', 'Temperature.png')

