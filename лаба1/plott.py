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
k,l -- костыли для распределений максвелла по проекциям скоросей
'''

def picture(file, n, Ox, Oy, Caption, name, a = 0, k = 0, l = 0):
    with open(file + '.csv'.format(1), 'r') as f:
        reader = csv.reader(f)
        rdr = [list(read) for read in reader][3::]
        
    fig, ax = plt.subplots()
    if a == 1:
        ax.plot( [np.log1p((float(r[1]))/(((float(r[0])**2)))) for r in rdr],[(float(r[0])**2) for r in rdr], '.')
    elif a == 2:
        for i in range(n):
            ax.plot([np.log1p(float(r[i+2])) for r in rdr], [(float(r[i + 5])**2) for r in rdr], 'r.-', color='#'+ str(int((i - 1)*(i - 2)*4.5))+ str(0) + str(i*(2 - i)*8) + str(0) + str(i*(i - 1)*4) + str(0))
    
    else:
        i = 0
        ax.plot(  [(8*float(r[0])) for r in rdr], [((float(r[1]) - 30000)/10000000 + 0.002) for r in rdr], 'r.-', color='#'+ str(int((i - 1)*(i - 2)*4.5))+ str(0) + str(i*(2 - i)*8) + str(0) + str(i*(i - 1)*4) + str(0))
        i = 1
        ax.plot([float(r[0]) for r in rdr], [((float(r[2]) )) for r in rdr], 'r.-', color='#'+ str(int((i - 1)*(i - 2)*4.5))+ str(0) + str(i*(2 - i)*8) + str(0) + str(i*(i - 1)*4) + str(0))
            
    ax.grid(linestyle='--', linewidth=0.5)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel(Oy)
    ax.set_xlabel(Ox)
    ax.set_title(Caption)
            
    plt.show()
    fig.savefig(name)



picture('Residuals',2, r'T', r'Невязки', 'Рассчет динамической пямяти', f"Residualss.png")
'''
picture('Velocities', 3, r'$V_i$', r'$N$', 'Распределение Максвелла по проекциям скоростей', 'Maxwell_pr.png', 0, 1, 5)
picture('Velocities', 3, r'$V_i^2$', r'$\ln{N}$', 'Линеаризованное распределение Максвелла по проекциям скоростей', 'Maxwell_pr_lin.png', 2, 5)
'''

#picture('Fluctuation', 1, r'$M$', r'$<E^2> - <E>^2$', 'Флуктуации энергии', 'Fluctuation.png')
