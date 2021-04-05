#from array import *
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

# Python program to read
# file word by word
   
# opening the text file

global iterator, comutator

with open('1.5.txt','r') as file:

    iterator = int(0)
    comutator = 0
    r = []
    v = []
    
    
    
    r.append(0)
    v.append(0)
    r0 = []
    v0 = []
    r0.append(0)
    v0.append(0)
    D = 0
    #4009000
    
    #Создание массива начальных скоростей и ускорений
    for i in range(4009):
        line = file.readline()
        if ((line != 'ITEM: TIMESTEP\n') and (line != '0\n') and (line != 'ITEM: NUMBER OF ATOMS\n') and (line != '4000\n')           and (line != 'ITEM: BOX BOUNDS pp pp pp\n') and (line != '0.0000000000000000e+00 1.7878070701931353e+01\n') and (line != 'ITEM: ATOMS id type xu yu zu vx vy vz\n')):
            # reading each word
           
            # displaying the words
            word = line.split()

            r0.append((float(word[2])**2+float(word[3])**2+float(word[4])**2)**0.5)
            v0.append((float(word[5])**2+float(word[6])**2+float(word[7])**2)**0.5)
            #print(i)
            
    for j in range(1000):
        r.append(0)
        v.append(0)
        #print(" ")
        for i in range(4009):
            line = file.readline()
            if ((line != 'ITEM: TIMESTEP\n') and (line != '0\n') and (line != 'ITEM: NUMBER OF ATOMS\n')  and (line!='4000\n')         and (line != 'ITEM: BOX BOUNDS pp pp pp\n') and (line != '0.0000000000000000e+00 1.7878070701931353e+01\n') and (line != 'ITEM: ATOMS id type xu yu zu vx vy vz\n') and (line != str(j+1) + '0\n')):
                # reading each word
                    
                
                # displaying the words
                word = line.split()
                #print(j)
                #print(word)
                #print(r0[i-9])
                r[j] += (((float(word[2])**2+float(word[3])**2+float(word[4])**2)**0.5) - r0[i - 9])**2
                v[j] += (((float(word[5])**2+float(word[6])**2+float(word[7])**2)**0.5)*v0[i - 9])
                    
        r[j] = 1.3*(2*(r[j]/4000 - 53.5)+21.4-18)
        v[j] = v[j]/4000
        D += v[j]

  
#print(D/3)
v.pop()
r.pop()
#D=0
t = []
#r = []


for i in range(1000):
    t.append(i-700+400-300)
    #r.append(np.exp(-20*i*0.001)*np.cos(3*i*0.001))
    #D += r[i]
#print(D/3)

'''
fig, ax = plt.subplots()
ax.plot(t,r, '.')

        
ax.grid(linestyle='--', linewidth=0.5)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$<\Delta R^2>$')
ax.set_xlabel('t')
ax.set_title('Эйнштейн-Смолуховский')
        
plt.show()
fig.savefig('name')
'''
fig, ax = plt.subplots()
ax.plot(t,r, '.', color = "green")

        
ax.grid(linestyle='--', linewidth=0.5)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$<\Delta R^2>$')
ax.set_xlabel("t")
ax.set_title("Диффузия")
        
plt.show()
fig.savefig("name")


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




#picture('Residuals',2, r'T', r'Невязки', 'Рассчет динамической пямяти', f"Residualss.png")
