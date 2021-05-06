import math
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

x = []
y = []
z = []
vx = []
vy = []
vz = []


with open('b.txt','r') as f:
    line = f.readline()
    for line in f:
        if ((line!= "ITEM: TIMESTEP\n")
        and (line!= "ITEM: NUMBER OF ATOMS\n")
        and (line!= "ITEM: BOX BOUNDS pp pp pp\n")
        and (line!= "0.0000000000000000e+00 1.5874010519681994e+01\n")
        and (line!= "ITEM: ATOMS id type xu yu zu\n")
        and (len(line) > 15)):
            
            word = line.split()
            
            x.append(float(word[2]))
            y.append(float(word[3]))
            z.append(float(word[4]))
            
r  = []
 
for i in range(4000):
    r.append(((x[4000*998 + i] - x[i])**2 + (y[4000*998 + i] - y[i])**2 + (z[4000*998 + i] - z[i])**2)**0.5)
        
r.sort()

n = - r[0] + r[3999]

N = []
k = []
now = 0

for i in range(1000):
    N.append(0)
    for j in range(4000):
        if (r[j] < now + n/20) and (r[j] > now - n/20):
            N[i] = N[i] + 1
    k.append(now)
    now = now + n/1000
        

print(n)

'''

fig, ax = plt.subplots()
ax.plot(k,N, '.')
        
ax.grid(linestyle='--', linewidth=0.5)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel('N')
ax.set_xlabel(r"$\sqrt{\Delta r^2}$")
ax.set_title('Распределение отклонений частиц')
        
plt.show()
fig.savefig('name')

'''

for i in range(999):
    N[i+1] = np.log1p(N[i+1]/(k[i+1]**2))
    k[i+1] = k[i+1]**2

N.pop()
k.pop()

fig, ax = plt.subplots()
ax.plot(k,N, '.')
        
ax.grid(linestyle='--', linewidth=0.5)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r"$\ln(\frac{N}{\Delta r^2})$")
ax.set_xlabel(r"$\Delta r^2$")
ax.set_title('Лин. распределение отклонений частиц')
        
plt.show()
fig.savefig('nam')
