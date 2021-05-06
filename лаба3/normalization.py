#from array import *
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


t = []
r0 = []
r1 = []
i = 0
j = 0
l = 0
flag = 0
E = []
T = []


with open('log.txt','r') as f:
    line = f.readline()
    
    for line in f:
        
            #print(float(word[1]))
        word = line.split()
        E.append(float(word[1]))
        T.append(float(word[3]))
        

    

rdf = []
r = []




fig, ax = plt.subplots()
ax.plot(E,T, '.', color = "green")

        
ax.grid(linestyle='--', linewidth=0.5)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel('E')
ax.set_xlabel('T')
ax.set_title('Зависимость Энергии от Температуры')
plt.show()

fig.savefig('5.png')

