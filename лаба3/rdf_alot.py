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
R = []
RDF = []


with open('rdf.txt','r') as f:
    line = f.readline()
    
    for line in f:
        if (len(line) > 15):
            #print(float(word[1]))
            word = line.split()
            R.append(float(word[1]))
            RDF.append(float(word[2]))
        

    

rdf = []
r = []

#print(l)

for i in range(200):
    for j in range(500):
        r.append(R[500*i+j])
        rdf.append(RDF[500*i+j])

    fig, ax = plt.subplots()
    ax.plot(r,rdf, '.', color = "green")

            
    ax.grid(linestyle='--', linewidth=0.5)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylabel('RDF')
    ax.set_xlabel('r')
    ax.set_title('Радиальная функция распределения')
            
    #plt.show()
    
    fig.savefig(str(i) +'.png')
    
    
    for j in range(500):
        r.pop()
        rdf.pop()


