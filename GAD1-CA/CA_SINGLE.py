#CA.py
from math import *
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import matplotlib.pylab as pl
import numpy as np
#import scipy.stats as stats
#from scipy.interpolate import make_interp_spline

def main():
    
    amino = int(input('Number of amino acids: '))
    #temperature = input('Please input the temperature (NUMBERS): ')
    temperature='CA-T240K.txt'
    CA = 0
    count = 0

    COM_X = 0
    COM_Y = 0
    COM_Z = 0
    v_com = [COM_X, COM_Y, COM_Z]
    
    data = []
    distance_skip = []
    with open(temperature) as f:
        for line in f:
            aline = line.strip().split()
            vector = aline[-3:]
            data.append(vector)
		
    CA = len(data)
    frame = int(CA/amino)
    file = open("r_CA.dat","w")
    file.write('FRAME  DISTANCE\n')
    for j in range(amino):
        file.write("Amino acid No%d\n" %(j+1))
        r_skip = []
        for i in range(frame):
            radius = calDistance(data[i*amino+j], v_com)
            file.write("%5d %8.5f\n" % ((i+1), radius))
            if i > 1/3*frame:
                count += 1
                r_skip.append(radius)
            #print(data[i*amino+j])
        
        distance_skip.append(r_skip)
    file.close()


    x = [0]*amino
    y = [0]*amino
    
    for i in range(amino):
        
        end_int = ceil(max(distance_skip[i])/0.05)
        begin_int = floor(min(distance_skip[i])/0.05)
        nbins= end_int-begin_int

        colors = mcd.XKCD_COLORS["xkcd:" + overlap[i]].upper()
        marker='o'
        if i == 0:
            marker='^'
        elif i == amino-2:
            colors = mcd.XKCD_COLORS["xkcd:" + overlap[i+3]].upper()
        elif i == amino-1:
            colors = mcd.XKCD_COLORS["xkcd:" + overlap[i+5]].upper()
            marker='D'

        #density = stats.gaussian_kde(distance[i])

        fig, ax = plt.subplots(tight_layout=True)
        y,x,_ = ax.hist(distance_skip[i],bins=nbins,range=[begin_int*0.05,end_int*0.05],ec='black',density=True)
        plt.title('Location of ${Cα}$ for Amino Acid No%d\n\n'%(i+1),fontweight ="bold")
        plt.legend(ncol=2,frameon=True, fontsize=5, numpoints=1)
        
        plt.xlabel('r$_{Cα}$ (nm)')
        plt.ylabel('P(r$_{Cα}$)')
        plt.savefig("hist%d.pdf" %(i+1), format="pdf")
        plt.close()

        #Writing Radius Distribution Function g(r)
        yy = y/(4*np.pi*(x[0:-1]+0.025)**2)
        print(x,y,type(y))
        if first3(GAD1[i]) in hydrophobic:
                yy = -yy

        plt.subplot(2, 1, 1)
        plt.xlim(0.8,3.3)
        plt.plot(x[0:-1]+0.025,yy,linewidth=0.3,marker=marker,markersize=1,color=colors)
    plt.axhline(y=0, c="k", linestyle="dashed", linewidth=1, zorder=0)
    #plt.legend(ncol=2,frameon=True, fontsize=5, numpoints=1)
    #plt.title('T240 to T260K',fontweight ="bold")
    plt.title('T240K  ',fontsize=8, loc='right', y=1.0, pad=-10)
    #plt.xlabel('r$_{Cα}$ (nm)')
    plt.ylabel('P(r$_{Cα}$)')
        #plt.plot(x,density(x))

#########################################################################
    temperature1='CA-T260K.txt'
    CA = 0
    count = 0

    COM_X = 0
    COM_Y = 0
    COM_Z = 0
    v_com = [COM_X, COM_Y, COM_Z]
    
    data = []
    distance_skip = []
    with open(temperature1) as f:
        for line in f:
            aline = line.strip().split()
            vector = aline[-3:]
            data.append(vector)
		
    CA = len(data)
    frame = int(CA/amino)
    file = open("r_CA.dat","w")
    file.write('FRAME  DISTANCE\n')
    for j in range(amino):
        file.write("Amino acid No%d\n" %(j+1))
        r_skip = []
        for i in range(frame):
            radius = calDistance(data[i*amino+j], v_com)
            file.write("%5d %8.5f\n" % ((i+1), radius))
            if i > 1/3*frame:
                count += 1
                r_skip.append(radius)
            #print(data[i*amino+j])
        
        distance_skip.append(r_skip)
    file.close()


    x = [0]*amino
    y = [0]*amino
    
    for i in range(amino):
        
        end_int = ceil(max(distance_skip[i])/0.05)
        begin_int = floor(min(distance_skip[i])/0.05)
        nbins= end_int-begin_int

        colors = mcd.XKCD_COLORS["xkcd:" + overlap[i]].upper()
        marker='o'
        if i == 0:
            marker='^'
        elif i == amino-2:
            colors = mcd.XKCD_COLORS["xkcd:" + overlap[i+3]].upper()
        elif i == amino-1:
            colors = mcd.XKCD_COLORS["xkcd:" + overlap[i+5]].upper()
            marker='D'

        #density = stats.gaussian_kde(distance[i])
        fig, ax = plt.subplots(tight_layout=True)
        y,x,_ = ax.hist(distance_skip[i],bins=nbins,range=[begin_int*0.05,end_int*0.05],ec='black',density=True)
        plt.title('Location of ${Cα}$ for Amino Acid No%d\n\n'%(i+1),fontweight ="bold")
        #plt.xlabel('r$_{Cα}$ (nm)')
        #plt.ylabel('P(r$_{Cα}$)')
        plt.legend(ncol=2,frameon=True, fontsize=5, numpoints=1)
        plt.savefig("hist%d.pdf" %(i+1), format="pdf")
        plt.close()

        print(x,y,type(y))
        #Writing Radius Distribution Function g(r)
        yy = y/(4*np.pi*(x[0:-1]+0.025)**2)
        if first3(GAD1[i]) in hydrophobic:
                yy = -yy

        plt.subplot(2, 1, 2)
        plt.xlim(0.8,3.3)
        plt.plot(x[0:-1]+0.025,yy,linewidth=0.3,marker=marker,markersize=1,color=colors,label='%s'%(GAD1[i]))
    plt.axhline(y=0, c="k", linestyle="dashed", linewidth=1, zorder=0)
    plt.legend(ncol=2,frameon=True, fontsize=5, numpoints=1)
    #plt.title('T260K',fontweight ="bold")
    plt.title('T260K  ',fontsize=8, loc='right', y=1.0, pad=-10)
    plt.xlabel('r$_{Cα}$ (nm)')
    plt.ylabel('P(r$_{Cα}$)')
        #plt.plot(x,density(x))

###############################################################################
    #plt.axhline(y=0, c="k", linestyle="dashed", linewidth=1, zorder=0) 
    #plt.xlabel('r$_{Cα}$ (nm)')
    #plt.ylabel('P(r$_{Cα}$)')
    #plt.xlim(0.7,3.3)
    #plt.legend(ncol=2,frameon=True, fontsize=5, numpoints=1)
    #plt.title('Location of ${Cα}$ at %s K\n'%(temperature),fontweight ="bold")
    plt.savefig("%s.pdf" %temperature[:-3], format="pdf")
    plt.show()
    plt.close()
        
def calDistance(v_amino, v_com):
    return sqrt((float(v_amino[0])-v_com[0])**2 + \
                       (float(v_amino[1])-v_com[1])**2 + (float(v_amino[2])-v_com[2])**2)

GAD1= ['Phe (1)','Ile (2)','His (3)','His (4)','Ile (5)','Ile (6)','Gly (7)','Trp (8)','Ile (9)',\
                      'Ser (10)','His (11)','Gly (12)','Val (13)','Arg (14)','Ala (15)','Ile (16)','His (17)',\
                      'Arg (18)','Ala (19)','Ile (20)','His (21)']
hydrophobic=['Ala','Val','Ile', 'Leu','Met', 'Phe', 'Tyr', 'Trp']
overlap = [name for name in mcd.CSS4_COLORS
           if "xkcd:" + name in mcd.XKCD_COLORS]
def first3(string):
    return string[:3]


main()
