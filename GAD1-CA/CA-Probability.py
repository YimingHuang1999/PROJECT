from math import *
import matplotlib.pyplot as plt
import numpy as np
#import scipy.stats as stats
#from scipy.interpolate import make_interp_spline

def main():
    
    #frame = 200
    #frame = int(input('Number of frames: '))
    #amino = 9
    amino = int(input('Number of amino acids: '))
    temperature = input('Please input the temperature (NUMBERS): ')
    CA = 0
    count = 0

    COM_X = 0
    COM_Y = 0
    COM_Z = 0
    v_com = [COM_X, COM_Y, COM_Z]
    
    data = []
    distance_skip = []
    with open('CA-T240K.txt') as f:
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

        #density = stats.gaussian_kde(distance[i])
        fig, ax = plt.subplots(tight_layout=True)
        y[i],x[i],_ = ax.hist(distance_skip[i],bins=nbins,range=[begin_int*0.05,end_int*0.05],ec='black',density=True)
        plt.title('Location of ${Cα}$ for Amino Acid No%d\n\n'%(i+1),fontweight ="bold")
        plt.xlabel('r$_{Cα}$ (nm)')
        plt.ylabel('P(r$_{Cα}$)')
        plt.savefig("hist_%d.png" %(i+1), format="png")
        plt.close()

        #print(x,y)
        plt.plot(x[i][0:-1]+0.025,y[i],'o-.',label='Amino Acid No%d'%(i+1))
        
        #plt.plot(x,density(x))
        #print(np.sum(y))
        #bin_centers = 0.5*(x[1:]+x[:-1])
        #X_Y_Spline = make_interp_spline(x, n)
        #X_ = np.linspace(x.min(), x.max(), 500)
        #Y_ = X_Y_Spline(X_)
        #plt.plot(X_,Y_)
    
    plt.xlabel('r$_{Cα}$ (nm)')
    plt.ylabel('P(r$_{Cα}$)')
    plt.ylim(0.0,5.0)
    plt.legend(frameon=True, fontsize=8.0, numpoints=1, loc='upper left')
    plt.title('Location of ${Cα}$ for %s K\n'%(temperature),fontweight ="bold")
    plt.savefig("T%sK_hist.eps" %temperature, format="eps")
    plt.show()
    plt.close()
        
def calDistance(v_amino, v_com):
    return sqrt((float(v_amino[0])-v_com[0])**2 + \
                       (float(v_amino[1])-v_com[1])**2 + (float(v_amino[2])-v_com[2])**2)

#sqrt((data[i*amino+j][3]-COM[i][2])**2 + (data[i*amino+j][4]-COM[i][3])**2 + \
        #(data[i*amino+j][5]-COM[i][4])**2)

aminos_names_list =  ['Phe (1)','Ile (2)','His (3)','His (4)','Ile (5)','Ile (6)','Gly (7)','Trp (8)','Ile (9)',\
                      'Sys (10)','His (11)','Gly (12)','Val (13)','Arg (14)','Ala (15)','Ile (16)','His (17)',\
                      'Arg (18)','Ala (19)','Ile (20)','His (21)'],
main()
