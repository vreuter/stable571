'''
Created on Dec 17, 2019


@author: kyang

These functions will allow us to proces the majiq deltapsi output as we would like.
'''
import math
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
sns.set()

def calc_delta_deltapsi(dir="/Users/kyang/bonini/stable571/majiq_analysis/results/",
                     comparisons_str="Ime4 Mcherry"):
    #dir is the directory where the majiq deltapsi.tsv data output from majiq tsv view is
    #comparisons_str is the space-delimited list of conditions for which delta-deltapsi of
    #HS vs control should be calculated for
    for c in comparisons_str.split(" "):
        LSVdict1 = dict()
        with open(dir+c+"_control_Input_"+c+"_control_m6A.deltapsi.tsv") as inF:
            for i, line in enumerate(inF):
                if i == 0:
                    continue
                else:
                    linelist = line.split("\t")
                    LSVdict1[linelist[1]] = linelist[3].split(";")
        LSVdict2 = dict()
        with open(dir+c + "_HS_Input_" + c + "_HS_m6A.deltapsi.tsv") as inF:
            for i, line in enumerate(inF):
                if i == 0:
                    continue
                else:
                    linelist = line.split("\t")
                    LSVdict2[linelist[1]] = linelist[3].split(";")
        intersect_keys = [x for x in LSVdict1.keys() if x in LSVdict2.keys()]
        print("There are "+str(len(intersect_keys))+" LSVs in commmon")
        #reorder from great
        deltadict = dict()
        for k in intersect_keys:
            #here this represents HS minus control
            deltadict[k] = [float(y) - float(x) for x,y in zip(LSVdict1[k],LSVdict2[k])]
        #output two files, one sorting delta deltapsi by min values - output most negative values at top
        #and another sorting delta deltapsi by max values - output most positive values at top
        deltalist = [[k,v] for k,v in deltadict.items()]
        sorted_min = sorted(deltalist,key=lambda x:min(x[1]))
        deltamin = [str(y[0])+"\t"+";".join([str(z) for z in y[1]])+"\n"
                    for y in sorted_min]
        sorted_max = sorted(deltalist,key=lambda x:max(x[1]),reverse=True)
        deltamax = [str(y[0])+"\t"+";".join([str(z) for z in y[1]])+"\n"
                    for y in sorted_max]
        with open(dir+c+"_deltamin.tsv",'w+') as inF:
            inF.write("LSV name\tDelta Deltapsi (HS - Control)\n")
            inF.write("".join(deltamin))
        with open(dir+c+"_deltamax.tsv",'w+') as inF:
            inF.write("LSV name\tDelta Deltapsi (HS - Control)\n")
            inF.write("".join(deltamax))
def normalize(dir="/Users/kyang/bonini/stable571/majiq_analysis/results/"):
    #dir is the directory where the the calc_delta_deltapsi results are
    #note this function uses deltamin but there's no reason deltamax can't be used instead
    exp_dict = dict()
    with open(dir+"Ime4"+"_deltamin.tsv") as inF:
        for i,line in enumerate(inF):
            if i == 0:
                continue
            if "\n" == line[-1]:
                line = line[:-1]
            linelist = line.split("\t")
            exp_dict[linelist[0]] = linelist[1].split(";")
    ctrl_dict = dict()
    with open(dir + "Mcherry" + "_deltamin.tsv") as inF:
        for i, line in enumerate(inF):
            if i == 0:
                continue
            if "\n" == line[-1]:
                line = line[:-1]
            linelist = line.split("\t")
            ctrl_dict[linelist[0]] = linelist[1].split(";")
    intersect_keys_temp = [x for x in exp_dict if x in ctrl_dict]
    print("There are " + str(len(intersect_keys_temp)) + " LSVs in commmon")
    exp_final = []
    ctrl_final = []
    temp_list = []
    # Remove redundant keys since some of the RefSeq genes are identical, just under a different NM number
    intersect_keys_id = []
    intersect_keys = []
    for k in intersect_keys_temp:
        if k.split(":")[-1] in intersect_keys_id:
            continue
        intersect_keys_id.append(k.split(":")[-1])
        intersect_keys.append(k)
    print("After removing redundant LSVs, there are "+str(len(intersect_keys))+" LSVs left over")
    for k in intersect_keys:
        normed_list = [float(x)-float(y) for x,y in zip(exp_dict[k],ctrl_dict[k])]
        temp_list.append([k,exp_dict[k],ctrl_dict[k],normed_list])

    #Now, sort the list by the max absolute value of normed list
    temp_list.sort(key=lambda x:max([abs(z) for z in x[3]]),reverse=True)
    out_list = ["\t".join([x[0].split(":")[0],x[0],";".join(x[1]),";".join(x[2]),";".join([str(z) for z in x[3]])])for x in temp_list]
    with open(dir+"Ime4_normalized.tsv",'w+') as inF:
        inF.write("Gene ID\tLSV ID\tIme4 delta deltapsi\tMcherry delta deltapsi\tIme4 delta deltapsi - Mcherry deltadeltapsi\n")
        inF.write("\n".join(out_list))
    #Make 2 plots: 1 is the scatterplot of max delta_deltapsi junction per LSV in Ime4 and its corresponding junction in
    #Mcherry; 2 is the histogram of fold changes
    x_final,y_final,z_final = [],[],[]
    for x in temp_list:
        x[1] = [float(z) for z in x[1]]
        x[2] = [float(z) for z in x[2]]
        index_max = np.argmax([abs(z) for z in x[1]])
        x_final.append(x[1][index_max])
        y_final.append(x[2][index_max])
        index_max = np.argmax([abs(z) for z in x[3]])
        z_final.append(x[3][index_max])
    fig,ax1= plt.subplots()
    scatterpd = pd.DataFrame(list(zip(x_final,y_final)),columns=["Ime4 delta deltapsi","Mcherry delta deltapsi"])
    sns.scatterplot(x="Ime4 delta deltapsi",y="Mcherry delta deltapsi",data=scatterpd,ax=ax1)
    plt.savefig(dir+"Fig_1.png",dpi=200)
    plt.show()
    fig, ax2 = plt.subplots()
    sns.distplot(z_final,ax=ax2)
    ax2.set(xlabel="Ime4 delta deltapsi - Mcherry delta deltapsi")
    plt.savefig(dir + "Fig_2.png", dpi=200)
    plt.show()