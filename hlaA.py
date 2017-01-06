from collections import defaultdict
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

input_file = 'A_nuc.txt'

def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt

seq = import_file(input_file)

def Set_Start(seq):
    Set_Start = []
    for i in range(len(seq)):
        if seq[i].startswith(' cDNA'):
            Set_Start.append(i)
    return Set_Start

#print set_start 
#-> [5, 3652, 7299, 10946, 14593, 18240, 21887, 25534, 29181, 32828, 36475, 40122, 43769, 47416, 51063, 54710]

#print len(set_start)
#-> total number of set = 16

def Set_Dict(seq):
    set_dict = {}
    set_start = Set_Start(seq)
    for n in range(len(set_start)):
        set_name = 's'+str(n+1)
        set_dict[set_name] = []
        if n in range(len(set_start)-1):
            for m in range(set_start[n], set_start[n+1]):
                if seq[m].startswith(' A*'):
                    set_dict[set_name].append(seq[m])
        else:
            for m in range(set_start[n], len(seq)):
                if seq[m].startswith(' A*'):
                    set_dict[set_name].append(seq[m])
    return set_dict

def Group_List(seq):
    group_list = []
    set_dict = Set_Dict(seq)
    for line in set_dict['s1']:
        group_name = line[1:5]
        if group_name not in group_list:
            group_list.append(group_name)
    return group_list

def Group_Dict(seq):
    group_dict = {}
    set_dict = Set_Dict(seq)
    for line in set_dict['s1']:
        group_name = line[1:5]
        if group_name not in group_dict.keys():
            group_dict[group_name] = 0
        if group_name in group_dict.keys():
            group_dict[group_name] += 1
    return group_dict

#print group_dict 
#-> {'A*24': 482, 'A*69': 4, 'A*29': 122, 'A*30': 134, 'A*34': 21, 'A*26': 176, 'A*25': 50, 'A*11': 335, 'A*23': 98, 'A*31': 141, 'A*32': 120, 'A*33': 157, 'A*80': 4, 'A*01': 297, 'A*66': 30, 'A*03': 347, 'A*02': 867, 'A*74': 33, 'A*43': 1, 'A*68': 220, 'A*36': 5}

#print len(group_dict.keys()) 
#-> total number of groups: 21 

#print sum(group_dict.values()) 
#-> total number of alleles: 3644

def Total_Allele_Dict(seq):
    set_dict = Set_Dict(seq)
    total_allele_dict = {}
    for s in sorted(set_dict.keys(), key=lambda x: int(x.replace('s', ''))):
        for allele in set_dict[s]:
            s_allele = allele.split()
            allele_name = s_allele[0]
            allele_seq = ' '.join(s_allele[1:len(s_allele)])
            if allele_name not in total_allele_dict.keys():
                total_allele_dict[allele_name] = ''
            if allele_name in total_allele_dict.keys():
                total_allele_dict[allele_name] += allele_seq
    return total_allele_dict

def Allele_Dict(seq):
    total_allele_dict = Total_Allele_Dict(seq)
    for allele in total_allele_dict.keys():
        total_allele_dict[allele] = total_allele_dict[allele].replace(' ', '')
    return total_allele_dict

# lenth of sequence without space: 1300

def Ref_Allele(seq):
    allele_dict = Allele_Dict(seq)
    for allele in allele_dict.keys():
        seq = allele_dict[allele]
        if '-' not in seq:
            ref_allele = allele
            ref_seq = allele_dict[allele]
    return ref_allele
    
def Score_Dict(seq):
    allele_dict = Allele_Dict(seq)
    ref_allele = Ref_Allele(seq)
    ref_seq = allele_dict[ref_allele]
    score_dict = {}
    for i in range(len(ref_seq)):
        position = str(i+1)
        score_dict[position] = 0
        for seq in allele_dict.values():
            if seq[i] not in '*-.':
                score_dict[position] += 1
    return score_dict

def Score_Graph(seq):
    score_dict = Score_Dict(seq)
    position = sorted(score_dict.keys(), key=lambda x: int(x))
    score = []
    for p in position:
        score.append(score_dict[p])
    print position
    print score
    plt.plot(position, score)
    return plt.savefig('hlaA.png')

print Score_Graph(seq)
