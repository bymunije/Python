#!/usr/bin/env python
import numpy as np 
import pandas as pd 
import sys
import os
import time
import re
# import re
# from re_module import re_module_search
reload(sys)
sys.setdefaultencoding('utf-8')

print('{}\t{}\n'.format(time.ctime(), 'start to analyze!'))


'''
add arguement
'''
import argparse
parser = argparse.ArgumentParser()  
parser.add_argument("-i", "--input", help="input file for analysis", required=True)
#parser.add_argument("-o", "--output", help="output file for analysis", required=True)  
args = parser.parse_args()  


SL_ReDisease_score_withHeader = pd.read_csv('{}'.format(args.input),sep='\t').values
# os.system('head -n 1 {} > {}'.format('SL_ReDisease_score_withHeader.txt','SL_ReDisease_score_withHeader_header.txt'))

ReDisease_MutationGene_Score = pd.read_csv('ReDisease_MutationGene_Score.txt',sep='\t').values

print(len(SL_ReDisease_score_withHeader))
print(len(ReDisease_MutationGene_Score))
# print('Step1 has been completed!')

SL_ReDisease_score = {}
ReDisease_MutationGene = {}

for i in range(len(SL_ReDisease_score_withHeader)):
	_key = str(SL_ReDisease_score_withHeader[i][0]).lower()+"#"+str(SL_ReDisease_score_withHeader[i][1]).lower()+"#"+str(SL_ReDisease_score_withHeader[i][2]).lower()+"#"+str(SL_ReDisease_score_withHeader[i][3]).lower()+"#"+str(SL_ReDisease_score_withHeader[i][4]).lower()+"#"+str(SL_ReDisease_score_withHeader[i][5]).lower()
	_value = np.array(SL_ReDisease_score_withHeader[i][0:],dtype = str).tolist()
	SL_ReDisease_score[_key] = _value

print(len(SL_ReDisease_score))
# print('Step2 has been completed!')

for i in range(len(ReDisease_MutationGene_Score)):
	_key = str(ReDisease_MutationGene_Score[i][0]).lower()+"#"+str(ReDisease_MutationGene_Score[i][4]).lower()+"#"+str(ReDisease_MutationGene_Score[i][1]).lower()
	_value = [str (i) for i in [str(ReDisease_MutationGene_Score[i][1])]+[str(ReDisease_MutationGene_Score[i][4])]+[str(ReDisease_MutationGene_Score[i][5])]]
	ReDisease_MutationGene[_key] = _value

print(len(ReDisease_MutationGene))
# print('Step3 has been completed!')
# print('Dictionary has been made! Start to search!')
#print(SL_ReDisease_score['10018#lymphoma;large-cell;follicular#bcl2'])
#print(ReDisease_MutationGene['10018#lymphoma;large-cell;follicular#bcl2'])
ReDisease_MutationGene_Score={}

ReDisease_MutationGene_index = {}
for _key in ReDisease_MutationGene.keys():
    id_3 = _key.split('#')[1]
    id_4 = _key.split('#')[2]
    new_key = id_3 + '#' + id_4
    _index = ReDisease_MutationGene.keys().index(_key)

    if ReDisease_MutationGene_index.has_key(new_key):
        ReDisease_MutationGene_index[new_key].append(_index)
    else:
        ReDisease_MutationGene_index[new_key] = [_index]
# print(ReDisease_MutationGene_index['lymphoma, follicular#mtor'])
def re_module_search(input1, input2):
    char = ''
    # for key1 in input1.keys():
    #     char += key1 + '\n'
    
    # n = re.findall('.*' + 'ki-1+ anaplastic large cell lymphoma' + '.*', char, re.MULTILINE)
    # if len(n) > 0:
    # 	print(n)
    # 	sys.exit(0)
    # else:
    # 	pass

    # sys.exit(0)

    for key2 in input2.keys():
        id_1 = key2.split('#')[2]
        id_2 = key2.split('#')[5]

        key3 = id_1 + '#' + id_2

        if ReDisease_MutationGene_index.has_key(key3):
            key3_index = ReDisease_MutationGene_index[key3]

            choosed_keys = np.array(input1.keys())[key3_index].tolist()

            for item in choosed_keys:
                char += str(item) + '\n'

            m = re.findall('.*' + id_1 + '#' + id_2 + '.*', char, re.MULTILINE)
            # m = re.findall('.*' + id_1 + '.*' + '#' + id_2 + '.*', char, re.MULTILINE)

            # if len(re.findall('.*' + id_1 + '#' + id_2 + '.*', char, re.MULTILINE)) > 0:
            #     m = re.findall('.*' + id_1 + '#' + id_2 + '.*', char, re.MULTILINE)
            # elif len(re.findall('.*' + id_1 + '#' + id_2 + '.*', char, re.MULTILINE)) == 0 and len(re.findall('.*' + id_1 + '.*' + '#' + id_2 + '.*', char, re.MULTILINE)) > 0:
            #     m = re.findall('.*' + id_1 + '.*' + '#' + id_2 + '.*', char, re.MULTILINE)
            # else:
            #     m = []

            # n = re.findall('.*anaplastic large cell lymphoma.*', key2)

            # if len(n) > 0:
            #     #print([key2]+[m[0]]+ReDisease_MutationGene[m[0]])
            #     print([n, id_1 + "#" + id_2, m])
            #     # print(m)
            #     #print(n)
            #     sys.exit(0)
            # else:
            #     pass

            #sys.exit(0)

            if len(m) > 0:            
                for i in range(len(m)):
                    id_3 = m[i].split('#')[1]
                    id_4 = m[i].split('#')[2]
                    if id_1 == id_3 and id_2 == id_4:
                        #print([key2,m[i]])
                        if ReDisease_MutationGene_Score.has_key(key2):
                            ReDisease_MutationGene_Score[key2].append(SL_ReDisease_score[key2] +ReDisease_MutationGene[m[i]])
                        else:
                            ReDisease_MutationGene_Score[key2] = [SL_ReDisease_score[key2] + ReDisease_MutationGene[m[i]]]
                    else:
                        ReDisease_MutationGene_Score[key2]= [SL_ReDisease_score[key2] +['&&']]
            else:
                ReDisease_MutationGene_Score[key2]= [SL_ReDisease_score[key2] +['%']]
                #print([key2,])
        else:
            #pass
            print([key2, key3])
        # print(ReDisease_MutationGene_Score.values()[0:5])
    return ReDisease_MutationGene_Score

re_module_search(ReDisease_MutationGene, SL_ReDisease_score)
print(len(ReDisease_MutationGene_Score))

with open('./ReDiseaseScore.txt','w') as f:
    f.write('MutationGene\tRepositionDisease\tReDiseaseScore\n')
f.close()

'''
SL_ReDisease_score_withHeader.txt
SL_ReDisease_MutationGene_score.txt
'''
prefix = args.input.split('.')[0]
with open('{}_MutationGene_score.txt'.format(prefix),'w') as g:
    for key in ReDisease_MutationGene_Score.keys():
        for j in range(len(ReDisease_MutationGene_Score[key])):
            print >>g, ('\t').join([str(i) for i in ReDisease_MutationGene_Score[key][j]])
        

os.system('paste {} {} | tr -d "\r" > {}'.format('SL_ReDisease_score_withHeader_header.txt','ReDiseaseScore.txt','add_ReDiseasescore_header.txt'))
# 'SL_ReDisease_MutationGene_score_withHeader.txt'
os.system('sort {} | uniq > {}'.format('{}_MutationGene_score.txt'.format(prefix), '{}_MutationGene_score.temp.txt'.format(prefix)))
os.system('cat {} {} > {}'.format('add_ReDiseasescore_header.txt','{}_MutationGene_score.temp.txt'.format(prefix),'{}_MutationGene_score_withHeader.txt'.format(prefix)))
os.system('rm {}'.format('{}_MutationGene_score.temp.txt'.format(prefix)))

print('{}\t{}\n'.format(time.ctime(), 'end to analyze!'))

		