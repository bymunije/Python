import numpy as np
import pandas as pd
import os
import time
import sys
import importlib
importlib.reload(sys)

print('{}\t{}\n'.format(time.ctime(),'Start to analyze!'))

# reading DR_233009 and SL_drug_matching
DR_233009 = pd.read_csv('./Drug_repositioning_233009.txt', sep='\t').values
SL_drug_matching_353 = pd.read_csv('./Hungarian_max_matching/SL_drug_matching_withHeader.txt', sep='\t').values

DR_233009_dict = {}
SL_drug_matching_353_dict = {}

for i in range(len(DR_233009)):
    _key = ('#').join(np.array(DR_233009[i][:], dtype=str))    
    _value = np.array(DR_233009[i][:], dtype=str).tolist()
    DR_233009_dict[_key] = _value
print(len(DR_233009_dict))

for i in range(len(SL_drug_matching_353)):
    _key = str(SL_drug_matching_353[i][0]) + '#' + str(
        SL_drug_matching_353[i][1]) + '#' + str(
            SL_drug_matching_353[i][3]) + '#' + str(SL_drug_matching_353[i][4])
    _value = np.array(SL_drug_matching_353[i][:], dtype=str).tolist()
    SL_drug_matching_353_dict[_key] = _value
print(len(SL_drug_matching_353_dict))

# gather kowndisease
gather_kowndisease_DR = {}
gather_kowndisease_list = []
for key in DR_233009_dict:
    new_key = str(DR_233009_dict[key][0]) + "#" + '*'.join(
        DR_233009_dict[key][2:9])
    new_value = [DR_233009_dict[key][0]] + np.array(DR_233009_dict[key][2:],
                                                    dtype=str).tolist()
    if new_key in gather_kowndisease_DR:
        if DR_233009_dict[key][1] in gather_kowndisease_DR[new_key]:
            pass
        else:
            gather_kowndisease_DR[new_key].append(DR_233009_dict[key][1])
            gather_kowndisease_DR[new_key] = gather_kowndisease_DR[new_key]
    else:
        gather_kowndisease_DR[new_key] = new_value + [DR_233009_dict[key][1]]
print(len(gather_kowndisease_DR))

gather_kowndisease_DR_modify = {}
for key in gather_kowndisease_DR:
    _value = gather_kowndisease_DR[key][0:18] + [
        ','.join(gather_kowndisease_DR[key][18:])
    ]
    _value[-1]=",".join((lambda x:(x.sort(),x)[1])(list(_value[-1].split(',')[:])))
    gather_kowndisease_DR_modify[key] = _value
with open('./maxMatching_drugCom/gather_kowndisease_DR_168402.txt', 'w') as gkd:
    for key in gather_kowndisease_DR_modify:
        print('\t'.join([str(i) for i in gather_kowndisease_DR_modify[key]]),file=gkd)

# selected_DR and gather Redisease
selected_DR = {}
for key in gather_kowndisease_DR_modify:
    _new_key = str(gather_kowndisease_DR_modify[key][4]) + "#" + str(
        gather_kowndisease_DR_modify[key][0]) + "#" + str(
            gather_kowndisease_DR_modify[key][3]) + "#" + str(
                gather_kowndisease_DR_modify[key][5])
    _new_value = [gather_kowndisease_DR_modify[key][0]] + np.array(
        gather_kowndisease_DR_modify[key][2:],
        dtype=str).tolist() + [gather_kowndisease_DR_modify[key][1]]
    if _new_key in SL_drug_matching_353_dict:
        if _new_key in selected_DR:
            if gather_kowndisease_DR_modify[key][1].lower() in [i.lower() for i in 
                    selected_DR[_new_key][17].split(',')[:]]:
                pass
            else:
                selected_DR[_new_key].append(
                    gather_kowndisease_DR_modify[key][1])
                selected_DR[_new_key] = selected_DR[_new_key]
        else:
            selected_DR[_new_key] = _new_value
    else:
        pass
print(len(selected_DR))

gather_Redisease_DR = {}
for key in selected_DR:
    _value = selected_DR[key][0:18] + [
        ','.join(selected_DR[key][18:])
    ]
    _value[-1]=",".join((lambda x:(x.sort(),x)[1])(list(_value[-1].split(',')[:])))
    num = len(_value[-1].split(',')[:])
    gather_Redisease_DR[key] = _value + [str(num)]
with open('./maxMatching_drugCom/gather_Redisease_DR_353.txt', 'w') as grd:
    for key in gather_Redisease_DR:
        print('\t'.join([str(i) for i in gather_Redisease_DR[key]]),file=grd)

# drug_combination
drug_combination_matching = {}
for key, value in gather_Redisease_DR.items():
    drug_combination_matching_list = []
    mutation_g = value[4]
    target_g = value[2]
    ReDisease = value[-1]
    for key_in, value_in in gather_Redisease_DR.items():
        mutation_g_in = value_in[4]
        target_g_in = value_in[2]
        ReDisease_in = value_in[-1]
        if key == key_in:
            pass
        elif mutation_g == mutation_g_in and target_g != target_g_in and ReDisease == ReDisease_in:
            if key in drug_combination_matching:
                drug_combination_matching_list.append(value_in)
                num = len(drug_combination_matching_list)
                drug_combination_matching[key] = gather_Redisease_DR[
                        key]  + [str(num)] + drug_combination_matching_list             
            else:
                drug_combination_matching[key] = gather_Redisease_DR[
                        key] + ['1'] + gather_Redisease_DR[key_in] 
        else:
            pass
with open('./maxMatching_drugCom/drug_combination_matching_DR.txt', 'w') as dcm:
    for key in drug_combination_matching:
        print('\t'.join([str(i) for i in drug_combination_matching[key]]), file=dcm)

# split_drug_combination
with open('./maxMatching_drugCom/drug_combination_matching_split_DR.txt', 'w') as dcms:
    for key in drug_combination_matching:       
        if drug_combination_matching[key][20] == '1':
            for i in drug_combination_matching[key][21:]:
                if isinstance(i,str):
                    print('\t'.join([
                        str(i) for i in drug_combination_matching[key][0:20] +
                        drug_combination_matching[key][21:]
                    ]),file = dcms)
                else:
                    for m in drug_combination_matching[key][21:]:
                        print('\t'.join([
                            str(i) for i in drug_combination_matching[key][0:20] +
                            m
                        ]),file = dcms)
        else:
            j = 21
            for m in range(int(drug_combination_matching[key][20])):
                print('\t'.join([
                    str(i) for i in drug_combination_matching[key][0:20] +
                    drug_combination_matching[key][j + m]
                ]),file = dcms)

#----------------
# after select_before_after.sh
#---------------

# # read file that drug_combinatiion_3689
# select_before_after_modify = pd.read_csv(
#     './maxMatching_drugCom/select_before_after_modify.txt', sep='\t', header=None).values

# select_before_after_modify_dict = {}

# for i in range(len(select_before_after_modify)):
#     _key = ('#').join(
#         np.array(select_before_after_modify[i][:], dtype=str))
#     _value = np.array(select_before_after_modify[i][:],
#                       dtype=str).tolist()
#     select_before_after_modify_dict[_key] = _value
# print(len(select_before_after_modify_dict))

# # one Redisease and drug_combination
# ReDisease = []
# with open('./maxMatching_drugCom/drug_combination_split.txt', 'w') as dcs:
#     for key in select_before_after_modify_dict:
#         ReDisease = select_before_after_modify_dict[key][-2].split(',')[:]
#         for Redis in ReDisease:
#             dcs.write('\t'.join([
#                 str(i)
#                 for i in select_before_after_modify_dict[key][0:18] + select_before_after_modify_dict[key][20:38] ] + [Redis+'\n']
#             ))

print('{}\t{}\n'.format(time.ctime(),'End to analyze!'))