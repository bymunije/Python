import matplotlib
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas as pd
from networkx.algorithms import bipartite
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import community
import os
import sys,time
from pulp import *
from networkx.algorithms.bipartite import biadjacency_matrix
import scipy
import importlib
importlib.reload(sys)

print('{}\t{}\n'.format(time.ctime(), 'start to analyze!'))

# draw bipartite considering weight
def draw_bipartite(G_each,sl_tuple_list,drug_tuple_list,sl_drug_tuple_list,file,drug_list_each,sl_list_each):
    # import warnings
    # warnings.filterwarnings("ignore", category=UserWarning)
    # picture size
    fig = plt.figure(1, figsize=(40, 40))
    pos_each = dict()
    # Update position for node from each group
    if len(drug_tuple_list) > len(sl_tuple_list)  :
        pos_each.update({
            n[0]: [1 - 0.5, i  * len(drug_tuple_list) / len(sl_tuple_list)]
            for i, n in enumerate(drug_tuple_list)
        })
        pos_each.update({
            n[0]: [2 - 0.5, i * (len(drug_tuple_list) / len(sl_tuple_list) + len(sl_tuple_list)/1.5)]
            for i, n in enumerate(sl_tuple_list)
        })
    else:
        pos_each.update({
            n[0]: [1 - 0.5, i  * (len(sl_tuple_list) / len(drug_tuple_list))]
            for i, n in enumerate(drug_tuple_list)
        })
        pos_each.update({
            n[0]: [2 - 0.5, i * (len(sl_tuple_list) / len(drug_tuple_list) + len(drug_tuple_list)/1.5)]
            for i, n in enumerate(sl_tuple_list)
        })
    # nodes
    drug_nodes = nx.draw_networkx_nodes(G_each,
                                        pos_each,
                                        nodelist=drug_list_each,
                                        node_color = '#d25959',
                                        alpha=0.95,
                                        node_size=150,
                                        node_shape='s')
    SL_nodes = nx.draw_networkx_nodes(G_each,
                                      pos_each,
                                      nodelist=sl_list_each,
                                      node_color = '#681313',
                                      alpha=0.95,
                                      node_size=400)
    # edges
    edgewidth = []
    for (u,v,d) in G_each.edges(data = True):
        edgewidth.append(d['weight'] * 2)
    # print(pos_each)
    nx.draw_networkx_edges(G_each,
                           pos_each,
                           edgelist=sl_drug_tuple_list,
                           width=edgewidth,
                           alpha=0.7,
                           edge_color='#5c7658')
    # Set edge color 
    drug_nodes.set_edgecolor('#e6d385')
    SL_nodes.set_edgecolor('#e6d385')
    # for labeling outside the node
    pos_each_labels = {}
    for key in pos_each.keys():
        x, y = pos_each[key]
        if pos_each[key][0] == 0.5:
            pos_each_labels[key] = (x+0.08, y)
        else:
            pos_each_labels[key] = (x-0.06 , y-0.05)
    # labels name
    labels_name = {}
    for node in G_each.nodes():
        labels_name["".join(node)] = "".join(node)
    nx.draw_networkx_labels(G_each,
                            pos_each_labels,
                            labels_name,
                            font_size=17)
    #                         font_family="monospace",
    #                         font_weight="medium")
    # bankground color
    fig.set_facecolor("#262626")
    plt.savefig("/home/byzhang/work_data/bipatite/bipartite_pult/picture/{}_bipartite.pdf".format(file))
    plt.close()

def bipartite_pulp_drugAntogonism(sl,drug,sl_drug,file):
    print(file)
	# select sl and max value
    # Sl and drug corresponding score more than one
    sl_drug_selected_max = {}
    for line in sl_drug:
        _key = '#'.join(line[0:2])
        if _key in sl_drug_selected_max:
            _value.append(line[2])
            _max = max(_value[2:])
            sl_drug_selected_max[_key] = _value[0:2] + [_max]
        else:
            _value = line
            sl_drug_selected_max[_key] = _value
    for key in sl_drug_selected_max:
        sl_drug_selected_max[key] = sl_drug_selected_max[key][0:2] + [
            round(sl_drug_selected_max[key][2],2)
        ]
    # print(list(sl_drug_selected_max.values()))

    # transform the file format and set attibute of nodes and edges considering weight
    sl_tuple_list = [tuple(i) for i in sl]
    drug_tuple_list = [tuple(i) for i in drug]
    sl_drug_tuple_list = [tuple(i) for i in list(sl_drug_selected_max.values())]
    #Construct graph
    G_each = nx.Graph()
    G_each.add_nodes_from(sl_tuple_list, bipartite=0)
    G_each.add_nodes_from(drug_tuple_list, bipartite=1)
    G_each.add_weighted_edges_from(sl_drug_tuple_list,weight='weight')

    # drug and sl list 、weight
    drug_list_each = []
    for i in drug:
        for j in i:
            drug_list_each.append(j)
    sl_list_each = []
    for i in sl:
        for j in i:
            sl_list_each.append(j)

    if len(sl_list_each) < 2 or len(drug_list_each) < 2:
        pass
    else:
        draw_bipartite(G_each,sl_tuple_list,drug_tuple_list,sl_drug_tuple_list,file,drug_list_each,sl_list_each)
        # try:
        #     # draw bipartite considering weight
        #    draw_bipartite(G_each,sl_tuple_list,drug_tuple_list,sl_drug_tuple_list,flie) 
        # except:
        #     with open('./Error.txt','a') as err:
        #         err.write('%s have error! \n' % (file))

        # convert bipartite graph to adjacency matrix considering weight
        G_each_biadjacency_matrix = biadjacency_matrix(
            G_each,
            row_order=drug_list_each,
            column_order=sl_list_each,
            weight='weight').toarray()
        G_each_df = pd.DataFrame(G_each_biadjacency_matrix,
                                 index=drug_list_each,
                                 columns=sl_list_each)
        # print(G_each_df)
        G_each_df.to_csv('/home/byzhang/work_data/bipatite/bipartite_pult/adjacency_matrix/{}_G_each_df.txt'.format(file),sep = '\t')

        # convert dataframe to list which is consist of  drug dictory
        drug_list = []
        # print(G_each_df.shape[0])
        for index,row in G_each_df.iteritems():    
            dict_name = {}
            for i in range(len(row.index)):
                dict_name[row.index[i]] = row.values[i]
            drug_list.append(dict_name)
        #     for i in range(1,len(G_each_df.iloc[1,:])):
        #         for j in range(len(G_each_df.iloc[:,i])):
        # print(drug_list)

        '''
        print all drug combination
        '''
        # interger programming
        # 创建问题实例
        prob = LpProblem("bipartite_drug_combine", LpMaximize)

        # 构建Lp变量字典，变量名以Ingr开头，如Ingr_CHICKEN，下界是0
        ingredient_vars = LpVariable.dicts("Ingr",
                                           drug_list_each,
                                           lowBound=0,
                                           upBound=1,
                                           cat=LpInteger)
        # # 添加目标方程
        prob += lpSum([[sl[i] * ingredient_vars[i] for i in drug_list_each]
                       for sl in drug_list])
        # prob += lpSum(
        #     [[0.1 * sl[i] * ingredient_vars[i] for i in drug_list_each]
        #      for sl in drug_list] +i
        #     [[0.9 * drug[i] * ingredient_vars_sl[i] for i in sl_list_each]
        #     for drug in sl_list])
        # #添加约束条件
        prob += lpSum([ingredient_vars[i]
                       for i in drug_list_each]) <= 2  # only select two drug
        # prob += lpSum([drug_list[1][i] * ingredient_vars[i]
        #                for i in drug_list_each]) >= 1  # Choose a drug to cover more sl
        # prob += lpSum([drug_list[1][i] * ingredient_vars[i]
        #                for i in drug_list_each]) <= 2  # Choose a drug to cover more sl
        # prob += lpSum([[ingredient_vars[i[0]] + ingredient_vars[i[1]]
        #                if i in drug_antagonism] for i in drug_list_each]) <= 2
        #求解
        prob.solve()
        #查看解的状态
        print("Status:", LpStatus[prob.status])
        print('max result:', value(prob.objective))
        Max_value = value(prob.objective)

        # print best varible
        def comb_drug():
        #'  get the combination of variables
            print('The best variable:\n')
            selected_com_drug = []
            for v in range(len(prob.variables())):
                if prob.variables()[v].varValue> 0:
                    selected_com_drug.append(prob.variables()[v].name)
            #         for drug in drug_list_each:
            #             if str(v.name) == str(ingredient_vars[drug]):
            #                 name = drug
            #         selected_com_drug.append(name)
                    print(prob.variables()[v].name, "=", prob.variables()[v].varValue)
            return(selected_com_drug)

        def get_key ():
        #' Get an identifiable key according to selected_com_drug
            selected_com_drug=comb_drug()
            normal_vars=dict()
            vv=list()
            for i in ingredient_vars:
                normal_vars[i] = str(ingredient_vars[i]) 
            for j in selected_com_drug:
                x = [k for k, v in normal_vars.items() if v == j] 
                vv.append(x[0])
            return vv
        new_st=get_key()

        # Loop output
        st_list=list()
        st_list.append(new_st)
        while(value(prob.objective)==Max_value):
            prob = LpProblem("bipartite_drug_combine", LpMaximize)

            # 构建Lp变量字典，变量名以Ingr开头，如Ingr_CHICKEN，下界是0
            ingredient_vars = LpVariable.dicts("Ingr",
                                               drug_list_each,
                                               lowBound=0,
                                               upBound=1,
                                               cat=LpInteger)
            # # 添加目标方程
            prob += lpSum([[sl[i] * ingredient_vars[i] for i in drug_list_each]
                               for sl in drug_list])
            # #添加约束条件
            prob += lpSum([ingredient_vars[i] for i in drug_list_each]) <= 2
            
                  # only select two drug
            for new_st in st_list:
                # exec("prob+=lpSum([ingredient_vars[i] for i in new_st]) <= 1")
                prob+=lpSum([ingredient_vars[i] for i in new_st]) <= 1
            prob.solve()
            #查看解的状态
            new_st=get_key()
            st_list.append(new_st) #all combination
            print('max result:',value(prob.objective))

        # choose to cover more sl
        def most_sl(st_list, G_each_df):
            sl_number_dict = dict()
            for combdrug in st_list:
                df_cd = G_each_df.loc[combdrug]
                df_above0 = df_cd.apply(lambda x: x.sum())
                sl_number = len([i for i in df_above0 if i > 0])
                sl_number_dict[str(combdrug)] = sl_number
            most_sl = max(sl_number_dict.values())
            comb = [m for m, n in sl_number_dict.items() if n == most_sl]
            print('Those drug combination can cover ', most_sl, 'kinds SL')
            for i in comb:
                print(i)
            return(comb)
        combdrug_opt = most_sl(st_list[0:len(st_list)-1], G_each_df)
        print(combdrug_opt)

        # drug antagonism file
        import ast
        drug_antagonism = pd.read_csv(
            '/home/byzhang/work_data/bipatite/bipartite_pult/Drug_antagonism/select_antagonism.csv',
            sep=',').values.tolist()
        # Antagonistic drug name changed to all lowercase
        all_to_lower = []
        for i in drug_antagonism:
            each_to_lower = []
            for j in i:
                each_to_lower.append(str(j).lower())
            all_to_lower.append(each_to_lower)
        # print(all_to_lower)
        # Determine whether the drug combination is in drug antagonism
        choose_no_antagonism = []
        choose_antagonism = []
        combdrug_opt_tolist = [ast.literal_eval(i.lower()) for i in combdrug_opt]
        # print(combdrug_opt_tolist)
        # for i in combdrug_opt:
        for i in combdrug_opt_tolist:
        #     print((i))
            if i in all_to_lower:
                choose_antagonism.append(i) 
            else:
                choose_no_antagonism.append(i)
        print('choose_no_antagonism:{}'.format(choose_no_antagonism))
        print('chooseed_antagonism:{}'.format(choose_antagonism))
        # str to list
        # choose_no_antagonism_tolist = []
        # choose_no_antagonism_tolist = [ast.literal_eval(i) for i in choose_no_antagonism]
        # choose_antagosim_tolist = []
        # print(choose_no_antagonism_tolist)

        # drug combiantion mapping sl
        all_mapping_sl = []
        for i in choose_no_antagonism:
        # for i in choose_no_antagonism_tolist:
            each_mapping_sl = []
            for j in i:  
                for key in sl_drug_selected_max:
        #             if j in sl_drug_selected_max[key]:
                    if j.capitalize() in sl_drug_selected_max[key]:
                        each_mapping_sl.append(sl_drug_selected_max[key])
            all_mapping_sl.append(each_mapping_sl)
        print(all_mapping_sl)

        # write file
        with open('/home/byzhang/work_data/bipatite/bipartite_pult/Drugcom/{}_drugcom.txt'.format(file),'w') as fdc:
        # with open('/home/byzhang/work_data/bipatite/bipartite_pult/Drugcom/VHL_drugcom.txt','w') as fdc:
            drug_com_dict = {}
            num = 0
            for i in choose_no_antagonism:
#           for i in choose_no_antagonism_tolist:
                _key = str(i[0].capitalize()) + "#" + str(i[1].capitalize())
#                _key = str(i[0]) + "#" + str(i[1])
                drug1 = []
                drug2 = []
                len_sl = []
                if num < len(all_mapping_sl):   
                    for j in all_mapping_sl[num]:
                        if j[1] in len_sl:
                            pass
                        else:
                            len_sl.append(j[1])
                        if j[0] == _key.split('#')[0]:
                            drug1.append(j[1])
                        elif j[0] == _key.split('#')[1]:
                            drug2.append(j[1])
                        else:
                            continue
                        _value = ';'.join(
                        [str(i) for i in drug1]) + '\t' + ';'.join([str(i) for i in drug2]) + '\t' + str(len(len_sl))
                        drug_com_dict[_key] = _value  
                    num += 1
                else:
                    pass
            print(list([drug_com_dict.items()]))
            for key in drug_com_dict:
                key_trans = '+'.join(key.split('#'))
                fdc.write('\t'.join([str(i) for i in [key_trans] + [file] + drug_com_dict[key].split('\t')]) + '\n')
    return()

# Cyclic reading of files considering weight
all_files = os.listdir(
   '/home/byzhang/work_data/bipatite/bipartite_pult/three_file_every_mulTarget_weight/')
for file_each in all_files:
    mut_three_files = os.listdir(
        '/home/byzhang/work_data/bipatite/bipartite_pult/three_file_every_mulTarget_weight/'
        + file_each + '/')
    for every_file in mut_three_files:
        if every_file.split('_')[1] == 'SL':
            sl_file = pd.read_csv(
                '/home/byzhang/work_data/bipatite/bipartite_pult/three_file_every_mulTarget_weight/'
                + file_each + '/' + every_file,header = None).values.tolist()
        elif every_file.split('_')[1] == 'drug':
            drug_file = pd.read_csv(
                '/home/byzhang/work_data/bipatite/bipartite_pult/three_file_every_mulTarget_weight/'
                + file_each + '/' + every_file,header = None).values.tolist()
        else:
            sl_drug_file = pd.read_csv(
                '/home/byzhang/work_data/bipatite/bipartite_pult/three_file_every_mulTarget_weight/'
                + file_each + '/' + every_file,
                sep='\t',header = None).values.tolist()
    bipartite_pulp_drugAntogonism(sl_file,drug_file,sl_drug_file,file_each)

print('{}\t{}\n'.format(time.ctime(), 'end to analyze!'))