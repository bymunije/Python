import timeit
from collections import deque
import networkx as nx
import numpy as np
import pandas as pd
from networkx.algorithms import bipartite
import os
import time
import sys
import importlib
importlib.reload(sys)

print('{}\t{}\n'.format(time.ctime(),'Start to analyze!'))

#==============================================================================
# 构建图和图矩阵
#==============================================================================

# read file about SL_score,drug_score and SL_drug relationship
SL_point = pd.read_csv('/home/byzhang/work_data/bipatite/SL_score_563.txt',
                       sep='\t').values.tolist()

drug_point = pd.read_csv(
    '/home/byzhang/work_data/bipatite/drug_score_max_withheader.txt',
    sep='\t').values.tolist()

SL_drug_edge_type = pd.read_csv(
    '/home/byzhang/work_data/bipatite/drug_SL_type_withHeader.txt',
    sep='\t').values.tolist()

# transform the file format and set attibute of nodes and edges
for i in range(len(SL_drug_edge_type)):
    SL_drug_edge_type[i][2] = {'type': SL_drug_edge_type[i][2]}

for i in range(len(SL_point)):
    SL_point[i][1] = {'SL_score': SL_point[i][1]}

for i in range(len(drug_point)):
    drug_point[i][1] = {'drug_score': drug_point[i][1]}

# SL_point_list = eval('[%s]' % repr(SL_point).replace('[','').replace(']',''))
SL_point_list = [tuple(i) for i in SL_point]
drug_point_list = [tuple(i) for i in drug_point]
SL_drug_edge_type_list = [tuple(i) for i in SL_drug_edge_type]

# Construct graph
G = nx.Graph()
G.add_nodes_from(SL_point_list, bipartite=0)
G.add_nodes_from(drug_point_list, bipartite=1)
G.add_edges_from(SL_drug_edge_type_list)

# graph matrix
G_matrix = nx.to_numpy_matrix(G)


#==============================================================================
# 匈牙利算法
#==============================================================================
class HungarianAlgorithm(object):
    def __init__(self, graph):
        """
        @graph:图的矩阵表示
        """
        self.graph = graph
        self.n = len(graph)

    def find(self, x):
        #         Matched = []
        for i in range(self.n):
            if int(np.array(
                    self.graph[x])[0].tolist()[i]) == 1 and (not self.used[i]):
                self.used[i] = 1  #放入交替路
                if self.match[i] == -1 or self.find(self.match[i]) == 1:
                    self.match[i] = x
                    self.match[x] = i
                    print(x + 1, '->', i + 1)
                    #                     Matched.append((x+1,'->',i+1))
                    return 1
        return 0

    def hungarian1(self):
        """递归形式
        """
        self.match = [-1] * self.n  #记录匹配情况
        self.used = [False] * self.n  #记录是否访问过
        m = 0
        for i in range(self.n):
            if self.match[i] == -1:
                self.used = [False] * self.n
                print('开始匹配:', i + 1)
                m += self.find(i)
#                 print(Matched)
        return m

    def hungarian2(self):
        """循环形式
        """
        #         matched = []
        match = [-1] * self.n  #记录匹配情况
        used = [-1] * self.n  #记录是否访问过
        Q = deque()  #设置队列
        ans = 0
        prev = [0] * self.n  #代表上一节点
        for i in range(self.n):
            if match[i] == -1:
                Q.clear()
                Q.append(i)
                prev[i] = -1  #设i为出发点
                flag = False  #未找到增广路
                while len(Q) > 0 and not flag:
                    u = Q.popleft()
                    for j in range(self.n):
                        if not flag and int(
                                np.array(self.graph[u])[0].tolist()
                            [j]) == 1 and used[j] != i:
                            used[j] = i
                            if match[j] != -1:
                                Q.append(match[j])
                                prev[match[j]] = u  #记录点的顺序
                            else:
                                flag = True
                                d = u
                                e = j
                                while (d != -1):  #将原匹配的边去掉加入原来不在匹配中的边
                                    t = match[d]
                                    match[d] = e
                                    match[e] = d
                                    d = prev[d]
                                    e = t


                                print('mathch:',match)
                                print('prev:',prev)
                                print('deque',Q)
                if match[i] != -1:  #新增匹配边
                    ans += 1
        return ans, match

# call function
h=HungarianAlgorithm(G_matrix)
print (h.hungarian1())

h=HungarianAlgorithm(G_matrix)
print (h.hungarian2())
matched = h.hungarian2()[1]

# edge of max matching 
matching_edge = []
count = 0
for i in range(563):    
    if i not in matched:
        pass
    else:
        _index_SL = matched.index(i)
        print(i+1,'->',(_index_SL+1))
        matching_edge.append(str(i+1) +'->'+ str(_index_SL+1))
        count += 1
print(count,matching_edge)

# drug and SL list
drug_list_all = []
for n in drug_point_list:
    drug_list_all.append(n[0])
SL_list_all = []
for n in SL_point_list:
    SL_list_all.append(n[0])

# index mapping symbol of SL and drug
with open('./Hungarian_max_matching/SL_drug_matching.txt', 'w') as m:
    for i in matching_edge:
        SL_index = int(i.split('->')[0])
        drug_index = int(i.split('->')[1])
        SL_matched = SL_list_all[SL_index - 1]
        A_gene = SL_matched.split('+')[0]
        B_gene = SL_matched.split('+')[1]
        drug_matched = drug_list_all[drug_index - 563 - 1]
        for i in SL_drug_edge_type:
            if SL_matched in i and drug_matched in i:
                if i[2]['type'] == 'B':
                    m.write('\t'.join(
                        [SL_matched, drug_matched, 'B', B_gene, A_gene]) +
                            '\n')
                else:
                    m.write('\t'.join(
                        [SL_matched, drug_matched, 'A', A_gene, B_gene]) +
                            '\n')

with open('./Hungarian_max_matching/SL_drug_matching_header.txt', 'w') as sldh:
    sldh.write('SL\tdrug\tType\tTarget_g\tMutation_g\n')

os.system('cat {} {} > {}'.format(
    './Hungarian_max_matching/SL_drug_matching_header.txt',
    './Hungarian_max_matching/SL_drug_matching.txt',
    './Hungarian_max_matching/SL_drug_matching_withHeader.txt'))

print('{}\t{}\n'.format(time.ctime(),'End to analyze!'))