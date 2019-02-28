# -*- coding: UTF-8 -*- 
#!/usr/bin/env python
import numpy as np 
import pandas as pd 
import os,re,sys,time
import subprocess

reload(sys)
sys.setdefaultencoding('utf-8')

print('{}\t{}\n'.format(time.ctime(),'Start to analyzed!'))

split_mut_wild = pd.read_csv('split_mut_wild_all_selected.txt',sep = '\t').values
os.system('head -n 1 {} > {}'.format('split_mut_wild_all_selected.txt','split_mut_wild_all_selected_header.txt'))
 

with open('BestCondidation_4034DR.txt') as BestCondidation_4034DR:
	count = 0
	Name = []
	Name_exist = []
	for line in BestCondidation_4034DR:
		count += 1
		try:
			line = line.replace('\n','')#去掉末尾换行符
			SL = line.split('\t')[5]
			Index_drug = line.split("\t")[20]
			drug = line.split("\t")[0]
			strings = 'grep {} split_mut_wild_all_selected.txt | grep "Cancer_Specific" | grep {} | wc -l'.format(SL,drug)
			f = os.popen(strings)
			data = f.readline()
			f.close()
			if int(data) > 0 :
				strings_g = 'grep {} split_mut_wild_all_selected.txt | grep {} '.format(SL,drug)
				g = os.popen(strings_g)
				data_g = g.read()
				fileName = 'split_mut_wild_{}_{}-{}.txt'.format(SL.strip(),drug.strip(), Index_drug.strip())
# 				fileName
				align1 = fileName.split("-")[0]
				# print(Name)
				if Name == []:
					Name.append(align1)
					with open('./split_every_drug_little/' + fileName,'w') as l:
							l.write('%s \n' % (data_g))
							os.system('cat {} {} > {}'.format("split_mut_wild_all_selected_header.txt",'./split_every_drug_little/' + fileName,'./split_every_drug_withHeader/split_mut_wild_{}_{}_{}_withHeader.txt'.format(SL.strip(),drug.strip(), Index_drug.strip())))
				else:
					if align1 in Name:
						Name_exist.append(fileName)
						print(Name_exist)
						with open('./Name_exist.txt','a') as b:
							b.write('{}\n'.format(fileName))
					else:
						with open('./split_every_drug_little/' + fileName,'a') as l:
							l.write('%s \n' % (data_g))
							os.system('cat {} {} > {}'.format("split_mut_wild_all_selected_header.txt",'./split_every_drug_little/' + fileName,'./split_every_drug_withHeader/split_mut_wild_{}_{}_{}_withHeader.txt'.format(SL.strip(),drug.strip(), Index_drug.strip())))
					Name.append(align1)
				g.close()						
			else:
				pass
		except:
			with open('./Error.txt','a') as h:
				h.write('%s:%s have error! \n' % (str(count),line))



print('{}\t{}\n'.format(time.ctime(),'Ended to analyzed!'))
