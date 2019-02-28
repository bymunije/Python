import numpy as np 
import pandas as pd 
import sys
import os
import time

part_number = 10

reload(sys)
sys.setdefaultencoding('utf-8')

print('{}\t{}\n'.format(time.ctime(), 'start to analyze!'))

SL_ReDisease_score_withHeader = pd.read_csv('SL_ReDisease_score_withHeader.txt',sep='\t').values
os.system('head -n 1 {} > {}'.format('SL_ReDisease_score_withHeader.txt','SL_ReDisease_score_withHeader_header.txt'))
print(len(SL_ReDisease_score_withHeader))

a = len(SL_ReDisease_score_withHeader)

for i in range(part_number):
	if i < part_number - 1:
		SL_ReDisease_score_withHeader_part = SL_ReDisease_score_withHeader[int(i*(a/10)):int((i+1)*(a/10))].tolist()
	else:
		SL_ReDisease_score_withHeader_part = SL_ReDisease_score_withHeader[int(i*(a/10)):].tolist()

	filename = "SL_ReDisease_score_withHeader_{}_temp.txt".format(i+1)
	with open(filename,'w') as m:
		for j in SL_ReDisease_score_withHeader_part:
			print >>m, ('\t').join([str(x) for x in j])

	os.system('cat {} {} > {}'.format('SL_ReDisease_score_withHeader_header.txt', filename, "SL_ReDisease_score_withHeader_{}.txt".format(i+1)))
	os.system('rm -r *temp.txt')

