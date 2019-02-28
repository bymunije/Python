from multiprocessing import Pool
import os

os.system('python2.7 divide_ten.py')

def sh(*args):
	os.system('python2.7 search.py -i {}'.format(*args))

pool = Pool(processes=10)
for i in range(1,11):
	infile = '{}_{}.txt'.format('SL_ReDisease_score_withHeader',i)
	pool.apply_async(sh, (infile,))
	#sh(infile)
pool.close()
pool.join()
