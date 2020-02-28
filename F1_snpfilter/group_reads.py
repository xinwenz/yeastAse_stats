import sys
import numpy as np

def reviseSt(x):
	newstr = ''
	i=0
	while i < len(x):
		if x[i].isdigit():
			if x[i+1].isdigit():
				skip = int(str(x[i])+ str(x[i+1]))
				i += skip + 1
			else:
				i += int(x[i])
			
			newstr = newstr[:-2] + "X"  # "+" and the '.' before it 
		elif x[i] == "^":
			i += 1
		elif x[i] == "$":
			i += 0
		else:
			newstr = newstr + x[i]
		i += 1 
	return newstr


def findrefPos(s):
    return [i for i, letter in enumerate(s) if (letter == ',' or letter == '.')]


f1 = open(sys.argv[1],'r')
gN_store = 0
for line in f1:
	tmpls = line.rstrip().split()
	gN = int(tmpls[0])
	cv = int(tmpls[3])
	
	## make sure the coverage matches positons and names 
	cv_st = reviseSt(tmpls[4])
	cv_name = np.array(tmpls[6].split(','))
	n_post = len(tmpls[5].split(',')) if tmpls[5] != '*' else 0
	n_name = len(cv_name) if cv_name[0] != '*' else 0
	if not (cv == n_post and cv == n_name) : 
		raise Exception('coverage not match')
	
	if gN != gN_store :
		if gN_store > 0:
			print gN_store,len(pool), ",".join(str(e) for e in pool)
		pool = set([])
		pool = pool |  set(cv_name[findrefPos(cv_st)])
	else:
		pool = pool | set(cv_name[findrefPos(cv_st)])
	#print cv == cv_st 
	#print '%s\t%s\t%s\t%s\t%s\t%s'%(gN,cv,cv_st,cv_st_n,cv_posit_n,cv_name_n)
	gN_store = gN
print gN, len(pool), ",".join(str(e) for e in pool)
f1.close()
