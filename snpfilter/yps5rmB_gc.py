import sys

f1 = open(sys.argv[1],'r')
f2 = open(sys.argv[2],'r')

ff1 = f1.readlines()
ff2 = f2.readlines()

f1.close()
f2.close()

i = 0 
j = 0 

while i < len(ff1) or j < len(ff2):
	if i >= len(ff1):
		p1 = [38000]
		p2 = ff2[j].split(" ")
	elif j >= len(ff2):
		p2 = [38000]
		p1 = ff1[i].split(" ")
	else:
		p1 = ff1[i].split(" ")
		p2 = ff2[j].split(" ")
	
	if int(p1[0]) == int(p2[0]):
		group = p1[0]
		comset = set(p1[2].split(",")) & set(p2[2].split(","))
		#print set(p1[2].split(","))
		com_n = len(comset)
		p1Uniq_n = int(p1[1]) - com_n
		p2Uniq_n = int(p2[1]) - com_n
		print group,p1Uniq_n,p2Uniq_n,com_n
		i += 1 
		j += 1
	elif int(p1[0]) > int(p2[0]):
		group = p2[0]
		print group,0,int(p2[1]),0
		j += 1 
	else:
		group = p1[0]
		print group,int(p1[1]),0,0
		i += 1 
