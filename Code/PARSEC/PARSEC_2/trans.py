import os

lines = open('list.txt','r').readlines()

for line in lines:
	cos = line.split()
	fl = cos[5].split('tmp/')[-1]

	os.system('cp '+fl+' '+cos[6])