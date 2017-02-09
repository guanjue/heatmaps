import numpy as np
def bin_row_col(inputname,row_bin,col_bin,outputname):
	data0=open(inputname,'r')
	data01_name=[]
	data01_num=[]
	data0.readline()
	for records in data0:
		tmp=[x.strip() for x in records.split('\t')]
		data01_name.append(tmp[0:2])
		
		if (len(tmp)-2)%row_bin != 0:
			data01_num.append(np.array(tmp[2:-((len(tmp)-2)%row_bin)])) ### if the row cannot be binned, remove end cells
		else:
			data01_num.append(np.array(tmp[2:]))

	data01_num=np.array(data01_num,dtype=float)
	print(data01_num.shape)
	if data01_num.shape[0]%col_bin != 0: ### if the col cannot be binned, remove end cells
		col_num=data01_num.shape[0]-data01_num.shape[0]%col_bin /col_bin 
		data01_num=data01_num[ 0:(data01_num.shape[0]-data01_num.shape[0]%col_bin),: ] 
	else:
		col_num=data01_num.shape[0]/col_bin
	row_num=data01_num.shape[1]/row_bin
	data01_num=data01_num.reshape(col_num,col_bin,row_num,row_bin).mean(axis=3).mean(axis=1)

	print(data01_num.shape)

	header=range(0,row_num)

	result=open(outputname,'w')
	result.write('YORF'+'\t'+'NAME'+'\t')

	for records in header:
		result.write(str(records)+'\t')
	result.write('\n')

	for name,num in zip(data01_name,data01_num):
		result.write(name[0]+'\t'+name[1]+'\t')
		for rec in num:
			result.write(str(rec)+'\t')
		result.write('\n')
	result.close()

############################################################################
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:r:c:o:")
	except getopt.GetoptError:
		print 'python bin_row_col.py -i <inputname> -r <row bin size> -c <column bin size> -o <outputname>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python bin_row_col.py -i <inputname> -r <row bin size> -c <column bin size> -o <outputname>'
			sys.exit()
		elif opt=="-i":
			inputname=str(arg.strip())
		elif opt=="-r":
			row_bin=int(arg.strip())
		elif opt=="-c":
			col_bin=int(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())
	bin_row_col(inputname,row_bin,col_bin,outputname)

if __name__=="__main__":
	main(sys.argv[1:])


