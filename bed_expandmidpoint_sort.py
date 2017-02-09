import numpy as np

def bed_expand_sort(file_name,upstream_expan,downstream_expan,roman2arabic_tf,outputname):
	roman2arabic = {'chrI':'chr1','chrII':'chr2','chrIII':'chr3','chrIV':'chr4','chrV':'chr5','chrVI':'chr6'
	,'chrVII':'chr7','chrVIII':'chr8','chrIX':'chr9','chrX':'chr10','chrXI':'chr11','chrXII':'chr12','chrXIII':'chr13'
	,'chrXIV':'chr14','chrXV':'chr15','chrXVI':'chr16',}

	data0=open(file_name,'r')
	data01=[]

	end_num = 0
	for records in data0:
		if len(records.split())!=0:
				tmp_inf = [x.strip() for x in records.split('\t')]
				if roman2arabic_tf =='T':
					#print(tmp_inf)
					if tmp_inf[0] in roman2arabic: ### remove header & convert roman to arabic number
						midpoint = (int(tmp_inf[1])+int(tmp_inf[2]) )/2   ### get the midpoint of a gene
						if tmp_inf[5] == '+' and midpoint - upstream_expan >=0:
							expand_inf = [roman2arabic[tmp_inf[0]] ] + [ midpoint-upstream_expan, midpoint+downstream_expan, tmp_inf[3], tmp_inf[4], tmp_inf[5], int(tmp_inf[2])-int(tmp_inf[1]) ]
						elif tmp_inf[5] == '-' and midpoint - downstream_expan >=0:
							expand_inf = [roman2arabic[tmp_inf[0]] ] + [ midpoint-downstream_expan, midpoint+upstream_expan, tmp_inf[3], tmp_inf[4], tmp_inf[5], int(tmp_inf[2])-int(tmp_inf[1]) ]
						else:
							end_num = end_num + 1
						data01.append(expand_inf)

				else: ### if the chr in already arabic number
					if tmp_inf[0][0:3] =='chr': ### remove header
						midpoint = (int(tmp_inf[1])+int(tmp_inf[2]) )/2
						if tmp_inf[5] == '+' and midpoint - upstream_expan >=0:
							expand_inf = [tmp_inf[0], midpoint-upstream_expan, midpoint+downstream_expan, tmp_inf[3], tmp_inf[4], tmp_inf[5], int(tmp_inf[2])-int(tmp_inf[1]) ]
						elif tmp_inf[5] == '-' and midpoint - downstream_expan >=0:
							expand_inf = [tmp_inf[0], midpoint-downstream_expan, midpoint+upstream_expan, tmp_inf[3], tmp_inf[4], tmp_inf[5], int(tmp_inf[2])-int(tmp_inf[1]) ]
						data01.append(expand_inf)

	print('number of genes expand the chrom: '+str(end_num) ) ### print the number of genes that expand outside the chrom
	#print(data01)
	data01 = np.array(data01)
	data01 = data01[data01[:,6].astype(int).argsort()] ### sort the genes by the length of the gene

	result=open(outputname,'w')
	for rec in data01:
		result.write(rec[0]+'\t'+rec[1]+'\t'+rec[2]+'\t'+rec[3]+'\t'+rec[4]+'\t'+rec[5]+'\t'+rec[6]+'\n')

	data0.close()
	result.close()

############################################################################
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hf:u:d:r:o:")
	except getopt.GetoptError:
		print 'one_conv_layer_net.py -s <image d2> -t <image d3> -f <image d4> -i <filter1 out> -e <filter2 out> -u <full connect1 out> -l <filter1 d1> -o <filter1 d2> -r <filter2 d1> -v <filter2 d1> -a iter_num -b batch_size -p train_pos_num -n train_neg_num -x max_pooling_1 -y max_pooling_2 -k visual_k-mer_mutation -m visual_Second_k-mer_mutation -d training_speed -j seq_num -q shape_num -w file_list_all -z file_list_random'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'one_conv_layer_net.py -s <image d2> -t <image d3> -f <image d4> -i <filter1 out> -e <filter2 out> -u <full connect1 out> -l <filter1 d1> -o <filter1 d2> -r <filter2 d1> -v <filter2 d1> -a iter_num -b batch_size -p train_pos_num -n train_neg_num -x max_pooling_1 -y max_pooling_2 -k visual_k-mer_mutation -m visual_Second_k-mer_mutation -d training_speed -j seq_num -q shape_num -w file_list_all -z file_list_random'
			sys.exit()
		elif opt=="-f":
			file_name=str(arg.strip())
		elif opt=="-u":
			upstream_expan=int(arg.strip())
		elif opt=="-d":
			downstream_expan=int(arg.strip())
		elif opt=="-r":
			roman2arabic_tf=str(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())
	bed_expand_sort(file_name, upstream_expan, downstream_expan, roman2arabic_tf,outputname)

if __name__=="__main__":
	main(sys.argv[1:])




