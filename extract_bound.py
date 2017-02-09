def extract_bound(boundgenelist,target_file,output_boundlist):
	data0=open(boundgenelist,'r')
	data01={}
	for rec in data0:
		gene=[x.strip() for x in rec.split('\t')][0]
		data01[gene]=' '

	data1=open(target_file,'r')
	result=open(output_boundlist,'w')
	for records in data1:
		tmp_data=[x.strip() for x in records.split('\t')]
		if tmp_data[8] in data01:
			result.write(tmp_data[0]+'\t'+str(int(tmp_data[3])-1)+'\t'+tmp_data[4]+'\t'+tmp_data[8]+'\t'+tmp_data[7]+'\t'+tmp_data[6]+'\n')
	result.close()
	data1.close()
	data0.close()

############################################################################
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hb:t:o:")
	except getopt.GetoptError:
		print 'one_conv_layer_net.py -s <image d2> -t <image d3> -f <image d4> -i <filter1 out> -e <filter2 out> -u <full connect1 out> -l <filter1 d1> -o <filter1 d2> -r <filter2 d1> -v <filter2 d1> -a iter_num -b batch_size -p train_pos_num -n train_neg_num -x max_pooling_1 -y max_pooling_2 -k visual_k-mer_mutation -m visual_Second_k-mer_mutation -d training_speed -j seq_num -q shape_num -w file_list_all -z file_list_random'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'one_conv_layer_net.py -s <image d2> -t <image d3> -f <image d4> -i <filter1 out> -e <filter2 out> -u <full connect1 out> -l <filter1 d1> -o <filter1 d2> -r <filter2 d1> -v <filter2 d1> -a iter_num -b batch_size -p train_pos_num -n train_neg_num -x max_pooling_1 -y max_pooling_2 -k visual_k-mer_mutation -m visual_Second_k-mer_mutation -d training_speed -j seq_num -q shape_num -w file_list_all -z file_list_random'
			sys.exit()
		elif opt=="-b":
			boundgenelist=str(arg.strip())
		elif opt=="-t":
			target_file=str(arg.strip())
		elif opt=="-o":
			output_boundlist=str(arg.strip())
	extract_bound(boundgenelist,target_file,output_boundlist)

if __name__=="__main__":
	main(sys.argv[1:])





