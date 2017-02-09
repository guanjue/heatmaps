def gene_group_split(target_filename,target_id_column,gene_group_list,gene_group_id_column,header):
	gene_group0=[]
	gene_group_file0=open(gene_group_list,'r')
	for records in gene_group_file0:
		gene_group0.append(records.split()[0])
	#gene_group0=['RP_137_genes_Xu_2009','SAGA_Xu_2009_ORF_Ts_V64','TFIID_Xu_2009_ORF_Ts_V64']
	target=open(target_filename,'r')
	target0=[]
	for records in target:
		target0.append(records.split())

	target0_genegroup = open(target_filename+'.genegroup','w')
	if header=='T':
		## write header
		for rec in target0[0]:
			target0_genegroup.write(rec+'\t')
		target0_genegroup.write('\n')

	for gene_group in gene_group0:
		### read gene group gene names
		#gene_group_file = open('/Volumes/MAC_Data/data/labs/pugh_lab/master_ref/nfr_based/gene_group/'+gene_group+'.genename.txt','r')
		gene_group_file = open(gene_group,'r')
		gene_group1={}
		for records in gene_group_file:
			gene_group1[records.split()[gene_group_id_column]]=' '
		gene_group_file.close()

		## write data
		if header=='T':
			start_row=1
		else:
			start_row=0
		for records in target0[start_row:]:
			if records[target_id_column] in gene_group1:
				for rec in records:
					target0_genegroup.write(rec+'\t')
				target0_genegroup.write('\n')
	target0_genegroup.close()


############################################################################
### python gene_group_split.py -t Xu_2009_ORF-Ts_V64_sort_pk_plus1.bed -a 3 -g gene_group_list.txt -b 0 -i F
import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"h:t:a:g:b:i:")
	except getopt.GetoptError:
		print 'python gene_group_split.py -t <target_filename> -a <target_id_column> -g <gene_group_list> -b <gene_group_id_column>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python gene_group_split.py -t <target_filename> -g <gene_group_list> -g <gene_group_list> -b <gene_group_id_column> -i <header T or F>'
			sys.exit()
		elif opt=="-t":
			target_filename=str(arg.strip())
		elif opt=="-a":
			target_id_column=int(arg.strip())
		elif opt=="-g":
			gene_group_list=str(arg.strip())
		elif opt=="-b":
			gene_group_id_column=int(arg.strip())
		elif opt=="-i":
			header=str(arg.strip())

	gene_group_split(target_filename,target_id_column,gene_group_list,gene_group_id_column,header)
if __name__=="__main__":
	main(sys.argv[1:])
