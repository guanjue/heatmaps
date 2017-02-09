#!/bin/bash
script_bin='/Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/pipeline_bin/'
while read info; do
	id=$(echo $info | awk -F ' ' -v OFS='\t' '{print $1}' )
	bam_file=$(echo $info | awk -F ' ' -v OFS='\t' '{print $2}' )
	echo $id'_experiment.events'
	echo $bam_file'.bam'

	rm -r $id'_gene_TSS_TES'
	mkdir $id'_gene_TSS_TES'
	cp $id'_experiment.events' $id'_gene_TSS_TES'
	cp $bam_file'.bam' $id'_gene_TSS_TES'
	cd $id'_gene_TSS_TES'
	#####################################
	echo 'get multiGPS pks'
	tail -n+8 $id'_experiment.events' | awk -F '\t' -v OFS='\t' '{print $1}' | awk -F ':' -v OFS='\t' '{print $1,$2,$2+1}' | sort -k1,1 -k2,2n > $id'_experiment.sort.midpoint'
	echo 'expand bed file and sort'
	python $script_bin'gff3_expand_sort.py' -f /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/Xu_2009_ORF-Ts_V64_NC_all.gff3 -u 500 -d 100 -r T -o Xu_2009_ORF-Ts_V64_sort.bed
	sort -k1,1 -k2,2n Xu_2009_ORF-Ts_V64_sort.bed > Xu_2009_ORF-Ts_V64_cosort.bed
	echo 'get bound genes'
	bedtools intersect -a Xu_2009_ORF-Ts_V64_cosort.bed -b $id'_experiment.sort.midpoint' -wa | sort -u > $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed'
	python $script_bin'bed_expand_sort.py' -f $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed' -u -500 -d -100 -r F -o $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed'	
	echo extract bound genes
	cat $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'_Xu_2009_ORF-Ts_V64.genelist.txt'
	echo 'expand bed file and sort'
	python $script_bin'bed_expandmidpoint_sort.py' -f $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' -u 2000 -d 2000 -r F -o $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F

	#echo 'generate bam.bai'
	samtools index $bam_file'.bam'
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed.genegroup' -s 6 -n 5 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	Rscript $script_bin'heatmap.R' $id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.tabular' $id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.png'
	### for all genes
	python $script_bin'gff3_expand_sort.py' -f /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/Xu_2009_ORF-Ts_V64_NC_all.gff3 -u 0 -d 0 -r T -o Xu_2009_ORF-Ts_V64_sort_OD.bed
	python $script_bin'bed_expandmidpoint_sort.py' -f Xu_2009_ORF-Ts_V64_sort_OD.bed -u 2000 -d 2000 -r F -o Xu_2009_ORF-Ts_V64_sort_midpoint_expand.bed
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t Xu_2009_ORF-Ts_V64_sort_midpoint_expand.bed -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c Xu_2009_ORF-Ts_V64_sort_midpoint_expand.bed.genegroup -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	python $script_bin'bin_row_col.py' -i 'Xu_2009_ORF-Ts_V64_sort_midpoint_expand_'$bam_file'_read1_combined.tabular' -r 25 -c 5 -o 'Xu_2009_ORF-Ts_V64_sort_midpoint_expand_'$bam_file'_read1_combined_binned.tabular'
	Rscript $script_bin'heatmap.R' 'Xu_2009_ORF-Ts_V64_sort_midpoint_expand_'$bam_file'_read1_combined_binned.tabular' 'Xu_2009_ORF-Ts_V64_midpointexpand_'$bam_file'_read1_combined.png'

	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.gene_len.png' -a 0 -c 20 -f png

	#####################################
	echo 'bound genes TSS'
	cat $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > $id'_Xu_2009_ORF-Ts_V64_TSS.bed'
	echo 'sort by dist to TSS'
	#cat $id'_experiment.sort.midpoint' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",".","+"}' > $id'_experiment.midpoint.pn.txt'
	#cat $id'_experiment.sort.midpoint' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",".","-"}' >> $id'_experiment.midpoint.pn.txt'
	#sort -k1,1 -k2,2n $id'_experiment.midpoint.pn.txt' > $id'_experiment.sort.midpoint.pn.txt'
	bedtools closest -a $id'_Xu_2009_ORF-Ts_V64_TSS.bed' -b $id'_experiment.sort.midpoint' -D a > $id'_Xu_2009_ORF-Ts_V64_TSS.bound.txt'
	sort -k10,10nr $id'_Xu_2009_ORF-Ts_V64_TSS.bound.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-750,$2+251,$4,$5,$6; else print $1,$2-250,$2+751,$4,$5,$6}' > $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed.genegroup' -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	Rscript $script_bin'heatmap.R' $id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' $id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.png'
	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.TSS_ref.png' -a 0 -c 20 -f png

	#####################################
	echo 'bound genes TES'
	cat $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$3-1,$3,$4,$5,$6; else print $1,$2,$2+1,$4,$5,$6}' | sort -k1,1 -k2,2n > $id'_Xu_2009_ORF-Ts_V64_TES.bed'
	echo 'sort by dist to TES'
	#cat $id'_experiment.sort.midpoint' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",".","+"}' > $id'_experiment.midpoint.pn.txt'
	#cat $id'_experiment.sort.midpoint' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",".","-"}' >> $id'_experiment.midpoint.pn.txt'
	#sort -k1,1 -k2,2n $id'_experiment.sort.midpoint.pn.txt' > $id'_experiment.sort.midpoint.pn.txt'
	bedtools closest -a $id'_Xu_2009_ORF-Ts_V64_TES.bed' -b $id'_experiment.sort.midpoint' -D a > $id'_Xu_2009_ORF-Ts_V64_TES.bound.txt'
	sort -k10,10n $id'_Xu_2009_ORF-Ts_V64_TES.bound.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6}' > $id'_Xu_2009_ORF-Ts_V64_TES.bound.sort.bed'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $id'_Xu_2009_ORF-Ts_V64_TES.bound.sort.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Xu_2009_ORF-Ts_V64_TES.bound.sort.bed.genegroup' -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	Rscript $script_bin'heatmap.R' $id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.tabular' $id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.png'
	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.TSS_ref.png' -a 0 -c 20 -f png

	#####################################
	##################
	### nucleosome part
	echo 'bound genes'
	python $script_bin'extract_bound.py' -b $id'_Xu_2009_ORF-Ts_V64.genelist.txt' -t /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/Yeast_plus_one_sacCer3.gff -o $id'_Yeast_plus_one_sacCer3.bed'
	echo 'sort by dist to TSS'
	bedtools closest -a $id'_Yeast_plus_one_sacCer3.bed' -b $id'_experiment.sort.midpoint' -D a > $id'_Yeast_plus_one_sacCer3.bound.txt'
	sort -k10,10nr $id'_Yeast_plus_one_sacCer3.bound.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-750,$2+251,$4,$5,$6; else print $1,$2-250,$2+751,$4,$5,$6}' > $id'_Yeast_plus_one_sacCer3.bound.sort.bed'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $id'_Yeast_plus_one_sacCer3.bound.sort.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Yeast_plus_one_sacCer3.bound.sort.bed.genegroup' -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	Rscript $script_bin'heatmap.R' $id'_Yeast_plus_one_sacCer3_'$bam_file'_read1_combined.tabular' $id'_Yeast_plus_one_sacCer3_'$bam_file'_read1_combined.png'
	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.TSS_ref.png' -a 0 -c 20 -f png

	cd ..
done < id_list.txt



##################
### Motif part
while read info; do
	id=$(echo $info | awk -F ' ' -v OFS='\t' '{print $1}' )
	motif_file=$(echo $info | awk -F ' ' -v OFS='\t' '{print $2}' )
	bam_file=$(echo $info | awk -F ' ' -v OFS='\t' '{print $3}' )
	echo $bam_file'.bam'
	echo $motif_file'.locations'
	cp $motif_file'.locations' $id'_gene_TSS_TES'

	cd $id'_gene_TSS_TES'
	#####################################
	tail -n+2 $motif_file'.locations' | awk -F '\t' -v OFS='\t' '{print $1}' | awk -F ':' -v OFS='\t' '{if ($3=="+") print $1,$2,$3; else print $1,$2,"_"}' | awk -F '-' -v OFS='\t' '{print $1,$2}' | awk -F '\t' -v OFS='\t' '{if ($4=="+") print $1,int(($2+$3)/2),int(($2+$3)/2)+1,".",".",$4; else print $1,int(($2+$3)/2),int(($2+$3)/2)+1,".",".","-"}' | sort -k1,1 -k2,2n > $motif_file'.bed'
	echo 'get bound genes'
	bedtools intersect -a Xu_2009_ORF-Ts_V64_cosort.bed -b $motif_file'.bed' -wa | sort -k1,1 -k2,2n -u > $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed'
	python $script_bin'bed_expand_sort.py' -f $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed' -u -500 -d -100 -r F -o $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed'
	echo 'sort by dist to TSS'
	cat $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.TSS.bed'
	bedtools closest -a $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.TSS.bed' -b $motif_file'.bed' -D a > $motif_file'_Xu_2009_ORF-Ts_V64_TSS.txt'
	sort -k13,13nr $motif_file'_Xu_2009_ORF-Ts_V64_TSS.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-750,$2+251,$4,$5,$6; else print $1,$2-250,$2+751,$4,$5,$6}' > $motif_file'_Xu_2009_ORF-Ts_V64_TSS.sort.bed'
	sort -k13,13nr $motif_file'_Xu_2009_ORF-Ts_V64_TSS.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-750,$2+251,$4,$5,$6; else print $1,$2-250,$2+751,$4,$5,$6,$13}' > $motif_file'_Xu_2009_ORF-Ts_V64_TSS.sort.bed.tmp'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $motif_file'_Xu_2009_ORF-Ts_V64_TSS.sort.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $motif_file'_Xu_2009_ORF-Ts_V64_TSS.sort.bed.genegroup' -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	Rscript $script_bin'heatmap.R' $motif_file'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' $motif_file'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.png'
	### TES
	echo 'sort by dist to TES'
	cat $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$3-1,$3,$4,$5,$6; else print $1,$2,$2+1,$4,$5,$6}' | sort -k1,1 -k2,2n > $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.TES.bed'
	bedtools closest -a $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.TES.bed' -b $motif_file'.bed' -D a > $motif_file'_Xu_2009_ORF-Ts_V64_TES.txt'
	sort -k13,13n $motif_file'_Xu_2009_ORF-Ts_V64_TES.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6}' > $motif_file'_Xu_2009_ORF-Ts_V64_TES.sort.bed'
	sort -k13,13n $motif_file'_Xu_2009_ORF-Ts_V64_TES.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+501,$4,$5,$6; else print $1,$2-500,$2+501,$4,$5,$6,$13}' > $motif_file'_Xu_2009_ORF-Ts_V64_TES.sort.bed.tmp'
	echo 'gene group split'
	python $script_bin'gene_group_split.py' -t $motif_file'_Xu_2009_ORF-Ts_V64_TES.sort.bed' -a 3 -g $script_bin'gene_group_list.txt' -b 0 -i F
	echo 'get heatmap table'
	java -jar $script_bin'TagPileup.jar' -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $motif_file'_Xu_2009_ORF-Ts_V64_TES.sort.bed.genegroup' -s 6 -n 1 -e true -r 0 -p false -a 1 -t 4 -w 0 -h true -m false 
	echo 'generate heatmap'
	Rscript $script_bin'heatmap.R' $motif_file'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.tabular' $motif_file'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.png'
	cd ..
	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$motif_file'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $motif_file'_Xu_2009_ORF-Ts_V64_TSS.png' -a 0 -c 20 -f png
done < id_list_motif.txt





