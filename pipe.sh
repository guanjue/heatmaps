#!/bin/bash
while read info; do
	id=$(echo $info | awk -F ' ' -v OFS='\t' '{print $1}' )
	bam_file=$(echo $info | awk -F ' ' -v OFS='\t' '{print $2}' )
	motif_file=$(echo $info | awk -F ' ' -v OFS='\t' '{print $3}' )
	echo $id'_experiment.events'
	echo $bam_file'.bam'
	echo $motif_file'.locations'
	#####################################
	echo 'get multiGPS pks'
	tail -n+8 $id'_experiment.events' | awk -F '\t' -v OFS='\t' '{print $1}' | awk -F ':' -v OFS='\t' '{print $1,$2,$2+1}' | sort -k1,1 -k2,2n > $id'_experiment.sort.midpoint'
	echo 'expand bed file and sort'
	python gff3_expand_sort.py -f Xu_2009_ORF-Ts_V64.gff3 -u 500 -d 100 -r T -o Xu_2009_ORF-Ts_V64_sort.bed
	sort -k1,1 -k2,2n Xu_2009_ORF-Ts_V64_sort.bed > Xu_2009_ORF-Ts_V64_cosort.bed
	echo 'get bound genes'
	bedtools intersect -a Xu_2009_ORF-Ts_V64_cosort.bed -b $id'_experiment.sort.midpoint' -wa | sort -u > $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed'
	python bed_expand_sort.py -f $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed' -u -500 -d -100 -r F -o $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed'
	echo extract bound genes
	cat $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{print $4}' > $id'_Xu_2009_ORF-Ts_V64.genelist.txt'
	echo 'expand bed file and sort'
	python bed_expandmidpoint_sort.py -f $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' -u 2000 -d 2000 -r F -o $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed'
	#echo 'generate bam.bai'
	samtools index $bam_file'.bam'
	echo 'get heatmap table'
	java -jar TagPileup.jar -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed' -s 6 -n 5 -e true -r 0 -p true -a 1 -t 4 -w 0 -h true -m true -o ./ -x $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.composite.txt'
	echo 'generate heatmap'
	java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.gene_len.png' -a 0 -c 20 -f png

	#####################################
	echo 'bound genes'
	cat $id'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > Xu_2009_ORF-Ts_V64_TSS.bed
	echo 'sort by dist to TSS'
	bedtools closest -a Xu_2009_ORF-Ts_V64_TSS.bed -b $id'_experiment.sort.midpoint' -D a > $id'_Xu_2009_ORF-Ts_V64_TSS.bound.txt'
	sort -k10,10nr $id'_Xu_2009_ORF-Ts_V64_TSS.bound.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+100,$4,$5,$6; else print $1,$2-100,$2+500,$4,$5,$6}' > $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed'
	echo 'get heatmap table'
	java -jar TagPileup.jar -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed' -s 6 -n 1 -e true -r 0 -p true -a 1 -t 4 -w 0 -h true -m true -o ./ -x $id'_Xu_2009_ORF-Ts_V64_TSS.bound.sort.composite.txt'
	echo 'generate heatmap'
	java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$id'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $id'_Xu_2009_ORF-Ts_V64.TSS_ref.png' -a 0 -c 20 -f png

done < id_list.txt
