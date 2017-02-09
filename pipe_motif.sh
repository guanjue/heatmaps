#!/bin/bash
while read info; do
	motif_file=$(echo $info | awk -F ' ' -v OFS='\t' '{print $1}' )
	bam_file=$(echo $info | awk -F ' ' -v OFS='\t' '{print $2}' )
	echo $bam_file'.bam'
	echo $motif_file'.locations'
	#####################################
	tail -n+2 $motif_file'.locations' | awk -F '\t' -v OFS='\t' '{print $1}' | awk -F ':' -v OFS='\t' '{if ($3=="+") print $1,$2,$3; else print $1,$2,"_"}' | awk -F '-' -v OFS='\t' '{print $1,$2}' | awk -F '\t' -v OFS='\t' '{if ($4=="+") print $1,int(($2+$3)/2),int(($2+$3)/2)+1,".",".",$4; else print $1,int(($2+$3)/2),int(($2+$3)/2)+1,".",".","-"}' | sort -k1,1 -k2,2n > $motif_file'.bed'
	echo 'expand bed file and sort'
	python gff3_expand_sort.py -f Xu_2009_ORF-Ts_V64.gff3 -u 500 -d 100 -r T -o Xu_2009_ORF-Ts_V64_sort.bed
	sort -k1,1 -k2,2n Xu_2009_ORF-Ts_V64_sort.bed > Xu_2009_ORF-Ts_V64_cosort.bed
	echo 'get bound genes'
	bedtools intersect -a Xu_2009_ORF-Ts_V64_cosort.bed -b $motif_file'.bed' -wa | sort -k1,1 -k2,2n -u > $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed'
	python bed_expand_sort.py -f $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed' -u -500 -d -100 -r F -o $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed'
	echo 'sort by dist to TSS'
	cat $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.TSS.bed'
	bedtools closest -a $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.TSS.bed' -b $motif_file'.bed' -D a > $motif_file'_Xu_2009_ORF-Ts_V64_TSS.txt'
	sort -k13,13nr $motif_file'_Xu_2009_ORF-Ts_V64_TSS.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+100,$4,$5,$6; else print $1,$2-100,$2+500,$4,$5,$6}' > $motif_file'_Xu_2009_ORF-Ts_V64_TSS.sort.bed'
	sort -k13,13nr $motif_file'_Xu_2009_ORF-Ts_V64_TSS.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+100,$4,$5,$6; else print $1,$2-100,$2+500,$4,$5,$6,$13}' > $motif_file'_Xu_2009_ORF-Ts_V64_TSS.sort.bed.tmp'
	echo 'get heatmap table'
	java -jar TagPileup.jar -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $motif_file'_Xu_2009_ORF-Ts_V64_TSS.sort.bed' -s 6 -n 1 -e true -r 0 -p true -a 1 -t 4 -w 0 -h true -m true -o ./ -x $motif_file'_Xu_2009_ORF-Ts_V64_TSS.sort.composite.txt'
	echo 'generate heatmap'
	Rscript heatmap.R $motif_file'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' $motif_file'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.png'
	### TES
	echo 'sort by dist to TES'
	cat $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$3-1,$3,$4,$5,$6; else print $1,$2,$2+1,$4,$5,$6}' | sort -k1,1 -k2,2n > $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.TES.bed'
	bedtools closest -a $motif_file'_Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.TES.bed' -b $motif_file'.bed' -D a > $motif_file'_Xu_2009_ORF-Ts_V64_TES.txt'
	sort -k13,13n $motif_file'_Xu_2009_ORF-Ts_V64_TES.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-100,$2+500,$4,$5,$6; else print $1,$2-500,$2+100,$4,$5,$6}' > $motif_file'_Xu_2009_ORF-Ts_V64_TES.sort.bed'
	sort -k13,13n $motif_file'_Xu_2009_ORF-Ts_V64_TES.txt' | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-100,$2+500,$4,$5,$6; else print $1,$2-500,$2+100,$4,$5,$6,$13}' > $motif_file'_Xu_2009_ORF-Ts_V64_TES.sort.bed.tmp'
	echo 'get heatmap table'
	java -jar TagPileup.jar -b $bam_file'.bam' -i $bam_file'.bam.bai' -c $motif_file'_Xu_2009_ORF-Ts_V64_TES.sort.bed' -s 6 -n 1 -e true -r 0 -p true -a 1 -t 4 -w 0 -h true -m true -o ./ -x $motif_file'_Xu_2009_ORF-Ts_V64_TES.sort.composite.txt'
	echo 'generate heatmap'
	Rscript heatmap.R $motif_file'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.tabular' $motif_file'_Xu_2009_ORF-Ts_V64_TES_'$bam_file'_read1_combined.png'

	#java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./$motif_file'_Xu_2009_ORF-Ts_V64_TSS_'$bam_file'_read1_combined.tabular' -x Dendrogram -- -o $motif_file'_Xu_2009_ORF-Ts_V64_TSS.png' -a 0 -c 20 -f png
done < id_list_motif.txt

