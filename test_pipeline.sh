for id in $(cat id_list.txt)
do
#####################################
echo 'get multiGPS pks'
tail -n+8 11766_experiment.events | awk -F '\t' -v OFS='\t' '{print $1}' | awk -F ':' -v OFS='\t' '{print $1,$2,$2+1}' | sort -k1,1 -k2,2n > 11766_experiment.sort.midpoint

echo 'expand bed file and sort'
python gff3_expand_sort.py -f Xu_2009_ORF-Ts_V64.gff3 -u 500 -d 100 -r T -o Xu_2009_ORF-Ts_V64_sort.bed
sort -k1,1 -k2,2n Xu_2009_ORF-Ts_V64_sort.bed > Xu_2009_ORF-Ts_V64_cosort.bed

echo 'get bound genes'
bedtools intersect -a Xu_2009_ORF-Ts_V64_cosort.bed -b 11766_experiment.sort.midpoint -wa | sort -u > Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed
python bed_expand_sort.py -f Xu_2009_ORF-Ts_V64.cosort.bound.uniq.bed -u -500 -d -100 -r F -o Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed

echo 'expand bed file and sort'
python bed_expandmidpoint_sort.py -f Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed -u 2000 -d 2000 -r F -o Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed
#echo 'generate bam.bai'
samtools index 11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216.bam

echo 'get heatmap table'
java -jar TagPileup.jar -b 11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216.bam -i 11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216.bam.bai -c Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.bed -s 6 -n 5 -e true -r 0 -p true -a 1 -t 4 -w 0 -h true -m true -o ./ -x Xu_2009_ORF-Ts_V64.cosort.bound.uniq.expand.composite.txt

echo 'generate heatmap'
java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./Xu_2009_ORF-Ts_V64_11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216_read1_combined.tabular -x Dendrogram -- -o gene_len.png -a 0 -c 20 -f png

java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./11766_Xu_2009_ORF-Ts_V64_11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216_read1_combined.tabular -x Dendrogram -- -o gene_len.png -a 0 -c 20 -f png

11766_Xu_2009_ORF-Ts_V64_11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216_read1_combined.tabular
#####################################
echo 'expand bed file and sort'
cat Xu_2009_ORF-Ts_V64.cosort.bound.uniq.noexpand.bed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > Xu_2009_ORF-Ts_V64_TSS.bed

echo 'sort by dist to TSS'
bedtools closest -a Xu_2009_ORF-Ts_V64_TSS.bed -b 11766_experiment.sort.midpoint -D a > Xu_2009_ORF-Ts_V64_TSS.bound.txt
sort -k10,10nr Xu_2009_ORF-Ts_V64_TSS.bound.txt | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+100,$4,$5,$6; else print $1,$2-100,$2+500,$4,$5,$6}' > Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed

echo 'get heatmap table'
java -jar TagPileup.jar -b 11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216.bam -i 11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216.bam.bai -c Xu_2009_ORF-Ts_V64_TSS.bound.sort.bed -s 6 -n 1 -e true -r 0 -p true -a 1 -t 4 -w 0 -h true -m true -o ./ -x Xu_2009_ORF-Ts_V64_TSS.bound.sort.composite.txt

echo 'generate heatmap'
java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./Xu_2009_ORF-Ts_V64_TSS_11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216_read1_combined.tabular -x Dendrogram -- -o TSS_ref.png -a 0 -c 20 -f png


#####################################
tail -n+2 transcripts_11766_Motif_0.locations | awk -F '\t' -v OFS='\t' '{print $1}' | awk -F ':' -v OFS='\t' '{if ($3=="+") print $1,$2,$3; else print $1,$2,"_"}' | awk -F '-' -v OFS='\t' '{print $1,$2}' | awk -F '\t' -v OFS='\t' '{if ($4=="+") print $1,$2,$3,".",".",$4; else print $1,$2,$3,".",".","-"}' | sort -k1,1 -k2,2n > transcripts_11766_Motif_0.bed

echo 'sort by dist to TSS'
bedtools closest -a Xu_2009_ORF-Ts_V64_TSS.bed -b transcripts_11766_Motif_0.bed -D a -s > Xu_2009_ORF-Ts_V64_TSS.transcripts_11766_Motif_0.txt
sort -k10,10nr Xu_2009_ORF-Ts_V64_TSS.transcripts_11766_Motif_0.txt | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2-500,$2+100,$4,$5,$6; else print $1,$2-100,$2+500,$4,$5,$6}' > Xu_2009_ORF-Ts_V64_TSS_transcripts_11766_Motif_0.sort.bed

echo 'get heatmap table'
java -jar TagPileup.jar -b 11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216.bam -i 11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216.bam.bai -c Xu_2009_ORF-Ts_V64_TSS_transcripts_11766_Motif_0.sort.bed -s 6 -n 1 -e true -r 0 -p true -a 1 -t 4 -w 0 -h true -m true -o ./ -x Xu_2009_ORF-Ts_V64_TSS_transcripts_11766_Motif_0.sort.composite.txt

echo 'generate heatmap'
java -jar /Volumes/MAC_Data/data/labs/pugh_lab/test_pipeline/TreeView-1.1.6r4-bin/TreeView.jar -r ./Xu_2009_ORF-Ts_V64_TSS_transcripts_11766_Motif_0_11766_Reb1_PughLab_63329_i5006_BY4741_None_YPD__XO_216_read1_combined.tabular -x Dendrogram -- -o TSS_transcripts_11766_Motif_0.png -a 0 -c 20 -f png


done
