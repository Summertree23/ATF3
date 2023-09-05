Rscript filter.counts.R sample.peak.remove.chry.chrun.random.xls 
Rscript cpm.normalized.R peak.filter.counts.xls comparison.xls 
Rscript atac.pca.right.order.R cpm_normalized.counts.rds comparison.xls
Rscript pearson.R cpm_normalized.counts.xls  pearson.pdf

for i in `cat comparison.group`;do Rscript deseq2.comparison.high.R peak.filter.counts.rds comparison.xls $i & done
for i in `cat comparison.group`;do Rscript deseq2.comparison.low.R peak.filter.counts.rds comparison.xls $i & done

for i in `ls *_*.xls|grep -v all|sed 's/.xls//g'`;do less $i.xls|cut -f 1|sed '1d'|cut -f 1|sed '1d'|sed 's/:/\t/g'|sed 's/-/\t/g'|grep -v "v1"|grep -v "v2" > $i.bed & done
for i in `ls *_*.bed`;do Rscript target_gene.100000tss.R $i & done 
for i in `ls *annotation.txt`;do le $i |cut -f 7 |sed '1d'|awk -F" " '{print $1}' > $i.line7 & done 
for i in `ls *.line7`;do Rscript /disk1/xilu/collaborate/niklas/4.differential_peaks/barplot2.R  $i ${i}.file & done
cat *file |sort|uniq|sed 's/.bedannotation.txt.line7//g'|grep -v region| sed '1iregion\tcluster\tnum' > accum.txt
Rscript accum.R

#TF enrichment
for i in `ls *_*.bed|sed 's/.bed//g'`;do findMotifsGenome.pl ${i}.bed hg38 ${i}_motif/ -nomotif -preparsedDir /disk1/xilu/glioblastoma/cell_line_51/compared_1_2_3 > $i.final.o & done

for i in `ls */knownResults.txt|sed 's/.txt$//g'`
do 
    cat ${i}.txt |cut -f 1,4,5,7 > ${i}.sort.txt & 
done
wait

for i in `ls */knownResults.txt|sed 's/.txt$//g'`
do
    cat ${i}.sort.txt | cut -f 1|awk -F"/" '{print $1}'|paste - ${i}.sort.txt |cut -f 1,3,4,5|sed 's/%//g' > ${i}.final.txt
done
wait

for i in `ls */knownResults.txt|sed 's/.txt$//g'`
do
    Rscript /disk1/xilu/pipeline/time_Series_analysis/tf.homer.rank.R  ${i}.final.txt ${i}.tf.rank.homer.pdf
done
wait
Rscript volcano.R  control_drug_all.xls  control_drug_all.volcano.pdf

perl add.peak.pl representative.gene.peak control_drug_all_2.xls > control_drug_all_3.xls
