#!/bin/bash

if [ "$#" -ne 5 ]
then
    echo "usage: $0 output_dir read1.fq read2.fq gene_list.bed clinvar.vcf.gz"
    exit 1
fi

BASE_DIR="/home/basespace"
TRIM=$BASE_DIR/ahcg_pipeline/lib/Trimmomatic-0.36/trimmomatic-0.36.jar
ADAPTER=$BASE_DIR/ahcg_pipeline/lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
BOWTIE=$BASE_DIR/ahcg_pipeline/lib/bowtie2-2.2.9/bowtie2
PICARD=$BASE_DIR/ahcg_pipeline/lib/picard.jar
GATK=$BASE_DIR/ahcg_pipeline/lib/GenomeAnalysisTK.jar

echo "Running Variant Pipeline..."

python3 $BASE_DIR/ahcg_pipeline/ahcg_pipeline.py \
    -t $TRIM \
    -b $BOWTIE \
    -p $PICARD \
    -g $GATK \
    -i $2 $3 \
    -w $BASE_DIR/ahcg_pipeline/bowtieIndex/hg19 \
    -d $BASE_DIR/ref/resources/dbsnp/dbsnp_138.hg19.vcf.gz \
    -r $BASE_DIR/ref/resources/genome/hg19.fa \
    -a $ADAPTER \
    -o $BASE_DIR/dcm/$1

echo "Finished Running... Now Recalibrating VCF"

BASE=$BASE_DIR
REF=$BASE/ref/resources/genome/hg19.fa
VARIANTS="$1/variants.vcf"
RECAL=$(dirname $VARIANTS)/output.recal
TRANCH=$(dirname $VARIANTS)/output.tranches
HAPMAP=$BASE/dcm/VQSR/hapmap_3.3.hg19.sites.vcf.gz
OMNI=$BASE/dcm/VQSR/1000G_omni2.5.hg19.sites.vcf.gz
PHASE=$BASE/dcm/VQSR/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
DBSNP=$BASE/ref/resources/dbsnp/dbsnp_138.hg19.vcf
java -Xmx4g -jar $GATK \
	-T VariantRecalibrator \
        -R $REF \
        -input $VARIANTS \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
        -resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
        -an DP \
        -an QD \
        -an FS \
        -an SOR \
        -an MQ \
        -an MQRankSum \
        -an ReadPosRankSum \
        -mode SNP \
        -recalFile $RECAL \
        -tranchesFile $TRANCH


suffix=_recal.vcf
outname="$(basename $VARIANTS | sed -e 's/\.vcf$//')$suffix"
RECALs=$(dirname $VARIANTS)/output.recal
TRANCHs=$(dirname $VARIANTS)/output.tranches
java -jar $GATK \
    -T ApplyRecalibration \
    -R $REF \
    -input $VARIANTS \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile $RECALs \
    --tranches_file $TRANCHs \
    -o $1/$outname

echo "Recalibration finished. Cross referencing to determine clinically related variants"

bedtools intersect -a $5 -b $4 -header > $1/clinvar_allfrombed.vcf
bedtools intersect -a $1/$outname -b $4 -header > $1/patient_dcm_final.vcf
bedtools intersect -b $1/patient_dcm_final.vcf -a $1/clinvar_allfrombed.vcf -header > $1/patient_intersect_clinvar.vcf

python3 parse_clnsig.py -i $1/patient_intersect_clinvar.vcf 2>&1 | tee $1/patient_simple_report.txt

cut -c 24- $1/patient_simple_report.txt

echo "Variants cross reference finished. Now performing steps to produce coverage plots..."

out="$(basename $VARIANTS | sed -e 's/\.bam$//').bga.bed"
final="$(basename $VARIANTS | sed -e 's/\.bam$//').join_final.bed"
depths="$(basename $VARIANTS | sed -e 's/\.bam$//').depths.bed"
bam="$(basename $VARIANTS | sed -e 's/\.fq$//')_final.bam"

samtools view -L $4 $1/$bam -b > $1/subset.bam
bedtools genomecov -ibam $1/subset.bam -bga > $1/$out
bedtools intersect -loj -F 0.10 -a $4 -b $1/$out -bed > $1/$final

#echo "awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$6,$7,$4,$8)}' $1/$final > $1/$depths"

for gene in $(cut -f4 $4 | sort -u | xargs)
	do
		grep $gene $1/$depths > $1/${gene}_raw.txt
		python cov.py $1/${gene}_raw.txt $1/${gene}.txt
		xvfb-run --server-args="-screen 0 1024x768x24" ./draw_depth.R $1/${gene}.txt
	done

convert $1/patient_simple_report.txt $1/*.png patient_report.pdf 
