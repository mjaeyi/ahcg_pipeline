##AHCG Pipeline Readme
Min Yi

  
````{sh}
********************************************************************************************************************************
*                                                                                                                              *
*                                         Download and install a Virtual Box                                                   *
*                                                                                                                              *
********************************************************************************************************************************
-https://www.virtualbox.org/wiki/Downloads



********************************************************************************************************************************
*                                                                                                                              *
*                                         Download and set up Basespace server                                                 *
*                                                                                                                              *
********************************************************************************************************************************
-https://da1s119xsxmu0.cloudfront.net/sites/developer/native/nativeappsvm/BaseSpace%20Native%20App%20VM%20(phix%20only)%20v9.ova
-Details on installing the server is found here:
   -https://developer.basespace.illumina.com/docs/content/documentation/native-apps/setup-dev-environment

-Once starting the server through Virtual Box, you can access the server through a shell (i.e. Putty)
  -Hostname: basespace@localhost
  -Password: basespace
  -Port: 2222



********************************************************************************************************************************
*                                                                                                                              *
*                                         Set up and run AHCG pipeline                                                         *
*                                                                                                                              *
********************************************************************************************************************************
-Clone the respository: https://github.com/shashidhar22/ahcg_pipeline
   -use command: git clone https://github.com/shashidhar22/ahcg_pipeline.git

-The requirements are:
   -Python3 - version 3.4.1
   -Trimmomatic - version 0.36
   -Bowtie2 - version 2.2.9
   -Picard - version 2.6.0
   -GATK - version 3.4
    **This can be manually installed or pulled from the above repository using command: git pull origin master**

-Download the reference genome
  -wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz
     **This resource file contains hg19 and the dbsnp vcf for hg19**
  -Extract via: tar -xvzf www.prism.gatech.edu/~sravishankar9/resources.tar.gz

-Build a BowTie Index
  -bowtie2-build -f /path/to/hg19.fa hg19
   -The above command designates that the reference is a fasta file with option "-f" and named as hg19

-To make installs of tools and dependencies easier, use Linuxbrew
  -ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"
  -echo 'export PATH="$HOME/.linuxbrew/bin:$PATH"' >>~/.bash_profile
   -Now use "brew install *toolname*" in order to easily install tools

-Download and install samtools
 -Manually
   -wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
   -cd .../samtools-1.3.1    # Within the unpacked release directory
   -./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.3.1
   -make all plugins-htslib
 -Using sudo
   -sudo apt-get install samtools

-Use samtools to build an index
  -samtools faidx /path/to/hg19.fa

-Install Java
 -Using sudo
   -sudo add-apt-repository ppa:webupd8team/java
   -sudo apt-get update
   -sudo apt-get install openjdk-8-jre
    **To run picard.jar you need at least java 8**

-Use Picard.jar to create a genome dictionary
  -java -jar /path/to/picard.jar CreateSequenceDictionary R=/path/to/hg19 O=hg19.dict 

-Download and prepare test data set
  -wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
  -wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
  -gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
  -gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
  -head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
  -head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq

-Run the ahcg_pipeline with the test data set
  -./ahcg_pipeline -h   **Shows a helpme**

   usage: VariantCaller [-h] [-t TRIM_PATH] [-b BOWTIE_PATH] [-p PICARD_PATH]
                        [-g GATK_PATH] [-i INPUT_PATH [INPUT_PATH ...]]
                        [-w INDEX_PATH] [-d DBSNP_PATH] [-r REF_PATH]
                        [-a ADAPTER_PATH] [-o OUT_PATH]

   optional arguments:
     -h, --help            show this help message and exit
     -t TRIM_PATH, --trimmomatic TRIM_PATH
                           Path to Trimmomatic
     -b BOWTIE_PATH, --bowtie BOWTIE_PATH
                           Path to Bowtie
     -p PICARD_PATH, --picard PICARD_PATH
                           Path to Picard
     -g GATK_PATH, --gatk GATK_PATH
                           Path to GATK
     -i INPUT_PATH [INPUT_PATH ...], --inputs INPUT_PATH [INPUT_PATH ...]
                           Path to Paired end reads
     -w INDEX_PATH, --index INDEX_PATH
                           Path to Reference bowtie index
     -d DBSNP_PATH, --dbsnp DBSNP_PATH
                           Path to dbSNP vcf file
     -r REF_PATH, --refrence REF_PATH
                           Path to Reference file
     -a ADAPTER_PATH, --adapter ADAPTER_PATH
                           Path to Adapter file
     -o OUT_PATH, --outpath OUT_PATH
                           Path to Ouput directory



********************************************************************************************************************************
*                                                                                                                              *
*                                         Setting up and updating GitHub                                                       *
*                                                                                                                              *
********************************************************************************************************************************
-First go and fork the original repository to personal github
-Then update the .git/config file and add personal github information
-Update .gitignore and place all folders that do not need to be pushed onto github
-Use command: git add *
  **This will add any updates and ready them to commit
-Set up global user.email and user.name (This is a one time step)
  -git config --global user.email "Email"
  -git config --global user.name "Name"
-Then commit the changes
  -git commit -m "Description of changes"
-Then push the changes onto github
  -git push origin master



********************************************************************************************************************************
*                                                                                                                              *
*                                         Finding sequence of exome using bedtools                                             *
*                                                                                                                              *
********************************************************************************************************************************
-Download the gene annotation file
  -wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
-Locate the BRCA1 genes
  -grep BRCA1 hg19_refGene.txt
-6 BRCA1 variants are found. Determine which variant is quantitatively the best
  -Search BRCA1 at https://dnasu.org/DNASU/Home.do
  -Select reference: HsCD00022337
  -Click on "Reference Sequence Alignment"
  -NM_007294 matched the reference alignment the best and was determined to be the best variant
-Extract the exome coordinates from hg19_refGene.txt
  -grep NM_007294 hg19_refGene.txt | cat > brca.txt
-Use bedconverter.py to convert the text file into a bed file
  -Usage: ./bedconverter.py -i brca.txt -o brca1.bed
-Install bedtools
  -brew install bedtools
-Use bedtools to get the sequences from the reference set
  -bedtools getfasta -s -fo results.fa -fi /path/to/hg19.fa -bed brca1.bed


********************************************************************************************************************************
*                                                                                                                              *
*                                         Extracting reads mapping to a region of interest                                     *
*                                                                                                                              *
********************************************************************************************************************************
-Dowload the NA12878 HiSeq Exome dataset from GIAB FTP site
  -ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/
   **There are 4 bam files you need to download**
-Use Samtools to subset the bam file to regions of interest (BRCA1)
  -samtools view -L <bed file> -b -o <output bam file> <input bam file>
-Use Bedtools to convert the bam file to a fastq file
  -bedtools bamtofastq -i <bam file> -fq <fastq r1> -fq2 <fastq r2>

********************************************************************************************************************************
*                                                                                                                              *
*                                         Variant Quality Score Recalibration (VQSR)                                           *
*                                                                                                                              *
********************************************************************************************************************************
java -jar GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R ~/ref/resources/genome/hg19.fa \
    -input ./variants.vcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./lib/VQSR/hapmap_3.3.hg19.sites.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 ./lib/VQSR/1000G_omni2.5.hg19.sites.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 ./lib/VQSR/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ~/ref/resources/dbsnp_138.hg19.vcf \
    -an DP \
    -an QD \
    -an FS \
    -an SOR \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches


java -jar GenomeAnalysisTK.jar \ 
    -T ApplyRecalibration \ 
    -R ~/ref/resources/genome/hg19.fa \ 
    -input ./variants.vcf \ 
    -mode SNP \ 
    --ts_filter_level 99.0 \ 
    -recalFile recalibrate_SNP.recal \ 
    -tranchesFile recalibrate_SNP.tranches \ 
    -o recalibrated_snps_raw_indels.vcf

********************************************************************************************************************************
*                                                                                                                              *
*                                                           DCM PIPELINE                                                       *
*                                                                                                                              *
********************************************************************************************************************************
To run this DCM pipeline, you must set up your directory as follows:

dcm/
	|---VQSR
		|---Omni.vcf
		|---Hapmap.vcf
		|---1000G.vcf
	|---cov.py
	|---draw_depth.R
	|---dcm_pipeline.sh
	|---parse_clnsig.py
ref/
	|---resources
		|---dbsnp
			|---dbsnp_138.hg19.vcf
		|---genome
			|---hg19.fa
ahcg_pipeline/
	|---lib
		|---GenomeAnalysisTK.jar

You must also these dependency files installed:
	--Bedtools
	--Samtools
	--Imagemagick
	--R
		--Rmodule "ggplot2"
	--Python 2.7+
	--Python 3+
		--Pymodule "click"
		--Pymodule "PyVCF"

The instructions to run the file:
	1) Be in the dcm/ directory
	2) bash dcm_pipeline.sh /name/of/output/dir/ /path/to/patient.bam /path/to/dcm_genelist.bed /path/to/clinvar.vcf
	3) Final patient report is generated in dcm/patient_report.pdf
