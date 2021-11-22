#!/bin/bash

./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene ./annovar/humandb/

echo "Downloading revel from annovar..."
./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar revel ./data/revel/
echo "Downloading revel done."

echo "Downloading cadd from annovar..."
./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cadd13 ./data/cadd/
echo "Downloading cadd done."

echo "Downloading dbnsfp from annovar..."
./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp41a ./data/dbnsfp/
echo "Downloading dbnsfp done."

echo "Downloading mcap from annovar..."
./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar mcap13 ./data/mcap/
echo "Downloading mcap done."

cd data
mkdir bigwig
echo "Downloading bigwig files"
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-DNase.pval.signal.bigwig ./data/bigwig
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K4me1.pval.signal.bigwig ./data/bigwig
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K4me3.pval.signal.bigwig ./data/bigwig
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K9me3.pval.signal.bigwig ./data/bigwig
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K27ac.pval.signal.bigwig ./data/bigwig
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K27me3.pval.signal.bigwig ./data/bigwig
wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E003-H3K36me3.pval.signal.bigwig ./data/bigwig
wget https://epigenomesportal.ca/tracks/Roadmap/hg19/28736.Roadmap.SRS004212_Combined_Libraries_424.WGB-Seq.signal.bigWig ./data/bigwig
wget https://www.encodeproject.org/files/ENCFF225MAO/@@download/ENCFF225MAO.bigWig ./data/bigwig
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/hg19.100way.phyloP100way.bw ./data/bigwig
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM923nnn/GSM923451/suppl/GSM923451_hg19_wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig ./data/bigwig
wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/E003_WGBS_FractionalMethylation.bigwig
echo "Downloading bigwig files done."


