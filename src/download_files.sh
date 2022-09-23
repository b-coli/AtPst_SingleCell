#!/bin/bash

startfile=20

## Define Metadata indexes in file metadata csv
sn_idx=$(awk -v RS=',' '/Sample_Name/{print NR; exit}' metadata/file_metadata.csv) ## Sample Name
acc_idx=$(awk -v RS=',' '/Accession/{print NR; exit}' metadata/file_metadata.csv) ## Accession ID
lane_idx=$(awk -v RS=',' '/Lane/{print NR; exit}' metadata/file_metadata.csv) ## Sequencing lane (technical)
smpl_idx=$(awk -v RS=',' '/Sample$/{print NR; exit}' metadata/file_metadata.csv) ## Sequencing sample ID (technical)
db_idx=$(awk -v RS=',' '/Database/{print NR; exit}' metadata/file_metadata.csv) ## Remote (or local) database
loc_idx=$(awk -v RS=',' '/Location/{print NR; exit}' metadata/file_metadata.csv) ## If file is local, location on filesystem
read_idx=$(awk -v RS=',' '/Read/{print NR; exit}' metadata/file_metadata.csv) ## Read (if paired)
seqtype_idx=$(awk -v RS=',' '/Sequence_Type/{print NR; exit}' metadata/file_metadata.csv) ## Read (if paired)
gen_idx=$(awk -v RS=',' '/Genomic_Read_No/{print NR; exit}' metadata/file_metadata.csv) ## Read (if paired)
bar_idx=$(awk -v RS=',' '/Barcode_Read_No/{print NR; exit}' metadata/file_metadata.csv) ## Read (if paired)

## Load metadata into a series of arrays
accessions=($(awk -F, -v col="$acc_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
samplenames=($(awk -F, -v col="$sn_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
lanes=($(awk -F, -v col="$lane_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
smpls=($(awk -F, -v col="$smpl_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
dbs=($(awk -F, -v col="$db_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
locs=($(awk -F, -v col="$loc_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
rds=($(awk -F, -v col="$read_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
seqtypes=($(awk -F, -v col="$seqtype_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
gens=($(awk -F, -v col="$gen_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))
bars=($(awk -F, -v col="$bar_idx" -v st="$startfile" '{if(NR > st) print $col}' metadata/file_metadata.csv))


## Count the number of filesets
len=${#accessions[@]}
i=0

## Generate the fastq directory
mkdir -p data/fastq/

## Loop over all file sets
while [ "$i" -lt "$len" ]; do
	## For each file set, define attributes
	accession=${accessions[$i]}
	samplename=${samplenames[$i]}
	lane=$(printf '%03i' "${lanes[$i]}") ## Format lane using illumina standard
	smpl=$(printf '%02i' "${smpls[$i]}") ## Format sample id using illumina standard
	db=${dbs[$i]}
	loc=${locs[$i]}
	rd=${rds[$i]}
	seqtype=${seqtypes[$i]}
	gen=${gens[$i]}
	bar=${bars[$i]}


	filebase="${samplename}_S${smpl}_L${lane}" ## Define "file base", the illumina-standard prefix for the file

	r1="${filebase}_R1_001.fastq.gz" ## Define Read 1 filename
	r2="${filebase}_R2_001.fastq.gz" ## Define Read 2 filename

	## Generate sample folder for storing fastq
	sample_folder="data/fastq/${samplename}"
	mkdir -p $sample_folder

	## If file is local, copy over to newly-created folder
	if [ "$db" == "Local" ]; then
		if [ "$rd" == "1" ]; then
			echo "$samplename"
			#cp ${loc} ${sample_folder}/${r1}
		fi
		if [ "$rd" == "2" ]; then
			echo "$samplename"
			#cp ${loc} ${sample_folder}/${r2}
		fi
	fi

	## If file is located in NCBI SRA, use fasterq-dump to download, then rename and move
	if [ "$db" == "NCBI_SRA" ]; then
		#singularity run docker://ncbi/sra-tools fasterq-dump -e 15 -S --include-technical $accession

		if [ "$seqtype" == "PE" ]; then
			echo "$samplename"
			#pigz -p 15 ${accession}_${bar}.fastq 
			#pigz -p 15 ${accession}_${gen}.fastq
			#mv ${accession}_${bar}.fastq.gz ${sample_folder}/$r1
			#mv ${accession}_${gen}.fastq.gz ${sample_folder}/$r2
		fi

		if [ "$seqtype" == "SE" ]; then
			echo "$samplename"
			#pigz -p 15 ${accession}_1.fastq
			#mv ${accession}_1.fastq.gz ${sample_folder}/$r1
		fi
	fi

	if [ "$db" == "NGDC" ]; then
		if [ "$rd" == "1" ]; then
			echo "${samplename}"
			wget ${loc} -O ${sample_folder}/${r1}
		fi

		if [ "$rd" == "2" ]; then
			echo "${samplename}"
			wget ${loc} -O ${sample_folder}/${r2}
		fi

	fi

        i=$(($i + 1))
done

## Download the Zhang data set

## Download the published Frommer dataset, to serve as a reference:
mkdir -p data/expression/FrommerPublished/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161332/suppl/GSE161332_barcodes.tsv.gz -O data/expression/KimPublished/barcodes.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161332/suppl/GSE161332_features.tsv.gz -O data/expression/KimPublished/features.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161332/suppl/GSE161332_matrix.mtx.gz -O data/expression/KimPublished/matrix.mtx.gz
wget https://www.arabidopsis.org/download_files/Public_Data_Releases/TAIR_Data_20200930/gene_aliases_20200930.txt.gz -O data/gene_aliases_20200930.txt.gz
