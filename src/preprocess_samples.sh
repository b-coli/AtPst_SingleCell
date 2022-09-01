#!/bin/bash

## Define Metadata indexes in file metadata csv
sn_idx=$(awk -v RS=',' '/Sample_Name/{print NR; exit}' metadata/sample_metadata.csv) ## Sample Name
seqtype_idx=$(awk -v RS=',' '/Sequence_Type/{print NR; exit}' metadata/sample_metadata.csv) ## Read (if paired)

## Load metadata into a series of arrays
samplenames=($(awk -F, -v col="$sn_idx" '{if(NR > 1) print $col}' metadata/sample_metadata.csv))
seqtypes=($(awk -F, -v col="$seqtype_idx" '{if(NR > 1) print $col}' metadata/sample_metadata.csv))

## Count the number of filesets
len=${#samplenames[@]}
i=0

## Loop over all file sets
while [ "$i" -lt "$len" ]; do
	## For each file set, define attributes
	samplename=${samplenames[$i]}
	seqtype=${seqtypes[$i]}

	if [ $seqtype == "RNAseq" ]; then
		./src/fastqc_trim_align_bulk.sh ${samplename}
	fi

	if [ $seqtype == "scRNAseq" ]; then
		./src/run_cellranger.sh ${samplename}
		./src/run_velocyto.sh ${samplename}
	fi
        i=$(($i + 1))
done
