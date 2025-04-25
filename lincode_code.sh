# Command used in Chapter 3 for stand-alone bioinformatic tools: core genome LIN code for H. influenzae

## 1) MSTclust: https://gitlab.pasteur.fr/GIPhy/MSTclust
## Program was run locally and through the University of Oxford Advanced Research Computing (ARC) computer cluster

PROGDIR="/dir/to/MSTclust/src"
infile="/dir/to/cgMLST_provile.tsv"

java -jar ${PROGDIR}/MSTclust -i ${infile} -o data  -l 1  -p 2- -e 0.024  -h  -m # maximum number missing core genome alleles in PubMLST is 25, therefore 25/1037 = 0.024

grep "data" data.graphml | awk -F '[<>]' '{print$3}' | sort -g | uniq | # searching for optimal clustering cut-off
	while read c ; do java -jar ${PROGDIR}/MSTclust.jar -i data.d -o optimal_out_${c} -c \
	$c -L 1037 -B 9 -t ; done 2>/dev/null

### choosing optimal cut-offs, stored in cutoff_arbitrary.txt

for f in $(cat cutoff_arbitrary.txt); # data perturbation and subsampling analyses
	do java -jar ${PROGDIR}/MSTclust.jar -i data.d -o clust_${f} -c $f -L 1037 -B 99 -R 1000 ;
done >> clust_replicate.out

### --- ###


## 2) pHierCC: https://github.com/zheminzhou/pHierCC
## Program was run locally and through the University of Oxford ARC computer cluster, in a personal conda environment

infile="/dir/to/cgMLST_provile.txt.gz"

pHierCC -p ${infile} -o output
HCCeval -p ${infile} -c output.HierCC.gz -o output_eval

### --- ###


## 3) FastANI: https://github.com/ParBLiSS/FastANI
## Program was run through the University of Oxford ARC computer cluster, in a personal conda environment

fastANI --ql q_many_hinfs.txt --rl ref_many_hinfs.txt -o hinf_to_hinfs_core.txt --matrix -t 16
	# q_many_hinfs.txt and ref_many_hinfs.txt each listed complete path to all fasta files (1 path per line)

### --- ###


## 4) RAxML: https://github.com/stamatak/standard-RAxML, and
## 5) ClonalFrameML: https://github.com/xavierdidelot/clonalframeml/wiki
## Program was run through the University of Oxford ARC computer cluster

module load RAxML/8.2.12-gompi-2020a-hybrid-avx2
module load ClonalFrameML/1.12-foss-2022a

infile="/dir/to/core_alignment.fasta"

raxmlHPC -m GTRGAMMA -p 619283 -T 48 -N 5 -s ${infile} -n lin_dataset1_core.nwk
ClonalFrameML /dir/to/RAxML_output/RAxML_bestTree.lin_dataset1_core.nwk ${infile} dataset1_CFML