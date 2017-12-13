#!/bin/bash

# Set default values
if [ ! -z "$CONDA_PREFIX" ] ; then
  frogs_dir="$CONDA_PREFIX"
  frogs_test=$frogs_dir/share/FROGS-*/test
else
  frogs_dir=".."
  frogs_test=$frogs_test/test
fi
nb_cpu=1
java_mem=4
out_dir=/tmp
delete_out_dir=false
run_pynast=false

# FUNCTION: display help message
function help(){
  printf "\n$0: a script to test FROGS distribution with sample data set. \n"
  printf "\n Optional arguments are:"
  printf "\n  -c <cores>         : nb of cores to use (default: 1)"
  printf "\n  -m <java_mem>      : nb of Gb to allocate to Java Runtime (default: 4)"
  printf "\n  -o <output-dir>    : working directory to use during testing (default: /tmp/)"
  printf "\n  -f <FROGS-home_dir>: home directory of FROGS (default: ".." if standalone; $CONDA_PREFIX if Conda environment)"
  printf "\n  -p                 : run Pynast step (default: false)"
  printf "\n  -d                 : delete output-dir when tests are completed (default: false)"
  printf "\n  -h                 : print this message" 
  printf "\n"
  exit 0
}

# Prepare arguments for processing
while getopts hdpc:m:o:f: opt
do
    case "$opt" in
      c)  nb_cpu="$OPTARG";;
      m)  java_mem="$OPTARG";;
      o)  out_dir="$OPTARG";;
      f)  frogs_dir="$OPTARG";;
      d)  delete_out_dir=true;;
      p)  run_pynast=true;;
      h)  help;;
      \?)	help;;
    esac
done
shift `expr $OPTIND - 1`

out_dir=${out_dir}/frogs-test-result

echo "Starting FROGS tests with:"
echo "  Nb cores         : $nb_cpu"
echo "  Java mem         : $java_mem"
echo "  Output dir       : $out_dir"
echo "  FROGS dir        : $frogs_dir"
echo "  FROGS test       : $frogs_test"
echo "  Run pynast       : $run_pynast"
echo "  Delete output dir: $delete_out_dir (will be done only if all tests are ok)"
echo ""

# Set ENV
if [ ! -v $CONDA_PREFIX ] ; then
  export PATH=$frogs_dir/libexec:$frogs_dir/app:$PATH
fi
export PATH=$frogs_dir/libexec:$frogs_dir/app:$PATH
export PYTHONPATH=$frogs_dir/lib:$PYTHONPATH

# Create output folder
if [ ! -d "$out_dir" ]
then
    mkdir $out_dir
fi


echo "Step preprocess `date`"

preprocess.py illumina \
 --min-amplicon-size 380 --max-amplicon-size 460 \
 --five-prim-primer GGCGVACGGGTGAGTAA --three-prim-primer GTGCCAGCNGCNGCGG \
  --R1-size 250 --R2-size 250 --expected-amplicon-size 420 \
 --input-archive $frogs_test/data/test_dataset.tar.gz \
 --output-dereplicated $out_dir/01-prepro.fasta \
 --output-count $out_dir/01-prepro.tsv \
 --summary $out_dir/01-prepro.html \
 --log-file $out_dir/01-prepro.log \
 --nb-cpus $nb_cpu --mismatch-rate 0.15
 
if [ $? -ne 0 ]
then
	echo "Error in preprocess" >&2
	exit 1;
fi


echo "Step clustering `date`"

clustering.py \
 --distance 3 \
 --denoising \
 --input-fasta $out_dir/01-prepro.fasta \
 --input-count $out_dir/01-prepro.tsv \
 --output-biom $out_dir/02-clustering.biom \
 --output-fasta $out_dir/02-clustering.fasta \
 --output-compo $out_dir/02-clustering_compo.tsv \
 --log-file $out_dir/02-clustering.log \
 --nb-cpus $nb_cpu

if [ $? -ne 0 ]
then
	echo "Error in clustering" >&2
	exit 1;
fi


echo "Step remove_chimera `date`"

remove_chimera.py \
 --input-fasta $out_dir/02-clustering.fasta \
 --input-biom $out_dir/02-clustering.biom \
 --non-chimera $out_dir/03-chimera.fasta \
 --out-abundance $out_dir/03-chimera.biom \
 --summary $out_dir/03-chimera.html \
 --log-file $out_dir/03-chimera.log \
 --nb-cpus $nb_cpu
 
if [ $? -ne 0 ]
then
	echo "Error in remove_chimera" >&2
	exit 1;
fi


echo "Step filters `date`"

filters.py \
 --min-abundance 0.00005 \
 --min-sample-presence 3\
 --input-biom $out_dir/03-chimera.biom \
 --input-fasta $out_dir/03-chimera.fasta \
 --output-fasta $out_dir/04-filters.fasta \
 --output-biom $out_dir/04-filters.biom \
 --excluded $out_dir/04-filters.excluded \
 --summary $out_dir/04-filters.html \
 --log-file $out_dir/04-filters.log 

if [ $? -ne 0 ]
then
	echo "Error in filters" >&2
	exit 1;
fi


echo "Step affiliation_OTU `date`"

affiliation_OTU.py \
 --reference $frogs_test/data/db.fasta \
 --input-fasta $out_dir/04-filters.fasta \
 --input-biom $out_dir/04-filters.biom \
 --output-biom $out_dir/04-affiliation.biom \
 --summary $out_dir/04-affiliation.html \
 --log-file $out_dir/04-affiliation.log \
 --nb-cpus $nb_cpu --java-mem $java_mem \
 --rdp

if [ $? -ne 0 ]
then
	echo "Error in affiliation_OTU" >&2
	exit 1;
fi


echo "Step clusters_stat `date`"

clusters_stat.py \
 --input-biom $out_dir/04-affiliation.biom \
 --output-file $out_dir/05-clustersStat.html \
 --log-file $out_dir/05-clustersStat.log

if [ $? -ne 0 ]
then
	echo "Error in clusters_stat" >&2
	exit 1;
fi


echo "Step affiliations_stat `date`"

affiliations_stat.py \
 --input-biom $out_dir/04-affiliation.biom \
 --output-file $out_dir/06-affiliationsStat.html \
 --log-file $out_dir/06-affiliationsStat.log \
 --tax-consensus-tag "blast_taxonomy" \
 --identity-tag "perc_identity" \
 --coverage-tag "perc_query_coverage" \
 --multiple-tag "blast_affiliations" \
 --rarefaction-ranks Family Genus Species

if [ $? -ne 0 ]
then
	echo "Error in affiliations_stat" >&2
	exit 1;
fi


echo "Step biom_to_tsv `date`"

biom_to_tsv.py \
 --input-biom $out_dir/04-affiliation.biom \
 --input-fasta $out_dir/04-filters.fasta \
 --output-tsv $out_dir/07-biom2tsv.tsv \
 --output-multi-affi $out_dir/07-biom2tsv.multi \
 --log-file $out_dir/07-biom2tsv.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_tsv" >&2
	exit 1;
fi

echo "Step biom_to_stdBiom `date`"


biom_to_stdBiom.py \
 --input-biom $out_dir/04-affiliation.biom \
 --output-biom $out_dir/08-affiliation_std.biom \
 --output-metadata $out_dir/08-affiliation_multihit.tsv \
 --log-file $out_dir/08-biom2stdbiom.log

if [ $? -ne 0 ]
then
	echo "Error in biom_to_stdBiom" >&2
	exit 1;
fi

echo "Step tsv_to_biom `date`"


tsv_to_biom.py \
 --input-tsv $out_dir/07-biom2tsv.tsv \
 --input-multi-affi $out_dir/07-biom2tsv.multi \
 --output-biom $out_dir/09-tsv2biom.biom \
 --output-fasta $out_dir/09-tsv2biom.fasta \
 --log-file $out_dir/09-tsv2biom.log 

if [ $? -ne 0 ]
then
	echo "Error in tsv_to_biom" >&2
	exit 1;
fi

echo $run_pynast

if $run_pynast ; then
	echo "Step tree : pynast `date`"
	
	tree.py \
	 --nb-cpus $nb_cpu  \
	 --input-otu $out_dir/04-filters.fasta \
	 --biomfile $out_dir/04-affiliation.biom \
	 --template-pynast $frogs_test/data/otus_pynast.fasta \
	 --out-tree $out_dir/10a-tree.nwk \
	 --html $out_dir/10a-tree.html \
	 --log-file $out_dir/10a-tree.log
	 
	if [ $? -ne 0 ]
	then
		echo "Error in tree : pynast" >&2
		exit 1;
	fi
fi

echo "Step tree : mafft `date`"

tree.py \
 --nb-cpus $nb_cpu \
 --input-otu $out_dir/04-filters.fasta \
 --biomfile $out_dir/04-affiliation.biom \
 --out-tree $out_dir/10b-tree.nwk \
 --html $out_dir/10b-tree.html \
 --log-file $out_dir/10b-tree.log

if [ $? -ne 0 ]
then
	echo "Error in tree : mafft" >&2
	exit 1;
fi

echo "Step r_import_data `date`"

r_import_data.py  \
 --biomfile $frogs_test/data/chaillou.biom \
 --samplefile $frogs_test/data/sample_metadata.tsv \
 --treefile $frogs_test/data/tree.nwk \
 --rdata $out_dir/11-phylo_import.Rdata --html $out_dir/11-phylo_import.html --log-file $out_dir/11-phylo_import.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_import_data " >&2
	exit 1;
fi

echo "Step r_composition `date`"

r_composition.py  \
 --varExp EnvType --taxaRank1 Kingdom --taxaSet1 Bacteria --taxaRank2 Phylum --numberOfTaxa 9 \
 --rdata $out_dir/11-phylo_import.Rdata \
 --html $out_dir/12-phylo_composition.html --log-file $out_dir/12-phylo_composition.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_composition " >&2
	exit 1;
fi

echo "Step r_alpha_diversity `date`"

r_alpha_diversity.py  \
 --varExp EnvType \
 --rdata $out_dir/11-phylo_import.Rdata --alpha-measures Observed Chao1 Shannon \
 --alpha-out $out_dir/13-phylo_alpha_div.tsv --html $out_dir/13-phylo_alpha_div.html --log-file $out_dir/13-phylo_alpha_div.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_alpha_diversity " >&2
	exit 1;
fi

echo "Step r_beta_diversity `date`"

r_beta_diversity.py  \
 --varExp EnvType --distance-methods cc,unifrac \
 --rdata $out_dir/11-phylo_import.Rdata \
 --matrix-outdir $out_dir --html $out_dir/14-phylo_beta_div.html --log-file $out_dir/14-phylo_beta_div.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_beta_diversity " >&2
	exit 1;
fi

echo "Step r_structure `date`"

r_structure.py  \
 --varExp EnvType --ordination-method MDS \
 --rdata $out_dir/11-phylo_import.Rdata --distance-matrix $out_dir/Unifrac.tsv \
 --html $out_dir/15-phylo_structure.html --log-file $out_dir/15-phylo_structure.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_structure " >&2
	exit 1;
fi

echo "Step r_clustering `date`"

r_clustering.py  \
 --varExp EnvType \
 --rdata $out_dir/11-phylo_import.Rdata --distance-matrix $out_dir/Unifrac.tsv \
 --html $out_dir/16-phylo_structure.html --log-file $out_dir/16-phylo_structure.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_clustering " >&2
	exit 1;
fi

echo "Step r_manova `date`"

r_manova.py  \
 --varExp EnvType \
 --rdata $out_dir/11-phylo_import.Rdata --distance-matrix $out_dir/Unifrac.tsv \
 --html $out_dir/17-phylo_structure.html --log-file $out_dir/17-phylo_structure.log

 
if [ $? -ne 0 ]
then
	echo "Error in r_manova " >&2
	exit 1;
fi

if [ $delete_out_dir ] ; then
  rm -rf $out_dir
fi

echo "Completed with success"
