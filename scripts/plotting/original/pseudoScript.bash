#!/bin/bash
export DF=2022-09-23
export DF_=20220923

export NR=~/Applications/Nextclade_CLI/Reference/sars-cov-2/
export NH=~/Applications/Nextclade_CLI/

##-- update nextclade_cli (https://github.com/nextstrain/nextclade/tree/master/packages/nextclade_cli)
cd $NH
##curl -fSL "https://github.com/nextstrain/nextclade/releases/download/1.11.0/nextclade-Linux-x86_64" -o "nextclade" && chmod +x nextclade
curl -fSL "https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu" -o "nextclade" && chmod +x nextclade

## update also reference and other files
cp $NR/qc.json $NR/qc.json.$DF
$NH/nextclade dataset get --name 'sars-cov-2' --output-dir $NR

cd $NR
diff -B4 qc.json qc.json.$DF
mv qc.json.$DF qc.json

##-- update pangolin (https://cov-lineages.org/pangolin_docs/updating.html
conda activate pangolin
cd /opt/pangolin/

## major release
git pull
conda update -n base -c defaults conda
conda env update -f environment.yml
pip install .

## minor release
pangolin --update
pangolin --update-data

pangolin -v
pangolin -pv
##pangolin -dv
conda deactivate

##-- update sc2rf 
cd ~/Work/Coronaseq/sc2rf
git pull

##-- download data from eltu
cd ~/Work/Coronaseq
scp -r eltu:cogcz/seq/df-${DF_} .
cd df-${DF_}

##-- run nextclade_cli
cd ~/Work/Coronaseq/df-${DF_}

$NH/nextclade --version
export NV=`$NH/nextclade --version | grep ^nextclade | cut -d" " -f2`

$NH/nextclade run \
	      --in-order \
              --input-dataset $NR \
	      --output-tsv datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv \
	      --output-insertions datafreeze-${DF}_12_weeks_good.nextclade.${NV}.insertions.csv \
	      --output-fasta datafreeze-${DF}_12_weeks_good.aligned.fasta \
	      datafreeze-${DF}_12_weeks_good.fasta.xz

##-- run pangolin
conda activate pangolin

pangolin -v
pangolin -pv
export PV=`pangolin -v  | grep ^pangolin   | cut -d" " -f2`
export PU=`pangolin -pv | grep ^pango-designation | cut -d" " -f2`

pangolin --analysis-mode "usher" ./datafreeze-${DF}_12_weeks_good.fasta.xz --threads 8

mv lineage_report.csv datafreeze-${DF}_12_weeks_good.pangolin-${PV}-UShER-${PU}.csv

conda deactivate

##-- some summary statistics
cut -d"," -f2 datafreeze-${DF}_12_weeks_good.nextclade.${NV}.insertions.csv | sort | uniq -c | sort -r -n
cut -f2 datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv | sort | uniq -c | sort -n -r | grep -v clade
cut -f5 datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv | sort | uniq -c
cut -f3 datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv | sort | uniq -c | sort -n -r | grep -v Nextclade_pango
cut -f2 -d"," datafreeze-${DF}_12_weeks_good.pangolin-${PV}-UShER-${PU}.csv      | sort | uniq -c | sort -n -r | grep -v lineage
cut -f2 -d"," datafreeze-${DF}_12_weeks_good.pangolin-${PV}-UShER-${PU}.csv      | sort | uniq -c | sort -n -r | grep -v lineage | wc -l

##-- compare nextclade-pango and UShER (columns are isolate, Nextclade-pango, UShER)
diff -y <(cut -f1,3 --output-delimiter="," datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv | sort -r) \
        <(cut -f1,2 -d"," datafreeze-${DF}_12_weeks_good.pangolin-${PV}-UShER-${PU}.csv | sort -r) | \
	sed -r 's/[[:space:]]+/,/g' | sed 's/,|,/,/g' | awk -F"," '{if (($1 == $3) || ($1 == "seqName")) {print $1 "\t" $2 "\t" $4}}' | \
	sed 's/seqName\tNextclade_pango\tlineage/Isolate\tNextclade_pango\tUShER/g' > lineageAssignmentChanges.tsv

##-- compare UShER and Nextclade 
diff -y <(cut -f1,2 -d"," datafreeze-${DF}_12_weeks_good.pangolin-${PV}-UShER-${PU}.csv | sort -r) \
        <(cut -f1,2       datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv | sort -r) | \
	sed -r 's/[[:space:]]+/,/g' | sed 's/,|,/,/g' | \
	awk -F"," '{if (($1 == $3) || ($1 == "taxon")) {print $1 "\t" $2 "\t" $4 "," $5}}' | \
	sed 's/taxon\tlineage\tclade,/Isolate\tUShER\tNextclade/g' | sed 's/"//g' > cladeAssignment.tsv

R CMD BATCH lineageAssignment.R
R CMD BATCH cladeAssignment.R

##-- collect sequences for recombinant reanalysis
## unassigned by pango
grep -i unassigned datafreeze-${DF}_12_weeks_good.pangolin-${PV}-UShER-${PU}.csv | cut -d"," -f1 | sort -u >  putativeRecombinants.tmp
## recombinant in nextclade
grep -i recombinant datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv | cut -f1 | sort -u >> putativeRecombinants.tmp

sort -u putativeRecombinants.tmp > putativeRecombinants.txt
rm putativeRecombinants.tmp

cat putativeRecombinants.txt | xargs -I % sh -c 'grep -A1 % datafreeze-${DF}_12_weeks_good.aligned.fasta' > putativeRecombinants.fasta

cd ../sc2rf/
conda deactivate
./sc2rf.py ../df-${DF_}/putativeRecombinants.fasta --csvfile ../df-${DF_}/putativeRecombinants.sc2rf.tmp
cd ../df-${DF_}

head -n  1 putativeRecombinants.sc2rf.tmp > datafreeze-${DF}_12_weeks_good.putativeRecombinants.sc2rf.csv
tail -n +2 putativeRecombinants.sc2rf.tmp | sort -t"," -k2 >> datafreeze-${DF}_12_weeks_good.putativeRecombinants.sc2rf.csv
rm putativeRecombinants.sc2rf.tmp

##tazbyk-joCpy5-jycqet
scp datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv eltu:cogcz/seq/df-${DF_}
scp datafreeze-${DF}_12_weeks_good.nextclade.${NV}.insertions.csv eltu:cogcz/seq/df-${DF_}
scp datafreeze-${DF}_12_weeks_good.pangolin-${PV}-UShER-${PU}.csv eltu:cogcz/seq/df-${DF_}
scp datafreeze-${DF}_12_weeks_good.putativeRecombinants.sc2rf.csv eltu:cogcz/seq/df-${DF_}

echo The nextclade and pangolin outputs are in the files datafreeze-${DF}_12_weeks_good.nextclade.${NV}.tsv, \
     datafreeze-${DF}_12_weeks_good.nextclade.${NV}.insertions.csv, \
     datafreeze-${DF}_12_weeks_good.pangolin-${PV}-UShER-${PU}.csv and \
     datafreeze-${DF}_12_weeks_good.putativeRecombinants.sc2rf.csv.

cd ..
scp -r df-${DF_} eltu:cogcz/kolarmi/
