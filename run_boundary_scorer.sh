#!/bin/bash                                                                                                                                                   
#SBATCH -o mp.%j.%N.out                                                                                                                                       
#SBATCH -e mp.%j.%N.err                                                                                                                                       
#SBATCH -J mp                                                                                                                                                 
#SBATCH --get-user-env                                                                                                                                        
#SBATCH --time=24:00:00                                                                                                                                       
#SBATCH -N 1 # number of nodes                                                                                                                                
#SBATCH -n 48                                                                                                                                                 
#SBATCH -p snowball                                                                                                                                              
#SBATCH --mem=96000                                                                                                                                           
#SBATCH --array=1-3

if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

#PATH=/home/$USER/bin/miniprot:$PATH
export PROC=miniprot
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

export KEY=Arabidopsis_thaliana
export ANNOT=~/runs_output/$KEY/miniprot/species_excluded/seq1.txt
config=~/boundary_scorer/config.txt

type=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
#sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

cd ~/boundary_scorer/miniprothint_${KEY}

mkdir tmp_${type} 
cd ./tmp_${type} 
echo "${KEY} : ${type}" >> output.txt

~/miniprot/miniprot -ut48 --aln ~/runs_output/$KEY/genome.fasta.masked ~/runs_output/$KEY/${type}.fa > ~/boundary_scorer/${KEY}_${type}.aln 2>> miniprot.stderr

~/miniprot-boundary-scorer/miniprot_boundary_scorer -o ~/boundary_scorer/${KEY}_${type}.gff -s ~/miniprot-boundary-scorer/blosum62.csv < ~/boundary_scorer/${KEY}_${type}.aln

~/miniprothint/miniprothint.py ~/boundary_scorer/${KEY}_${type}.gff --workdir ~/boundary_scorer/miniprothint_${KEY}/tmp_${type} --topNperSeed 10 --minScoreFraction 0.5


cat ~/boundary_scorer/miniprothint_${KEY}/tmp_${type}/miniprot_representatives.gff | awk '/intron/{print $3,$4,$5,$7}' | sed 's/ /_/g' > seq2_${type}.txt
~/overlap/overlapStat.pl $ANNOT seq2_${type}.txt > ~/boundary_scorer/miniprothint_${KEY}/repres_results_${type}.txt
rm seq2_${type}.txt 

echo "Done for ${KEY}:${type}" >> output.txt
#rm tmp_*



#~/GALBA/scripts/aln2hints.pl --in=miniprot_representatives.gtf --out=hints.gff --prg=miniprot


