#!/bin/bash
#$ -N dbcan4
#$ -cwd
#$ -pe mpi 4
#$ -j y
#$ -l h_vmem=8G

export PATH=/home/install/miniconda3/envs/run_dbcan4/bin:$PATH

source activate run_dbcan4

run_dbcan *.fna --dbCANFile /tres/DB/dbCAN4/dbCAN.txt  -c cluster --cgc_dis 2 --cgc_sig_genes all --db_dir /tres/DB/dbCAN4 --pul /tres/DB/dbCAN4/PUL.faa --odbcanpul prok --use_signalP all --out_dir output_dbcan4/

conda deactivate
