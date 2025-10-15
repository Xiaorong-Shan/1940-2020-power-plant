#!/bin/sh

#SBATCH --job-name=04hyads_to_pm25
#SBATCH --partition=normal
#SBATCH --constraint=intel
#SBATCH --output=/scratch/%u/logs/%x-%A_%a.out
#SBATCH --error=/scratch/%u/logs/%x-%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xshan2@gmu.edu

# 内存说明：若你想申请“总内存 50G”，用 --mem=50G；若坚持每核 50G，请保留 --mem-per-cpu=50G 并把 -c 设小
#SBATCH --mem=50G
#SBATCH --time=01-06:00
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --array=0-5   # 0→1940, 1→1950, ..., 5→1990


# ✅ 你要求的两个模块（必须保留）
module load gnu9
module load r-disperseR/0.1.0


mkdir -p /scratch/$USER/logs
cd /scratch/xshan2/R_Code/disperseR

Rscript  /scratch/xshan2/R_Code/disperseR/hyads_1940_1990_to_pm25.R
