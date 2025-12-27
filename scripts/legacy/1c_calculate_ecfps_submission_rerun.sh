#!/bin/bash
#

#SBATCH -J enamine-real--ecfp
#SBATCH --chdir=/aloy/home/acomajuncosa/Ersilia/ready-to-screen-enamine-real
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-89%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --output=/aloy/scratch/acomajuncosa/Ersilia/ready-to-screen-enamine-real_rerun/%x_%a.out
#SBATCH -p sbnb_cpu_sphr,sbnb_cpu_zen3

# Loads default environment configuration
export SINGULARITYENV_LD_LIBRARY_PATH=$LD_LIBRARY_PATH #/.singularity.d/libs
export SINGULARITY_BINDPATH="/home/sbnb:/aloy/home,/data/sbnb/data:/aloy/data,/data/sbnb/scratch:/aloy/scratch"

# Load cuda libraries
export LD_LIBRARY_PATH=/apps/manual/software/CUDA/11.6.1/lib64:/apps/manual/software/CUDA/11.6.1/targets/x86_64-linux/lib:/apps/manual/software/CUDA/11.6.1/extras/CUPTI/lib64/:/apps/manual/software/CUDA/11.6.1/nvvm/lib64/:$LD_LIBRARY_PATH

alpha=(88 80 65 57 31 990 987 979 977 976 974 971 967 965 961 959 957 955 950 946 944 943 936 932 930 925 917 914 913 912 911 906 905 904 901 900 898 896 894 886 885 880 871 869 868 854 853 850 846 845 844 843 841 839 832 831 825 824 823 822 819 816 812 807 798 796 787 785 778 775 774 771 770 768 767 764 763 762 754 752 748 743 740 735 733 729 727 726)

singularity exec --cleanenv /apps/singularity/ood_images/docker_irb_intel-optimized-tensorflow-avx512-2.13-pip-conda-jupyter-v6.sif ./scripts/1c_calculate_ecfps.sh ${alpha[$SLURM_ARRAY_TASK_ID]}