#!/bin/bash
#

#SBATCH -J enamine-real--ecfp
#SBATCH --chdir=/aloy/home/acomajuncosa/Ersilia/ready-to-screen-enamine-real
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --array=0-157%1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --output=/aloy/scratch/acomajuncosa/Ersilia/ready-to-screen-enamine-real_v3/%x_%a.out
#SBATCH -p sbnb_cpu_sphr,sbnb_cpu_zen3

# Loads default environment configuration
export SINGULARITYENV_LD_LIBRARY_PATH=$LD_LIBRARY_PATH #/.singularity.d/libs
export SINGULARITY_BINDPATH="/home/sbnb:/aloy/home,/data/sbnb/data:/aloy/data,/data/sbnb/scratch:/aloy/scratch"

# Load cuda libraries
export LD_LIBRARY_PATH=/apps/manual/software/CUDA/11.6.1/lib64:/apps/manual/software/CUDA/11.6.1/targets/x86_64-linux/lib:/apps/manual/software/CUDA/11.6.1/extras/CUPTI/lib64/:/apps/manual/software/CUDA/11.6.1/nvvm/lib64/:$LD_LIBRARY_PATH

alpha=(543 838 839 840 841 842 843 844 845 846 847 848 849 850 851 852 853 854 855 856 857 858 859 860 861 862 863 864 865 866 867 868 869 870 871 872 873 874 875 876 877 878 879 880 881 882 883 884 885 886 887 888 889 890 891 892 893 894 895 896 897 898 899 900 901 902 903 904 905 906 907 908 909 910 911 912 913 914 915 916 917 918 919 920 921 922 923 924 925 926 927 928 929 930 931 932 933 934 935 936 937 938 939 940 941 942 943 944 945 946 947 948 949 950 951 952 953 954 955 956 957 958 959 960 961 962 963 964 965 966 967 968 969 970 971 972 973 974 975 976 977 978 979 980 981 982 983 984 985 986 987 988 989 990 991 992 993)

singularity exec --cleanenv /apps/singularity/ood_images/docker_irb_intel-optimized-tensorflow-avx512-2.13-pip-conda-jupyter-v6.sif ./scripts/run.sh ${alpha[$SLURM_ARRAY_TASK_ID]}
