#!/bin/bash
#
# Jeff Eilbott, 2018, jeilbott@surveybott.com

#SBATCH --job-name=vsr_bott
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=6000
#SBATCH --time=6:00:00

# process inputs
while getopts "u:d:b:e:g:l:" OPTION
do
  case $OPTION in
    u)
      UPLOAD="$OPTARG"
      ;;
    d)
      DATA="$OPTARG"
      ;;
    b)
      BASE="$OPTARG"
      ;;
    e)
      EMAIL="$OPTARG"
      ;;
    g)
      GROUP="$OPTARG"
      ;;
    l)
      LOG="$OPTARG"
      ;;
  esac
done
if [ ! -d "$DATA" ]; then
  mkdir -p "$DATA"
fi
if [ ! -d "$UPLOAD" ] || [ ! -d "$BASE" ]; then
  echo  "usage: ABA_slurm.sh -u <upload dir> -d <data dir> -b <base dir> -e [e-mail addresses]"
  exit 1
fi
source $HOME/.bashrc

# setup log output
LOGFILE="$LOG/vsr_bott_$(date '+%Y%m%d_%H%M%S').log"

# run ABA_CTRL_batch
module load MATLAB/2018b
matlab -nodisplay -r "addpath(genpath('$BASE')); ABA_CTRL_batch('$DATA','upload','$UPLOAD','logFile','$LOGFILE','combine',true,'indiv',true,'overwrite',false); exit"

# send e-mail
if [ -n "$EMAIL" ];then
  mail -s "New VSR Processed" -r "VSR Bott<info@surveybott.com>" "$EMAIL" < $LOGFILE
  while [ -n "$(ps | grep sendmail)" ]; do
    sleep 2
  done
  #rm $SUMMARY
fi

# fix permissions
ANALYSIS=$(echo $DATA | sed 's/data/analysis/')
ERROR=$(echo $DATA | sed 's/data/error/')
chmod -R 750 $DATA
setfacl -R -m u:ejh22:rwx $DATA
setfacl -R -m u:je45:rwx $DATA
chmod -R 750 $ANALYSIS
setfacl -R -m u:ejh22:rwx $ANALYSIS
setfacl -R -m u:je45:rwx $ANALYSIS
chmod -R 770 $ERROR/*
setfacl -R -m u:ejh22:rwx $ERROR
setfacl -R -m u:je45:rwx $ERROR
chmod -R 750 $LOG
setfacl -R -m u:ejh22:rwx $LOG
setfacl -R -m u:je45:rwx $LOG
if [ -n "$GROUP" ]; then
   chown -R :$GROUP $DATA
   chown -R :$GROUP $ANALYSIS
   chown -R :$GROUP $ERROR
   chown -R :$GROUP $LOG
fi
