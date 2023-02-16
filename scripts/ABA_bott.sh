#!/bin/bash
#
# Jeff Eilbott, 2017, jeilbott@surveybott.com

# inputs
PROJECT="/ysm-gpfs/project/ejh22/bioinfo/startle"
UPLOAD="$PROJECT/upload"
DATA="$PROJECT/data"
LOG="$PROJECT/log"
CORES=12
QUEUE="general"
EMAIL="jeilbott@surveybott.com ellen.hoffman@yale.edu"
BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE="$(dirname $BASE)"
GROUP="hoffman"
# process inputs
while getopts "u:d:n:q:e:b:g:l:" OPTION
do
  case $OPTION in
    u)
      UPLOAD="$OPTARG"
      ;;
    d)
      DATA="$OPTARG"
      ;;
    n)
      CORES="$OPTARG"
      ;;
    q)
      QUEUE="$OPTARG"
      ;;
    e)
      EMAIL="$OPTARG"
      ;;
    b)
      BASE="$OPTARG"
      ;;
    g)
      GROUP="$OPTARG"
      ;;
    l)
      LOG="$OPTARG"
      ;;
  esac
done
if [ "$GROUP" != "$(id -gn)" ]; then
  newgrp $GROUP
fi
if [ ! -d "$LOG" ]; then
  mkdir -p $LOG/slurm
fi

sbatch --partition="$QUEUE" --cpus-per-task="$CORES" --output="$LOG/slurm/vsr_slurm_%j.log" $BASE/scripts/ABA_slurm.sh -u "$UPLOAD" -d "$DATA" -b "$BASE" -e "$EMAIL" -g "$GROUP" -l "$LOG" > /dev/null
