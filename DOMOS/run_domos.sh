#!/usr/bin/env bash
## created on 2025-03-05

#### Execute DOMOS

exec 9>"/dev/shm/$(basename $0).lock"
if ! flock -n 9  ; then
    echo "Another instance of $0 is running";
    exit 1
fi

ldir="$HOME/OREO/operation/LOGs/$(basename "$0")"
mkdir -p "$ldir"
LOG_FILE="$ldir/$(basename "$0")_$(date +%F_%R).log"
ERR_FILE="$ldir/$(basename "$0")_$(date +%F_%R).err"
exec  > >(tee -i "${LOG_FILE}")
exec 2> >(tee -i "${ERR_FILE}" >&2)

: "${ID:=$(hostname)}"
SCRIPT="$(basename "$0")"

info() { echo ; echo "$(date +'%F %T') ::${SCRIPT}::${ID}:: $* ::" ; echo ; }

echo ""
echo "****    $(date +"%F %T") $USER@$HOSTNAME    ****"
echo ""

## ignore errors
set +e



info "Activate Conda environment"
source /home/folder/miniconda/bin/activate && \
  conda activate oreo ||                         \
  { echo "Failed to activate environment"; exit 1; }


export PYTHONUNBUFFERED=1

info "Get raw data from ERA5"
./Step_00_get_ERA5_data.py


info "Do regrid on ERA5 data"
./Step_01_regrid_ERA5.py


info "Deactivate Conda environment"
conda deactivate

info "END in $SECONDS seconds"
exit 0
