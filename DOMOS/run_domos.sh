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



info "Get raw data from ERA5"
./Step_00_get_ERA5_data.py



info "END in $SECONDS seconds"
exit 0
