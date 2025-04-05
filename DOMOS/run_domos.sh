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

notification() {
  curl --silent --insecure --data chat_id=6849911952 --data disable_notification=false --data protect_content=false --data disable_web_page_preview=false --data text="$1
$2" "https://api.telegram.org/bot7814434886:AAFQXk24RajNIwCNIT37DI38MSMqtKd0Cgw/sendMessage" >/dev/null; }

end_status() {
  STATUS=$1
  NAME=$2
  REST=($@)
  if [ $STATUS == 0 ]; then
    notification "end: $ID $NAME" "Status: $STATUS  ${REST[@]:2:${#REST[@]}}"
  else
    notification "FAILED: $ID $NAME !!!" "Status: $STATUS  ${REST[@]:2:${#REST[@]}}"
  fi
}

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
notification "start: $ID Step_00_get_ERA5_data.py"
./Step_00_get_ERA5_data.py
end_status $? "Step_00_get_ERA5_data.py"


info "Regrid of ERA5 data"
notification "start: $ID Step_01_regrid_ERA5.py"
./Step_01_regrid_ERA5.py
end_status $? "Step_01_regrid_ERA5.py"


info "Merge ERA with LIVAS"
notification "start: $ID Step_02_read_LIVAS.py"
./Step_02_read_LIVAS.py
end_status $? "Step_02_read_LIVAS.py"


info "Deactivate Conda environment"
conda deactivate

info "END in $SECONDS seconds"
exit 0
