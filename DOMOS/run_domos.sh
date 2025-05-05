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
PROFILE="${1:-$(hostname)}"

info() { echo ; echo "$(date +'%F %T') ::${SCRIPT}::${ID}:: $* ::" ; echo ; }

notification() {
  curl --silent --insecure --data chat_id=6849911952 --data disable_notification=false --data protect_content=false --data disable_web_page_preview=false --data text="$1
$2" "https://api.telegram.org/bot7814434886:AAFQXk24RajNIwCNIT37DI38MSMqtKd0Cgw/sendMessage" >/dev/null; }

end_status() {
  STATUS=$1
  NAME="$(basename $2)"
  REST=($@)
  if [ $STATUS == 0 ]; then
    notification "$ID OK: $NAME $PROFILE" "$STATUS  ${REST[@]:2:${#REST[@]}}"
  else
    notification "$ID FAILED!!!: $NAME $PROFILE" "$STATUS  ${REST[@]:2:${#REST[@]}}"
  fi
}

start_status() {
  NAME="$(basename $1)"
  REST=($@)
  notification "$ID start: $NAME $PROFILE"
}

# Python command in conda
CONDA="conda run -n oreo --live-stream python"

## ignore errors
set +e
## python immediate output
export PYTHONUNBUFFERED=1

## make a commit of some important components before run
find . -type f \( \
  -iname 'Step*.py' \
  \) -print0 |
  xargs -t -0 git add

find ../run_profiles -type f \( \
  -iname '*.yaml' \
  \) -print0 |
  xargs -t -0 git add

## commit and push
git commit -uno -a -m "Commit before execution $(date +'%F %R')"
git push -f
git push -f --tag
git maintenance run --auto


exit

##  MAIN SEQUENCE  -------------------------------------------------------------
echo ""
echo "****    $(date +"%F %T") $USER@$HOSTNAME    ****"
echo ""


info "Get raw data from ERA5"
ascript="./Step_00_get_ERA5_data.py"
start_status "$ascript"
$CONDA "$ascript"  -p "$PROFILE"
end_status $? "$ascript"

info "Regrid of ERA5 data"
ascript="./Step_01_regrid_ERA5_parallel.py"
start_status "$ascript"
$CONDA "$ascript"  -p "$PROFILE"
end_status $? "$ascript"

info "Merge ERA with LIVAS"
ascript="./Step_02_read_LIVAS.py"
start_status "$ascript"
$CONDA "$ascript"  -p "$PROFILE"
end_status $? "$ascript"


info "END in $SECONDS seconds"
exit 0
