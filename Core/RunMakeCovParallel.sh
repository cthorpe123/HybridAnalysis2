#!/bin/bash
#
# Run MakeCov.C and/or MakeCovDetvar.C as several parallel ROOT instances,
# each one processing a chunk of the variables.
#
# Usage: Core/RunMakeCovParallel.sh [options]
#   -n N     number of variables per instance (default 3)
#   -j J     maximum number of instances running at once (default: number of cores)
#   -m MODE  which macros to run: cov, detvar or both (default both)
#   -v LIST  comma separated list of variables (default: var_names in Funcs/BranchList.h)
#   -h       show this help
#
# Each instance runs inside the uboone apptainer container (unless this script is
# already running inside it) with setup.sh sourced from the repo root.
# Logs go to logs/MakeCov_<timestamp>/.

set -u

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
APPTAINER_LAUNCHER=/cvmfs/uboone.opensciencegrid.org/bin/shell_apptainer.sh

vars_per_job=3
max_jobs=$(nproc)
mode=both
var_list=""

usage(){ sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while getopts "n:j:m:v:h" opt; do
  case $opt in
    n) vars_per_job=$OPTARG ;;
    j) max_jobs=$OPTARG ;;
    m) mode=$OPTARG ;;
    v) var_list=$OPTARG ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

if ! [[ $vars_per_job =~ ^[1-9][0-9]*$ ]]; then echo "-n must be a positive integer"; exit 1; fi
if ! [[ $max_jobs =~ ^[1-9][0-9]*$ ]]; then echo "-j must be a positive integer"; exit 1; fi

case $mode in
  cov)    macros=(MakeCov) ;;
  detvar) macros=(MakeCovDetvar) ;;
  both)   macros=(MakeCov MakeCovDetvar) ;;
  *) echo "-m must be cov, detvar or both"; exit 1 ;;
esac

# Variable list - by default read var_names from BranchList.h so it stays in sync
if [ -n "$var_list" ]; then
  IFS=',' read -r -a vars <<< "$var_list"
else
  mapfile -t vars < <(
    awk '/const std::vector<std::string> var_names *= *\{/{f=1;next} f&&/\};/{exit} f' \
      "$REPO_DIR/Funcs/BranchList.h" | grep -o '"[^"]*"' | tr -d '"')
fi
if [ ${#vars[@]} -eq 0 ]; then echo "No variables found"; exit 1; fi

# Split the variables into chunks of vars_per_job
chunks=()
for ((i=0; i<${#vars[@]}; i+=vars_per_job)); do
  chunk=$(IFS=','; echo "${vars[*]:i:vars_per_job}")
  chunks+=("$chunk")
done

log_dir="$REPO_DIR/logs/MakeCov_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$log_dir"

# Detect whether we are already inside the container
in_container=false
if [ -n "${APPTAINER_CONTAINER:-}" ] || [ -n "${SINGULARITY_CONTAINER:-}" ]; then in_container=true; fi

echo "Variables (${#vars[@]}): ${vars[*]}"
echo "Macros: ${macros[*]}"
echo "${#chunks[@]} chunk(s) of up to $vars_per_job variable(s) per macro, max $max_jobs instance(s) at once"
echo "Logs: $log_dir"

# Write a self-contained job script and launch it in the background
launch(){
  local macro=$1 chunk=$2 tag=$3
  local job_script="$log_dir/$tag.sh"
  local job_log="$log_dir/$tag.log"
  local dict_dir="$log_dir/dict_$tag"

  cat > "$job_script" <<EOF
#!/bin/bash
cd "$REPO_DIR"
source setup.sh > /dev/null 2>&1
# Separate dictionary dir per instance so parallel jobs don't compile into the same place
export DIS_DICT_DIR="$dict_dir"
mkdir -p "\$DIS_DICT_DIR"
root -l -b -q 'Core/$macro.C("$chunk")'
echo \$? > "$log_dir/$tag.status"
EOF
  chmod +x "$job_script"

  if $in_container; then
    bash "$job_script" > "$job_log" 2>&1 &
  else
    printf 'bash %q\nexit\n' "$job_script" | "$APPTAINER_LAUNCHER" > "$job_log" 2>&1 &
  fi
  echo "  started $tag: $macro($chunk)"
}

# Block until fewer than max_jobs background jobs are running
# (polls rather than using "wait -n", which needs bash >= 4.3 - the SL7 container has 4.2)
wait_for_slot(){
  while [ "$(jobs -rp | wc -l)" -ge "$max_jobs" ]; do
    sleep 5
  done
}

tags=()
for macro in "${macros[@]}"; do
  for i in "${!chunks[@]}"; do
    wait_for_slot
    tag="${macro}_$i"
    launch "$macro" "${chunks[$i]}" "$tag"
    tags+=("$tag:${chunks[$i]}")
  done
done
wait

# Summary - a job counts as failed if ROOT exited non-zero, never reported, or the log shows an error
echo
echo "Summary:"
n_fail=0
for entry in "${tags[@]}"; do
  tag=${entry%%:*}
  chunk=${entry#*:}
  status_file="$log_dir/$tag.status"
  log="$log_dir/$tag.log"
  status=$(cat "$status_file" 2>/dev/null || echo "none")
  if [ "$status" != "0" ] || grep -qE "terminate called|Error in <|\*\*\* Break|segmentation violation|error: " "$log"; then
    echo "  FAILED  $tag ($chunk)  exit=$status  log: $log"
    n_fail=$((n_fail+1))
  else
    echo "  ok      $tag ($chunk)"
  fi
done

echo
if [ $n_fail -gt 0 ]; then
  echo "$n_fail job(s) failed"
  exit 1
fi
echo "All jobs finished successfully"
