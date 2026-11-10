#!/bin/bash
#
# Run the full analysis chain (the .C macros in Core, FFTest and Generators) in
# parallel batches of variables. Each batch runs every stage of the pipeline in
# order for its own variables, so batches never touch each other's files.
#
# Usage: ./RunAnalysisParallel.sh [options]
#   -n N      number of variables per batch (default 3)
#   -j J      maximum number of batches running at once (default: number of cores)
#   -v LIST   comma separated list of variables (default: var_names in Funcs/BranchList.h + Enu)
#   -s LIST   only run these stages, e.g. "MakeCov,MakeCovDetvar" or "Core/MakePrediction"
#   -x LIST   skip these stages, e.g. "MakeBinning,MakeCovNuWro"
#   -k        keep going with a batch's later stages after one of its stages fails
#   -l        list the stages in the order they run, then exit
#   -h        show this help
#
# Stages can be named by macro (MakePrediction matches both Core/ and Generators/)
# or by dir/macro. Each macro runs from its own directory, inside the uboone
# apptainer container (unless already inside it) with setup.sh sourced.
# Logs go to logs/RunAnalysis_<timestamp>/batch_<i>/.

set -u

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APPTAINER_LAUNCHER=/cvmfs/uboone.opensciencegrid.org/bin/shell_apptainer.sh

# Pipeline, in dependency order. Each entry is <dir>/<macro> relative to the repo.
ALL_STAGES=(
  # Binning templates and the main histograms
  Core/MakeBinning
  Core/MakeCov
  Core/MakeCovDetvar
  FFTest/MakeCovNuWro           # appends NuWro fake data to Histograms.root
  # Generator predictions and forward folding
  Generators/MakeGeneratorXSec
  Generators/MakeFoldingIngredients
  Generators/Recipe1
  Generators/Recipe2
  Generators/Recipe3
  Generators/Recipe4
  Generators/Recipe5
  Generators/Recipe6
  Generators/Recipe7
  Generators/ShapeOnly/Recipe1_ShapeOnly
  Generators/ShapeOnly/Recipe2_ShapeOnly
  Generators/ShapeOnly/Recipe3_ShapeOnly
  Generators/CompareRecipes
  # Plots
  Core/MakePrediction
  Core/SysTest
  Core/CheckDetvars
  FFTest/Truth
  FFTest/SpecTruthCVRes
  FFTest/CVTruthSpecRes
  FFTest/CheckEff
  Generators/MakePrediction
  Generators/GetEff
  Generators/GetPurity
  Generators/PlotFoldingIngredients
  Generators/UnderstandingFlux
)

vars_per_batch=3
max_jobs=$(nproc)
var_list=""
only_list=""
skip_list=""
keep_going=false

usage(){ sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while getopts "n:j:v:s:x:klh" opt; do
  case $opt in
    n) vars_per_batch=$OPTARG ;;
    j) max_jobs=$OPTARG ;;
    v) var_list=$OPTARG ;;
    s) only_list=$OPTARG ;;
    x) skip_list=$OPTARG ;;
    k) keep_going=true ;;
    l) printf '%s\n' "${ALL_STAGES[@]}"; exit 0 ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

if ! [[ $vars_per_batch =~ ^[1-9][0-9]*$ ]]; then echo "-n must be a positive integer"; exit 1; fi
if ! [[ $max_jobs =~ ^[1-9][0-9]*$ ]]; then echo "-j must be a positive integer"; exit 1; fi

# Does stage (dir/macro) match an entry in a comma separated list of names?
stage_in_list(){
  local stage=$1 list=$2 name
  IFS=',' read -r -a names <<< "$list"
  for name in "${names[@]}"; do
    [ "$name" = "$stage" ] || [ "$name" = "$(basename "$stage")" ] && return 0
  done
  return 1
}

# Check every name given to -s / -x is a real stage, so typos don't silently do nothing
check_names(){
  local list=$1 flag=$2 name stage found
  IFS=',' read -r -a names <<< "$list"
  for name in "${names[@]}"; do
    found=false
    for stage in "${ALL_STAGES[@]}"; do
      if [ "$name" = "$stage" ] || [ "$name" = "$(basename "$stage")" ]; then found=true; fi
    done
    if ! $found; then echo "$flag: unknown stage '$name' (use -l to list stages)"; exit 1; fi
  done
}
[ -n "$only_list" ] && check_names "$only_list" -s
[ -n "$skip_list" ] && check_names "$skip_list" -x

stages=()
for stage in "${ALL_STAGES[@]}"; do
  if [ -n "$only_list" ] && ! stage_in_list "$stage" "$only_list"; then continue; fi
  if [ -n "$skip_list" ] && stage_in_list "$stage" "$skip_list"; then continue; fi
  stages+=("$stage")
done
if [ ${#stages[@]} -eq 0 ]; then echo "No stages selected"; exit 1; fi

# Variable list - by default read var_names from BranchList.h so it stays in sync
if [ -n "$var_list" ]; then
  IFS=',' read -r -a vars <<< "$var_list"
else
  mapfile -t vars < <(
    awk '/const std::vector<std::string> var_names *= *\{/{f=1;next} f&&/\};/{exit} f' \
      "$REPO_DIR/Funcs/BranchList.h" | grep -o '"[^"]*"' | tr -d '"')
  vars+=(Enu)
fi
if [ ${#vars[@]} -eq 0 ]; then echo "No variables found"; exit 1; fi

# Split the variables into batches
batches=()
for ((i=0; i<${#vars[@]}; i+=vars_per_batch)); do
  batches+=("$(IFS=','; echo "${vars[*]:i:vars_per_batch}")")
done

log_dir="$REPO_DIR/logs/RunAnalysis_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$log_dir"

in_container=false
if [ -n "${APPTAINER_CONTAINER:-}" ] || [ -n "${SINGULARITY_CONTAINER:-}" ]; then in_container=true; fi

echo "Variables (${#vars[@]}): ${vars[*]}"
echo "Stages (${#stages[@]}): ${stages[*]}"
echo "${#batches[@]} batch(es) of up to $vars_per_batch variable(s), max $max_jobs at once"
echo "Logs: $log_dir"

# Write a job script that runs every selected stage for one batch, in order,
# then launch it in the background. Each stage writes <NN>_<macro>.status
# containing "ok", "failed <exit>" or "skipped" (after an earlier failure).
launch(){
  local i_b=$1 chunk=$2
  local b_dir="$log_dir/batch_$i_b"
  mkdir -p "$b_dir"
  echo "$chunk" > "$b_dir/variables"
  local job_script="$b_dir/job.sh"

  {
    echo '#!/bin/bash'
    echo "cd \"$REPO_DIR\""
    echo 'source setup.sh > /dev/null 2>&1'
    echo '# Separate dictionary dir per batch so parallel jobs do not compile into the same place'
    echo "export DIS_DICT_DIR=\"$b_dir/dict\""
    echo 'mkdir -p "$DIS_DICT_DIR"'
    echo "keep_going=$keep_going"
    echo 'failed=false'
    local n=0 stage tag
    for stage in "${stages[@]}"; do
      n=$((n+1))
      tag=$(printf '%02d_%s' $n "$(basename "$stage")")
      cat <<EOF
if \$failed && ! \$keep_going; then
  echo skipped > "$b_dir/$tag.status"
else
  echo "\$(date +%T) $stage"
  cd "$REPO_DIR/$(dirname "$stage")"
  root -l -b -q '$(basename "$stage").C("$chunk")' > "$b_dir/$tag.log" 2>&1
  rc=\$?
  # ROOT exits 0 on compile errors, so check the log as well as the exit code
  if [ \$rc -ne 0 ] || grep -qE "error: |\*\*\* Break|segmentation violation|terminate called" "$b_dir/$tag.log"; then
    echo "failed \$rc" > "$b_dir/$tag.status"; failed=true
  else
    echo ok > "$b_dir/$tag.status"
  fi
fi
EOF
    done
  } > "$job_script"
  chmod +x "$job_script"

  if $in_container; then
    bash "$job_script" > "$b_dir/job.log" 2>&1 &
  else
    printf 'bash %q\nexit\n' "$job_script" | "$APPTAINER_LAUNCHER" > "$b_dir/job.log" 2>&1 &
  fi
  echo "  started batch_$i_b: $chunk"
}

# Block until fewer than max_jobs batches are running
# (polls rather than using "wait -n", which needs bash >= 4.3 - the SL7 container has 4.2)
wait_for_slot(){
  while [ "$(jobs -rp | wc -l)" -ge "$max_jobs" ]; do
    sleep 5
  done
}

for i in "${!batches[@]}"; do
  wait_for_slot
  launch "$i" "${batches[$i]}"
done
wait

# Summary
echo
echo "Summary:"
n_fail=0
for i in "${!batches[@]}"; do
  b_dir="$log_dir/batch_$i"
  echo "  batch_$i (${batches[$i]}):"
  n=0
  for stage in "${stages[@]}"; do
    n=$((n+1))
    tag=$(printf '%02d_%s' $n "$(basename "$stage")")
    status=$(cat "$b_dir/$tag.status" 2>/dev/null || echo "not run")
    # Variables skipped because their input files were missing (printed by FilesExist)
    skipped=""
    if [ -f "$b_dir/$tag.log" ]; then
      skipped=$(grep "^Skipping " "$b_dir/$tag.log" | cut -d' ' -f2 | sort -u | tr '\n' ' ')
    fi
    case $status in
      ok)      line="ok     " ;;
      skipped) line="skipped" ;;
      *)       line="FAILED "; n_fail=$((n_fail+1)) ;;
    esac
    msg="    $line $stage"
    [ "$status" != "ok" ] && [ "$status" != "skipped" ] && msg="$msg ($status)  log: $b_dir/$tag.log"
    [ -n "$skipped" ] && msg="$msg  [missing inputs: $skipped]"
    echo "$msg"
  done
done

echo
if [ $n_fail -gt 0 ]; then
  echo "$n_fail stage(s) failed"
  exit 1
fi
echo "All stages finished successfully"
