# shellcheck shell=bash
# GENERATED COPY of biowulf-snakemake-launcher v1.0.0 (279854e); do not edit here.
# Shared Snakemake launcher for Biowulf (SLURM + environment modules).
# Sourced by each pipeline's run_<name> wrapper, which then calls: main "$@"
# Author: Erica Pehrsson, PhD; Genetics Branch, NCI. Based on the CCBR snakemake template.
# Copyright (c) 2026 Erica Pehrsson, Genetics Branch, NCI; Copyright (c) 2021 CCBR. MIT License (see LICENSE in biowulf-snakemake-launcher).
# Pipelines get a copy at bin/launcher.sh via sync_launcher.sh; edit the copy in biowulf-snakemake-launcher only.

# Filled in by sync_launcher.sh in the copied file.
LAUNCHER_VERSION_STAMP="1.0.0 (279854e)"

set -eo pipefail

: "${PIPELINE_HOME:?PIPELINE_HOME must be set by the wrapper before sourcing launcher.sh}"
: "${PIPELINE_NAME:?PIPELINE_NAME must be set by the wrapper before sourcing launcher.sh}"

PIPELINE_DESCRIPTION="${PIPELINE_DESCRIPTION:-run the ${PIPELINE_NAME} pipeline}"
SBATCH_JOB_NAME="${SBATCH_JOB_NAME:-$PIPELINE_NAME}"
CLUSTER_CONFIG_SRC="${CLUSTER_CONFIG_SRC:-resources/cluster.yaml}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"
PYTHON_VERSION="${PYTHON_VERSION:-python/3.9}"
SNAKEMAKE_VERSION="${SNAKEMAKE_VERSION:-snakemake/7.19.1}"
SBATCH_MEM="${SBATCH_MEM:-10g}"
SBATCH_TIME="${SBATCH_TIME:-96:00:00}"
SBATCH_PARTITION="${SBATCH_PARTITION:-norm}"
SBATCH_CPUS="${SBATCH_CPUS:-2}"
PIPELINE_SBATCH_SETUP="${PIPELINE_SBATCH_SETUP:-}"
PIPELINE_USAGE_EXTRA="${PIPELINE_USAGE_EXTRA:-}"

LAUNCHER_HOME="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
SCRIPTNAME="$(basename "$0")"
SNAKEFILE="${PIPELINE_HOME}/workflow/Snakefile"
# copied from PIPELINE_HOME/config/ into WORKDIR/config/ at init
TEMPLATE_FILES="config.yaml samples.tsv"
ESSENTIAL_FILES="config.yaml samples.tsv cluster.yaml"
CLUSTER_SUBMIT='sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}'

RUNMODE=""
WORKDIR=""
SAMPLESHEET=""

function _read_version() {
  if [[ -f "$1/VERSION" ]]; then
    tr -d '[:space:]' < "$1/VERSION"
  else
    echo "unknown"
  fi
}

PIPELINE_VERSION="$(_read_version "$PIPELINE_HOME")"
LAUNCHER_VERSION="${LAUNCHER_VERSION_STAMP:-$(_read_version "$LAUNCHER_HOME")}"
GIT_COMMIT_TAG="$(git -C "$PIPELINE_HOME" rev-parse --short HEAD 2>/dev/null)" || GIT_COMMIT_TAG="n/a"
GIT_TAG="$(git -C "$PIPELINE_HOME" describe --tags --exact-match HEAD 2>/dev/null)" || GIT_TAG=""
GIT_COMMIT_TAG="${GIT_COMMIT_TAG}${GIT_TAG:+ (${GIT_TAG})}"

function usage() {
  cat << EOF
${SCRIPTNAME} (${PIPELINE_NAME} v${PIPELINE_VERSION})
  --> ${PIPELINE_DESCRIPTION}

USAGE:
  bash ${SCRIPTNAME} -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR> [options]

Required arguments:
  -m/--runmode=<RUNMODE>   One of:
                             init   : create WORKDIR and copy config files into WORKDIR/config/
                             dryrun : Snakemake dry run (writes WORKDIR/dryrun.log)
                             unlock : unlock WORKDIR after an interrupted run
                             run    : submit the pipeline to SLURM
  -w/--workdir=<WORKDIR>   Output folder (absolute or relative path)

Options:
  -s/--samplesheet=<PATH>  (init only) sample sheet to copy to WORKDIR/config/samples.tsv
                           instead of the template
  -h/--help                show this message
  -v/--version             show pipeline and launcher versions
EOF
  if [[ -n "$PIPELINE_USAGE_EXTRA" ]]; then
    echo ""
    echo "Pipeline-specific options:"
    echo "$PIPELINE_USAGE_EXTRA"
  fi
}

function err() {
  {
    echo ""
    echo "##### ERROR #####"
    echo "  $*"
    echo ""
    usage
  } >&2
  exit 1
}

function _timestamp() {
  date -r "$1" +%Y%m%d%H%M%S
}

function init() {
  local f
  if [[ -e "$WORKDIR" ]]; then err "Folder $WORKDIR already exists!"; fi
  for f in $TEMPLATE_FILES; do
    if [[ "$f" == "samples.tsv" && -n "$SAMPLESHEET" ]]; then continue; fi
    if [[ ! -f "${PIPELINE_HOME}/config/${f}" ]]; then err "Template ${PIPELINE_HOME}/config/${f} not found"; fi
  done
  if [[ ! -f "${PIPELINE_HOME}/${CLUSTER_CONFIG_SRC}" ]]; then
    err "Cluster config ${PIPELINE_HOME}/${CLUSTER_CONFIG_SRC} not found"
  fi

  mkdir -p "${WORKDIR}/config" "${WORKDIR}/logs" "${WORKDIR}/stats"

  for f in $TEMPLATE_FILES; do
    if [[ "$f" == "samples.tsv" && -n "$SAMPLESHEET" ]]; then
      echo "Copying sample sheet: ${SAMPLESHEET}"
      cp "$SAMPLESHEET" "${WORKDIR}/config/samples.tsv"
    else
      echo "Copying template: config/${f}"
      sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" \
          -e "s/WORKDIR/${WORKDIR//\//\\/}/g" \
          "${PIPELINE_HOME}/config/${f}" > "${WORKDIR}/config/${f}"
    fi
  done
  echo "Copying cluster config: ${CLUSTER_CONFIG_SRC}"
  cp "${PIPELINE_HOME}/${CLUSTER_CONFIG_SRC}" "${WORKDIR}/config/cluster.yaml"

  if declare -F pipeline_init_hook > /dev/null; then pipeline_init_hook; fi

  echo "Done initializing ${WORKDIR}"
  echo "Edit ${WORKDIR}/config/config.yaml and ${WORKDIR}/config/samples.tsv, then run with --runmode=dryrun"
}

function check_essential_files() {
  local f
  if [[ ! -d "$WORKDIR" ]]; then err "Folder $WORKDIR does not exist! Run with --runmode=init first."; fi
  for f in $ESSENTIAL_FILES; do
    if [[ ! -f "${WORKDIR}/config/${f}" ]]; then
      err "'${f}' not found in ${WORKDIR}/config/ ... initialize first!"
    fi
  done
}

function load_modules() {
  if ! declare -F module > /dev/null; then
    err "The 'module' command is not available; this launcher must be run on Biowulf."
  fi
  module purge
  module load "$PYTHON_VERSION" || err "Failed to load module $PYTHON_VERSION"
  module load "$SNAKEMAKE_VERSION" || err "Failed to load module $SNAKEMAKE_VERSION"
}

function _snakemake_args() {
  SNAKEMAKE_ARGS=(
    -s "$SNAKEFILE"
    --directory "$WORKDIR"
    --use-envmodules
    --printshellcmds
    --latency-wait "$LATENCY_WAIT"
    --configfile "${WORKDIR}/config/config.yaml"
    --cluster-config "${WORKDIR}/config/cluster.yaml"
    --cluster "$CLUSTER_SUBMIT"
    -j 500
    --rerun-incomplete
    --keep-going
  )
}

function dryrun() {
  check_essential_files
  load_modules
  if [[ -f "${WORKDIR}/dryrun.log" ]]; then
    mkdir -p "${WORKDIR}/logs/dryrun"
    mv "${WORKDIR}/dryrun.log" "${WORKDIR}/logs/dryrun/dryrun.$(_timestamp "${WORKDIR}/dryrun.log").log"
  fi
  _snakemake_args
  snakemake --dry-run "${SNAKEMAKE_ARGS[@]}" 2>&1 | tee "${WORKDIR}/dryrun.log"
}

function unlock() {
  check_essential_files
  load_modules
  _snakemake_args
  snakemake --unlock "${SNAKEMAKE_ARGS[@]}"
}

function _shell_quote() {
  if [[ "$1" =~ ^[A-Za-z0-9_./=:+-]+$ ]]; then
    printf '%s' "$1"
  else
    printf "'%s'" "${1//\'/\'\\\'\'}"
  fi
}

function runslurm() {
  local sbatch_file="${WORKDIR}/submit_script.sbatch" runinfo="${WORKDIR}/runinfo.yaml"
  local jobid ts f i quoted=""
  check_essential_files
  load_modules

  # archive the previous run
  if [[ -f "${WORKDIR}/snakemake.log" ]]; then
    ts=$(_timestamp "${WORKDIR}/snakemake.log")
    mv "${WORKDIR}/snakemake.log" "${WORKDIR}/stats/snakemake.${ts}.log"
    if [[ -f "${WORKDIR}/snakemake.stats" ]]; then
      mv "${WORKDIR}/snakemake.stats" "${WORKDIR}/stats/snakemake.${ts}.stats"
    fi
  fi
  for f in "${WORKDIR}"/slurm-*.out; do
    if [[ -f "$f" ]]; then mv "$f" "${WORKDIR}/logs/"; fi
  done
  if [[ -f "$runinfo" ]]; then mv "$runinfo" "${runinfo}.$(_timestamp "$runinfo")"; fi

  cat > "$runinfo" << EOF
pipeline: "${PIPELINE_NAME}"
version: "${PIPELINE_VERSION}"
pipeline_dir: "${PIPELINE_HOME}"
git_commit: "${GIT_COMMIT_TAG}"
launcher_version: "${LAUNCHER_VERSION}"
python: "$(python --version 2>&1)"
snakemake: "$(snakemake --version 2>&1)"
login: "$(whoami)"
date: "$(date)"
EOF

  _snakemake_args
  # pair each option with its value on one line
  for ((i = 0; i < ${#SNAKEMAKE_ARGS[@]}; i++)); do
    quoted+=" \\"$'\n'"  $(_shell_quote "${SNAKEMAKE_ARGS[i]}")"
    if [[ $((i + 1)) -lt ${#SNAKEMAKE_ARGS[@]} && "${SNAKEMAKE_ARGS[i+1]}" != -* ]]; then
      i=$((i + 1))
      quoted+=" $(_shell_quote "${SNAKEMAKE_ARGS[i]}")"
    fi
  done

  cat > "$sbatch_file" << EOF
#!/bin/bash
#SBATCH --job-name="${SBATCH_JOB_NAME}"
#SBATCH --mem=${SBATCH_MEM}
#SBATCH --partition="${SBATCH_PARTITION}"
#SBATCH --time=${SBATCH_TIME}
#SBATCH --cpus-per-task=${SBATCH_CPUS}

${PIPELINE_SBATCH_SETUP}

module purge
module load ${PYTHON_VERSION}
module load ${SNAKEMAKE_VERSION}

cd $(_shell_quote "$WORKDIR")

snakemake${quoted} \\
  --stats $(_shell_quote "${WORKDIR}/snakemake.stats") \\
  2>&1 | tee $(_shell_quote "${WORKDIR}/snakemake.log")
exit "\${PIPESTATUS[0]}"
EOF

  if [[ "${LAUNCHER_NO_SUBMIT:-0}" == "1" ]]; then
    echo "LAUNCHER_NO_SUBMIT=1: wrote ${sbatch_file} (not submitted)"
    return 0
  fi
  jobid=$(cd "$WORKDIR" && sbatch --parsable "$sbatch_file") || err "sbatch submission failed"
  echo "Submitted ${PIPELINE_NAME} as SLURM job ${jobid}"
}

function main() {
  local i
  if [[ $# -eq 0 ]]; then usage; exit 1; fi

  for i in "$@"; do
    case "$i" in
      -m=*|--runmode=*) RUNMODE="${i#*=}" ;;
      -w=*|--workdir=*) WORKDIR="${i#*=}" ;;
      -s=*|--samplesheet=*) SAMPLESHEET="${i#*=}" ;;
      -h|--help) usage; exit 0 ;;
      -v|--version)
        echo "${PIPELINE_NAME} v${PIPELINE_VERSION} (git: ${GIT_COMMIT_TAG}); launcher v${LAUNCHER_VERSION}"
        exit 0
        ;;
      *)
        if declare -F pipeline_parse_arg > /dev/null && pipeline_parse_arg "$i"; then continue; fi
        err "Unknown argument: $i"
        ;;
    esac
  done

  if [[ -z "$RUNMODE" ]]; then err "--runmode is required"; fi
  if [[ -z "$WORKDIR" ]]; then err "--workdir is required"; fi
  WORKDIR="$(readlink -m "$WORKDIR")"

  if [[ -n "$SAMPLESHEET" ]]; then
    if [[ "$RUNMODE" != "init" ]]; then err "--samplesheet can only be used with --runmode=init"; fi
    if [[ ! -f "$SAMPLESHEET" ]]; then err "Sample sheet $SAMPLESHEET not found"; fi
    SAMPLESHEET="$(readlink -f "$SAMPLESHEET")"
  fi

  echo "Pipeline: ${PIPELINE_NAME} v${PIPELINE_VERSION}"
  echo "Pipeline Dir: ${PIPELINE_HOME}"
  echo "Git Commit/Tag: ${GIT_COMMIT_TAG}"
  echo "Launcher: v${LAUNCHER_VERSION}"
  echo "Working Dir: ${WORKDIR}"

  case "$RUNMODE" in
    init) init ;;
    dryrun) dryrun ;;
    unlock) unlock ;;
    run) runslurm ;;
    *) err "Unknown runmode \"$RUNMODE\"" ;;
  esac
}
