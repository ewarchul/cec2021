#!/bin/bash

MODE=$1
CECPATH=/cec2020
SUITES=(
  "basic" "shift" "rot" "bias" "bias-rot" "bias-shift-rot" "bias-shift" "shift-rot"
)

######################################
#      Reproduce CEC experiments
# Globals:
#   -MODE reproduction mode:
#     single benchmark || all of them
#         (main function)
######################################

run_experiment () {
  mkdir -p -m=777 $CECPATH/data/reproduction
  cd $CECPATH
  if [[ $MODE == "all" ]]; then
    exec_all
  elif [[ $MODE == "2013" ]]; then
    exec_cec 2013
  elif [[ $MODE == "2017" ]]; then
    exec_cec 2017
  elif [[ $MODE == "2021" ]]; then
    exec_cec 2021
  else
    echo "Usage: bash exp_exec.sh [CEC_VERSION | all]";
    printf "Args:";
    printf "\n CEC_VERSION - reproduce one of the available CEC variant: 2013, 2017, or 2021 \n";
    printf "\n all - reproduce all experiments \n";
  fi
}

######################################
#              Run CEC 
# Locals:
#    -year 2013, 2017, or 2021 
#
######################################

exec_cec () {
  local year=$1;
  if [[ $year == "2013" || $year == "2017" ]]; then
    Rscript -e "library(cecb)"; 
    Rscript -e "cecb::run_benchmark('$CECPATH/configs/reproduction/csa/cec$year.yml')"; 
    Rscript -e "cecb::run_benchmark('$CECPATH/configs/reproduction/msr/cec$year.yml')"; 
    Rscript -e "cecb::run_benchmark('$CECPATH/configs/reproduction/ppmf/cec$year.yml')"; 
  elif [[ $year == "2021" ]]; then
    Rscript -e "library(cecb)"; 
    Rscript -e "cecb::run_benchmark('$CECPATH/configs/reproduction/csa/cec2021.yml')"; 
    Rscript -e "cecb::run_benchmark('$CECPATH/configs/reproduction/msr/cec2021.yml')"; 
    exec_ppmf_exps
  fi
 
}

######################################
#      Run CEC 2021 for PPMF 
# Globals: 
#   - SUITES 
# Helper function for exec_cec() 
# function. It iterates through global
# variable SUITES which contains
# CEC2021 suites names.
# 
######################################

exec_ppmf_exps () {
  for suite in ${SUITES[@]}; do
    Rscript -e "library(cecb)"; 
    Rscript -e "cecb::run_benchmark('$CECPATH/configs/reproduction/ppmf/$suite.yml')"; 
  done
}

######################################
#      Reproduce all experiments 
#
######################################

exec_all () {
  exec_cec 2013
  exec_cec 2017
  exec_cec 2021
}

run_experiment
