#!/bin/bash

MODE=$1

######################################
#      Reproduce CEC experiments
# Globals:
#   -MODE reproduction mode:
#     single benchmark || all of them
#
######################################

run_experiment () {
  mkdir /cec2020/data/reproduction
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
  Rscript -e "library(cecb)"; 
  Rscript -e "cecb::run_benchmark('configs/csa/cec$1.yml')"; 
  Rscript -e "cecb::run_benchmark('configs/msr/cec$1.yml')"; 
  Rscript -e "cecb::run_benchmark('configs/ppmf/cec$1.yml')"; 
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
