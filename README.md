# CEC 2021

## About 

This repository contains results of numerical experiments performed on CEC 2021 for the paper:

### *The step-size adaptation procedure utilizing fitness of the population midpoint*

It also contains the source code of used optimization algorithms and scripts to reproduce conducted experiments.  

## Content 

* `R/` stores implementation of used algorithms in **R** 
* `data/` stores obtained data from numerical experiments:
    - `cec13` has data from CEC 2013
    - `cec17` has data from CEC 2017
    - `cec21` has data from CEC 2021
* `time-complexity` contains **R** script to compute complexity of considered algorithm in the paper
* `configs` stores benchmark configuration files written in `YAML`.

### Algorithms

Each algorithm is implemented in **R** programming language and based on the source code of the CRAN package [{cmaesr}](https://cran.r-project.org/web/packages/cmaesr/index.html).

### Data 

Recorded data from experiments, i.e. error values (`F(x) - F(x*)`), is stored in `m/` and `M/` directories. The first directory contains the error values recorded with a budget step defined as in CEC 2021 specification:

```r
d^(0:15/5 - 3) * maxFES
```

where `d` stands for the dimension and `maxFES` is budget for objective function evaluation.

The data from the directory `M/` stores the error values recorded in the same manner as in CEC 2017:

```r
c(0.01, 0.02, 0.03, 0.05, seq(0.1, 1, 0.1)) * maxFES
```

### Benchmark configs

Experiments were performed using [{cecb}](https://github.com/ewarchul/cecb), i.e. the **R** package written by me; mainly for the self-purposes. The package executes CECs experiments using the `YAML` data-serialization language. See the documentation for further details. 

Also, we used [{cecs}](https://github.com/ewarchul/cecs) package which provides an interface from **R** to the benchmark functions implemented in **C**.


## Experiment reproduction

If you want to reproduce our results you have to have installed docker in version 19.03.8 on your machine and run the commands written below in the root directory of the repository:

```bash
docker build -t cec_experiments . 
docker run cec_experiments --cec="[CEC_VERSION | all]"
```

Those commands will build `docker` image and execute the **bash** script which takes option `cec`:

- set `--cec` equals to `CEC_VERSION` which takes values from {2013, 2017, 2021} to reproduce experiments only for the one CEC version 
- set `--cec` equals to `all` to reproduce all experiments presented in the paper.

## Contact 

Feel free to contact me if you have any suggestions, questions, etc.: [ewarchul@gmail.com](mailto:ewarchul@gmail.com?subject=[CEC2021])
