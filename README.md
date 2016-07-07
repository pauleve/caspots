# Caspots - Boolean network inference from time series data with perturbations

## Requirements

- python 2.7
- caspo (pip install caspo)

From http://potassco.sourceforge.net:
- clasp >= 3.1
- gringo >= 4.4
- gringo python module (shipped with gringo)
On Ubuntu, installation the gringo package should satisfy those dependencies.

- NuSMV >= 2.5 - http://nusmv.fbk.eu/NuSMV/download/getting-v2.html

## Installation

	pip install https://github.com/pauleve/caspots.git

This will install necessary dependencies, except the gringo python module and
NuSMV.

A docker image is also available

	docker pull pauleve/caspots

You can then create an alias to use `caspots` command:

	alias caspots='docker run --volume "$PWD":/wd --workdir /wd pauleve/caspots'


## Usage

In the following, we assume that
* PKN.sif is the SIF description of the PKN delimiting the domain of BNs, e.g.:
    benchmarks/1/pkn1_cmpr.sif
* DATASET.csv is the MIDAS description of the multiplex dataset, e.g.,
    benchmarks/1/dataset1_cmpr_bn_1.csv
* RESULTS.csv is a CSV description of a set of Boolean Networks, as outputted by
  our python scripts.
* python is the python interpreter in version 2.7.X. On some systems, you should
  use python2.


To identify all Boolean Networks, call

	caspots identify PKN.sif DATASET.csv RESULTS.csv

By default, the identification will return the subset-minimal BNs.
Add --family all to compute _all_ the BNs.
Add --family mincard to compute the cardinal-minimal BNs.

The option --true-positives invokes a model-checker (NuSMV) to ensure that only true
positive BNs are returned. The true positive rate is then displayed.
If the PKN is not compatible with the data, the estimated difference of MSE with
minimal MSE is displayed.

The minimal estimated MSE is obtained with

	caspots mse PKN.sif DATASET.csv

The option --check-exacts invokes a model-checker (NuSMV) until it finds a BN
and a trace with the estimated MSE: in such a case, the displayed MSE is the
actual minimal MSE of the PKN with respect to the dataset.



