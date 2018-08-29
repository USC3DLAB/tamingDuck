#!/bin/bash

echo Cleaning...
rm *sol *log *lp
rm -r DDD DDS SDS

echo Initiating test...

echo Running DDD...
mkdir DDD
./run ../datasets/ ./sd/ ./ 3d-nrel118-feb-curtail -setting deterministic deterministic deterministic > log.log
mv *sol *log *lp DDD

echo Running DDS...
mkdir DDS
./run ../datasets/ ./sd/ ./ 3d-nrel118-feb-curtail -setting deterministic deterministic stochastic > log.log
mv *sol *log *lp DDS

echo Running SDS...
mkdir SDS
./run ../datasets/ ./sd/ ./ 3d-nrel118-feb-curtail -setting stochastic deterministic stochastic > log.log
mv *sol *log *lp SDS

echo Completed.

