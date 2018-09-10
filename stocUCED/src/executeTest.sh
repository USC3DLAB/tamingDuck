#!/bin/bash

echo Cleaning...
rm *sol *log *lp
rm -r DDD* DDS* SDS*

echo Initiating test...

UCResName=10.0
EDResName=2.5
declare -a renNames=("1.00" "2.00" "3.00")
declare -a settingNames=("DDD" "DDS" "SDS")
declare -a inputNames=("deterministic deterministic deterministic" "deterministic deterministic stochastic" "stochastic deterministic stochastic")

length=${#inputNames[@]}
for ((i=0;i<$length;i++)); do
	settingName=${settingNames[$i]};
	inputName=${inputNames[$i]}
	for renName in ${renNames[@]}; do
		fileName=${settingName}_Ren${renName}_Res${UCResName}_${EDResName}

		sed -i "s/^renewableCoef.*/renewableCoef "$renName"/g" ../datasets/runParameters.txt

		echo Running $fileName...
		#echo $inputName---$renName---$fileName
		mkdir $fileName

		./run ../datasets/ ./sd/ ./ 3d-nrel118-feb-curtail -setting $inputName > log.log
		mv *sol *log *lp $fileName
	done
done

echo Completed.

