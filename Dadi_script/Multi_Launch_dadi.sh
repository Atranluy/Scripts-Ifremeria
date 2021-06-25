#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=01:00:00
#$ -V

##### script from Adrien Tran Lu Y
######### change header depending of type of HPC you will use (grid cluster engine or slurm)
######### im using pyenv environnement where dadi is installed but you can call the dadi version from module load or something else depending on your HPC
######### dadi v2.1 from https://gitlab.mbb.univ-montp2.fr/khalid/dadi/-/tree/master  see there to install dadi and pyenv to install it within a specific python3 environnement


FS=Spectrum_FS_nosingleton ###### Spectrum file


#### set name of each model do you want to launch 
#### it will create a single bash file for each run (1 to 11) for each model and run it independently 
#### it check before if the iteration is already done within the folder (current folder /IM/IM_1.txt) for exemple. 
for model in SC2N2m AM2N2m;  
do

for i in `seq 1 11`;
do
echo ${i}

if [ -d "$model" ]; 
then
	if [ -f "${model}/${model}_${i}_.txt" ]; 
	then
		echo "already done"
	else
		printf "
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=480:00:00
#$ -V
source activate dadi_sim
python3 call_all_model_NA_change.py $model $FS $i " > ${model}_${i}_.sh 
		qsub ${model}_${i}_.sh 
	fi
else
printf "
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=480:00:00
#$ -V
source activate dadi_sim
python3 call_all_model_NA_change.py $model $FS $i " > ${model}_${i}_.sh 
qsub ${model}_${i}_.sh 
fi
done
done

