#$ -N RunAll 

#$ -S /bin/bash

# Archivo de salida
#$ -o .../output/output_all_$JOBID.txt

# Archivo de salida de error
#$ -e .../error/error_all_$JOBID.txt

# JOB_ID=$SGE_TASK_ID

echo " starting job at: `date`  "
echo " executing job: " $JOBID  " at nodo"  bin/hostname 
echo "------------ Nodo ------------------"

/bin/hostname 

echo "---------- Current JOB_dir-------------"

pwd

echo "------------- Trabajo -----------------"

# echo $JOBID

# Job execution

cd .../Flow/restore_dca
source config_1.sh

root -b -l <<-EOF
	gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
	mpdloadlibs()
	.L restore_dca.c+
	.L .../Flow/get_dca/get_dca.cxx+
	.L .../Flow/get_dca/get_fit.cxx+
	.L .../Flow/get_dca/MakeFitDCA.cxx+
	.L .../Flow/get_centrality/get_multiplicity.cxx+
	restore_dca(".../mpddstold$JOBID.root",".../mpddstnew$JOBID.root")
	get_dca(".../mpddstnew$JOBID.root",".../getdca$JOBID.root")
	get_fit(".../getdca$JOBID.root",".../getfit$JOBID.root")
	MakeFitDCA(".../getfit$JOBID.root",".../makefitdca$JOBID.root")
	get_multiplicity(".../mpddstnew$JOBID.root",".../getmultiplicity$JOBID.root",".../makefitdca$JOBID.root")
	.q
EOF
rm .../mpddstnew$JOBID.root
echo " Ending job at: `date`"
