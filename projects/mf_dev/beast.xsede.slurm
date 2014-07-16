#!/bin/csh
#-------------------------------------------------------
#SBATCH -J Beast          # Job name
#SBATCH -n 4              # Total number of tasks
#SBATCH -N 2              # Total number of nodes
#SBATCH -p normal         # Queue name
#SBATCH -o AP4603.out     # Name of stdout file
#SBATCH -t 10:00:00       # Run time (hh:mm:ss)
#SBATCH -A TG-AST130057   # Account Assignment
#------------------------------------------------------

echo StartTime: `date`

module load launcher
setenv EXECUTABLE     $TACC_LAUNCHER_DIR/init_launcher 
setenv CONTROL_FILE   beast.in
setenv WORKDIR        $WORK/phat/beast

# 
# # Variable description:
# #
# #  EXECUTABLE          = full path to the job launcher executable
# #  CONTROL_FILE        = text input file which specifies
# #                        executable for each process
# #                        (should be located in WORKDIR)
# #  PHI_CONTROL_FILE    = text input file which specifies executable
# #                        for each process to be run on Intel Xeon Phi
# #                        (should be located in PHI_WORKDIR)
# #  WORKDIR             = location of working directory
# #  PHI_WORKDIR         = location of working directory for Intel Xeon Phi jobs
# #  TACC_LAUNCHER_SCHED = scheduling method for lines in CONTROL_FILE
# #                        options (k=process, n=num. lines, p=num. procs):
# #                          - interleaved (default): 
# #                              process k executes every k+nth line
# #                          - block:
# #                              process k executes lines [ k(n/p)+1 , (k+1)(n/p) ]
# #                          - dynamic:
# #                              process k executes first available unclaimed line
# #
# #      <------ End Setup Parameters ------>
# #--------------------------------------------------------
# #--------------------------------------------------------
#

if ( ! -e $WORKDIR ) then
        echo " "
	echo "Error: unable to change to working directory."
	echo "       $WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $EXECUTABLE ) then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
endif

if ( ! -f $WORKDIR/$CONTROL_FILE ) then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
endif

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

echo " EXECUTING: $TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE"
$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE

echo EndTime: `date`

echo " "
echo " Job Complete"
echo " "
