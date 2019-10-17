#!/usr/bin/env python

def write_sbatch_file(
        file_name,
        job_command,
        directory,
        modules=['module load anaconda3','source activate bdev'],
        job_name='beast',
        stdout_file=None,
        egress=False,
        queue='LM',
        run_time='1:00:00',
        mem='128GB',
        array=None,
):
    """
    Write out a file to submit to slurm using sbatch:
    > sbatch [file_name]

    More info here:
    https://portal.tacc.utexas.edu/archives/stampede#slurm-job-control

    Parameters
    ----------
    file_name : string
        file in which to save the job submission stuff

    job_command : string or list of strings
        the command(s) to execute, or file(s) that have the command(s)

    directory : string
        the directory that slurm should assume you're working in

    modules :  string or list of strings (default=['module load anaconda3','source activate bdev']
        modules to load before running job

    job_name : string (default='beast')
        name to give the slurm job (shows in squeue)

    stdout_file : string (default=None)
        If set, any output to the terminal will go into this file

    egress : boolean (default=False)
        If your job will need connections to external sites (e.g., for downloading
        stellar tracks), set this to True.

    queue : string (default='LM')
        the queue to submit your job to
        'RM' = bridges regular
        'LM' = bridges large

    run_time : string (default='1:00:00')
        Maximum run time (hh:mm:ss).  If your job is shorter, it's fine.

    mem : string (default='128GB')
        For bridges large, the memory to allocate to the job.  If your job uses
        less memory than this, the full amount will still be charged.

    array : list of two ints (default=None)
        If set, #SBATCH --array=[0]-[1] will be included.  In this case, make
        sure to use "${SLURM_ARRAY_TASK_ID}" somewhere in the job command.

    """

    with open(file_name,'w') as f:

        f.write('#!/bin/bash \n')
        f.write('\n')

        f.write('#SBATCH -J '+job_name+'\n')

        if stdout_file is not None:
            f.write('#SBATCH -o '+stdout_file+'\n')

        if egress == True:
            f.write('#SBATCH -C EGRESS \n')

        f.write('#SBATCH -p '+queue+'\n')
        f.write('#SBATCH -t '+run_time+'\n')
        f.write('#SBATCH --mem '+mem+'\n')

        if array is not None:
            f.write('#SBATCH --array={0}-{1}\n'.format(array[0], array[1]))

        f.write('\n')

        f.write('# move to appropriate directory\n')
        f.write('cd '+directory+'\n')
        f.write('\n')
        
        f.write('# Load any necessary modules.\n')
        f.write('# Loading modules in the script ensures a consistent environment.\n')
        if type(modules) == str:
            f.write(modules+'\n')
        elif type(modules) == list:
            for item in modules:
                f.write(item+'\n')
        f.write('\n')
        
        f.write('# Launch a job\n')
        if type(job_command) == str:
            f.write(job_command+'\n')
        elif type(job_command) == list:
            for item in job_command:
                f.write(item+'\n')
                
        f.write('\n')
