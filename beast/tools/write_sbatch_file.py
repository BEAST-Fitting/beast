#!/usr/bin/env python

import argparse


def write_sbatch_file(
    file_name,
    job_command,
    directory,
    modules=["module load anaconda3", "source activate bdev"],
    job_name="beast",
    stdout_file=None,
    queue="EM",
    nodes="24",
    run_time="1:00:00",
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
        file in which to save the job submission settings

    job_command : string or list of strings
        the command(s) to execute, or file(s) that have the command(s)

    directory : string
        the directory that slurm should assume you're working in

    modules : string or list of strings (default=['module load anaconda3','source activate bdev']
        modules to load before running job

    job_name : string (default='beast')
        name to give the slurm job (shows in squeue)

    stdout_file : string (default=None)
        If set, any output to the terminal will go into this file

    queue : string (default='EM')
        the queue to submit your job to
        'RM' = bridges2 regular
        'EM' = bridges2 extreme

    nodes : string (default = "24")
        the number of nodes (min for EM on Bridges2 is 24)

    run_time : string (default='1:00:00')
        Maximum run time (hh:mm:ss).  If your job is shorter, it's fine.

    array : list of two ints (default=None)
        If set, #SBATCH --array=[0]-[1] will be included. In this case, make
        sure to use "${SLURM_ARRAY_TASK_ID}" somewhere in the job command.

    """

    with open(file_name, "w") as f:

        f.write("#!/bin/bash\n")
        f.write("\n")

        f.write("#SBATCH -J " + job_name + "\n")

        if stdout_file is not None:
            f.write("#SBATCH -o " + stdout_file + "\n")

        f.write("#SBATCH -p " + queue + "\n")

        f.write("#SBATCH -n " + nodes + "\n")

        f.write("#SBATCH -t " + run_time + "\n")

        if array is not None:
            f.write("#SBATCH --array={0}-{1}\n".format(array[0], array[1]))

        f.write("\n")

        f.write("# move to appropriate directory\n")
        f.write("cd " + directory + "\n")
        f.write("\n")

        f.write("# Load any necessary modules.\n")
        f.write("# Loading modules in the script ensures a consistent environment.\n")
        if isinstance(modules, str):
            f.write(modules + "\n")
        elif isinstance(modules, list):
            for item in modules:
                f.write(item + "\n")
        f.write("\n")

        f.write("# Launch a job\n")
        if isinstance(job_command, str):
            f.write(job_command + "\n")
        elif isinstance(job_command, list):
            for item in job_command:
                f.write(item + "\n")
        f.write("\n")


if __name__ == "__main__":  # pragma: no cover
    # commandline parser
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "file_name",
        type=str,
        help="file in which to save the job submission settings",
    )
    parser.add_argument(
        "job_command",
        type=str,
        help="""the command(s) to execute, or file(s) that have the command(s),
                separated by newlines ('\n') as needed""",
    )
    parser.add_argument(
        "directory",
        type=str,
        help="the directory that slurm should assume you're working in",
    )
    parser.add_argument(
        "--modules",
        type=str,
        nargs="+",
        default=["module load anaconda3", "source activate bdev"],
        help="modules to load before running job",
    )
    parser.add_argument(
        "--job_name",
        type=str,
        default="beast",
        help="name to give the slurm job (shows in squeue)",
    )
    parser.add_argument(
        "--stdout_file",
        type=str,
        default=None,
        help="If set, any output to the terminal will go into this file",
    )
    parser.add_argument(
        "--queue",
        type=str,
        default="EM",
        help="the queue to submit your job to",
    )
    parser.add_argument(
        "--run_time",
        type=str,
        default="1:00:00",
        help="Maximum run time (hh:mm:ss).  If your job is shorter, it's fine.",
    )
    parser.add_argument(
        "--nodes",
        type=str,
        default="24",
        help="""Number of nodes to allocate to a job""",
    )
    parser.add_argument(
        "--array",
        type=int,
        nargs=2,
        default=None,
        help="""If set, #SBATCH --array=[0]-[1] will be included.  In this case,
                make sure to use "${SLURM_ARRAY_TASK_ID}" somewhere in the job
                command.""",
    )

    args = parser.parse_args()

    write_sbatch_file(
        args.file_name,
        args.job_command.split(r"\n"),
        args.directory,
        modules=args.modules,
        job_name=args.job_name,
        stdout_file=args.stdout_file,
        queue=args.queue,
        run_time=args.run_time,
        nodes=args.nodes,
        array=args.array,
    )

    # print help if no arguments
    if not any(vars(args).values()):
        parser.print_help()
