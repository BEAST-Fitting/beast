from beast.tools.write_sbatch_file import write_sbatch_file

def test_sbatch_file():
    write_sbatch_file("create_LMC_mastergrid.script", "./mastergrid_LMC/model_batch_jobs/create_physicsmodel_\"${SLURM_ARRAY_TASK_ID}\".job",
    "/pylon5/as5pi7p/lhagen", modules=["module load anaconda3", "source activate bdev"], job_name="LMCgrid", egress=True, queue="LM",
    stdout_file="/pylon5/as5pi7p/lhagen/mastergrid_LMC/model_batch_jobs/logs/%A_%a.out", run_time="35:00:00", mem="570GB", array=[1,9])

    file = open("create_LMC_mastergrid.script")
    content = file.read()
    expected = "#!/bin/bash\n\n#SBATCH -J LMCgrid\n#SBATCH -o /pylon5/as5pi7p/lhagen/mastergrid_LMC/model_batch_jobs/logs/%A_%a.out\n#SBATCH -C EGRESS\n#SBATCH -p LM\n#SBATCH -t 35:00:00\n#SBATCH --mem 570GB\n#SBATCH --array=1-9\n\n# move to appropriate directory\ncd /pylon5/as5pi7p/lhagen\n\n# Load any necessary modules.\n# Loading modules in the script ensures a consistent environment.\nmodule load anaconda3\nsource activate bdev\n\n# Launch a job\n./mastergrid_LMC/model_batch_jobs/create_physicsmodel_\"${SLURM_ARRAY_TASK_ID}\".job\n\n"

    assert content == expected, "The created sbatch file does not have the expected output."
