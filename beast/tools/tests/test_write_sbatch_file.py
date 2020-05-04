from tempfile import NamedTemporaryFile
from beast.tools.write_sbatch_file import write_sbatch_file

def test_sbatch_file():
    temp_file = NamedTemporaryFile(suffix=".script")
    write_sbatch_file(temp_file.name, "./mastergrid_LMC/model_batch_jobs/create_physicsmodel_\"${SLURM_ARRAY_TASK_ID}\".job",
    "/pylon5/as5pi7p/lhagen", modules=["module load anaconda3", "source activate bdev"], job_name="LMCgrid", egress=True, queue="LM",
    stdout_file="/pylon5/as5pi7p/lhagen/mastergrid_LMC/model_batch_jobs/logs/%A_%a.out", run_time="35:00:00", mem="570GB", array=[1,9])

    file = open(temp_file.name)
    content = file.read()

    expected = ("""#!/bin/bash

#SBATCH -J LMCgrid
#SBATCH -o /pylon5/as5pi7p/lhagen/mastergrid_LMC/model_batch_jobs/logs/%A_%a.out
#SBATCH -C EGRESS
#SBATCH -p LM
#SBATCH -t 35:00:00
#SBATCH --mem 570GB
#SBATCH --array=1-9

# move to appropriate directory
cd /pylon5/as5pi7p/lhagen

# Load any necessary modules.
# Loading modules in the script ensures a consistent environment.
module load anaconda3
source activate bdev

# Launch a job
./mastergrid_LMC/model_batch_jobs/create_physicsmodel_\"${SLURM_ARRAY_TASK_ID}\".job

""")

    assert content == expected, "The created sbatch file does not have the expected output."
