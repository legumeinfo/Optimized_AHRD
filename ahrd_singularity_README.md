Create/modify nextflow.config:
This file must be in the directory the pipeline is launched from. Bind mounts must be provided in “runOptions” for the database/s directory, as well as any other directories involved in the command (including the working directory). Memory overrides are important to ensure resources for the computationally demanding processes. 

Note: The example .config file in this repository is configured for Ceres. In an environment using Singularity, the name of that section should be changed from "apptainer" to "singularity" in the config, and the flag "-with-apptainer" in the run command changed to "-with-singularity". Additionally, different options under "process" in the .config may not be necessary. On dal, for example, only the line
env.PYTHONPATH = '/usr/local/lib64/python3.9/site-packages'
is required there.

Load modules (names also dependent on environment):
```
module load nextflow
module load apptainer
```
This will load the apptainer/1.3.1 and nextflow/24.10.1 defaults. 

Example command:
```
nextflow run /project/legume_project/common_data/ahrd/ahrd_singularity.nf -with-apptainer /project/legume_project/singularity_images/AHRD.sif --chunksize 500 --threads_per_process 4 --total_threads 200 --input_fasta pissa.Cameor.gnm1.ann1.7SZR.protein_primary.faa --outdir pissa_test --databases databases.csv --gaf /reference/data/Uniprot/2025-07-07/goa_uniprot_all.gaf.gz
```
Arguments:
	
--input_fasta (required)
The proteins to analyze; the .faa extension is not required, though diamond will exit with an error message if a nucleotide file is provided.

--databases (required)
3 column .csv file including database name, full path to database, and diamond index file, respectively. Any number of desired databases >= 1 is acceptable, though the names must match one of the following strings: uniref90, swissprot, trembl, medtr_lis, glyma_refseq
Uncompressed and gzipped .fasta files are both acceptable.

--gaf (required for GO terms)

--chunksize <500>
	number of records (- 1) each “chunk” of the input_fasta will include. 

--total_threads <128>
	
--threads_per_process <4>
	total_threads / threads_per_process will equal the processes created, not accounting for the controller.

--outdir <working directory>
	Directory name for output folder. Relative path is acceptable.

Notes:

•	The pipeline may be restarted with progress retained after alignment and/or interproscan processes have completed on all .fasta “chunks” and been concatenated, as well as when the ahrd_config.yml is generated.

•	Because the pipeline will not overwrite these files (or the .fasta chunks), a new out_dir should be provided for different runs (when the ‘resume’ functionality is not required/desired), or the contents of the out_dir deleted prior. There is a warning message to this effect when reusing an existing out_dir.

•	The work directory “work” persists as a general nextflow feature. It is advised to delete this periodically in the absence of special cause to retain it. 

•	For reference: An input of ~89k Pisum sativum proteins completes in ~11 hours with total_threads = 400 and threads_per_process = 4. Expect the AHRD portion to run for 1-2 hours. Completion time will vary with queue load. 


