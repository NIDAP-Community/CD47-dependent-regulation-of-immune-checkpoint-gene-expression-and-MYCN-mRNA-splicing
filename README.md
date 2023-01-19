# CD47-dependent-regulation-of-immune-checkpoint-gene-expression-and-MYCN-mRNA-splicing

This code accompanies the paper entitled: CD47-dependent regulation of immune checkpoint gene expression and MYCN mRNA splicing in murine CD8 and Jurkat T cells


To reproduce these results, follow these steps:

1.  Clone this GitHub repo (i.e. the page you are on):
    * ```git clone https://github.com/NIDAP-Community/CD47-dependent-regulation-of-immune-checkpoint-gene-expression-and-MYCN-mRNA-splicing.git```

2.  The input files for this pipeline will be available upon request. Please reach out to the authors before continue to following steps

3.  Install docker and build the docker container:
    * Navigate to the cloned repository directory. 
    * Move to the ./Docker_file/ directory of this repo

4.  Build the container:
    * ```docker build --tag ccbr1188_r3_final .```

5.  Navidate to the cloned repository directory, Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container:
    * ```docker run -ti -v $(pwd)/src:/tmp ccbr1188_r3_final```

6.  Activate the conda environment in the docker container:
    * ```conda activate ccbr1188_r3_final```

5.  Run the following code
    * ```cd /tmp```
    * ```bash run_pipeline.sh```

