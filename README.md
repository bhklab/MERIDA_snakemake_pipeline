# MERIDA Snakemake pipeline

## Snakemake

This pipeline leverages the `snakemake` Python package for workflow management. 
As a result the pipeline and its dependencies are easily
installed from this repository, allowing quick setup, configuration and 
deployment.

For more information on Snakemake, please see: 
https://snakemake.readthedocs.io/en/stable/.

## Requirements

Dependency management for this pipeline is handled via `conda` for Python 
and `renv` for R. To get started with setup you can install
miniconda3 using the instructions available here: 
https://docs.conda.io/en/latest/miniconda.html.

Alternatively you can install it directly from CRAN
as described here: https://cran.r-project.org/.

## Setting Up Your Software Environment

The first step to deploying an analysis pipeline is to install the various
software packages it depends on. We have included the `envs/merida.yaml` and 
`renv.lock` files here to easily accomplish this.

All commands should be executed from the top level directory of this
repository.

### IBM ILOG CPLEX

This pipeline requires the use of IBM ILOG CPLEX logical solver, which
is closed source software. An academic license can be obtained from IBM
to use this software for free. Follow the instructions on 
[Compute Canada](https://docs.computecanada.ca/wiki/CPLEX/en).

We recommend also setting the CPLEX_HOME environmental variable to make it
easy for compiled code to find your installation. You can set it with:

**NOTE:** replace /path/to/CPLEX_StudioXYZ with your CPLEX path and version 
before running! A back-up of your .bashrc will be saved in ~/.bashrc.bak if
something goes wrong with this script.

```
# Delete existing the line containg CPLEX_HOME in-place, if it exists
sed -i.bak '/CPLEX_HOME=.*/d' ~/.bashrc
# Add new CPLEX_HOME
echo "export CPLEX_HOME=/path/to/CPLEX_StudioXYZ" >> ~/.bashrc
# Reload .bashrc
source ~/.bashrc
# Check that the variable exists
echo $CPLEX_HOME
```

This will set an environmental variable in your ~/.bashrc file, then reload
your enviroment variables from .bashrc and print the path. If the command 
doesn't print your the path you configured it has failed and may need to be 
fixed manually.

If you are on Windows, the environmental variable should be created 
automatically by the CPLEX installer. Please note this pipeline has not been 
tested on Windows and may not work correctly.

### Python and System Dependencies

Conda can be used to install all Python and most OS system dependencies
using:

`conda env create --file envs/merida.yml`

This will take some time to run as it gathers and installs the correct
package versions. The environent it creates should be called `merida`.

If it is not automatically activated after installation please run 
`conda activate merida` before proceeding to the next step.

### R Dependencies

The `renv` package can be used to install all R dependencies (both CRAN and
Bioconductor). R version 4.1 and `renv` are included as dependencies in the 
`merida.yml` file and should be installed automatically when setting up your 
conda environment. If R is not installed, you can install it via conda using 
the command: `conda -c conda-forge r-base==4.1`.

To initialize this project with renv run:

`Rscript -e 'library(renv); renv::init()'`

If you wish to isolate the R dependencies from your Conda environment R 
libraries, you can use this command instead:

`Rscript -e 'library(renv); renv::isolate(); renv::init(bare=TRUE)'`

If intialization doesn't trigger dependency installation, you can do so manually using:

`Rscript -e 'renv::restore()'`

For more information on renv and how it can be used to manage dependencies in
your project, please see: https://rstudio.github.io/renv/articles/renv.html.

If the `renv` commands fail, you may need to install `renv` manually with
`Rscript -e 'install.packages("renv")'` then retry the above commands.

## Configuring the Pipeline

This pipeline assumes the following directory structure:

```
.
├── envs
├── src
├── bin
├── metadata
├── rawdata
├── procdata
├── results
└── scripts
```

Please at minimum create the `rawdata` and `metadata` directories. 
The remainder will be created automatically during pipeline execution.

### config.yaml

This file hold the relevant pipeline documentation. Here you can specify the paths
to all the parameters for your current pipeline use case. Documentation is provided
in the `config.yaml` file on what each field should contain.

## Using the Pipeline

Make sure you have set all values in `config.yml` before trying to run the
pipeline! This file contains all user configuration necessary to get the
pipleine running with your specific data and project structure.

### Compiling MERIDA

`snakemake --cores 2 download_and_compile_merida`

### Downloading the Data

`snakemake --cores`

### Extracting Feature Matrices

`snakemake --cores`

### Configuring Hyperparemeter Search Space

`snakemake --cores`
