import os
import re
from glob import glob
import pandas as pd
from collections import OrderedDict

configfile: "config.yml"

# -- 0.0 Parse config.yml files into paths and parameters of the pipeline
# For compliling merida
src_path = config["src_path"]
merida_repo = config["merida_repo"]
bin_path = config["bin_path"]
cplex_path = config["cplex_path"]

# Paths to input and output data
dataset = config["dataset"]
rawdata = config["rawdata"]
metadata = config["metadata"]
procdata = config["procdata"]
results = config["results"]
# add trailing slash for MERIDA_ILP compatibility
out_dir = os.path.join(results, "")

# Analysis name and paths to feature and response files
analysis_name = config["analysis_name"]
feature_matrix = os.path.join(procdata, config["feature_matrix"])
response_vector = os.path.join(procdata, config["response_vector"])

# Input data dimensions for file names
samples = config["num_samples"]
features = config["num_features"]

# Cross validation and a priori feature selection
cv = config["cross_validation"]

preselected_features = config["preselected_features"]
if preselected_features is not None:
    preselected_features = os.path.join(procdata, preselected_features)

# -- 0.1 Make MERIDA input files for to grid search hyperparameters
threshold = config["threshold"]
M_range = config["M"]
v_function = config["v"]

M_list = list(range(M_range["start"], M_range["stop"] + 1, M_range["step"]))

# -- 0.2 Build input and output file names
input_files = glob(f"{procdata}/{dataset}_*_Input_Matrix.txt")
drug_names = [re.sub(f"{procdata}/{dataset}_|_Input_Matrix.txt", "", f) for f in input_files]

## FIXME:: Remove this before a real run! Just a hack to run less drugs
drug_names = ["Sorafenib"]

nFeature = config["nFeature"]
nSolution = config["nSolution"]
merida_output = expand(
    f"{results}/{dataset}_{{drug}}.mRMRe_nsol{{nSolution}}_nfeat{{nFeature}}.MERIDA_M{{m}}_v{{v}}_cv{{cv}}.tar.gz",
    nSolution=nSolution, nFeature=nFeature, drug=drug_names, m=M_list, v=v_function, cv=cv
    )

# -- 1. Rule to gather results from other rules
rule all:
    input:
        merida_output


# -- 2.0 Project tool installation
rule download_and_compile_merida:
    params:
        repo=merida_repo,
        cplex=cplex_path
    output:
        src=directory(os.path.join(src_path)),
        bin=os.path.join(bin_path, "MERIDA_ILP")
    shell:
        """
        # Clone the repo and rename to to src_path
        git clone {params.repo} {output.src}
        # Make a build directory
        cd {output.src}
        mkdir build
        cd build
        # Make the MERIDA binary and add it to $HOME/bin
        cmake -DCPLEX_PATH={params.cplex} ..
        make
        ln -s `pwd`/MERIDA_ILP {output.bin}
        """


# NOTE: nSolution is number of solutions to ensemble and will output
#   nSolution x nFeature total features in the output matrix, not 3 10-feature
#   output matrices
rule mrmr_feature_selection:
    input:
        f"{procdata}/{dataset}_{{drug}}_Input_Matrix.txt"
    params:
        nSolution=nSolution,
        nFeature=nFeature
    output:
        # Inside an expand statment, a wildcard must be specified with {{
        # Inside an f-string, the first { is removed
        # Therefore: { is interpolated, {{ is literal {, {{{{ is literal {{
        expand(
            f"{procdata}/{analysis_name}/{dataset}_{{{{drug}}}}.mRMRe_nsol{{nSolution}}_nfeat{{nFeature}}.txt",
            nSolution=nSolution, nFeature=nFeature)
    script:
        "scripts/mRMReFeatureSelection.R"



# -- 4.0 Build configuration files required as MERIDA_ILP input
rule build_merida_input_files:
    input:
        feature_matrix=expand(
            f"{procdata}/{analysis_name}/{dataset}_{{{{drug}}}}.mRMRe_nsol{{nSolution}}_nfeat{{nFeature}}.txt",
            nSolution=nSolution, nFeature=nFeature
        ),
        response_vector=f"{procdata}/{dataset}_{{drug}}_AAC_file.txt"
    params:
        M_list=M_list,
        v_function=v_function,
        threshold=threshold,
        cv=cv,
        out_dir=directory(out_dir),
        container=config["storage_container"]
    log:
        "logs/build_merida_input_files_{drug}.log"
    output:
        expand(f"{procdata}/{analysis_name}/{dataset}_{{{{drug}}}}.mRMRe_nsol{{nSolution}}_nfeat{{nFeature}}.MERIDA_M{{m}}_v{{v}}_cv{{cv}}.txt",
            nSolution=nSolution, nFeature=nFeature, m=M_list, v=v_function, cv=cv)
    run:
        print(os.listdir())
        conf = OrderedDict([
            ("File", None),
            ("Directory", None),
            ("M1", None),
            ("M2", None),
            ("L1", None),
            ("L2", None),
            ("WeightFunction", None),
            ("ErrorFunction", None),
            ("Threshold", None),
            ("IC50ValueFile", None),
            ("CVMode", None),
            ("alpha", None),
            ("preselected_features", None)
        ])
        grid = [(m, v) for m in params.M_list for v in params.v_function]
        for fl, p in zip(output, grid):
            # NOTE: MERIDA requires tab delimitation or will fail to read config
            conf["File"] = f"File\t{input.feature_matrix}"
            result_dir = os.path.join(
                params.container,
                params.out_dir,
                os.path.splitext(os.path.basename(fl))[0]
            )
            conf["Directory"] = f"Directory\t{result_dir}/"
            os.makedirs(result_dir, mode=0o774, exist_ok=True)
            conf["M1"] = f"M1\t{p[0]}"
            conf["M2"] = f"M2\t{p[0]}"
            conf["WeightFunction"] = f"WeightFunction\t{p[1]}"
            conf["IC50ValueFile"] = f"IC50ValueFile\t{input.response_vector}"
            conf["Threshold"] = f"Threshold\t{params.threshold}"
            # Add additional configuration when doing cross validation
            if isinstance(params.cv, int) and params.cv > 0:
                conf["L1"] = f"L1\t0"
                conf["L2"] = f"L2\t0"
                conf["ErrorFunction"] = f"ErrorFunction\tmisclassification"
                conf["CVMode"] = f"CVMode\tfold"
                conf["alpha"] = f"alpha\t0"
            # Add a priori feature selection, if specified in config.yml
            # if input.preselected is not "" and os.path.exists(input.preselected):
            #     conf["preselected_features"] = \
            #         f"preselected_features\t{input.preselected}"
            # Write the configuration files to the procdata directory
            conf_file = "\n".join([
                *[s for s in conf.values() if s is not None],
                ""
            ])
            with open(fl, "w+") as f:
                f.write(conf_file)


# -- 5.0 Run MERIDA_ILP using config files
rule run_merida:
    input:
        feature_matrix=f"{procdata}/{analysis_name}/{dataset}_{{drug}}.mRMRe_nsol{{nSolution}}_nfeat{{nFeature}}.txt",
        response_vector=f"{procdata}/{dataset}_{{drug}}_AAC_file.txt",
        file_name=f"{procdata}/{analysis_name}/{dataset}_{{drug}}.mRMRe_nsol{{nSolution}}_nfeat{{nFeature}}.MERIDA_M{{m}}_v{{v}}_cv{{cv}}.txt"
    log:
        f"logs/{analysis_name}_{dataset}_{{drug}}.mRMRe_nsol{{nSolution}}_nfeat{{nFeature}}.MERIDA_M{{m}}_v{{v}}_cv{{cv}}.log"
    params:
        cv=cv,
        out_dir=os.path.join(config["storage_container"], out_dir),
        runtime=config["runtime_per_job"],
        cpu=config["cpu_per_job"],
        mem=config["mem_per_job"],
        slurm_output=config["slurm_output"],
        features=features,
        samples=samples
    singularity:
        "bhklabbatch.azurecr.io/aks/aks-snakemake-cplex-merida-r"
    output:
        f"{results}/{dataset}_{{drug}}.mRMRe_nsol{{nSolution}}_nfeat{{nFeature}}.MERIDA_M{{m}}_v{{v}}_cv{{cv}}.tar.gz"
    shell:
        # Double curly brace is literal curly brace inside f-string, single
        #  curly brace is interpolated python variable
        f'''
        echo $(whoami) >> "{{log}}"
        echo $PATH >> "{{log}}"

        if [ {{params.cv}} = "no" ]
        then
            MERIDA_ILP "{{input.file_name}}" {{params.cv}} &>> "{{log}}"
        else
            MERIDA_ILP "{{input.file_name}}" yes {{params.cv}} &>> "{{log}}"
        fi

        ls -lah -R >> "{{log}}"
        tar -cvf "{{output}}" "{{params.out_dir}}/{dataset}_{{wildcards.drug}}.mRMRe_nsol{{wildcards.nSolution}}_nfeat{{wildcards.nFeature}}.MERIDA_M{{wildcards.m}}_v{{wildcards.v}}_cv{{wildcards.cv}}"
        '''