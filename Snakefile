import os
import re
from glob import glob
import pandas as pd


configfile: "config.yml"


# -- 0.0 Configuration variables
src_path = config["src_path"]
merida_repo = config["merida_repo"]
bin_path = config["bin_path"]
cplex_path = config["cplex_path"]

rawdata = config["rawdata"]
metadata = config["metadata"]
procdata = config["procdata"]
results = config["results"]
# add trailing slash for MERIDA_ILP compatibility
out_dir = os.path.join(results, "")

analysis_name = config["analysis_name"]
feature_matrix = os.path.join(procdata, config["feature_matrix"])
response_vector = os.path.join(procdata, config["response_vector"])
cv = config["cross_validation"]


# -- 0.1 Make MERIDA input files for to grid search hyperparameterss
threshold = config["threshold"]
M_range = config["M"]
v_function = config["v"]

M_list = list(range(M_range["start"], M_range["stop"] + 1, M_range["step"]))
param_grid = [(m, v) for m in M_list for v in v_function]

# -- 0.2 Build input and output file names
merida_files = [
    f"procdata/merida_input/{analysis_name}_M{m}_v{v}_cv{cv}_{os.path.basename(feature_matrix)}" 
        for m, v in param_grid
]
merida_results = [
    f"results/{analysis_name}_M{m}_v{v}_cv{cv}_{os.path.basename(feature_matrix)}" 
        for m, v in param_grid
]

# -- 0.2 Dimensions of input matrix
feature_df = pd.read_csv(feature_matrix, sep="")
features, samples = feature_df.shape


# -- 1. Rule to gather results from other rules
rule all:
    input:
        merida_results


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


# -- 3.0 Download the Project Data
## TODO:: Add download and preprocessing scripts
# pset_list = config["pset_list"]
# rule download_psets:
#     input:
#         os.path.join(bin_path, "MERIDA_ILP")
#     params:
#         psets=


# -- 4.0 Build configuration files required as MERIDA_ILP input
rule build_merida_input_files:
    params:
        grid=param_grid,
        feature_matrix=feature_matrix,
        response_vector=response_vector,
        out_dir=out_dir,
        threshold=threshold,
        cv=cv
    output:
        merida_files
    run:
        for fl, p in zip(output, params.grid):
            # NOTE: MERIDA requires tab delimitation or will fail to read config
            config = (
                f"File\t{params.feature_matrix}\n"
                f"Directory\t{params.out_dir}\n"
                f"M1\t{p[0]}\n"
                f"M2\t{p[0]}\n"
                f"WeightFunction\t{p[1]}\n"
                f"IC50ValueFile\t{params.response_vector}\n"
                f"Threshold\t{params.threshold}\n" # fails without trailing \n
            )
            # Add additional configuration when doing cross validation
            if isinstance(params.cv, int) and params.cv > 0:
                config = (
                    f"{config}"
                    f"L1\t0\n"
                    f"L2\t0\n"
                    f"ErrorFunction\tmisclassification\n"
                    f"CVMode\tfold\n"
                    f"alpha\t0\n"
                )  
            with open(fl, "w+") as f:
                f.write(config)


# -- 5.0 Run MERIDA_ILP using config files
rule run_merida:
    input:
        "procdata/merida_input/roche_M{m}_v{v}_{file}.txt"
    params:
        cv=cv,
        jobname="roche_M{m}_v{v}_{file}.job",
        runtime="2:0:0",
        features=2816,
        samples=9
    resources:
        time="6:0:0"
    output:
        f"{results}/roche_M{{m}}_v{{v}}_{{file}}.txt"
    shell:
        # Double curly brace is literal curly brance inside f-string, single
        #  curly brace is interpolated variable
        f"""
        source ~/.bashrc
        conda activate merida
        if [ param.cv = "no" ]
        then
            MERIDA_ILP {{input}} {{params.cv}} > \
                results/roche_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}.log
        else
            MERIDA_ILP {{input}} yes {{params.cv}} > \
                results/roche_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}.log
        mv {results}/Results_M_{{m}}_Feat_{{params.features}}Sample_{{params.samples}}.txt \
            {{output}} || echo "failed" > {{output}}
        """
