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
param_grid = [(m, v) for m in M_list for v in v_function]

# -- 0.2 Build input and output file names
merida_files = [
    f"{procdata}/merida_input/{analysis_name}_M{m}_v{v}_cv{cv}_{os.path.basename(feature_matrix)}" 
        for m, v in param_grid
]
merida_results = [
    f"{results}/{analysis_name}_M{m}_v{v}_cv{cv}_{os.path.splitext(os.path.basename(feature_matrix))[0]}/{analysis_name}_M{m}_v{v}_cv{cv}_{os.path.basename(feature_matrix)}" 
        for m, v in param_grid
]


# -- 0.2 Dimensions of input matrix
feature_df = pd.read_csv(feature_matrix, sep=" ")
samples, features = feature_df.shape


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
        cv=cv,
        preselected=preselected_features
    output:
        merida_files
    run:
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
        for fl, p in zip(output, params.grid):
            # NOTE: MERIDA requires tab delimitation or will fail to read config
            conf["File"] = f"File\t{params.feature_matrix}"
            result_dir = os.path.join(
                params.out_dir,
                os.path.splitext(os.path.basename(fl))[0]
            )
            conf["Directory"] = f"Directory\t{result_dir}/"
            os.makedirs(result_dir, mode=0o774, exist_ok=True)
            conf["M1"] = f"M1\t{p[0]}"
            conf["M2"] = f"M2\t{p[0]}"
            conf["WeightFunction"] = f"WeightFunction\t{p[1]}"
            conf["IC50ValueFile"] = f"IC50ValueFile\t{params.response_vector}"
            conf["Threshold"] = f"Threshold\t{params.threshold}"
            # Add additional configuration when doing cross validation
            if isinstance(params.cv, int) and params.cv > 0:
                conf["L1"] = f"L1\t0"
                conf["L2"] = f"L2\t0"
                conf["ErrorFunction"] = f"ErrorFunction\tmisclassification"
                conf["CVMode"] = f"CVMode\tfold"
                conf["alpha"] = f"alpha\t0"
            # Add a priori feature selection, if specified in config.yml
            if params.preselected is not None and os.path.exists(params.preselected):
                conf["preselected_features"] = \
                    f"preselected_features\t{params.preselected}"
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
        f"procdata/merida_input/{analysis_name}_M{{m}}_v{{v}}_{{file}}.txt"
    log:
        f"logs/{analysis_name}_M{{m}}_v{{v}}_{{file}}.log"
    params:
        cv=cv,
        jobname=f"{analysis_name}_M{{m}}_v{{v}}_{{file}}.job",
        runtime=config["runtime_per_job"],
        cpu=config["cpu_per_job"],
        mem=config["mem_per_job"],
        features=features,
        samples=samples
    resources:
        time="6:0:0"
    output:
        f"{results}/{analysis_name}_M{{m}}_v{{v}}_{{file}}/{analysis_name}_M{{m}}_v{{v}}_{{file}}.txt"
    shell:
        # Double curly brace is literal curly brance inside f-string, single
        #  curly brace is interpolated variable
        f"""
        set +u
        source ~/.bashrc
        conda activate merida

        MERIDA_ILP {{input}} yes {{params.cv}} &> logs/merida_errors.log
        """

# Run MERIDA_ILP
# MERIDA_ILP {{input}} yes {{params.cv}}

# > {results}/{analysis_name}_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}/{analysis_name}_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}.log
        
#if [ {{params.cv}} = "no" ]
#then
#    echo "No VCV"
#    MERIDA_ILP {{input}} {{params.cv}} > {results}/{analysis_name}_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}/{analysis_name}_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}.log || touch {{output}}
#els
#fi
# Fix file names after execution
# mv {results}/Results_M_{{wildcards.m}}_Feat_{{params.features}}Sample_{{params.samples}}.txt {{output}} || echo "failed" > {{output}}
# mv ./Solution.sol {results}/{analysis_name}_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}/{analysis_name}_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}.solution.sol || echo "solution failed"
# mv ./Feasible.txt \
#    {results}/{analysis_name}_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}/{analysis_name}_M{{wildcards.m}}_v{{wildcards.v}}_{{wildcards.file}}.feasible.txt || echo "feasible failed"