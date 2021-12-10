import os
import re
from glob import glob

configfile: "config.yml"

# -- 0.0 Configuration variables
src_path = config["src_path"]
merida_repo = config["merida_repo"]
bin_path = config["bin_path"]
cplex_path = config["cplex_path"]

# -- 0.1 Make MERIDA input files for to grid search hyperparameterss
analysis_name = config["analysis_name"]
feature_matrix = config["feature_matrix"]
response_vector = config["response_vector"]
out_dir = config["out_dir"]
threshold = config["threshold"]

M_range = config["M"]
v_function = config["v"]

M_list = list(range(M_range["start"], M_range["stop"] + 1, M_range["step"]))
param_grid = [(m, v) for m in M_list for v in v_function]

merida_files = [
    f"procdata/merida_input/{analysis_name}_M{m}_v{v}_{os.path.basename(feature_matrix)}" 
        for m, v in param_grid
]

## -- 3.0 Run MERIDA_ILP
merida_results = [
    f"results/{analysis_name}_M{m}_v{v}_{os.path.basename(feature_matrix)}" 
        for m, v in param_grid
]


rule all:
    input:
        merida_results

# -- 0.1 Project tool installation
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

# -- 1.0 Download the Project Data
## TODO:: Add download and preprocessing scripts
# pset_list = config["pset_list"]
# rule download_psets:
#     input:
#         os.path.join(bin_path, "MERIDA_ILP")
#     params:
#         psets=


rule build_merida_input_files:
    params:
        grid=param_grid,
        feature_matrix=feature_matrix,
        response_vector=response_vector,
        out_dir=out_dir,
        threshold=threshold
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
                f"Threshold\t{params.threshold}"
            )
            with open(fl, "w+") as f:
                f.write(config)


rule run_merida:
    input:
        "procdata/merida_input/roche_{file}.txt"
    params:
        "no"
    output:
        "results/roche_{file}.txt"
    shell:
        """
        MERIDA_ILP {input} {params} || touch {output}
        """
