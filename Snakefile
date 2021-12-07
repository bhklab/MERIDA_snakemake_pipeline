import os
import re
from glob import glob

configfile: "config.yml"

# -- 0.0 Configuration variables
src_path = config["src_path"]
merida_repo = config["merida_repo"]
bin_path = config["bin_path"]
cplex_path = config["cplex_path"]

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