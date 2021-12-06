import os
import re
from glob import glob
from datatable import dt, f, fread, join, update

configfile: config.yml

# -- 0.0 Configuration variables
src_path = config["src_path"]
merida_repo = config["merida_repo"]
bin_path = config["bin_path"]

# -- 0.1 Project tool installation
rule download_and_install_merida:
    params:
        repo=merida_repo
    output:
        src=os.path.join(src_path),
        bin=os.path.join(bin_path)
    shell:
        """
        git clone {params.repo}
        mv {params.repo} {params.src}
        cd {params.src}
        mkdir build
        cd build
        cmake -DEIGEN_PATH=$CONDA_PREFIX/cplex -DCPLEX_PATH=$CONDA_PREFIX/lib/python3.8/site-packages/cplex
        """

