biimport os
import re
from glob import glob
from datatable import dt, f, fread, join, update

configfile: config.yml

# -- 0.0 Configuration variables


# -- 0.1 Project tool installation
rule download_merida:
    params:
        target=
    shell:
        

rule install_merida:
    params:
        target=
    shell:

