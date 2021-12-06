import os
import re
from glob import glob
from datatable import dt, f, fread, join, update

configfile: config.yaml

# -- 0. Configuration
rule download_merida:


# -- 1. Export 