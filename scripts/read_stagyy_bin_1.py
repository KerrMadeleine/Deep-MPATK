import numpy as np
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import csv
import os
import sys
import stagyypythonmodule as spm
from cmcrameri import cm

import StagModelClass as smc 
from modelparameters import STAGYY_OUTPUT_FOLDER, MESH_DIR, STAGYY_MODS, STAGYY_MODS_YDIM, STAGYY_MODS_ZDIM, STAGYY_MODS_NMAX, STAGYY_MODS_TCMB, STAGYY_MODS_TM
from modelparameters import TSURF, CORE_RADIUS, DOMAIN_THICKNESS


