"""
gb_network
"""
from __future__ import absolute_import
from __future__ import print_function

try:
    import argparse
except ImportError:
    raise ImportError("require argparse")

#numerical / data packages
try:
    import numpy as np
    np.set_printoptions(threshold=10)
except ImportError:
    raise ImportError("require numpy")
try:
    import pandas as pd
except ImportError:
    raise ImportError("require pandas")
try:
    import pytraj as pt
except ImportError:
    raise ImportError("require pytraj")
try:
    import scipy
    import scipy.sparse
    import scipy as sp
except ImportError:
    raise ImportError("require scipy")

#utilities
import os
import sys
import gc
import copy
import time
import re
import collections
import subprocess
import glob
import sqlite3
from itertools import islice
from itertools import combinations
from sqlite3 import Error

#Others
import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import networkx as nx

#self defined functions
import functions.correlation_data_utilities as corr_utils
import functions.database_mod as db_m

# Handle versioneer
from _version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
