#!/usr/bin/python

import pandas as pd
import glob
import sys 
import subprocess
import yaml
import os
import numpy as np
import os
from copy import deepcopy

#segment cells
subprocess.call(['/S3segmenter/run_S3segmenterWrapperDocker.sh','/opt/mcr/v94/','/config/S3segmenter_config.yml'])

