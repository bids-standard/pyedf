# -*- coding: utf-8 -*-
# Copyright (C) 2018 Phillip Alday <phillip.alday@mpi.nl>
# License: BSD (3-clause)
"""Tests for writing EDF files"""

from nose.tools import assert_dict_equal, assert_raises

from config_vars import data_path
from os.path import join