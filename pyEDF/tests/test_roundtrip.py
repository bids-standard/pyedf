# -*- coding: utf-8 -*-
# Copyright (C) 2018 Phillip Alday <phillip.alday@mpi.nl>
# License: BSD (3-clause)
"""Round Trip Data IO tests.

    Note that these tests are *slow*. The entirety of each file is effectively
    read in 3 times and written to disk once.
"""

from nose.tools import assert_dict_equal, assert_true, assert_sequence_equal
import numpy as np

from config_vars import data_path
from os.path import join

import pyEDF

try:
    # Python 2
    from os import tempnam, unlink
except ImportError:
    # Python 3
    from os import unlink
    from tempfile import NamedTemporaryFile

    tempnam = lambda: NamedTemporaryFile(delete=False).name

def _roundtrip(edf_file):
    edfin = pyEDF.EDFReader(edf_file)

    fout = tempnam()
    edfout = pyEDF.EDFWriter(fout)

    header = edfin.readHeader()
    edfout.writeHeader(header)

    meas_info = header[0]
    for i in range(int(meas_info['n_records'])):
        edfout.writeBlock(edfin.readBlock(i))

    edfin.close()
    edfout.close()

    original = pyEDF.EDFReader(edf_file)
    copy = pyEDF.EDFReader(fout)

    copy_header = copy.readHeader()

    assert_dict_equal(header[0], copy_header[0])
    for key in header[1]:
        assert_sequence_equal(list(header[1][key]), list(copy_header[1][key]))

    for i in range(int(meas_info['n_records'])):
        # although this is floating point, it's not really numerics.
        # it's just comparing copies of the same data, so exact equality
        # should be a doable goal.
        for ch_orig, ch_copy in zip(original.readBlock(i), copy.readBlock(i)):
            assert_sequence_equal(list(ch_orig), list(ch_copy))

    unlink(fout)

def test_roundtrip_0601_s():
    '''Roundtrip of file 0601_s.edf'''
    # this file seems to have bad physical_min values
    _roundtrip(join(data_path, '0601_s.edf'))

def test_roundtrip_composition1_0s_to_1892s_fs20_15channels_tap127():
    '''Roundtrip of file composition1_0s_to_1892s_fs20_15channels_tap127.edf'''
    _roundtrip(join(data_path, 'composition1_0s_to_1892s_fs20_15channels_tap127.edf'))


def test_roundtrip_NY394_VisualLoc_R1():
    '''Roundtrip of file NY394_VisualLoc_R1.edf'''
    _roundtrip(join(data_path, 'NY394_VisualLoc_R1.edf'))


def test_roundtrip_shhs1_200001():
    '''Roundtrip of file shhs1-200001.edf'''
    _roundtrip(join(data_path, 'shhs1-200001.edf'))

def test_roundtrip_testAlphaIR20170321_0():
    '''Roundtrip of file testAlphaIR20170321-0.edf'''
    _roundtrip(join(data_path, 'testAlphaIR20170321-0.edf'))

def test_roundtrip_test_generator():
    '''Roundtrip of file test_generator.edf'''
    _roundtrip(join(data_path, 'test_generator.edf'))

