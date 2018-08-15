# -*- coding: utf-8 -*-
# Copyright (C) 2018 Phillip Alday <phillip.alday@mpi.nl>
# License: BSD (3-clause)
"""Round Trip Data IO tests."""

from nose.tools import assert_dict_equal, assert_true

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
    assert_dict_equal(header[1], copy_header[1])

    for i in range(int(meas_info['n_records'])):
        # although this is floating point, it's not really numerics.
        # it's just comparing copies of the same data, so exact equality
        # should be a doable goal.
        assert_equal(original.readBlock(), copy.readBlock())

    unlink(fout)

def test_roundtrip_0601_s():
    '''Roundtrip of file 0601_s.edf'''
    _roundtrip(join(data_path, '0601_s.edf'))

def test_roundtrip_composition1_0s_to_1892s_fs20_15channels_tap127():
    '''Roundtrip of file composition1_0s_to_1892s_fs20_15channels_tap127.edf'''
    _roundtrip(join(data_path, 'composition1_0s_to_1892s_fs20_15channels_tap127'))


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

