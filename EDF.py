#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" pyedf is a python package to read from and write EEG data to European Data 
    Format files. Since EDF is such a widely used format, there exist multiple 
    Python implementations for reading and writing EDF files. However, most of 
    these Python modules consist of wrappers around the C-code implementation, 
    which makes installation more cumbersome and reduces portability. This 
    implementation is in pure python with limited dependencies on external 
    packages while having support for Python 2.7 and 3.x.

    Note: the EDF header is represented as a tuple of (meas_info, chan_info)
        meas_info should have the following: [
            'record_length', 'file_ver', 'hour', 'subject_id',
            'recording_id', 'n_records', 'month', 'subtype',
            'second', 'nchan', 'data_size', 'data_offset', 
            'lowpass', 'year', 'highpass', 'day', 'minute']
        chan_info should have the following: [
            'physical_min', 'transducers', 'physical_max', 
            'digital_max', 'ch_names', 'n_samps', 'units', 'digital_min']
    
    The EDF standard was introduced in 1992. The extension of EDF with 
    annotations was first described in 1998 and more formalized with the EDF+ 
    standard that was published in 2003. To learn more about both standards and 
    implementation details, check out https://www.edfplus.info/index.html
"""
__version__ = "0.2.0"
__author__ = "Robert Oostenveld <r.oostenveld@gmail.com> | Sandeepan B <bsandeepan95.work@gmail.com>"
__copyright__ = "https://bids.neuroimaging.io/"
__credits__ = ["various"]
__license__ = "BSD 3-Clause License"
__maintainer__ = "https://github.com/bids-standard"
__status__ = "Production"

from copy import deepcopy
from math import ceil, floor
from struct import pack, unpack
import calendar
import datetime
import numpy as np
import os
import re
from warnings import warn


def padtrim(buf, num):
    num -= len(buf)
    if num >= 0:
        # pad the input to the specified length
        return (str(buf) + ' ' * num)
    else:
        # trim the input to the specified length
        return buf[0:num]

def writebyte(file, content, encoding='ascii'):
    try:
        # Py3 onwards bytes and strings are separate data format
        content = bytes(content, encoding)
    except TypeError :
        # Py2.7 Support
        content = content.encode(encoding)
    except Exception as e:
        print(type(e))
        print(str(e))
        print(
            "If you see this message, please go to " + 
            "https://github.com/bids-standard/pyedf/issues" + 
            " and open an issue there regarding this. Thank you!")
    finally:
        file.write(content)

def set_offset(chan_obj):
    # return a tuple of offset value and calibrate value
    physical_range = chan_obj['physical_max'] - chan_obj['physical_min']
    digital_range = chan_obj['digital_max'] - chan_obj['digital_min']
    calibrate = physical_range / digital_range
    calibrated_dig_min = calibrate * chan_obj['digital_min']
    offset = chan_obj['physical_min'] - calibrated_dig_min
    return (offset, calibrate)


class EDFWriter():
    def __init__(self, fname=None):
        self._initial_state()
        if fname:
            self.open(fname)
    
    def _initial_state(self):
        self.fname = None
        self.meas_info = None
        self.chan_info = None
        self.calibrate = None
        self.offset    = None
        self.n_records = 0
    
    def open(self, fname):
        with open(fname, 'wb') as fid:
            assert fid.tell() == 0
        self.fname = fname

    def close(self):
        # it is still needed to update the number of records in the header
        # this requires copying the whole file content
        meas_info = self.meas_info
        chan_info = self.chan_info

        # update the n_records value in the file
        tempname = self.fname + '.bak'
        os.rename(self.fname, tempname)

        with open(tempname, 'rb') as fid1:
            assert fid1.tell() == 0
            with open(self.fname, 'wb') as fid2:
                assert fid2.tell() == 0
                fid2.write(fid1.read(236))

                # skip this part
                fid1.read(8)
                # but write this instead
                writebyte(fid2, padtrim(str(self.n_records), 8))

                writebyte(fid2, str(fid1.read(meas_info['data_offset'] - 236 - 8)))
                
                blocksize = np.sum(chan_info['n_samps']) * meas_info['data_size']
                
                for block in range(self.n_records):
                    writebyte(fid2, str(fid1.read(blocksize)))
        
        os.remove(tempname)
        self._initial_state()

    def write_header(self, header):
        meas_info = header[0]
        chan_info = header[1]
        meas_size = 256
        chan_size = 256 * meas_info['nchan']

        with open(self.fname, 'wb') as fid:
            assert fid.tell() == 0

            # fill in the missing or incomplete information
            if not 'subject_id' in meas_info:
                meas_info['subject_id'] = ''
            
            if not 'recording_id' in meas_info:
                meas_info['recording_id'] = ''
            
            if not 'subtype' in meas_info:
                meas_info['subtype'] = 'edf'
            
            nchan = meas_info['nchan']
            
            if ((not 'ch_names' in chan_info) or (len(chan_info['ch_names']) < nchan)):
                chan_info['ch_names'] = [str(i) for i in range(nchan)]
            
            if ((not 'transducers' in chan_info) or (len(chan_info['transducers']) < nchan)):
                chan_info['transducers'] = ['' for i in range(nchan)]
            
            if ((not 'units' in chan_info) or (
                len(chan_info['units']) < nchan)):
                chan_info['units'] = ['' for i in range(nchan)]

            if meas_info['subtype'] in ('24BIT', 'bdf'):
                meas_info['data_size'] = 3  # 24-bit (3 byte) integers
            else:
                meas_info['data_size'] = 2  # 16-bit (2 byte) integers

            writebyte(fid, (padtrim('0', 8)))
            writebyte(fid, (padtrim(meas_info['subject_id'], 80)))
            writebyte(fid, (padtrim(meas_info['recording_id'], 80)))

            writebyte(fid, (padtrim('{:0>2d}.{:0>2d}.{:0>2d}'.format(
                meas_info['day'], meas_info['month'], meas_info['year']), 8)))
            
            writebyte(fid, (padtrim('{:0>2d}.{:0>2d}.{:0>2d}'.format(
                meas_info['hour'], meas_info['minute'], meas_info['second']), 8)))
            
            writebyte(fid, (padtrim(str(meas_size + chan_size), 8)))
            writebyte(fid, (' ' * 44))

            # the final n_records should be inserted on byte 236
            writebyte(fid, (padtrim(str(-1), 8)))  
            writebyte(fid, (padtrim(str(meas_info['record_length']), 8)))
            writebyte(fid, (padtrim(str(meas_info['nchan']), 4)))

            # ensure that these are all np arrays rather than lists
            for key in ['physical_min', 'transducers', 'physical_max', 
                'digital_max', 'ch_names', 'n_samps', 'units', 'digital_min']:
                chan_info[key] = np.asarray(chan_info[key])

            for i in range(meas_info['nchan']):
                writebyte(fid, (padtrim(chan_info['ch_names'][i], 16)))

            for i in range(meas_info['nchan']):
                writebyte(fid, (padtrim(chan_info['transducers'][i], 80)))
            
            for i in range(meas_info['nchan']):
                writebyte(fid, (padtrim(chan_info['units'][i], 8)))
            
            for i in range(meas_info['nchan']):
                writebyte(fid, (padtrim(str(chan_info['physical_min'][i]), 8)))
            
            for i in range(meas_info['nchan']):
                writebyte(fid, (padtrim(str(chan_info['physical_max'][i]), 8)))
            
            for i in range(meas_info['nchan']):
                writebyte(fid, (padtrim(str(int(chan_info['digital_min'][i])), 8)))
            
            for i in range(meas_info['nchan']):
                writebyte(fid, (padtrim(str(int(chan_info['digital_max'][i])), 8)))
            
            for i in range(meas_info['nchan']):
                writebyte(fid, (' ' * 80)) # prefiltering
            
            for i in range(meas_info['nchan']):
                writebyte(fid, (padtrim(str(chan_info['n_samps'][i]), 8)))
            
            for i in range(meas_info['nchan']):
                writebyte(fid, (' ' * 32)) # reserved
            
            meas_info['data_offset'] = fid.tell()

        self.meas_info = meas_info
        self.chan_info = chan_info
        self.offset, self.calibrate = set_offset(chan_info)

        channels = list(range(meas_info['nchan']))
        for ch in channels:
            if self.calibrate[ch]<0:
              self.calibrate[ch] = 1
              self.offset[ch]    = 0

    def write_block(self, data):
        meas_info = self.meas_info
        chan_info = self.chan_info
        with open(self.fname, 'ab') as fid:
            assert fid.tell() > 0
            for i in range(meas_info['nchan']):
                raw = deepcopy(data[i])

                assert len(raw) == chan_info['n_samps'][i]
                if min(raw) < chan_info['physical_min'][i]:
                    warn('Value exceeds physical_min: ' + str(min(raw)))
                if max(raw) > chan_info['physical_max'][i]:
                    warn('Value exceeds physical_max: '+ str(max(raw)))

                raw -= self.offset[i]  # FIXME I am not sure about the order of calibrate and offset
                raw /= self.calibrate[i]

                raw = np.asarray(raw, dtype=np.int16)
                buf = [pack('h', x) for x in raw]
                for val in buf:
                    fid.write(val)
            self.n_records += 1


class EDFReader():
    def __init__(self, fname=None):
        self._initial_state()
        if fname:
            self.open(fname)

    def _initial_state(self):
        self.fname = None
        self.meas_info = None
        self.chan_info = None
        self.calibrate = None
        self.offset    = None

    def open(self, fname):
        with open(fname, 'rb') as fid:
            assert fid.tell() == 0
        self.fname = fname
        self.read_header()
        return self.meas_info, self.chan_info

    def close(self):
        self._initial_state()

    def read_header(self):
        meas_info = {}
        chan_info = {}
        with open(self.fname, 'rb') as fid:
            assert fid.tell() == 0

            meas_info['file_ver']     = fid.read(8).strip().decode()
            meas_info['subject_id']   = fid.read(80).strip().decode()
            meas_info['recording_id'] = fid.read(80).strip().decode()

            day, month, year     = [int(x) for x in re.findall('(\d+)', fid.read(8).decode())]
            hour, minute, second = [int(x) for x in re.findall('(\d+)', fid.read(8).decode())]

            meas_info['day'] = day
            meas_info['month'] = month
            meas_info['year'] = year
            meas_info['hour'] = hour
            meas_info['minute'] = minute
            meas_info['second'] = second

            date = datetime.datetime(year + 2000, month, day, hour, minute, second)

            meas_info['meas_date'] = calendar.timegm(date.utctimetuple())
            meas_info['data_offset'] = header_nbyte = int(fid.read(8).decode())

            subtype = fid.read(44).strip().decode()[:5]
            if len(subtype) > 0:
                meas_info['subtype'] = subtype
            else:
                meas_info['subtype'] = os.path.splitext(self.fname)[1][1:].lower()

            if meas_info['subtype'] in ('24BIT', 'bdf'):
                meas_info['data_size'] = 3  # 24-bit (3 byte) integers
            else:
                meas_info['data_size'] = 2  # 16-bit (2 byte) integers

            meas_info['n_records'] = int(fid.read(8).decode())

            # record length in seconds
            record_length = float(fid.read(8).decode())
            if record_length == 0:
                meas_info['record_length'] = record_length = 1.
                warn("Header measurement information is incorrect for recordlength. Default record length set to 1.")
            else:
                meas_info['record_length'] = record_length

            meas_info['nchan'] = int(fid.read(4).decode())
            channels = list(range(meas_info['nchan']))

            chan_info['ch_names']     = [fid.read(16).strip().decode() for ch in channels]
            chan_info['transducers']  = [fid.read(80).strip().decode() for ch in channels]
            chan_info['units']        = [fid.read(8).strip().decode() for ch in channels]
            chan_info['physical_min'] = np.array([float(fid.read(8).decode()) for ch in channels])
            chan_info['physical_max'] = np.array([float(fid.read(8).decode()) for ch in channels])
            chan_info['digital_min']  = np.array([float(fid.read(8).decode()) for ch in channels])
            chan_info['digital_max']  = np.array([float(fid.read(8).decode()) for ch in channels])

            prefiltering = [fid.read(80).strip().decode() for ch in channels][:-1]
            highpass = np.ravel([re.findall('HP:\s+(\w+)', filt) for filt in prefiltering])
            lowpass = np.ravel([re.findall('LP:\s+(\w+)', filt) for filt in prefiltering])
                
            high_pass_default = 0.
            if highpass.size == 0:
                meas_info['highpass'] = high_pass_default
            elif all(highpass):
                if highpass[0] == 'NaN':
                    meas_info['highpass'] = high_pass_default
                elif highpass[0] == 'DC':
                    meas_info['highpass'] = 0.
                else:
                    meas_info['highpass'] = float(highpass[0])
            else:
                meas_info['highpass'] = float(np.max(highpass))
                warn("Channels contain different highpass filters. Highest filter setting will be stored.")

            if lowpass.size == 0:
                meas_info['lowpass'] = None
            elif all(lowpass):
                if lowpass[0] == 'NaN':
                    meas_info['lowpass'] = None
                else:
                    meas_info['lowpass'] = float(lowpass[0])
            else:
                meas_info['lowpass'] = float(np.min(lowpass))
                warn('%s' % ("Channels contain different lowpass filters. Lowest filter setting will be stored."))
            # number of samples per record
            chan_info['n_samps'] = n_samps = np.array([int(fid.read(8).decode()) for ch in channels])

            fid.read(32 *meas_info['nchan']).decode()  # reserved
            assert fid.tell() == header_nbyte

            if meas_info['n_records'] == -1:
                # happens if n_records is not updated at the end of recording
                file_size = os.path.getsize(self.fname)
                total_samples = (file_size - meas_info['data_offset']) / meas_info['data_size']
                meas_info['n_records'] = total_samples / sum(n_samps)

        self.offset, self.calibrate = set_offset(chan_info)

        for ch in channels:
            if self.calibrate[ch] < 0:
              self.calibrate[ch] = 1
              self.offset[ch]    = 0

        self.meas_info = meas_info
        self.chan_info = chan_info
        return (meas_info, chan_info)

    def read_block(self, block):
        assert block >= 0
        meas_info = self.meas_info
        chan_info = self.chan_info
        data = []
        with open(self.fname, 'rb') as fid:
            assert fid.tell() == 0
            blocksize = np.sum(chan_info['n_samps']) * meas_info['data_size']
            fid.seek(meas_info['data_offset'] + block * blocksize)
            
            for i in range(meas_info['nchan']):
                buf = fid.read(chan_info['n_samps'][i]*meas_info['data_size'])
                raw = np.asarray(unpack('<{}h'.format(chan_info['n_samps'][i]), buf), dtype=np.float32)

                raw *= self.calibrate[i]
                raw += self.offset[i]  # FIXME I am not sure about the order of calibrate and offset
                data.append(raw)
        return data

    def read_samples(self, channel, begsample, endsample):
        meas_info = self.meas_info
        chan_info = self.chan_info
        n_samps = chan_info['n_samps'][channel]
        begblock = int(floor((begsample) / n_samps))
        endblock = int(floor((endsample) / n_samps))
        data = self.read_block(begblock)[channel]
        for block in range(begblock+1, endblock+1):
            data = np.append(data, self.read_block(block)[channel])
        begsample -= begblock*n_samps
        endsample -= begblock*n_samps
        return data[begsample:(endsample + 1)]

    def get_signal_text_labels(self):
        return [str(x) for x in self.chan_info['ch_names']]

    def get_n_signals(self):
        return self.meas_info['nchan']

    def get_signal_freqs(self):
        return self.chan_info['n_samps'] / self.meas_info['record_length']

    def get_n_samples(self):
        return self.chan_info['n_samps'] * self.meas_info['n_records']

    def read_signal(self, chanindx):
        begsample = 0
        endsample = (self.chan_info['n_samps'][chanindx] * self.meas_info['n_records']) - 1
        return self.read_samples(chanindx, begsample, endsample)


if __name__ == "__main__":
    input_str = (
        "Type 1 to only read from test EDF file.\n" +
        "Type 2 to only write from test EDF file.\n" +
        "Type 3 to read and write from test EDF file.\n" +
        "Press Enter after typing your choice: ")
    user_input = int(input(input_str))

    # Edit the filename variable to test out EDFReader and EDFWriter class 
    # with different files.
    filename = 'test_generator_2.edf'
    file_in = EDFReader()
    file_in.open(filename)
    header = file_in.read_header()

    if user_input in (1, 3):
        print("Following are data blocks from the EDF file.\n")
        print(file_in.read_samples(0, 0, 0))
        print('\n')
        print(file_in.read_samples(0, 0, 128))
        print('\n')
        print(header)

    if user_input in (2, 3):
        file_out = EDFWriter()
        file_out.open('copy of ' + filename)
        file_out.write_header(header)
        
        meas_info = header[0]
        for i in range(meas_info['n_records']):
            data = file_in.read_block(i)
            file_out.write_block(data)
        file_out.close()

    file_in.close()
