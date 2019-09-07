import EDF

__author__ = "Sandeepan Bhattacharyya <bsandeepan95.work@gmail.com>"
__version__ = "0.1.1"

# Edit the name of files to test out EDFWriter class
file_in = EDF.EDFReader()
file_in.open('test_generator_2.edf')

file_out = EDF.EDFWriter()
file_out.open('test_generator_2 copy.edf')

header = file_in.readHeader()
print(header)
file_out.writeHeader(header)

meas_info = header[0]
for i in range(meas_info['n_records']):
    data = file_in.readBlock(i)
    file_out.writeBlock(data)

file_in.close()
file_out.close()