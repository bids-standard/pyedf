import EDF

__author__ = "Sandeepan Bhattacharyya <bsandeepan95.work@gmail.com>"
__version__ = "0.1.1"

# Edit the name of files to test out EDFReader and EDFWriter class
input_str = (
    "Type 1 to only read from test EDF file.\n" +
    "Type 2 to only write from test EDF file.\n" +
    "Type 3 to read and write from test EDF file.\n" +
    "Press Enter after typing your choice: ")
user_input = int(input(input_str))

filename = 'test_generator_2.edf'
file_in = EDF.EDFReader()
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
    file_out = EDF.EDFWriter()
    file_out.open('copy of ' + filename)
    file_out.write_header(header)
    
    meas_info = header[0]
    for i in range(meas_info['n_records']):
        data = file_in.read_block(i)
        file_out.write_block(data)
    file_out.close()

file_in.close()
