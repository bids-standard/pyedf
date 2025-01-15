# pyedf

Python package to read from and write EEG data to European Data Format files.

An experimental Python 3 compatible version is available from the `py3k` branch.

The EDF standard was introduced in [1992](https://doi.org/10.1016/0013-4694(92)90009-7). The extension of EDF with annotations was first described in [1998](https://doi.org/10.1016/S0013-4694(98)00029-7) and more formalized with the EDF+ standard that was published in [2003](https://doi.org/10.1016/S1388-2457(03)00123-8). The [EDF](http://www.edfplus.info) website describes both standards and discusses implementation details.

Since EDF is such a widely used format, there exist multiple implementations for reading and writing EDF files. Most Python implementations consist of wrappers around the C-code implementation, which makes installation more cumbersome and reduces portability. This implementation is in pure Python with limited dependencies on external packages.

## See also
  - http://www.edfplus.info
  - https://github.com/holgern/pyedflib
  - https://github.com/the-siesta-group/edfio
  - https://github.com/MNE-tools/MNE-python
