# pyedf

pyedf is a python package to read from and write EEG data to European Data Format files.

Since EDF is such a widely used format, there exist multiple Python implementations for reading and writing EDF files. However, most of these Python modules consist of wrappers around the C-code implementation `edflib`, which makes installation more cumbersome and reduces portability. This implementation is in pure python with limited dependencies on external packages while having support for Python 2.7 and 3.x.

The EDF standard was introduced in [1992](https://doi.org/10.1016/0013-4694(92)90009-7). The extension of EDF with annotations was first described in [1998](https://doi.org/10.1016/S0013-4694(98)00029-7) and more formalized with the EDF+ standard that was published in [2003](https://doi.org/10.1016/S1388-2457(03)00123-8). The [EDF](http://www.edfplus.info) website describes both standards and discusses implementation details.

## See also
  - http://www.edfplus.info
  - https://github.com/holgern/pyedflib
  - https://github.com/MNE-tools/MNE-python
