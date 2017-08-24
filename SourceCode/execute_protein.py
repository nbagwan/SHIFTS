__author__ = 'nbagwan'

## need py2exe as a parameter


# Used successfully in Python2.5 with matplotlib 0.91.2 and PyQt4 (and Qt 4.3.3)
from distutils.core import setup
import py2exe, sys, os
import scipy
import matplotlib
sys.setrecursionlimit(5000)
# opts = {
#     'py2exe': { "includes" : ["matplotlib.backends",  "matplotlib.backends.backend_qt4agg",
#                                "matplotlib.figure","pylab", "numpy",
#                                "matplotlib.backends.backend_tkagg", "scipy", "Tkinter", "ttk", r"scipy.sparse.csgraph._validation",
#                           r"scipy.special._ufuncs_cxx"]
#               }
#        }
#
#
# # for console program use 'console = [{"script" : "scriptname.py"}]
# setup(console=[{"script" : "Guimethod.py"}])

setup(
    data_files = matplotlib.get_py2exe_datafiles(),
    console=['SHIFTS.py'],
    options={
        'py2exe': {
            r'includes': ["matplotlib.backends",  "matplotlib.backends.backend_qt4agg",
                                "matplotlib.figure","pylab", "numpy",
                               "matplotlib.backends.backend_tkagg", 'scipy', 'scipy.integrate', 'scipy.special.*','scipy.linalg.*', "Tkinter", "ttk", r'scipy.sparse.csgraph._validation',
                          r'scipy.special._ufuncs_cxx'],
            'excludes': ['zmq.libzmq'],
            'dll_excludes': ['libzmq.pyd']
        }
    }
)
