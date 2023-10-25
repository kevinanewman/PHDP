"""
    Post Horiba Data Processing

"""

code_version = "0.0.0"
print('loading phdp version %s' % code_version)

import os, sys

# if 'darwin' in sys.platform:
#     os.environ['QT_MAC_WANTS_LAYER'] = '1'  # for pyqtgraph on MacOS

import traceback

try:
    import time

    import pandas as pd
    # from warnings import simplefilter
    # simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
    pd.set_option('chained_assignment', 'raise')
    # from pandas.api.types import is_numeric_dtype

    import numpy as np
    np.seterr(all='raise')

    import copy

    from common import phdp_globals
    from common.phdp_types import *
    from common.phdp_functions import *
    from common import file_io, phdp_log

    import CFR1065

    from tkinter import filedialog

    class PHDPSettings(PHDPBase):
        """
        Define the settings required for a simulation session

        """
        def __init__(self):
            """
            Create an OMEGASessionSettings object with default settings used for testing and development.

            The primary way to create an OMEGASessionSettings object is via the batch process.

            """
            import time

            self.session_name = 'PHDP'
            self.session_unique_name = 'PHDP'

            self.verbose = True

            path = os.path.dirname(os.path.abspath(__file__)) + os.sep
            self.output_folder_base = path + 'out' + os.sep
            self.output_folder = self.output_folder_base
            self.phdp_path = path
            self.logfile_prefix = 'phdplog_'
            self.logfilename = ''
            self.timestamp_str = time.strftime('%Y%m%d_%H%M%S')

            self.horiba_file = None  # path + '/test_inputs/HD02.202309061402.00097.GHGTRNS.Tn.xlsxm'

            self.start_time = 0
            self.end_time = 0

            self.run_profiler = False
            self.multiprocessing = False and not self.run_profiler and not getattr(sys, 'frozen', False)

            self.notification_destination = None
            self.notification_email = None
            self.notification_password = None
            self.encoding = {  # input file encoding (decoding) methods
                'HD02': 'cp1252',  # observed file encoding from test data
            }

            self.output_encoding = 'utf-8'  # cp1253 supports greek letters on Windows machines
            self.chemical_balance_convergence_tolerance = 0.01


except:
    print("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
    os._exit(-1)