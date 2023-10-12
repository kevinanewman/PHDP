"""

PHDP top level code

----

**CODE**

"""

import sys, os

path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(path, '..'))  # picks up omega_model sub-packages

from phdp import *


def init_phdp(runtime_options):
    """
    Initialize PHDP

    Args:
        runtime_options:

    Returns:
        List of init errors, else empty list on success

    """
    phdp_globals.options = runtime_options

    phdp_log.init_logfile()

    phdp_log.logwrite("Initializing PHDP %s:" % code_version)

    phdp_globals.test_data = dict()
    phdp_globals.validated_data = dict()

    init_fail = []

    return init_fail


def get_unitized_columns(filename, sheet_name=None, ignore_units=['Text']):
    """

    Args:
        filename:
        sheet_name:
        ignore_units:

    Returns:

    """
    if sheet_name:
        columns = pd.read_excel(filename, header=None, nrows=1, sheet_name=sheet_name)
        units = pd.read_excel(filename, header=None, skiprows=1, nrows=1, sheet_name=sheet_name)
    else:
        columns = pd.read_csv(filename, header=None, nrows=1, encoding_errors='ignore')
        units = pd.read_csv(filename, header=None, skiprows=1, nrows=1, encoding_errors='ignore')

    unitized_columns = []

    for col, unit in zip(columns.values[0], units.values[0]):
        if unit not in ignore_units:
            unitized_columns.append('%s_%s' % (col, unit))
        else:
            unitized_columns.append('%s' % col)

    if phdp_globals.options.verbose:
        for idx, uc in enumerate(unitized_columns):
            phdp_log.logwrite('%s: %s' % (get_letters_from_index(idx).ljust(2), uc))
        phdp_log.logwrite('')

    return unitized_columns


def load_data():
    """
        Load test data into phdp_globals.test_data dict

    """
    # read the Horiba file
    phdp_log.logwrite('reading %s...' % phdp_globals.options.horiba_file)
    unitized_columns = get_unitized_columns(phdp_globals.options.horiba_file, sheet_name='ContinuousData1')
    phdp_globals.test_data['horiba_data'] = \
        pd.read_excel(phdp_globals.options.horiba_file, names=unitized_columns, header=1, skiprows=0,
                      sheet_name='ContinuousData1')
    # read in the rest of the files just in case
    os.chdir(file_io.get_filepath(phdp_globals.options.horiba_file))
    input_files = sorted([f for f in os.listdir() if f.endswith('.csv')])
    for input_file in input_files:
        file_name = input_file.rsplit('.', 2)[-2]
        if file_name != 'Processing':
            phdp_log.logwrite('reading %s...' % input_file)
            unitized_columns = get_unitized_columns(input_file)
            phdp_globals.test_data[file_name] = \
                pd.read_csv(input_file, names=unitized_columns, encoding_errors='ignore', header=1, skiprows=0)


def run_phdp(runtime_options):
    """

    Args:
        runtime_options:

    Returns:

    """

    runtime_options.start_time = time.time()

    try:
        init_fail = init_phdp(runtime_options.copy())

        if not init_fail:
            if phdp_globals.options.horiba_file is None:
                phdp_globals.options.horiba_file = \
                    filedialog.askopenfilename(title='Select Horiba Tn File', filetypes=[('xlsm', '*.xlsm')])

            horiba_filename = file_io.get_filename(phdp_globals.options.horiba_file)

            test_site, test_datetime, test_num, test_type = horiba_filename.replace('.Tn', '').split('.')

            phdp_log.logwrite('\nProcessing test %s (%s) from %s...\n' % (test_num, test_type, test_site))

            load_data()



    except:
        phdp_log.logwrite("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        print("### Check PHDP log for error messages ###")
        phdp_log.end_logfile("\nSession Fail")


if __name__ == "__main__":
    try:
        run_phdp(PHDPSettings())
    except:
        print("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        os._exit(-1)
