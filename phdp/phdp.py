"""

PHDP top level code

----

**CODE**

"""

import sys, os

import numpy as np
import pandas as pd
import math

import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(path, '..'))  # picks up omega_model sub-packages

from phdp import *
from constants import constants
import test_sites
from phdp_reports import *

from iteration_alignment_and_drift_correction import *
from calculations_stats_validation import *


def init_logfile():
    phdp_log.init_logfile()
    phdp_log.logwrite("Initializing PHDP %s:" % code_version)


def init_phdp(runtime_options):
    """
    Initialize PHDP

    Args:
        runtime_options:

    Returns:
        List of init errors, else empty list on success

    """
    phdp_globals.options = runtime_options

    phdp_globals.test_data = dict()
    phdp_globals.validated_data = dict()

    init_fail = []

    return init_fail


def get_unitized_columns(filename, sheet_name=None, ignore_units=('Text', ''), encoding='utf-8', units_nrows=1):
    """
    Combine column labels and units row into a single combined string to identify the column.

    Args:
        filename (str): name of the file to read
        sheet_name (str): for reading a particular Excel sheet, or ``None`` for CSV files
        ignore_units (list of str): unit values to ignore, default ``Text``
        encoding (str): file encoding (decoding) method name
        units_nrows (int): number of units rows, if any

    Returns:
        List of combined column headers and units as in ``['ColumnName_units', ...]``

    """
    if sheet_name:
        columns = pd.read_excel(filename, header=None, nrows=1, sheet_name=sheet_name)
        if units_nrows > 0:
            units = pd.read_excel(filename, header=None, skiprows=1, nrows=units_nrows, sheet_name=sheet_name)
        else:
            units = pd.DataFrame({'units': [''] * columns.shape[1]}).transpose()
    else:
        columns = pd.read_csv(filename, header=None, nrows=1,  encoding=encoding,
                              encoding_errors='strict')
        if units_nrows > 0:
            units = pd.read_csv(filename, header=None, skiprows=1, nrows=units_nrows, encoding=encoding,
                                encoding_errors='strict')
        else:
            units = pd.DataFrame({'units': [''] * columns.shape[1]}).transpose()

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


def load_data(test_filename_prefix, test_site):
    """
        Load test data into phdp_globals.test_data dict

    Args:
        test_filename_prefix (str): test filename prefix, e.g. 'HD05.202402141518.00033.FTP'
        test_site (str): the name of the site where the test was performed, e.g. 'HD02'

    """
    required_file_names = ['BagData', 'BagDriftCheck', 'ContinuousData', 'CycleDefinition', 'CycleDefinition',
                           'DriftCheck', 'EmsCalResults', 'EmsComponents', 'EngineData', 'Header', 'MapResults',
                           'ModalTestData', 'ModeValidationResults', 'TestDetails', 'TestParameters',
                           'drift_corrected_BagData', 'Workstation', 'PreTest', 'CFR1065EMS', 'CFR1065CVS',
                           'CFR1065PM']

    optional_file_names = ['CVSDLSFlows', 'CVSDLSSampleResults', 'CFR1065CVS', 'CFR1065PM']

    required_file_names.extend(optional_file_names)  # add optional files

    input_files = sorted([f for f in os.listdir() if f.startswith(test_filename_prefix) and f.endswith('.csv')])
    for input_file in input_files:
        file_name = input_file.rsplit('.', 2)[-2]

        if file_name in required_file_names and file_name != 'Processing':
            phdp_log.logwrite('reading %s...' % input_file)
            if 'tad' not in file_name:
                unitized_columns = get_unitized_columns(input_file, encoding=phdp_globals.options.encoding[test_site])
                phdp_globals.test_data[file_name] = \
                    pd.read_csv(input_file, names=unitized_columns,
                                encoding=phdp_globals.options.encoding[test_site], encoding_errors='strict',
                                header=1, skiprows=0)
            else:
                get_unitized_columns(input_file, units_nrows=0, encoding=phdp_globals.options.output_encoding)  # dump tad columns to logfile for reference
                phdp_globals.test_data[file_name] = (
                    pd.read_csv(input_file, encoding=phdp_globals.options.output_encoding))

    for optional_file in optional_file_names:
        if optional_file not in phdp_globals.test_data:
            phdp_globals.test_data[optional_file] = None


def proportionality_check(ref, meas, skip_secs=5, test_type=None):
    """
    Calculate proportionality check percentage, the ratio of standard error to the mean

    Args:
        ref (list): A list of reference values.
        meas (list): A list of corresponding measured values.
        skip_secs (int): The number of seconds to skip at the start of the data, if any
        test_type (str): test type, i.e. 'transient' or 'modal'
        mode_numbers (dataframe): CVSDLSFlows mode number data for modal test processing

    Returns:
        Proportionality percent, the ratio of standard error to the mean

    """
    from statistics import linear_regression

    if test_type == 'transient' and skip_secs > 0:
        # skip the first few seconds of the data
        ref_skip = \
            ref.iloc[int(skip_secs / constants['MeasurementPeriod_s']):]
        meas_skip = \
            meas.iloc[int(skip_secs / constants['MeasurementPeriod_s']):]
    else:
        ref_skip = ref
        meas_skip = meas

    samples_per_second = int(1 / constants['MeasurementPeriod_s'])
    sample_length_1Hz = int(len(ref_skip) / samples_per_second) * samples_per_second

    # truncate data to an even multiple of the 1Hz sample period
    ref_skip_trunc = ref_skip.iloc[0: sample_length_1Hz]
    meas_skip_trunc = meas_skip.iloc[0: sample_length_1Hz]

    # average the data in 1Hz intervals
    ref_skip_1Hz = ref_skip_trunc.values.reshape(
        int(len(ref_skip_trunc) / samples_per_second), samples_per_second).mean(1)

    meas_skip_1Hz = meas_skip_trunc.values.reshape(
        int(len(meas_skip_trunc) / samples_per_second), samples_per_second).mean(1)

    # calculate a linear fit of the measured data, with a zero intercept
    slope, intercept = linear_regression(ref_skip_1Hz, meas_skip_1Hz,
                                         proportional=True)
    # calculate proportionality percent
    SEE = calc_STEYX(ref_skip_1Hz, meas_skip_1Hz, slope, intercept)

    meas_skip_1Hz_mean = meas_skip_1Hz.mean()

    proportionality_pct = SEE / meas_skip_1Hz_mean * 100

    return proportionality_pct


mode_max_filter_pressure_drop_kPa = 0
mode_max_proportionality_pct = 0


def get_pm_measurement(prompt):
    """
    This function prompts the user to enter a PM (particulate matter) measurement in mg.

    The input is validated to ensure it's a numeric value and must be confirmed by the user.

    If pm_measurement_mg is None or an invalid number is entered, the function will ask for a new input until
    a valid one is provided.

    Args:
        prompt (str): A string message to display to the user requesting the PM measurement.
        pm_measurement_mg (str): The initial input by the user for the PM measurement. Default is None.

    Returns:
        float: The validated PM measurement in mg.

    """
    valid_value = False
    while not valid_value:
        try:
            pm_measurement_mg = (
                float(input(prompt)))
            verify = input('Verify value %.4f (Y/n)' % pm_measurement_mg) or 'Y'
            if verify.strip().lower() == 'y':
                valid_value = True
        except:
            print('Invalid value entered, input must be numeric')

    return pm_measurement_mg


def particulate_matter_calculations(emissions_cycle_number, test_type, calc_mode, drift_corrected_time_aligned_data,
                                    validation_results):
    """

    Args:
        emissions_cycle_number (int): emissions cycle number (or mode for modal tests)
        test_type (str): test type, i.e. 'transient' or 'modal'
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        drift_corrected_time_aligned_data (DataFrame): drift-corrected time-aligned data

    Returns:

    """
    global mode_max_filter_pressure_drop_kPa, mode_max_proportionality_pct

    skip_secs = 5

    pm_sample_pts = phdp_globals.test_data['ContinuousData']['PMSampling_Logical'] == True

    drift_corrected_time_aligned_data['Kh'] = 9.953 * drift_corrected_time_aligned_data[
        'xH2Oint_mol/mol'] + 0.832

    if calc_mode != 'dilute-bag' and any(pm_sample_pts):
        CVSDLSFlows = phdp_globals.test_data['CVSDLSFlows']
        ContinuousData = phdp_globals.test_data['ContinuousData']

        if test_type == 'transient':
            pm_sample_start_index = ContinuousData.loc[pm_sample_pts, :].index[0]
            pm_sample_end_index = ContinuousData.loc[pm_sample_pts, :].index[-1]
            pm_sample_start_time = ContinuousData['Time_Date'].loc[pm_sample_start_index]
            pm_sample_end_time = ContinuousData['Time_Date'].loc[pm_sample_end_index]

            cvsdls_sample_start = CVSDLSFlows['Time_Date'] == pm_sample_start_time
            cvsdls_sample_end = CVSDLSFlows['Time_Date'] == pm_sample_end_time
            cvsdls_sample_start_index = CVSDLSFlows.loc[cvsdls_sample_start, :].index[0]
            cvsdls_sample_end_index = CVSDLSFlows.loc[cvsdls_sample_end, :].index[0]

            cvs_mass_flow_kgps = CVSDLSFlows['CVSMassFlow_kg/s'].loc[
                                 cvsdls_sample_start_index:cvsdls_sample_end_index+1]

            transfer_mass_flow_kgps = CVSDLSFlows['TransferMassFlow_g/s'].loc[
                                      cvsdls_sample_start_index:cvsdls_sample_end_index+1] / 1000
        else:  # modal
            cvs_mass_flow_kgps = CVSDLSFlows['CVSMassFlow_kg/s']
            transfer_mass_flow_kgps = CVSDLSFlows['TransferMassFlow_g/s'] / 1000

        if test_type == 'transient' and emissions_cycle_number not in validation_results['PM_results']:
            validation_results['PM_results'][emissions_cycle_number] = dict()

            cvs_mass_kg = cvs_mass_flow_kgps.sum() * constants['MeasurementPeriod_s']
            transfer_mass_kg = transfer_mass_flow_kgps.sum() * constants['MeasurementPeriod_s']

            dilution_factor = (cvs_mass_kg + transfer_mass_kg) / transfer_mass_kg

            pre_test_pm_measurement_mg = \
                get_pm_measurement('Enter pre-test PM test filter mass (mg)')

            post_test_pm_measurement_mg = \
                get_pm_measurement('Enter post-test PM test filter mass (mg)')

            pm_net_filter_mass_mg = post_test_pm_measurement_mg - pre_test_pm_measurement_mg
            validation_results['PM_results'][emissions_cycle_number]['pm_mass_g'] = \
                pm_net_filter_mass_mg * dilution_factor / 1000

            skip_index_count = int(skip_secs / constants['MeasurementPeriod_s'])

            validation_results['PM_results'][emissions_cycle_number]['Proportionality_pct'] = \
                proportionality_check(cvs_mass_flow_kgps.iloc[skip_index_count:],
                                      transfer_mass_flow_kgps.iloc[skip_index_count:], skip_secs=0,
                                      test_type=test_type)

            validation_results['PM_results'][emissions_cycle_number]['CVS Dilution Ratio'] = \
                (drift_corrected_time_aligned_data['CVSFlow_mol/s'].iloc[skip_index_count:] /
                 drift_corrected_time_aligned_data['nexh_mol/s'].iloc[skip_index_count:])

            cvsdls_sample_start = CVSDLSFlows['Time_Date'] == drift_corrected_time_aligned_data['Time_Date'].iloc[0]
            cvsdls_sample_end = CVSDLSFlows['Time_Date'] == drift_corrected_time_aligned_data['Time_Date'].iloc[-1]
            cvsdls_sample_start_index = CVSDLSFlows.loc[cvsdls_sample_start, :].index[0] + skip_index_count
            cvsdls_sample_end_index = CVSDLSFlows.loc[cvsdls_sample_end, :].index[0]

            validation_results['PM_results'][emissions_cycle_number]['TransferMassFlow_g/s'] = \
                CVSDLSFlows['TransferMassFlow_g/s'][
                cvsdls_sample_start_index:cvsdls_sample_end_index + 1].values

            validation_results['PM_results'][emissions_cycle_number]['FilterMassFlow_g/s'] = \
                CVSDLSFlows['FilterMassFlow_g/s'][
                cvsdls_sample_start_index:cvsdls_sample_end_index + 1].values

            validation_results['PM_results'][emissions_cycle_number]['FilterPressureDrop_kPa'] = \
                abs(CVSDLSFlows['FilterPressureDrop_kPa'][
                cvsdls_sample_start_index:cvsdls_sample_end_index + 1].values)

            validation_results['PM_results'][emissions_cycle_number]['PSU Dilution Ratio'] = \
                (validation_results['PM_results'][emissions_cycle_number]['FilterMassFlow_g/s'] /
                 validation_results['PM_results'][emissions_cycle_number]['TransferMassFlow_g/s'])

            validation_results['PM_results'][emissions_cycle_number]['Overall Dilution Ratio'] = \
                (validation_results['PM_results'][emissions_cycle_number]['PSU Dilution Ratio'] *
                 validation_results['PM_results'][emissions_cycle_number]['CVS Dilution Ratio'])

        elif test_type == 'modal':
            mode_number = emissions_cycle_number

            if mode_number not in validation_results['PM_results']:
                mode_numbers = CVSDLSFlows['ModeNumber_Integer']
                mode_indices = mode_numbers.loc[mode_numbers == mode_number].index[skip_secs:]

                pre_test_pm_measurement_mg = \
                    get_pm_measurement('Enter pre-test PM test filter mass (mg) for mode %d' % mode_number)

                post_test_pm_measurement_mg = \
                    get_pm_measurement('Enter post-test PM test filter mass (mg) for mode %d' % mode_number)

                cvs_mass_kg = cvs_mass_flow_kgps.loc[mode_indices].sum() * constants['MeasurementPeriod_s']
                transfer_mass_kg = transfer_mass_flow_kgps.loc[mode_indices].sum() * constants['MeasurementPeriod_s']

                dilution_factor = (cvs_mass_kg + transfer_mass_kg) / transfer_mass_kg

                pm_net_filter_mass_mg = post_test_pm_measurement_mg - pre_test_pm_measurement_mg

                mode_pressure_drop_kPa = abs(CVSDLSFlows['FilterPressureDrop_kPa'].loc[mode_indices[-1]] -
                                          CVSDLSFlows['FilterPressureDrop_kPa'].loc[mode_indices[0]])

                if mode_pressure_drop_kPa > mode_max_filter_pressure_drop_kPa:
                    mode_max_filter_pressure_drop_kPa = mode_pressure_drop_kPa

                mode_proportionality_pct = \
                    proportionality_check(cvs_mass_flow_kgps.loc[mode_indices],
                                          transfer_mass_flow_kgps.loc[mode_indices], skip_secs=0,
                                          test_type=test_type)

                if mode_proportionality_pct > mode_max_proportionality_pct:
                    mode_max_proportionality_pct = mode_proportionality_pct

                validation_results['PM_results'][mode_number] = dict()

                # store max pressure drop and proportionality percent
                validation_results['PM_results'][1]['FilterPressureDrop_kPa'] = \
                    [0, mode_max_filter_pressure_drop_kPa]

                validation_results['PM_results'][1]['Proportionality_pct'] = \
                    pd.Series(mode_max_proportionality_pct)

                # store various mode values
                validation_results['PM_results'][mode_number]['pm_mass_g'] = (
                        pm_net_filter_mass_mg * dilution_factor / 1000)

                validation_results['PM_results'][mode_number]['TransferMassFlow_g/s'] = \
                    CVSDLSFlows['TransferMassFlow_g/s'].loc[mode_indices]

                validation_results['PM_results'][mode_number]['FilterMassFlow_g/s'] = \
                    CVSDLSFlows['FilterMassFlow_g/s'].loc[mode_indices]

                validation_results['PM_results'][mode_number]['CVS Dilution Ratio'] = \
                    (drift_corrected_time_aligned_data['CVSFlow_mol/s'] /
                     drift_corrected_time_aligned_data['nexh_mol/s'])

                validation_results['PM_results'][mode_number]['PSU Dilution Ratio'] = \
                    (validation_results['PM_results'][mode_number]['FilterMassFlow_g/s'] /
                     validation_results['PM_results'][mode_number]['TransferMassFlow_g/s'])

                validation_results['PM_results'][mode_number]['Overall Dilution Ratio'] = \
                    (validation_results['PM_results'][mode_number]['PSU Dilution Ratio'] /
                    validation_results['PM_results'][mode_number]['CVS Dilution Ratio'].item())

                validation_results['PM_results'][mode_number]['PSU Dilution Ratio'] = \
                    validation_results['PM_results'][mode_number]['PSU Dilution Ratio'].min()

                validation_results['PM_results'][mode_number]['Overall Dilution Ratio'] = \
                    validation_results['PM_results'][mode_number]['Overall Dilution Ratio'].min()


def run_phdp(runtime_options):
    """
    Run Post-Horiba Data Processing and generate output summary and report files

    Args:
        runtime_options (PHDPSettings): processor settings

    Returns:
        dict of emissions results

    """
    runtime_options.start_time = time.time()

    try:
        init_fail = init_phdp(runtime_options.copy())

        if not init_fail:
            if phdp_globals.options.horiba_file is None:
                phdp_globals.options.horiba_file = \
                    filedialog.askopenfilename(title='Select any test file', filetypes=[('csv', '*.csv')])

            horiba_filename = file_io.get_filename(phdp_globals.options.horiba_file)

            os.chdir(file_io.get_filepath(phdp_globals.options.horiba_file))

            phdp_globals.options.output_folder = (file_io.get_filepath(phdp_globals.options.horiba_file) + os.sep +
                                                  horiba_filename.rsplit('.', 1)[0] + '.PHDP' + os.sep)

            phdp_globals.options.output_folder = phdp_globals.options.output_folder

            init_logfile()

            test_site, test_datetime, test_num, test_name, _ = horiba_filename.replace('.Tn', '').split('.')

            sampled_crank = False
            min_mode_number = 0

            if test_name in ('FTP', 'LLC'):
                test_type = 'transient'
                sampled_crank = True
            elif test_name == 'RMC':
                test_type = 'transient'
                min_mode_number = 1
            elif test_name in ('GHGTRNS', 'GHGCRSE'):
                test_type = 'transient'
            else:
                test_type = 'modal'

            if test_name in ('GHGTRNS', 'GHGCRSE'):
                vehicle_test = True
            else:
                vehicle_test = False

            phdp_log.logwrite('\nProcessing test %s (%s) from %s...\n' % (test_num, test_type, test_site))

            test_filename_prefix = horiba_filename.rsplit('.', 1)[0]

            # load data
            load_data(test_filename_prefix, test_site)

            emissions_cycles = \
                [ecn for ecn in phdp_globals.test_data['ContinuousData']['EmissionsCycleNumber_Integer'].unique()
                 if ecn >= 1]

            if test_type == 'modal':
                modes = phdp_globals.test_data['ModalTestData']['ModeNumber_Integer'].values
            else:
                modes = None

            constants['MeasurementPeriod_s'] = (
                phdp_globals.test_data['TestParameters']['ContinuousLoggerPeriod_s'].item())

            test_valid, validation_results = (
                validate_data(test_name, test_type, horiba_filename.rsplit('.', 1)[0], emissions_cycles, modes,
                              do_plots=False))

            if test_type == 'modal':
                emissions_cycles = modes

            if not test_valid:
                # Warning('\n!!! Test Validation Failed !!!')
                phdp_log.logwrite('\n!!! Test Validation Failed !!!')

            emissions_available = phdp_globals.test_data['BagData']['RbSpanValue_ppm'].max() > 0

            if emissions_available:

                report_output_prefix = test_filename_prefix + '-'
                report_filename = (phdp_globals.options.output_folder + report_output_prefix + 'report.xlsx')
                generate_pre_test_check_report(report_filename, test_datetime)
                generate_cycle_validation_report(report_filename, validation_results)

                if [p for p in phdp_globals.test_data['EmsComponents']['ParameterName'] if 'raw' in p.lower()]:
                    calc_modes = ['dilute', 'raw', 'dilute-bag']  # NOTE: 'dilute' mode must be first for reports
                else:
                    calc_modes = ['dilute', 'dilute-bag']

                if test_type == 'modal' or test_name == 'LLC':
                    # no bag data for LLC tests, duration is too long for bagging
                    calc_modes.remove('dilute-bag')

                for calc_mode in calc_modes:
                    results = \
                        {'1036_calculations': [], 'tad': [], 'tadsummary': [], 'dctad': [], 'dctadsummary': []}

                    for ecn in emissions_cycles:
                        print('\nProcessing %s %s %d ...' % (test_type, calc_mode, ecn))

                        # pull in raw data and time align as necessary
                        if test_type == 'transient':
                            emissions_cycle_number = ecn
                            time_aligned_data = time_align_continuous_data(test_site, vehicle_test, sampled_crank,
                                                                           emissions_cycle_number, min_mode_number,
                                                                           validation_results)
                        else:  # modal test
                            from test_sites import site_info

                            test_sites.init_site_info(test_site, test_type)

                            emissions_cycle_number = 1
                            mode_number = ecn

                            time_aligned_data = phdp_globals.test_data['ModalTestData'].loc[
                                phdp_globals.test_data['ModalTestData']['ModeNumber_Integer'] == mode_number].copy()

                            non_numeric_columns = [c for c in time_aligned_data.columns
                                                   if pd.api.types.is_object_dtype(time_aligned_data[c])]

                            time_aligned_data = time_aligned_data.drop(non_numeric_columns, axis=1)
                            time_aligned_data['tqShaft_Nm'] = time_aligned_data['tqShaft_Avg_Nm']
                            time_aligned_data['spDyno_rev/min'] = time_aligned_data['spDyno_Avg_rev/min']
                            # time_aligned_data['qmFuel_g/h'] = time_aligned_data['qmFuel_Avg_g/h']
                            time_aligned_data['tIntakeAir_°C'] = time_aligned_data['tIntakeAir_Avg_°C']
                            time_aligned_data['tCellDewPt_°C'] = time_aligned_data['tCellDewPt_Avg_°C']
                            time_aligned_data['pCellAmbient_kPa'] = time_aligned_data['pCellAmbient_Avg_kPa']

                            phdp_globals.test_data['ContinuousData']['time_s'] = (
                                    phdp_globals.test_data['ContinuousData']['Time_Date'] * 24 * 3600)

                            stabilization_pts = (
                                    (phdp_globals.test_data['ContinuousData']['ModeNumber_Integer'] == mode_number) &
                                    (phdp_globals.test_data['ContinuousData']['ModeStatus'] == 'Stabilizing'))

                            mode_pts = (
                                    (phdp_globals.test_data['ContinuousData']['ModeNumber_Integer'] == mode_number) &
                                    (phdp_globals.test_data['ContinuousData']['InModeLog_Logical'] == True))

                            continuous_data = phdp_globals.test_data['ContinuousData'].loc[mode_pts]

                            constants['SamplePeriod_s'] = len(continuous_data) * constants['MeasurementPeriod_s']

                            time_aligned_data['elapsed_time_s'] = constants['SamplePeriod_s']

                            time_aligned_data['stabilization_time_s'] = \
                                sum(stabilization_pts) * constants['MeasurementPeriod_s']

                            time_aligned_data['spDyno_Avg_rev/min_mean'] = continuous_data['spDyno_Avg_rev/min'].mean()
                            time_aligned_data['spDyno_Avg_rev/min_min'] = continuous_data['spDyno_Avg_rev/min'].min()
                            time_aligned_data['spDyno_Avg_rev/min_max'] = continuous_data['spDyno_Avg_rev/min'].max()

                            time_aligned_data['tqShaft_Avg_Nm_mean'] = continuous_data['tqShaft_Avg_Nm'].mean()
                            time_aligned_data['tqShaft_Avg_Nm_min'] = continuous_data['tqShaft_Avg_Nm'].min()
                            time_aligned_data['tqShaft_Avg_Nm_max'] = continuous_data['tqShaft_Avg_Nm'].max()

                            for optional_signal in site_info['optional_modal_signals']:
                                if optional_signal not in time_aligned_data:
                                    time_aligned_data[optional_signal] = (
                                        site_info)['optional_modal_signals'][optional_signal]

                            time_aligned_data = time_aligned_data.dropna(axis=1).reset_index(drop=True)

                        # add calculated values
                        time_aligned_data = pre_chemical_balance_calculations(time_aligned_data, calc_mode, test_type)

                        # chemical balance iteration to calculate xDil/Exh_mol/mol, xH2Oexh_mol/mol and xCcombdry_mol/mol
                        iterate_chemical_balance(time_aligned_data, calc_mode, emissions_cycle_number)

                        post_chemical_balance_calculations(time_aligned_data, calc_mode)

                        time_aligned_data_summary = (
                            calc_summary_results(time_aligned_data, calc_mode, emissions_cycle_number))

                        drift_corrected_time_aligned_data = time_aligned_data.copy()

                        # drift-correct concentrations
                        for signal_name in \
                                [col for col in drift_corrected_time_aligned_data.columns if col.startswith('con')]:
                            drift_correct_continuous_data(drift_corrected_time_aligned_data, signal_name,
                                                          emissions_cycle_number, test_name)

                        # drift-correct bag values
                        phdp_globals.test_data['drift_corrected_BagData'] = phdp_globals.test_data['BagData'].copy()
                        for idx in phdp_globals.test_data['drift_corrected_BagData'].index:
                            if phdp_globals.test_data['drift_corrected_BagData'].loc[idx, 'RbComponent'] != 'NMHC':
                                drift_correct_bag_data(phdp_globals.test_data['drift_corrected_BagData'], idx)

                        # add calculated values
                        drift_corrected_time_aligned_data = (
                            pre_chemical_balance_calculations(drift_corrected_time_aligned_data, calc_mode, test_type))

                        # chemical balance iteration to calculate xDil/Exh_mol/mol, xH2Oexh_mol/mol and xCcombdry_mol/mol
                        iterate_chemical_balance(drift_corrected_time_aligned_data, calc_mode, emissions_cycle_number,
                                                 drift_corrected=True)

                        post_chemical_balance_calculations(drift_corrected_time_aligned_data, calc_mode)

                        drift_corrected_time_aligned_data_summary = (
                            calc_summary_results(drift_corrected_time_aligned_data, calc_mode, emissions_cycle_number,
                                                 drift_corrected=True))

                        calculations_1036 = calc_1036_results(calc_mode, drift_corrected_time_aligned_data,
                                                              drift_corrected_time_aligned_data_summary,
                                                              emissions_cycle_number, test_type, vehicle_test)

                        if test_type == 'transient':
                            time_aligned_data_summary['EmissionsCycleNumber_Integer'] = emissions_cycle_number

                            if calc_mode != 'dilute-bag':
                                drift_corrected_time_aligned_data_summary['BagFillProportionality'] = (
                                    proportionality_check(time_aligned_data['CVSFlow_mol/s'],
                                                          time_aligned_data['BagFillFlow_Avg_mol/s']))

                            drift_corrected_time_aligned_data_summary['EmissionsCycleNumber_Integer'] = (
                                emissions_cycle_number)

                            calculations_1036['EmissionsCycleNumber_Integer'] = emissions_cycle_number
                        else:
                            time_aligned_data_summary['ModeNumber_Integer'] = mode_number

                            drift_corrected_time_aligned_data_summary['ModeNumber_Integer'] = mode_number

                            calculations_1036['ModeNumber_Integer'] = mode_number

                        # PM measurements and related calculations:
                        particulate_matter_calculations(ecn, test_type, calc_mode, drift_corrected_time_aligned_data,
                                                        validation_results)

                        results['tad'].append(time_aligned_data)
                        results['dctad'].append(drift_corrected_time_aligned_data)
                        results['tadsummary'].append(time_aligned_data_summary)
                        results['dctadsummary'].append(drift_corrected_time_aligned_data_summary)
                        results['1036_calculations'].append(calculations_1036)

                    print('\nSaving results...\n')

                    testdata_output_prefix = horiba_filename.rsplit('.', 1)[0] + '-%s-' % calc_mode

                    if test_type == 'modal':
                        index_name = 'ModeNumber_Integer'
                    else:
                        index_name = 'EmissionsCycleNumber_Integer'
                        phdp_globals.test_data['drift_corrected_BagData'].to_csv(
                            phdp_globals.options.output_folder + testdata_output_prefix + 'dcbagdata.csv', index=False,
                            encoding=phdp_globals.options.output_encoding, errors='replace')

                    pd.concat(results['tad']).set_index(index_name).to_csv(
                        phdp_globals.options.output_folder + testdata_output_prefix + 'tad.csv',
                        encoding=phdp_globals.options.output_encoding, errors='replace')

                    pd.concat(results['dctad']).set_index(index_name).to_csv(
                        phdp_globals.options.output_folder + testdata_output_prefix + 'dctad.csv',
                        encoding=phdp_globals.options.output_encoding, errors='replace')

                    pd.concat([pd.DataFrame(ts) for ts in results['tadsummary']]).set_index(index_name).to_csv(
                        phdp_globals.options.output_folder + testdata_output_prefix + 'tadsummary.csv',
                        encoding=phdp_globals.options.output_encoding, errors='replace')

                    # transfer PM results into dctadsummary
                    for idx, mode in enumerate(validation_results['PM_results'].keys()):
                        results['dctadsummary'][idx]['mPM_g'] = validation_results['PM_results'][mode]['pm_mass_g']

                    pd.concat([pd.DataFrame(ts) for ts in results['dctadsummary']]).set_index(index_name).to_csv(
                        phdp_globals.options.output_folder + testdata_output_prefix + 'dctadsummary.csv',
                        encoding=phdp_globals.options.output_encoding, errors='replace')

                    pd.concat([pd.DataFrame(ts) for ts in results['1036_calculations']]).set_index(index_name).to_csv(
                        phdp_globals.options.output_folder + testdata_output_prefix + '1036_calculations.csv',
                        encoding=phdp_globals.options.output_encoding, errors='replace')

                    if calc_mode == 'dilute':
                        generate_general_report(report_filename, calc_mode, results, validation_results,
                                                test_type, test_datetime, test_site)

                    if test_type == 'transient':
                        generate_transient_report(report_filename, calc_mode, results, validation_results,
                                                  test_datetime, test_type, test_num, test_site, vehicle_test)
                    else:
                        generate_modal_report(report_filename, calc_mode, results, validation_results,
                                              test_datetime, test_type, test_num, test_site)

                    generate_driftcheck_report(report_filename, results, test_type, test_name)

                print('done!')

                return results
            else:
                phdp_log.logwrite('No emissions data available, skipping emissions calcs.')
    except:
        phdp_log.logwrite("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        print("### Check PHDP log for error messages ###")
        phdp_log.end_logfile("\nSession Fail")


if __name__ == "__main__":
    try:
        results = run_phdp(PHDPSettings())
    except:
        print("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        os._exit(-1)
