"""

PHDP top level code

----

**CODE**

"""

import sys, os

import numpy as np
import pandas as pd
from scipy.stats import linregress
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


def time_align_continuous_data(test_site, vehicle_test, sampled_crank, emissions_cycle_number, min_mode_number,
                               validation_results):
    """
    Time-align continuous data, add Vehicle Moving from cycle definition

    Args:
        test_site (str): the name of the site where the test was performed, e.g. 'HD02'
        vehicle_test (bool): ``True`` if test has an associated vehicle speed trace
        sampled_crank (bool): ``True`` if test has sampled engine cranking and startup
        emissions_cycle_number (int): emissions cycle number to process
        min_mode_number (int): the starting mode number of the data
        validation_results (dict): dict of cycle validation data

    Returns:
        Dataframe of time-aligned data

    """
    from test_sites import site_info

    test_sites.init_site_info(test_site, 'transient')

    SamplePeriod_s = phdp_globals.test_data['TestParameters']['ContinuousLoggerPeriod_s'].item()

    constants['SamplePeriod_s'] = SamplePeriod_s  # nominal sample period for purposes of time-aligning data

    time_aligned_data = pd.DataFrame(index=phdp_globals.test_data['ContinuousData'].index)

    cycle_shift = validation_results['regression_results'][emissions_cycle_number]['time_shift']

    site_info['signals_and_delays']['ContinuousData']['spDyno_Avg_rev/min'] = cycle_shift
    site_info['signals_and_delays']['ContinuousData']['tqShaft_Avg_Nm'] = cycle_shift

    for source in site_info['signals_and_delays'].keys():
        for signal in site_info['signals_and_delays'][source]:
            delay_s = site_info['signals_and_delays'][source][signal]
            if delay_s is not None:
                delay_samples = round(delay_s / SamplePeriod_s)
                time_aligned_data = pd.concat([time_aligned_data,
                                               pd.DataFrame({signal: phdp_globals.test_data[source][signal]
                                                            .iloc[delay_samples:].values})], axis=1)

    time_aligned_data = \
        (time_aligned_data[time_aligned_data['EmissionsCycleNumber_Integer'] == emissions_cycle_number].
         reset_index(drop=True))

    if vehicle_test:
        # add vehicle moving flag
        test_cycle_definition = \
            phdp_globals.test_data['CycleDefinition'][phdp_globals.test_data['CycleDefinition']
                                               ['EmissionsCycleNumber_Integer'] == emissions_cycle_number]

        vehicle_moving_int = test_cycle_definition['VehicleMoving_Logical']

        time_aligned_data = pd.concat([time_aligned_data,
                                       pd.DataFrame({'VehicleMoving_Logical': vehicle_moving_int.values})], axis=1)

        time_aligned_data['VehicleMoving_Logical'] = 1 * time_aligned_data['VehicleMoving_Logical'].fillna(False)

    if sampled_crank:
        test_start_index = time_aligned_data[time_aligned_data['EngDynoMode'] == 'Starting'].index[0]
    else:
        test_start_index = time_aligned_data[time_aligned_data['ModeNumber_Integer'] >= min_mode_number].index[0]

    time_aligned_data = time_aligned_data.drop('EngDynoMode', axis=1)  # can't have strings in the tad

    test_end_index = time_aligned_data[time_aligned_data['ModeNumber_Integer'] == -1].index[0]

    time_aligned_data['time_s'] = time_aligned_data['Time_Date'] * 24 * 3600

    time_aligned_data['elapsed_time_s'] = (
            time_aligned_data['time_s'].iloc[test_end_index] - time_aligned_data['time_s'].iloc[test_start_index])

    time_aligned_data = time_aligned_data.iloc[test_start_index:test_end_index]

    time_aligned_data = time_aligned_data.fillna(0)  # TODO: for now, need to confirm...

    return time_aligned_data


def pre_chemical_balance_calculations(time_aligned_data, calc_mode, test_type):
    """
    Calculate values required for chemical balance iteration

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        test_type (str): 'transient' or 'modal'

    Returns:
        Nothing, updates time_aligned_data

    """
    from constants import constants

    if calc_mode == 'dilute-bag':
        if len(time_aligned_data) > 1:
            time_aligned_data['Power_kW'] = (
                np.maximum(0, time_aligned_data['tqShaft_Avg_Nm'] * time_aligned_data['spDyno_Avg_rev/min'] / 9548.8))

            # calculate average and total values
            time_aligned_data_avg = time_aligned_data.apply(lambda x: [x.mean()])
            time_aligned_data_sum = time_aligned_data.apply(lambda x: [x.sum()])

            time_aligned_data = (
                pd.DataFrame({'elapsed_time_s': [time_aligned_data['elapsed_time_s'].iloc[0]],
                              'EmissionsCycleNumber_Integer': time_aligned_data['EmissionsCycleNumber_Integer']
                             .iloc[0]}))

            time_aligned_data['Power_kW'] = time_aligned_data_avg['Power_kW']

            time_aligned_data['qmIntakeAir_Avg_kg'] = (
                    time_aligned_data_sum['qmIntakeAir_Avg_kg/h'] / 3600 * constants['SamplePeriod_s'])

            time_aligned_data['qmFuel_Avg_g'] = (
                    time_aligned_data_sum['qmFuel_Avg_g/h'] / 3600 * constants['SamplePeriod_s'])

            time_aligned_data['DEFMassFlowRate_Avg_g'] = (
                    time_aligned_data_sum['DEFMassFlowRate_Avg_g/h'] / 3600 * constants['SamplePeriod_s'])

            time_aligned_data['IntakeAirPress_Avg_kPa'] = time_aligned_data_avg['IntakeAirPress_Avg_kPa']
            time_aligned_data['tCellDewPt_Avg_°C'] = time_aligned_data_avg['tCellDewPt_Avg_°C']

            time_aligned_data['CVSDilAirTemp_Avg_°C'] = time_aligned_data_avg['CVSDilAirTemp_Avg_°C']

            time_aligned_data_sum['BagFillFlow_Avg_m³/s'] = time_aligned_data_sum['BagFillFlow_Avg_l/min'] / 60000

            if 'CVSDilAirRH_Avg_%' in time_aligned_data_avg:
                time_aligned_data['CVSDilAirRH_Avg_%'] = time_aligned_data_avg['CVSDilAirRH_Avg_%']

            if 'CVSFlow_Avg_m³/s' in time_aligned_data_sum:
                time_aligned_data['CVSFlow_mol'] = ((
                        (time_aligned_data_sum['BagFillFlow_Avg_m³/s'] + time_aligned_data_sum['CVSFlow_Avg_m³/s'])
                        / 0.024055 * constants['SamplePeriod_s'] +
                        phdp_globals.test_data['TestParameters']['DiluteSampleMolarFlow_mol/s'] *
                        time_aligned_data['elapsed_time_s']))

            if 'CVSDilAirDPTemp_Avg_°C' in time_aligned_data_avg:
                time_aligned_data['CVSDilAirDPTemp_Avg_°C'] = time_aligned_data_avg['CVSDilAirDPTemp_Avg_°C']

            if 'CVSMolarFlow_Avg_mol/s' in time_aligned_data_avg:
                time_aligned_data['CVSFlow_mol'] = (
                        (time_aligned_data_sum['CVSMolarFlow_Avg_mol/s'] +
                         time_aligned_data_sum['BagFillFlow_Avg_m³/s'] / 0.024055) * constants['SamplePeriod_s'] +
                        phdp_globals.test_data['TestParameters']['DiluteSampleMolarFlow_mol/s'] *
                        time_aligned_data['elapsed_time_s'])

            time_aligned_data['tIntakeAir_Avg_°C'] = time_aligned_data_avg['tIntakeAir_Avg_°C']

            time_aligned_data['pCellAmbient_Avg_kPa'] = time_aligned_data_avg['pCellAmbient_Avg_kPa']

        unit_rate = ''
    else:
        unit_rate = '/s'

        time_aligned_data['Power_kW'] = (
            np.maximum(0, time_aligned_data['tqShaft_Avg_Nm'] * time_aligned_data['spDyno_Avg_rev/min'] / 9548.8))

        time_aligned_data['BagFillFlow_Avg_m³/s'] = time_aligned_data['BagFillFlow_Avg_l/min'] / 60000

        time_aligned_data['BagFillFlow_Avg_mol/s'] = time_aligned_data['BagFillFlow_Avg_m³/s'] / 0.024055

        if test_type == 'transient':
            if 'CVSMolarFlow_Avg_mol/s' not in time_aligned_data:
                if 'DiluteSampleVolumeFlow_m³/s' in phdp_globals.test_data['TestParameters']:
                    time_aligned_data['CVSFlow_mol/s'] = \
                        (time_aligned_data['CVSFlow_Avg_m³/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] +
                        phdp_globals.test_data['TestParameters']['DiluteSampleVolumeFlow_m³/s'].item()) / 0.024055
                else:
                    time_aligned_data['CVSFlow_mol/s'] = \
                        (time_aligned_data['CVSFlow_Avg_m³/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] +
                        phdp_globals.test_data['TestParameters']['DiluteSampleVolumeFlow_l/s'].item() / 1000) / 0.024055
            else:
                # TODO: need to verify if this is correct for HD05 transient (FTP) test
                time_aligned_data['CVSFlow_mol/s'] = (
                        time_aligned_data['CVSMolarFlow_Avg_mol/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] /
                        0.024055 + phdp_globals.test_data['TestParameters']['DiluteSampleMolarFlow_mol/s'].item())
        else:
            if 'CVSMolarFlow_Avg_mol/s' not in time_aligned_data:
                time_aligned_data['qvCVS_Avg_m³/s'] = time_aligned_data['qvCVS_Avg_m³/min'] / 60

                time_aligned_data['CVSFlow_mol/s'] = (
                    (time_aligned_data['qvCVS_Avg_m³/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] +
                     phdp_globals.test_data['TestParameters']['DiluteSampleVolumeFlow_l/s'] / 1000) / 0.024055)
            else:
                time_aligned_data['CVSFlow_mol/s'] = (
                        time_aligned_data['CVSMolarFlow_Avg_mol/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] /
                        0.024055 + phdp_globals.test_data['TestParameters']['DiluteSampleMolarFlow_mol/s'])

    time_aligned_data['Tsat_K'] = time_aligned_data['CVSDilAirTemp_Avg_°C'] + 273.15

    if 'CVSDilAirRH_Avg_%' in time_aligned_data:
        # calculations based on relative humidity
        # 1065.645-1:
        time_aligned_data['pH2Odilsat_kPa'] = \
            CFR1065.vapor_pressure_of_water_kPa(time_aligned_data['Tsat_K'])
        time_aligned_data['pH2Odilscal_Pa'] = time_aligned_data['pH2Odilsat_kPa'] * \
                                              time_aligned_data['CVSDilAirRH_Avg_%'] / 100 * 1000
        # 1065.645-5
        time_aligned_data['Tdewdil_K'] = CFR1065.dewpoint_temp_K(time_aligned_data['pH2Odilscal_Pa'])
        time_aligned_data['Tdewdil_°C'] = time_aligned_data['Tdewdil_K'] - 273.15
        time_aligned_data['CVSDilAirDPTemp_°C'] = time_aligned_data['Tdewdil_°C']
        time_aligned_data['pH2Odil_kPa'] = CFR1065.vapor_pressure_of_water_kPa(time_aligned_data['Tdewdil_K'])

        # 1065.645-3:
        if test_type == 'transient':
            time_aligned_data['xH2Odil_mol/mol'] = \
                time_aligned_data['pH2Odil_kPa'] / time_aligned_data['pCellAmbient_Avg_kPa']
        else:
            time_aligned_data['xH2Odil_mol/mol'] = (
                    time_aligned_data['CVSDilAirRH_Avg_%'] / 100 *
                    time_aligned_data['pH2Odilsat_kPa'] / time_aligned_data['pCellAmbient_Avg_kPa'])
    else:
        # calculations based on dew point
        time_aligned_data['Tdewdil_K'] = time_aligned_data['CVSDilAirDPTemp_Avg_°C'] + 273.15
        # 1065.645-1:
        time_aligned_data['pH2Odilsat_kPa'] = CFR1065.vapor_pressure_of_water_kPa(time_aligned_data['Tdewdil_K'])
        if test_type == 'transient':
            # TODO: verify this for transient non-RH-based calculations
            # 1065.645-3:
            time_aligned_data['xH2Odil_mol/mol'] = (time_aligned_data['pH2Odilsat_kPa'] /
                                                    time_aligned_data['pCellAmbient_Avg_kPa'])
        else:
            # 1065.645-3:
            time_aligned_data['xH2Odil_mol/mol'] = (time_aligned_data['pH2Odilsat_kPa'] /
                                                    time_aligned_data['pCellAmbient_Avg_kPa'])

    from constants import constants, update_constants
    update_constants()  # update constants that rely on test fuel properties, etc

    if unit_rate == '/s':
        greek_units = 'g/h'
    else:
        greek_units = 'g'

    # 1065.655-20:
    time_aligned_data['alpha'] = \
        CFR1065.alpha(time_aligned_data['qmFuel_Avg_%s' % greek_units],
                      time_aligned_data['DEFMassFlowRate_Avg_%s' % greek_units], units=greek_units)

    # 1065.655-21:
    time_aligned_data['beta'] = \
        CFR1065.beta(time_aligned_data['qmFuel_Avg_%s' % greek_units],
                     time_aligned_data['DEFMassFlowRate_Avg_%s' % greek_units], units=greek_units)

    # 1065.655-23:
    time_aligned_data['delta'] = \
        CFR1065.delta(time_aligned_data['qmFuel_Avg_%s' % greek_units],
                      time_aligned_data['DEFMassFlowRate_Avg_%s' % greek_units], units=greek_units)

    time_aligned_data['gamma'] = 0  # no gamma for now

    time_aligned_data['Tint_K'] = time_aligned_data['tIntakeAir_Avg_°C'] + 273.15
    time_aligned_data['Tdewint_°C'] = time_aligned_data['tCellDewPt_Avg_°C']

    # 1065.645-1:
    time_aligned_data['Tdewint_K'] = time_aligned_data['Tdewint_°C'] + 273.15

    # 1065.645-1:
    time_aligned_data['pH2Oamb_kPa'] = \
        CFR1065.vapor_pressure_of_water_kPa(time_aligned_data['Tdewint_K'])

    # 1065.645-3:
    time_aligned_data['xH2Oint_mol/mol'] = \
        time_aligned_data['pH2Oamb_kPa'] / time_aligned_data['pCellAmbient_Avg_kPa']

    if calc_mode == 'raw':
        time_aligned_data['xH2Odil_mol/mol'] = time_aligned_data['xH2Oint_mol/mol']

    return time_aligned_data


def iterate_chemical_balance(time_aligned_data, calc_mode, emissions_cycle_number, drift_corrected=False):
    """
    Iterate the chemical balance equations

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        emissions_cycle_number (int): emissions cycle number to process
        drift_corrected (bool): use drift-correct background if ``True``

    Returns:
        Nothing, updates time_aligned_data

    """
    if not drift_corrected:
        BagData = phdp_globals.test_data['BagData']
    else:
        BagData = phdp_globals.test_data['drift_corrected_BagData']

    if calc_mode == 'raw':
        ctype = 'conRaw'
        COrange = 'H'
    else:
        ctype = 'con'
        COrange = 'L'

    if calc_mode == 'dilute-bag':
        # grab bag sample data for calcs
        time_aligned_data['conCO2_Avg_%vol'] = (
                BagData.loc[(BagData['RbComponent'] == 'CO2') &
                            (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number),
                'RbSmpConc_ppm'].item() / 10000)

        time_aligned_data['conLCO_Avg_ppm'] = (
            BagData.loc[(BagData['RbComponent'] == 'CO') &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number),
            'RbSmpConc_ppm'].item())

        time_aligned_data['conNOX_Avg_ppm'] = (
            BagData.loc[(BagData['RbComponent'] == 'NOX') &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number),
            'RbSmpConc_ppm'].item())

        time_aligned_data['conN2O_Avg_ppm'] = (
            BagData.loc[(BagData['RbComponent'] == 'N2O') &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number),
            'RbSmpConc_ppm'].item())

        time_aligned_data['conTHC_Avg_ppmC'] = (
            BagData.loc[(BagData['RbComponent'] == 'THC') &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number),
            'RbSmpConc_ppm'].item())

        time_aligned_data['conCH4_Avg_ppm'] = (
            BagData.loc[(BagData['RbComponent'] == 'CH4') &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number),
            'RbSmpConc_ppm'].item())

    time_aligned_data['xDil/Exh_mol/mol'] = 0.8
    time_aligned_data['xH2Oexh_mol/mol'] = 2 * time_aligned_data['xH2Oint_mol/mol']
    time_aligned_data['xCcombdry_mol/mol'] = \
        time_aligned_data['%sCO2_Avg_%%vol' % ctype] / 100 + (time_aligned_data['%sTHC_Avg_ppmC' % ctype] +
                                                      time_aligned_data['%s%sCO_Avg_ppm' % (ctype, COrange)]) / 1e6
    time_aligned_data['xH2dry_μmol/mol'] = 0
    time_aligned_data['xint/exhdry_mol/mol'] = 0

    if calc_mode == 'dilute-bag':
        # special treatment for bag calculations:
        phdp_globals.test_data['EmsComponents']['ResidualH2O_%vol'] = time_aligned_data['xH2Oexh_mol/mol'].item() * 100

    time_aligned_data_prior = pd.DataFrame()
    converged = False
    iteration = 0
    while not converged:
        # 1065.655-11
        time_aligned_data['xH2Ointdry_mol/mol'] = \
            time_aligned_data['xH2Oint_mol/mol'] / (1 - time_aligned_data['xH2Oint_mol/mol'])

        # 1065.655-13
        time_aligned_data['xH2Odildry_mol/mol'] = \
            time_aligned_data['xH2Odil_mol/mol'] / (1 - time_aligned_data['xH2Odil_mol/mol'])

        # 1065.655-14
        residual_H2O_pctvol = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['InputName'] == '%s%sCO' % (ctype, COrange),
            'ResidualH2O_%vol'].item()

        time_aligned_data['xCOdry_μmol/mol'] = \
            time_aligned_data['%s%sCO_Avg_ppm' % (ctype, COrange)] / (1 - residual_H2O_pctvol / 100)

        # 1065.655-15
        residual_H2O_pctvol = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['InputName'] == 'conCO2', 'ResidualH2O_%vol'].item()

        time_aligned_data['xCO2dry_%'] = \
            time_aligned_data['%sCO2_Avg_%%vol' % ctype] / (1 - residual_H2O_pctvol / 100)

        # 1065.655-16
        residual_H2O_pctvol = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['InputName'] == 'conNOX', 'ResidualH2O_%vol'].item()
        time_aligned_data['xNOdry_μmol/mol'] = \
            time_aligned_data['%sNOX_Avg_ppm' % ctype] * 0.75 / (1 - residual_H2O_pctvol / 100)

        # 1065.655-17
        time_aligned_data['xNO2dry_μmol/mol'] = \
            time_aligned_data['%sNOX_Avg_ppm' % ctype] * 0.25 / (1 - residual_H2O_pctvol / 100)

        # 1065.655(c)(1)
        ambient_CO2_conc_ppm = \
            BagData.loc[(BagData['RbComponent'] == 'CO2') &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number), 'RbAmbConc_ppm'].item()

        # 1065.655(c)(1)
        time_aligned_data['xCO2intdry_μmol/mol'] = \
            ambient_CO2_conc_ppm / (1 - time_aligned_data['xH2Oint_mol/mol'])

        # 1065.655(c)(1)
        time_aligned_data['xCO2dildry_μmol/mol'] = \
            ambient_CO2_conc_ppm / (1 - time_aligned_data['xH2Odil_mol/mol'])

        # 1065.655-10
        time_aligned_data['xCO2int_μmol/mol'] = \
            time_aligned_data['xCO2intdry_μmol/mol'] / (1 + time_aligned_data['xH2Ointdry_mol/mol'])

        # 1065.655-12
        time_aligned_data['xCO2dil_μmol/mol'] = \
            time_aligned_data['xCO2dildry_μmol/mol'] / (1 + time_aligned_data['xH2Odildry_mol/mol'])

        # 1065.655-9
        time_aligned_data['xO2int_%'] = (0.20982 - (time_aligned_data['xCO2intdry_μmol/mol'] / 1e6)) / \
                                        (1 + time_aligned_data['xH2Ointdry_mol/mol'])

        # 1065.655-6
        time_aligned_data['xdil/exhdry_mol/mol'] = \
            time_aligned_data['xDil/Exh_mol/mol'] / (1 - time_aligned_data['xH2Oexh_mol/mol'])

        # 1065.655-18
        time_aligned_data['xTHCdry_μmol/mol'] = \
            time_aligned_data['%sTHC_Avg_ppmC' % ctype] / (1 - time_aligned_data['xH2Oexh_mol/mol'])

        time_aligned_data['xraw/exhdry_mol/mol'] = CFR1065.rawexhdry(time_aligned_data)

        time_aligned_data['xH2Oexhdry_mol/mol'] = CFR1065.xH2Oexhdry(time_aligned_data)

        time_aligned_data['xH2dry_μmol/mol'] = CFR1065.xH2exhdry(time_aligned_data)

        time_aligned_data['xint/exhdry_mol/mol'] = CFR1065.xintexhdry(time_aligned_data)

        # 1065.655-1
        time_aligned_data['xDil/Exh_mol/mol'] = \
            1 - time_aligned_data['xraw/exhdry_mol/mol'] / (1 + time_aligned_data['xH2Oexhdry_mol/mol'])

        # 1065.655-2
        time_aligned_data['xH2Oexh_mol/mol'] = \
            time_aligned_data['xH2Oexhdry_mol/mol'] / (1 + time_aligned_data['xH2Oexhdry_mol/mol'])

        # 1065.655-3
        time_aligned_data['xCcombdry_mol/mol'] = CFR1065.xccombdry(time_aligned_data)

        if calc_mode == 'dilute-bag':
            # special treatment for bag calculations:
            phdp_globals.test_data['EmsComponents']['ResidualH2O_%vol'] = (
                    time_aligned_data['xH2Oexh_mol/mol'].item() * 100)

        if iteration == 0:
            converged = False
        else:
            converged = ((time_aligned_data - time_aligned_data_prior).abs() <=
                         phdp_globals.options.chemical_balance_convergence_tolerance *
                         time_aligned_data.abs()).all().all() or iteration == 99

        time_aligned_data_prior = time_aligned_data.copy()

        if iteration > 10:
            print(emissions_cycle_number, iteration, time_aligned_data['xCcombdry_mol/mol'].iloc[0])

        iteration = iteration + 1

    print('Converged in %d iterations' % iteration)

    if iteration == 100:
        print('-- Possible failed convergence, convergence iteration limit exceeded --')

    if iteration > 10:
        print()


def post_chemical_balance_calculations(time_aligned_data, calc_mode):
    """
    Calculate values after chemical balance iteration

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.

    Returns:
        Nothing, updates time_aligned_data

    """
    from constants import constants

    # CFR 1065.640-9
    time_aligned_data['Mmix_intake_g/mol'] = \
        constants['Mair_g/mol'] * (1 - time_aligned_data['xH2Oint_mol/mol']) + \
        constants['MH2O_g/mol'] * time_aligned_data['xH2Oint_mol/mol']

    if calc_mode == 'dilute-bag':
        unit_rate = ''
    else:
        time_aligned_data['nint_mol/s'] = \
            time_aligned_data['qmIntakeAir_Avg_kg/h'] / 3.6 / time_aligned_data['Mmix_intake_g/mol']

        unit_rate = '/s'

    if calc_mode == 'raw':
        ctype = 'conRaw'
        COrange = 'H'
        ems_prefix = 'Raw'
        testparam_prefix = 'Raw'

        # CFR 1065.655-24
        time_aligned_data['nexh_mol/s'] = \
            (time_aligned_data['nint_mol/s'] /
             (1 + (time_aligned_data['xint/exhdry_mol/mol'] - time_aligned_data['xraw/exhdry_mol/mol']) /
              (1 + time_aligned_data['xH2Oexhdry_mol/mol'])))

        flow_mol = time_aligned_data['nexh_mol/s']
    else:
        ctype = 'con'
        COrange = 'L'
        ems_prefix = 'Dil'
        testparam_prefix = 'Dilute'

        if calc_mode == 'dilute-bag':
            flow_mol = time_aligned_data['CVSFlow_mol']

            time_aligned_data['nint_mol'] = \
                time_aligned_data['qmIntakeAir_Avg_kg'] * 1000 / time_aligned_data['Mmix_intake_g/mol']

            time_aligned_data['nexh_mol'] = (time_aligned_data['xraw/exhdry_mol/mol'] - time_aligned_data[
                'xint/exhdry_mol/mol']) * flow_mol + time_aligned_data['nint_mol']
        else:
            flow_mol = time_aligned_data['CVSFlow_mol/s']

            # CFR 1065.655-26
            time_aligned_data['nexh_mol/s'] = \
                (time_aligned_data['xraw/exhdry_mol/mol'] - time_aligned_data['xint/exhdry_mol/mol']) * \
                (1 - time_aligned_data['xH2Oexh_mol/mol']) * flow_mol + time_aligned_data['nint_mol/s']

    if calc_mode == 'dilute-bag':
        time_aligned_data['xCO2exh_%mol'] = time_aligned_data['conCO2_Avg_%vol']
        time_aligned_data['xTHCexh_μmol/mol'] = time_aligned_data['conTHC_Avg_ppmC']
        time_aligned_data['xCH4exh_μmol/mol'] = time_aligned_data['conCH4_Avg_ppm']
        # CFR 1065.660-05
        time_aligned_data['xNMHCexh_μmol/mol'] = (
                time_aligned_data['xTHCexh_μmol/mol'] -
                phdp_globals.test_data['TestParameters']['BagRFCH4_Fraction'] * time_aligned_data['xCH4exh_μmol/mol'])
        time_aligned_data['xCOexh_μmol/mol'] = time_aligned_data['conLCO_Avg_ppm']
        time_aligned_data['xNOexh_μmol/mol'] = time_aligned_data['conNOX_Avg_ppm'] * 0.75
        time_aligned_data['xNO2exh_μmol/mol'] = time_aligned_data['conNOX_Avg_ppm'] * 0.25
    else:
        # CFR 1065.659-1
        xH2OCO2dilmeas = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['ParameterName'] == '%sCO2_System' %
            ems_prefix]['ResidualH2O_%vol'].item()

        time_aligned_data['xCO2exh_%mol'] = time_aligned_data['%sCO2_Avg_%%vol' % ctype] * \
                                            ((1 - time_aligned_data['xH2Oexh_mol/mol']) / (1 - xH2OCO2dilmeas / 100))
        # CFR 1065.659-1
        time_aligned_data['xTHCexh_μmol/mol'] = \
            time_aligned_data['%sTHC_Avg_ppmC' % ctype] - \
            phdp_globals.test_data['TestParameters']['Initial%sTHC_ppmC' % ems_prefix].item()

        # CFR 1065.660-9
        RFPFC2H6_Fraction = phdp_globals.test_data['TestParameters']['%sRFPFC2H6_Fraction' % testparam_prefix].item()
        RFCH4_Fraction = phdp_globals.test_data['TestParameters']['%sRFCH4_Fraction' % testparam_prefix].item()
        time_aligned_data['xCH4exh_μmol/mol'] = \
            (time_aligned_data['%sCH4cutter_Avg_ppmC' % ctype] -
             time_aligned_data['%sTHC_Avg_ppmC' % ctype] * RFPFC2H6_Fraction) / \
            (1 - RFPFC2H6_Fraction * RFCH4_Fraction)

        # CFR 1065.660-4
        time_aligned_data['xNMHCexh_μmol/mol'] = \
            (time_aligned_data['xTHCexh_μmol/mol'] -
             time_aligned_data['xCH4exh_μmol/mol'] * RFCH4_Fraction) / \
            (1 - RFPFC2H6_Fraction * RFCH4_Fraction)

        # CFR 1065.659-1
        xH2OCOdilmeas = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['ParameterName'] == '%sCO_System' %
            ems_prefix]['ResidualH2O_%vol'].item()

        time_aligned_data['xCOexh_μmol/mol'] = \
            time_aligned_data['%s%sCO_Avg_ppm' % (ctype, COrange)] * \
            ((1 - time_aligned_data['xH2Oexh_mol/mol']) / (1 - xH2OCOdilmeas / 100))

        # CFR 1065.659-1
        xH2ONOxdilmeas = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['ParameterName'] == '%sNOx_System' %
            ems_prefix]['ResidualH2O_%vol'].item()

        time_aligned_data['xNOexh_μmol/mol'] = \
            time_aligned_data['%sNOX_Avg_ppm' % ctype] * 0.75 * \
            ((1 - time_aligned_data['xH2Oexh_mol/mol']) / (1 - xH2ONOxdilmeas / 100))

        # CFR 1065.659-1
        time_aligned_data['xNO2exh_μmol/mol'] = \
            time_aligned_data['%sNOX_Avg_ppm' % ctype] * 0.25 * \
            ((1 - time_aligned_data['xH2Oexh_mol/mol']) / (1 - xH2ONOxdilmeas / 100))

    # CFR 1065.670-1
    time_aligned_data['xNOxcorrected_μmol/mol'] = \
        (time_aligned_data['xNOexh_μmol/mol'] + time_aligned_data['xNO2exh_μmol/mol']) * \
        (9.953 * time_aligned_data['xH2Oint_mol/mol'] + 0.832)

    # CFR 1065.650-5
    time_aligned_data['mCO2_g%s' % unit_rate] = \
        time_aligned_data['xCO2exh_%mol'] / 100 * constants['MCO2_g/mol'] * flow_mol

    # CFR 1065.650-5
    time_aligned_data['mCO_g%s' % unit_rate] = \
        time_aligned_data['xCOexh_μmol/mol'] / 1e6 * constants['MCO_g/mol'] * flow_mol

    # CFR 1065.650-5
    time_aligned_data['mCH4_g%s' % unit_rate] = \
        time_aligned_data['xCH4exh_μmol/mol'] / 1e6 * constants['MCH4_g/mol'] * flow_mol

    # CFR 1065.650-5
    time_aligned_data['mTHC_g%s' % unit_rate] = \
        time_aligned_data['xTHCexh_μmol/mol'] / 1e6 * constants['MTHC_g/mol'] * flow_mol

    # CFR 1065.650-5
    time_aligned_data['mNMHC_g%s' % unit_rate] = \
        time_aligned_data['xNMHCexh_μmol/mol'] / 1e6 * constants['MNMHC_g/mol'] * flow_mol

    # CFR 1065.650-5
    time_aligned_data['mNOxuncorrected_g%s' % unit_rate] = \
        (time_aligned_data['xNOexh_μmol/mol'] + time_aligned_data['xNO2exh_μmol/mol']) / 1e6 * \
        constants['MNOx_g/mol'] * flow_mol

    # CFR 1065.650-5
    time_aligned_data['mNOxcorrected_g%s' % unit_rate] = \
        (time_aligned_data['xNOxcorrected_μmol/mol']) / 1e6 * \
        constants['MNOx_g/mol'] * flow_mol

    # CFR 1065.650-5
    if '%sN2O_Avg_ppm' % ctype in time_aligned_data:
        time_aligned_data['mN2O_g%s' % unit_rate] = \
            time_aligned_data['%sN2O_Avg_ppm' % ctype] / 1e6 * constants['MN2O_g/mol'] * flow_mol
    else:
        time_aligned_data['mN2O_g%s' % unit_rate] = 0

    # CFR 1065.667(C)
    time_aligned_data['ndil_mol%s' % unit_rate] = (
            time_aligned_data['CVSFlow_mol%s' % unit_rate] - time_aligned_data['nexh_mol%s' % unit_rate])


def drift_correct_continuous_data(time_aligned_data, signal_name, emissions_cycle_number, test_name):
    """
    Perform drift correction for continuous data, per CFR1065.672-1

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        signal_name (str): signal name, e.g. 'conRawCO2_Avg_%vol'
        emissions_cycle_number (int): emissions cycle number to process
        test_name (str): test type name, e.g. 'FTP', 'RMC', etc

    Returns:
        Nothing, updates dataframe with drift corrected signal.

    """
    xrefzero = 0

    EmsCalResults = phdp_globals.test_data['EmsCalResults']

    EmsCalResults_ecns = EmsCalResults['EmissionsCycleNumber_Integer'].unique()

    if len(EmsCalResults_ecns) == 1:
        emscal_emissions_cycle_number = EmsCalResults_ecns[0]
    else:
        emscal_emissions_cycle_number = emissions_cycle_number

    component, driftline, scale_factor = handle_emscal_driftline(signal_name)

    xpre_data = (
        EmsCalResults)[(EmsCalResults['DriftComponent'] == component) & (EmsCalResults['DriftLine'] == driftline) &
                       (EmsCalResults['EmissionsCycleNumber_Integer'] == emscal_emissions_cycle_number)]

    if not xpre_data.empty:
        xrefspan = xpre_data['DriftSpanValue_ppm'].item()
        xprezero = xpre_data['DriftZero2Measured_ppm'].item()
        xprespan = xpre_data['DriftSpanMeasured_ppm'].item()

        DriftCheck = phdp_globals.test_data['DriftCheck']

        DriftCheck_ecns = DriftCheck['EmissionsCycleNumber_Integer'].unique()

        if len(DriftCheck_ecns) == 1:
            drift_check_emissions_cycle_number = DriftCheck_ecns[0]
        else:
            drift_check_emissions_cycle_number = emissions_cycle_number

        xpost_data = (
            DriftCheck)[(DriftCheck['DriftComponent'] == component) & (DriftCheck['DriftLine'] == driftline) &
                        (DriftCheck['EmissionsCycleNumber_Integer'] == drift_check_emissions_cycle_number)]

        xpostzero = xpost_data['DriftZeroMeasured_ppm'].item()
        xpostspan = xpost_data['DriftSpanMeasured_ppm'].item()

        time_aligned_data[signal_name] = xrefzero + (xrefspan - xrefzero) * (
                2 * time_aligned_data[signal_name] * scale_factor - (xprezero + xpostzero)) / (
                (xprespan + xpostspan) - (xprezero + xpostzero)) / scale_factor
    else:
        phdp_log.logwrite('*** Warning, no zero span zero data available for signal "%s", driftline "%s" ***' %
                          (signal_name, driftline))


def drift_correct_bag_data(bag_data, idx):
    """
    Perform drift correction for bag data, per CFR1065.672-1

    Args:
        bag_data (dataframe): emissions bag data
        idx (int): row index into bag_data

    Returns:
        Nothing, updates bag data with drift corrected signal.

    """

    xrefzero = 0

    component = bag_data.loc[idx, 'RbComponent']

    xrefspan = bag_data.loc[idx, 'RbSpanValue_ppm']
    xprezero = bag_data.loc[idx, 'RbZero2CalMeasured_ppm']
    xprespan = bag_data.loc[idx, 'RbSpanCalMeasured_ppm']

    DriftCheck = phdp_globals.test_data['BagDriftCheck']
    xpost_data = DriftCheck[(DriftCheck['DriftComponent'] == component) &
                            (DriftCheck['EmissionsCycleNumber_Integer'] ==
                             bag_data.loc[idx, 'EmissionsCycleNumber_Integer'])]

    xpostzero = xpost_data['DriftZeroMeasured_ppm'].item()
    xpostspan = xpost_data['DriftSpanMeasured_ppm'].item()

    bag_data.loc[idx, 'RbSmpConc_ppm'] = xrefzero + (xrefspan - xrefzero) * (
            2 * bag_data.loc[idx, 'RbSmpConc_ppm'] - (xprezero + xpostzero)) / (
            (xprespan + xpostspan) - (xprezero + xpostzero))

    bag_data.loc[idx, 'RbAmbConc_ppm'] = xrefzero + (xrefspan - xrefzero) * (
            2 * bag_data.loc[idx, 'RbAmbConc_ppm'] - (xprezero + xpostzero)) / (
            (xprespan + xpostspan) - (xprezero + xpostzero))


def calc_summary_results(time_aligned_data, calc_mode, emissions_cycle_number, drift_corrected=False):
    """
    Calculate summary results for the given time aligned data

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        emissions_cycle_number (int): emissions cycle number to process
        drift_corrected (bool): use drift-correct background if ``True``

    Returns:
        Summary results in a pandas Series

    """
    from constants import constants

    if calc_mode == 'raw':
        ctype = 'conRaw'
        testparam_prefix = 'Raw'
    elif calc_mode == 'dilute':
        ctype = 'con'
        testparam_prefix = 'Dilute'
    else:  # 'dilute-bag'
        ctype = 'con'
        testparam_prefix = 'Bag'

    # calculate summary values
    summary_results = pd.DataFrame(index=[0])
    summary_results['avg_xCO2exh_%mol'] = time_aligned_data['xCO2exh_%mol'].mean()
    summary_results['avg_xCOexh_μmol/mol'] = time_aligned_data['xCOexh_μmol/mol'].mean()
    summary_results['avg_xNOxcorrected_μmol/mol'] = time_aligned_data['xNOxcorrected_μmol/mol'].mean()
    summary_results['avg_xTHCexh_μmol/mol'] = time_aligned_data['xTHCexh_μmol/mol'].mean()
    summary_results['avg_xCH4exh_μmol/mol'] = time_aligned_data['xCH4exh_μmol/mol'].mean()
    if '%sN2O_Avg_ppm' % ctype in time_aligned_data:
        summary_results['avg_xN2O_μmol/mol'] = time_aligned_data['%sN2O_Avg_ppm' % ctype].mean()
    else:
        summary_results['avg_xN2O_μmol/mol'] = 0
    summary_results['avg_xNMHCexh_μmol/mol'] = time_aligned_data['xNMHCexh_μmol/mol'].mean()

    if calc_mode == 'dilute-bag':
        SamplePeriod_s = 1
        summary_results['cycle_work_kWh'] = (
                (time_aligned_data['Power_kW'] > 0) * time_aligned_data['Power_kW'] *
                time_aligned_data['elapsed_time_s'] / 3600)
        unit_rate = ''
    else:
        SamplePeriod_s = constants['SamplePeriod_s']
        summary_results['cycle_work_kWh'] = (
                ((time_aligned_data['Power_kW'] > 0) * time_aligned_data['Power_kW'] *
                 SamplePeriod_s).sum() / 3600)
        unit_rate = '/s'

    summary_results['total_dilute_flow_mol'] = time_aligned_data['ndil_mol%s' % unit_rate].sum() * SamplePeriod_s

    # calculate sample mass grams
    summary_results['mCO2_g'] = time_aligned_data['mCO2_g%s' % unit_rate].sum() * SamplePeriod_s
    summary_results['mCO_g'] = time_aligned_data['mCO_g%s' % unit_rate].sum() * SamplePeriod_s
    summary_results['mNOx_g'] = time_aligned_data['mNOxcorrected_g%s' % unit_rate].sum() * SamplePeriod_s
    summary_results['mTHC_g'] = time_aligned_data['mTHC_g%s' % unit_rate].sum() * SamplePeriod_s
    summary_results['mCH4_g'] = time_aligned_data['mCH4_g%s' % unit_rate].sum() * SamplePeriod_s
    summary_results['mN2O_g'] = time_aligned_data['mN2O_g%s' % unit_rate].sum() * SamplePeriod_s
    summary_results['mNMHC_g'] = time_aligned_data['mNMHC_g%s' % unit_rate].sum() * SamplePeriod_s
    summary_results['mNMHC_g+mNOx_g'] = summary_results['mNMHC_g'] + summary_results['mNOx_g']

    # calculate background mass grams
    from constants import constants

    if not drift_corrected:
        BagData = phdp_globals.test_data['BagData']
    else:
        BagData = phdp_globals.test_data['drift_corrected_BagData']

    if 'dilute' in calc_mode:
        # perform background emissions for dilute calcs
        for background_component in ['CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O']:
            background_conc = (
                BagData.loc[(BagData['RbComponent'] == str.upper(background_component)) &
                            (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number),
                'RbAmbConc_ppm'].item())
            summary_results['m%sbkgrnd_g' % background_component] = (
                    summary_results['total_dilute_flow_mol'] *
                    constants['M%s_g/mol' % background_component] *
                    background_conc / 10 ** 6)

        # handle NMHC mass grams
        THC_background_conc = (
            BagData.loc[(BagData['RbComponent'] == 'THC') &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number), 'RbAmbConc_ppm'].item())

        CH4_background_conc = (
            BagData.loc[(BagData['RbComponent'] == 'CH4') &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number), 'RbAmbConc_ppm'].item())
        DiluteRFCH4_Fraction = phdp_globals.test_data['TestParameters']['%sRFCH4_Fraction' % testparam_prefix].item()

        summary_results['mNMHCbkgrnd_g'] = \
            (summary_results['total_dilute_flow_mol'] * constants['MNMHC_g/mol'] *
             (THC_background_conc - CH4_background_conc * DiluteRFCH4_Fraction) / 10 ** 6)

        summary_results['mNMHCbkgrnd_g+mNOxbkgrnd_g'] = (
                summary_results['mNMHCbkgrnd_g'] + summary_results['mNOxbkgrnd_g'])

        # calculate net mass grams
        for component in ['CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O', 'NMHC']:
            summary_results['m%snet_g' % component] = \
                (summary_results['m%s_g' % component] -
                 summary_results['m%sbkgrnd_g' % component])
    else:
        # no background emissions for raw calcs
        for component in ['CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O', 'NMHC']:
            summary_results['m%sbkgrnd_g' % component] = None
        summary_results['mNMHCbkgrnd_g+mNOxbkgrnd_g'] = None

        for component in ['CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O', 'NMHC']:
            summary_results['m%snet_g' % component] = summary_results['m%s_g' % component]

    # handle NMHC + NOx net mass grams
    summary_results['mNMHCnet_g+mNOxnet_g'] = + summary_results['mNMHCnet_g'] + summary_results['mNOxnet_g']

    # calculate brake-specific net mass
    for component in ['CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O', 'NMHC']:
        if summary_results['cycle_work_kWh'].item() > 0:
            summary_results['m%snet_g/kWh' % component] = (
                    summary_results['m%snet_g' % component] /
                    summary_results['cycle_work_kWh'])
        else:
            summary_results['m%snet_g/kWh' % component] = 0

    # handle brake-specific NMHC + NOx net mass grams
    if summary_results['cycle_work_kWh'].item() > 0:
        summary_results['mNMHCnet_g/kWh+mNOxnet_g/kWh'] = (
                summary_results['mNMHCnet_g/kWh'] + summary_results['mNOxnet_g/kWh'])
    else:
        summary_results['mNMHCnet_g/kWh+mNOxnet_g/kWh'] = 0

    summary_results['elapsed_time_s'] = time_aligned_data['elapsed_time_s'].iloc[0]

    if 'stabilization_time_s' in time_aligned_data:
        summary_results['stabilization_time_s'] = time_aligned_data['stabilization_time_s'].iloc[0]

    return summary_results


def calc_1036_results(calc_mode, drift_corrected_time_aligned_data, drift_corrected_time_aligned_data_summary_results,
                      emissions_cycle_number, test_type, vehicle_test):
    """
    Calculate 1036 summary and carbon balance error check results

    Args:
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        drift_corrected_time_aligned_data (dataframe): drift-corrected, time-aligned continuous test data
        drift_corrected_time_aligned_data_summary_results (series): drift correct summary results
        emissions_cycle_number (int): emissions cycle number to process
        test_type (str): 'transient' or 'modal'
        vehicle_test (bool): ``True`` if test has an associated vehicle speed trace

    Returns:
        Pandas series of 1036 calculation results

    """

    from constants import constants

    calculations_1036 = pd.DataFrame(index=[0])

    SamplePeriod_s = constants['SamplePeriod_s']

    # rho_DEF_g_per_ml = 1.09  # REFERENCE / UNUSED ?

    EmfuelCmeas = phdp_globals.test_data['EngineData']['FuelLowHeatingValue_MJ/kg'].item()

    # CFR 1065.655-19
    wCmeas = constants['wcFuel']

    calculations_1036['CO2 Energy_Corr g/kWh'] = (
            drift_corrected_time_aligned_data_summary_results['mCO2net_g/kWh'] * EmfuelCmeas /
            constants['EmfuelCref	MJ/kg'] / wCmeas)

    BagData = phdp_globals.test_data['drift_corrected_BagData']

    # CFR 1036.550-1
    xCO2int_corr = BagData.loc[(BagData['RbComponent'] == 'CO2') & (
            BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number), 'RbAmbConc_ppm'].item()

    # CFR 1036.535-2
    if 'DEFMassFlowRate_Avg_g' in drift_corrected_time_aligned_data:
        mDEF_g = drift_corrected_time_aligned_data['DEFMassFlowRate_Avg_g']
    else:
        mDEF_g = drift_corrected_time_aligned_data[
                 'DEFMassFlowRate_Avg_g/h'].sum() * SamplePeriod_s / 3600

    mdot_avg_CO2DEF = (mDEF_g * constants['MCO2_g/mol'] * constants['wCH4N2O_Mass Fraction of urea in DEF'] /
                       constants['MCH4N2O_g/mol'])

    if calc_mode == 'raw':
        flow_mol = drift_corrected_time_aligned_data['nexh_mol/s']
        # CFR 1036.535-4
        mfuel_g = drift_corrected_time_aligned_data['qmFuel_Avg_g/h'].sum() * SamplePeriod_s / 3600
        # CFR 1065.643-1
        nint_mol = drift_corrected_time_aligned_data['nint_mol/s'].sum() * SamplePeriod_s
    elif calc_mode == 'dilute':
        flow_mol = drift_corrected_time_aligned_data['CVSFlow_mol/s']
        # CFR 1036.535-4
        mfuel_g = drift_corrected_time_aligned_data['qmFuel_Avg_g/h'].sum() * SamplePeriod_s / 3600
        # CFR 1065.643-1
        nint_mol = drift_corrected_time_aligned_data['nint_mol/s'].sum() * SamplePeriod_s
    else:  # 'dilute-bag'
        flow_mol = drift_corrected_time_aligned_data['CVSFlow_mol']
        # CFR 1036.535-4
        mfuel_g = drift_corrected_time_aligned_data['qmFuel_Avg_g']
        nint_mol = drift_corrected_time_aligned_data['nint_mol']

    if calc_mode == 'dilute-bag':
        # CFR 1036.540-6
        mfuel_cycle = constants['MC_g/mol'] / wCmeas * (drift_corrected_time_aligned_data['xCcombdry_mol/mol'] /
                                                        (1 + drift_corrected_time_aligned_data['xH2Oexhdry_mol/mol'])
                                                        * drift_corrected_time_aligned_data['CVSFlow_mol'] -
                                                        mdot_avg_CO2DEF / constants['MCO2_g/mol'])
    else:
        # CFR 1036.540-5
        mfuel_term = (flow_mol * drift_corrected_time_aligned_data['xCcombdry_mol/mol'] /
                      (1 + drift_corrected_time_aligned_data['xH2Oexhdry_mol/mol']) * SamplePeriod_s)

        mfuel_cycle = constants['MC_g/mol'] / wCmeas * (sum(mfuel_term) - mdot_avg_CO2DEF / constants['MCO2_g/mol'])

    calculations_1036['mfuelcor_meas'] = mfuel_g * EmfuelCmeas / constants['EmfuelCref	MJ/kg'] / constants['wCref']

    # CFR 1036.535-4
    calculations_1036['mfuelcor_dil'] = mfuel_cycle * EmfuelCmeas / constants['EmfuelCref	MJ/kg'] / constants['wCref']

    if test_type == 'modal':
        # units g/s instead of g:
        calculations_1036['mfuelcor_meas'] = calculations_1036['mfuelcor_meas'] / SamplePeriod_s
        calculations_1036['mfuelcor_dil'] = calculations_1036['mfuelcor_dil'] / SamplePeriod_s

    TestDetails = phdp_globals.test_data['TestDetails']

    if test_type == 'transient' and vehicle_test and calc_mode != 'dilute-bag':
        simulation_average_vehicle_speed_mps = \
            TestDetails[TestDetails['EmissionsCycleNumber_Integer'] == emissions_cycle_number][
                'CycleAverageVehicleSpeed_m/s'].item()

        calculations_1036['simulation_average_vehicle_speed_mps'] = simulation_average_vehicle_speed_mps

        calculations_1036['CycleAverageEngineWork_kWh'] = sum(
            drift_corrected_time_aligned_data['Power_kW'] * drift_corrected_time_aligned_data[
                'VehicleMoving_Logical'] * (
                    drift_corrected_time_aligned_data['Power_kW'] > 0)) * SamplePeriod_s / 3600

        calculations_1036['CycleAverageIdleSpeed_rpm'] = (
            (drift_corrected_time_aligned_data['spDyno_Avg_rev/min'][
                drift_corrected_time_aligned_data['VehicleMoving_Logical'] == 0]).mean())

        calculations_1036['CycleAverageTorque_Nm'] = (
            (drift_corrected_time_aligned_data['tqShaft_Avg_Nm'][
                drift_corrected_time_aligned_data['VehicleMoving_Logical'] == 0]).mean())

        calculations_1036['EngineToVehicleSpeedRatio_rev/mi'] = (
                (drift_corrected_time_aligned_data['spDyno_Avg_rev/min'][drift_corrected_time_aligned_data[
                                                                         'VehicleMoving_Logical'] == 1]).mean()
                / 60 / simulation_average_vehicle_speed_mps)

    # Carbon Balance Error Check calculations:
    # CFR 1065.643-6
    calculations_1036['mCexh_g'] = (
            constants['MC_g/mol'] *
            (drift_corrected_time_aligned_data_summary_results['mCO2net_g'] / constants['MCO2_g/mol'] +
             drift_corrected_time_aligned_data_summary_results['mCOnet_g'] / constants['MCO_g/mol'] +
             drift_corrected_time_aligned_data_summary_results['mTHCnet_g'] / constants['MTHC_g/mol']))

    calculations_1036['mCair_g'] = constants['MC_g/mol'] * nint_mol * xCO2int_corr / 10 ** 6

    # CFR 1065.643-9
    calculations_1036['mCfluidj_g'] = (
            constants['wcFuel'] * mfuel_g + constants['wcDef'] * mDEF_g)

    # CFR 1065.643-7
    calculations_1036['eaC_g'] = (
            calculations_1036['mCexh_g'] - calculations_1036['mCair_g'] -
            calculations_1036['mCfluidj_g'])

    # CFR 1065.643-8
    calculations_1036['erC_rel_err_%'] = (
            100 * calculations_1036['eaC_g'] /
            (calculations_1036['mCair_g'] + calculations_1036['mCfluidj_g']))

    # CFR 1065.643-8
    t_s = drift_corrected_time_aligned_data['elapsed_time_s'].iloc[0]

    calculations_1036['eaCrate_g/h'] = calculations_1036['eaC_g'] / t_s * 3600

    # limit from CFR 1065.543 (b)(2)(iii)
    if abs(calculations_1036['erC_rel_err_%'].item()) > 2.0:
        calculations_1036['erC_rel_err_%_check'] = 'FAIL'
    else:
        calculations_1036['erC_rel_err_%_check'] = 'Pass'

    # limit from CFR 1065.543-1
    calculations_1036['eaC_g_limit'] = (
        ASTM_round(phdp_globals.test_data['MapResults']['EngPeakPower_kW'] * 0.007, 3))

    if abs(calculations_1036['eaC_g'].item()) > calculations_1036['eaC_g_limit'].item():
        calculations_1036['eaC_g_check'] = 'FAIL'
    else:
        calculations_1036['eaC_g_check'] = 'Pass'

    # limit from CFR 1065.543-1
    calculations_1036['eaCrate_g/h_limit'] = (
        ASTM_round(phdp_globals.test_data['MapResults']['EngPeakPower_kW'] * 0.31, 3))

    if abs(calculations_1036['eaCrate_g/h'].item()) > calculations_1036['eaCrate_g/h_limit'].item():
        calculations_1036['eaCrate_g/h_check'] = 'FAIL'
    else:
        calculations_1036['eaCrate_g/h_check'] = 'Pass'

    return calculations_1036


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


def calc_STEYX(ref, meas, slope, intercept):
    """
    This function calculates the standard error of prediction (STEYX) given a reference dataset (ref),
    measured dataset (meas), a slope and an intercept. The function first calculates predicted values for
    the measurement based on the slope and intercept. Then, it computes the difference between the actual
    measurements and the predicted values. This difference is squared, summed up and divided by `(n - 2)`,
    where `n` is the length of the reference dataset. Finally, the result is square rooted to get the standard
    error of prediction.

    Arguments:
        ref (list): A list of reference values.
        meas (list): A list of corresponding measured values.
        slope (float): The slope value for the linear regression.
        intercept (float): The intercept value for the linear regression.

    Returns:
        std_error (float): The standard error of prediction.

    """
    n = len(ref)
    y_pred = slope * ref + intercept
    STEYX = (((meas - y_pred) ** 2).sum() / (n - 2)) ** 0.5

    return STEYX


def calc_stats(ref, meas):
    """
    Performs a linear regression analysis on two input data sets and returns the slope, intercept, correlation
    coefficient, R-squared value, and Standard Error of the Estimate.

    Args:
        ref (list or array like): A list or array like containing x-values.
        meas (list or array like): A list or array like containing corresponding y-values.

    Returns:
        dict: A dictionary containing the slope, intercept, R-squared value and Standard Error of the Estimate.

    """
    try:
        slope, intercept, r_value, _, _ = linregress(ref, meas)

        STEYX = calc_STEYX(ref, meas, slope, intercept)
    except:
        slope = intercept = STEYX = None
        r_value = 0

    return {'slope': slope, 'intercept': intercept, 'R2': r_value**2, 'SEE': STEYX}


def validate_data(test_name, test_type, output_prefix, emissions_cycles, modes=None, do_plots=False):
    """
    Validate measured data against reference cycle data

    Args:
        test_name (str): test type name, e.g. 'FTP', 'RMC', etc
        test_type (str): test type, i.e. 'transient' or 'modal'
        output_prefix (str): the prefix to be added to the filename of the generated report file.
        emissions_cycles (list): list of emission cycle numbers
        modes (list): list of emissions modes for modal tests
        do_plots (bool): generate regression plots if ``True``

    Returns:
        ``True`` if measured data passes validation regression test, else ``False``

    """
    # set up cycle data
    reference = dict()

    if test_name == 'GHGTRNS':
        cycledef_sample_period_s = 0.1
    else:
        cycledef_sample_period_s = 1.0

    cycle_valid = [False] * len(emissions_cycles)

    best_validation = {'lowest_fail_count': [math.inf] * len(emissions_cycles),
                       'description': [''] * len(emissions_cycles),
                       'validation_data': [None] * len(emissions_cycles),
                       'regression_results': dict(),
                       'PM_results': dict()}

    continuous_data = phdp_globals.test_data['ContinuousData'].copy()

    for ecn in emissions_cycles:
        linear_regression_fail = False

        if test_type == 'transient':
            if test_name == 'RMC':
                # generate 1 Hz RMC cycle reference data
                cycle_definition = phdp_globals.test_data['CycleDefinition']
                time_index = list(range(1, cycle_definition['CycleTime_s'].iloc[-1] + 1))

                speed_ramp_targets = list(cycle_definition['SpeedDemand_rpm'][1:])
                speed_ramp_targets.extend([cycle_definition['SpeedDemand_rpm'].iloc[-1]])

                cycledef_time_s = [1] + sorted(list(cycle_definition['CycleTime_s']) + list(
                    cycle_definition['CycleTime_s'] + cycle_definition['RampTime_s'].iloc[-1]))

                cycledef_speed_rpm = ([cycle_definition['SpeedDemand_rpm'].iloc[0]] +
                                      [val for pair in zip(cycle_definition['SpeedDemand_rpm'],
                                                           speed_ramp_targets) for val in pair])

                torque_ramp_targets = list(cycle_definition['TorqueDemand_Nm'][1:])
                torque_ramp_targets.extend([cycle_definition['TorqueDemand_Nm'].iloc[-1]])

                cycledef_torque_Nm = ([cycle_definition['TorqueDemand_Nm'].iloc[0]] +
                                      [val for pair in zip(cycle_definition['TorqueDemand_Nm'],
                                                           torque_ramp_targets) for val in pair])

                reference['speed_rpm'] = pd.Series(np.interp(time_index, cycledef_time_s, cycledef_speed_rpm))
                reference['torque_Nm'] = pd.Series(np.interp(time_index, cycledef_time_s, cycledef_torque_Nm))
            elif test_name == 'GHGTRNS':
                reference['speed_rpm'] = phdp_globals.test_data['CycleDefinition']['SpeedDemand_rpm']. \
                    loc[phdp_globals.test_data['CycleDefinition']['EmissionsCycleNumber_Integer'] == ecn]

                reference['torque_Nm'] = phdp_globals.test_data['CycleDefinition']['TorqueDemand_Nm']. \
                    loc[phdp_globals.test_data['CycleDefinition']['EmissionsCycleNumber_Integer'] == ecn]
            else:
                reference['speed_rpm'] = phdp_globals.test_data['CycleDefinition']['SpeedDemand_rpm']
                reference['torque_Nm'] = phdp_globals.test_data['CycleDefinition']['TorqueDemand_Nm']
        else:
            in_mode_pts = continuous_data['InModeLog_Logical'] == True

            reference['speed_rpm'] = pd.Series([np.nan] * len(continuous_data))
            reference['torque_Nm'] = pd.Series([np.nan] * len(continuous_data))

            for mode in modes:
                mode_pts = in_mode_pts & (continuous_data['ModeNumber_Integer'] == mode)
                reference['speed_rpm'].loc[mode_pts] = phdp_globals.test_data['CycleDefinition']['SpeedDemand_rpm'].iloc[mode-1]
                reference['torque_Nm'].loc[mode_pts] = phdp_globals.test_data['CycleDefinition']['TorqueDemand_Nm'].iloc[mode-1]

            reference['speed_rpm'] = reference['speed_rpm'].loc[in_mode_pts]
            reference['torque_Nm'] = reference['torque_Nm'].loc[in_mode_pts]

            continuous_data = continuous_data.loc[in_mode_pts].reset_index()

        reference['power_kW'] = reference['speed_rpm'] * reference['torque_Nm'] / 9548.8

        # set up limits
        warm_idle_rpm = phdp_globals.test_data['MapResults']['EngLowIdleSpeed_rev/min'].item()
        T_max_mapped_Nm = phdp_globals.test_data['MapResults']['EngMaxTorque_Nm'].item()
        P_max_mapped_kW = phdp_globals.test_data['MapResults']['EngPeakPower_kW'].item()
        maximum_test_speed_rpm = phdp_globals.test_data['MapResults']['EngNrefCFR1065_rev/min'].item()

        limits = dict({'speed_rpm': dict(), 'torque_Nm': dict(), 'power_kW': dict()})
        limits['speed_rpm']['slope'] = (0.95, 1.03)
        limits['speed_rpm']['intercept'] = (-0.1 * warm_idle_rpm, 0.1 * warm_idle_rpm)
        limits['speed_rpm']['R2'] = (0.97, np.inf)
        limits['speed_rpm']['SEE'] = (-np.inf, 0.05 * maximum_test_speed_rpm)

        limits['torque_Nm']['slope'] = (0.83, 1.03)
        limits['torque_Nm']['intercept'] = (-0.02 * T_max_mapped_Nm, 0.02 * T_max_mapped_Nm)
        limits['torque_Nm']['R2'] = (0.85, np.inf)
        limits['torque_Nm']['SEE'] = (-np.inf, 0.1 * T_max_mapped_Nm)

        limits['power_kW']['slope'] = (0.83, 1.03)
        limits['power_kW']['intercept'] = (-0.02 * P_max_mapped_kW, 0.02 * P_max_mapped_kW)
        limits['power_kW']['R2'] = (0.91, np.inf)
        limits['power_kW']['SEE'] = (-np.inf, 0.1 * P_max_mapped_kW)

        SamplePeriod_s = phdp_globals.test_data['TestParameters']['ContinuousLoggerPeriod_s'].item()

        # step by 10s when reference is 1Hz, step by 1s when data and reference are 10Hz:
        sample_index_stepsize = int(1 / (SamplePeriod_s / cycledef_sample_period_s))

        regression_results = dict()

        if test_type == 'transient':
            shift_range = range(0, int(1.0/SamplePeriod_s + 1))
        else:  # no need to shift modal data
            shift_range = [0]

        for shift in shift_range:
            time_shift = shift * SamplePeriod_s

            may_omit_options = list(range(0, 2**5))

            for may_omit_option in may_omit_options:
                may_omit_str = format(may_omit_option, '05b')

                validation_data = dict()

                if test_type == 'transient':
                    if test_name == 'GHGTRNS':
                        start_condition = (continuous_data['EmissionsCycleNumber_Integer'] == ecn)
                    else:
                        start_condition = ((continuous_data['ModeNumber_Integer'] == 1) &
                                           (continuous_data['EmissionsCycleNumber_Integer'] == ecn))

                    if test_name == 'RMC':
                        start_condition_index = 9
                    elif test_name == 'GHGTRNS':
                        start_condition_index = 0
                    else:
                        start_condition_index = -1

                    end_condition = ((continuous_data['ModeNumber_Integer'] == -1) &
                                     (continuous_data['EmissionsCycleNumber_Integer'] == ecn))
                    end_condition_index = 0
                else:
                    start_condition = [True] * len(continuous_data)
                    start_condition_index = 0

                    end_condition = [True] * len(continuous_data)
                    end_condition_index = -1

                start_index = (
                    continuous_data.loc[start_condition, :].index)[start_condition_index]
                max_index = continuous_data.index[-1]
                end_index = (
                    min(max_index, continuous_data.loc[end_condition, :]
                        .index[end_condition_index]))

                # select throttle data, but don't shift it:
                validation_data['measured_throttle_pct'] = (
                    continuous_data['pctThrottle_Avg_%'].loc[
                    start_index:end_index:sample_index_stepsize].values)

                if len(validation_data['measured_throttle_pct']) < len(reference['speed_rpm']):
                    # append last available data point if test data is cut short
                    validation_data['measured_throttle_pct'] = np.append(validation_data['measured_throttle_pct'],
                              continuous_data['pctThrottle_Avg_%'].iloc[-1])

                # allow speed and torque shift:
                start_index += shift
                max_index = continuous_data.index[-1]
                end_index = min(max_index, end_index + shift)

                validation_data['mode_number'] = (
                    continuous_data['ModeNumber_Integer'].loc[
                    start_index:end_index:sample_index_stepsize].values)

                if len(validation_data['mode_number']) < len(reference['speed_rpm']):
                    # append last available data point if test data is cut short
                    validation_data['mode_number'] = np.append(validation_data['mode_number'],
                              continuous_data['ModeNumber_Integer'].iloc[-1])

                validation_data['reference_speed_rpm'] = reference['speed_rpm']

                validation_data['measured_speed_rpm'] = (
                    continuous_data['spDyno_Avg_rev/min'].loc[
                    start_index:end_index:sample_index_stepsize].values)

                if len(validation_data['measured_speed_rpm']) < len(reference['speed_rpm']):
                    # append last available data point if test data is cut short
                    validation_data['measured_speed_rpm'] = np.append(validation_data['measured_speed_rpm'],
                              continuous_data['spDyno_Avg_rev/min'].iloc[-1])

                validation_data['reference_torque_Nm'] = reference['torque_Nm']

                validation_data['measured_torque_Nm'] = (
                    continuous_data['tqShaft_Avg_Nm'].loc[
                    start_index:end_index:sample_index_stepsize].values)

                if len(validation_data['measured_torque_Nm']) < len(reference['speed_rpm']):
                    # append last available data point if test data is cut short
                    validation_data['measured_torque_Nm'] = np.append(validation_data['measured_torque_Nm'],
                              continuous_data['tqShaft_Avg_Nm'].iloc[-1])

                validation_data['reference_power_kW'] = reference['power_kW']

                validation_data['measured_power_kW'] = (
                        validation_data['measured_speed_rpm'] * validation_data['measured_torque_Nm'] / 9548.8)

                # from CFR, torque/power check "Motoring point":
                motoring_at_min_demand = ((validation_data['measured_throttle_pct'] <= 1.0) &
                            (reference['torque_Nm'] < 0))
                validation_data['motoring_at_min_demand'] = motoring_at_min_demand

                # from CFR, speed/power check:
                idling_at_min_demand = ((1 & may_omit_option) != 0) & (
                        (validation_data['measured_throttle_pct'] <= 1.0) &
                        (reference['speed_rpm'] == warm_idle_rpm) & (reference['torque_Nm'] == 0) &
                        ((reference['torque_Nm'] - 0.02 * T_max_mapped_Nm) < validation_data['measured_torque_Nm']) &
                        (validation_data['measured_torque_Nm'] < (reference['torque_Nm'] + 0.02 * T_max_mapped_Nm))
                )
                validation_data['idling_at_min_demand'] = idling_at_min_demand

                # from CFR, torque/power check "Torque/Power No Load, Torque > Reference":
                somewhat_above_torque_at_min_demand = ((2 & may_omit_option) != 0) & (
                        (validation_data['measured_throttle_pct'] <= 1.0) &
                        (validation_data['measured_torque_Nm'] > reference['torque_Nm']) &
                        ~((validation_data['measured_speed_rpm'] > reference['speed_rpm'] * 1.02) &
                          (validation_data['measured_torque_Nm'] > (reference['torque_Nm'] + 0.02 * T_max_mapped_Nm)))
                )
                validation_data['somewhat_above_torque_at_min_demand'] = somewhat_above_torque_at_min_demand

                # from CFR, speed/power check "Speed/Power No Load, Speed > Reference"?:
                somewhat_above_speed_at_min_demand = ((4 & may_omit_option) != 0) & (
                        (validation_data['measured_throttle_pct'] <= 1.0) &
                        (validation_data['measured_speed_rpm'] > reference['speed_rpm']) &
                         ~somewhat_above_torque_at_min_demand &
                        ~((validation_data['measured_speed_rpm'] > reference['speed_rpm'] * 1.02) &
                          (validation_data['measured_torque_Nm'] > (reference['torque_Nm'] + 0.02 * T_max_mapped_Nm)))
                )
                validation_data['somewhat_above_speed_at_min_demand'] = somewhat_above_speed_at_min_demand

                # from CFR, torque/power check "Torque/Power Full Load, Torque < Reference":
                somewhat_below_torque_at_max_demand = ((8 & may_omit_option) != 0) & (
                        (validation_data['measured_throttle_pct'] >= 99.0) &
                        (validation_data['measured_torque_Nm'] < reference['torque_Nm']) &
                        ~((validation_data['measured_speed_rpm'] < reference['speed_rpm'] * 0.98) &
                          (validation_data['measured_torque_Nm'] < (reference['torque_Nm'] - 0.02 * T_max_mapped_Nm)))
                )
                validation_data['somewhat_below_torque_at_max_demand'] = somewhat_below_torque_at_max_demand

                # from CFR, speed/power check "Speed/Power Full Load, Speed < Reference":
                somewhat_below_speed_at_max_demand = ((16 & may_omit_option) != 0) & (
                        (validation_data['measured_throttle_pct'] >= 99.0) &
                        (validation_data['measured_speed_rpm'] < reference['speed_rpm']) &
                         ~somewhat_below_torque_at_max_demand &
                        ~((validation_data['measured_speed_rpm'] < reference['speed_rpm'] * 0.98) &
                          (validation_data['measured_torque_Nm'] < (reference['torque_Nm'] - 0.02 * T_max_mapped_Nm)))
                )
                validation_data['somewhat_below_speed_at_max_demand'] = somewhat_below_speed_at_max_demand

                validation_data['speed_rpm_omittable'] = (
                        idling_at_min_demand | somewhat_above_speed_at_min_demand | somewhat_below_speed_at_max_demand)

                validation_data['torque_Nm_omittable'] = (
                        motoring_at_min_demand | somewhat_above_torque_at_min_demand |
                        somewhat_below_torque_at_max_demand)

                validation_data['power_kW_omittable'] = (
                        validation_data['speed_rpm_omittable'] | validation_data['torque_Nm_omittable'])

                pass_fail = dict()
                regression_results['%d-%.1f-%s' % (ecn, time_shift, may_omit_str)] = dict()
                fail_count = 0
                for stp in ['speed_rpm', 'torque_Nm', 'power_kW']:
                    # print(stp)
                    ref = reference[stp].loc[~validation_data['%s_omittable' % stp]]
                    meas = validation_data['measured_%s' % stp][~validation_data['%s_omittable' % stp]]

                    non_ref = reference[stp].loc[validation_data['%s_omittable' % stp]]
                    non_meas = validation_data['measured_%s' % stp][validation_data['%s_omittable' % stp]]

                    stats = calc_stats(ref, meas)

                    if stats['slope'] is None:
                        pass_fail[stp] = 'FAIL'
                        linear_regression_fail = True
                    else:
                        pass_fail[stp] = all([pass_fail_range(stats[k], limits[stp][k]) == 'pass' for k in stats])

                    if do_plots:
                        fig, ax1 = plt.subplots()
                        ax1.plot(ref, meas, 'b.', label='regression points')
                        ax1.plot(non_ref, non_meas, 'rx', label='omitted')
                        ax1.grid(True, which='both')
                        ax1.set_xlabel('reference', fontsize=9)
                        ax1.set_ylabel('measured', fontsize=9)
                        title_str = '%.1f-%s-%s-%s-regression_shift.png' % (time_shift, may_omit_str, stp,
                                                                            {True: 'pass', False: 'FAIL'}
                                                                            [pass_fail[stp]])
                        ax1.set_title(title_str)
                        ax1.legend()
                        fig.savefig(phdp_globals.options.output_folder + output_prefix + '-' + title_str)
                        plt.close(fig)

                    fail_count += sum([int(pass_fail_range(stats[k], limits[stp][k]) == 'FAIL') for k in stats])

                    regression_results['%d-%.1f-%s' % (ecn, time_shift, may_omit_str)].update({
                        '%s_Slope' % stp: stats['slope'],
                        '%s_Slope_limit_min' % stp: limits[stp]['slope'][0],
                        '%s_Slope_limit_max' % stp: limits[stp]['slope'][1],
                        '%s_SlopeOK' % stp: pass_fail_range(stats['slope'], limits[stp]['slope']) == 'pass',
                        '%s_Intercept' % stp: stats['intercept'],
                        '%s_Intercept_limit_min' % stp: limits[stp]['intercept'][0],
                        '%s_Intercept_limit_max' % stp: limits[stp]['intercept'][1],
                        '%s_InterceptOK' % stp: pass_fail_range(stats['intercept'], limits[stp]['intercept']) == 'pass',
                        '%s_Rsq' % stp: stats['R2'],
                        '%s_Rsq_limit_min' % stp: limits[stp]['R2'][0],
                        '%s_Rsq_limit_max' % stp: limits[stp]['R2'][1],
                        '%s_RsqOK' % stp: pass_fail_range(stats['R2'], limits[stp]['R2']) == 'pass',
                        '%s_StdErr' % stp: stats['SEE'],
                        '%s_StdErr_limit_min' % stp: limits[stp]['SEE'][0],
                        '%s_StdErr_limit_max' % stp: limits[stp]['SEE'][1],
                        '%s_StdErrOK' % stp: pass_fail_range(stats['SEE'], limits[stp]['SEE']) == 'pass',
                        '%s_Points' % stp: len(meas),
                        '%s_CheckFailCount' % stp: fail_count,
                        'time_shift': time_shift,
                        'Emissions Cycle Number': ecn,
                        'descriptor': '%d-%.1f-%s' % (ecn, time_shift, may_omit_str),
                    })

                if not cycle_valid[ecn-1]:
                    cycle_valid[ecn-1] = all([pass_fail[k] is True for k in pass_fail])

                if all([pass_fail[k] is True for k in pass_fail]):
                    print('$$$ Cycle %d %.1f %s time_shift PASS $$$\n' % (ecn, time_shift, may_omit_str))

                if fail_count < best_validation['lowest_fail_count'][ecn-1]:
                    best_validation['lowest_fail_count'][ecn-1] = fail_count
                    best_validation['description'][ecn-1] = '%d-%.1f-%s' % (ecn, time_shift, may_omit_str)
                    best_validation['validation_data'][ecn-1] = validation_data
                    best_validation['regression_results'][ecn] = regression_results['%d-%.1f-%s' %
                                                                                    (ecn, time_shift, may_omit_str)]

        df = pd.DataFrame(regression_results).transpose()
        df.to_csv(phdp_globals.options.output_folder + output_prefix +
                  '-cycle-%d-regression_shift.csv' % ecn, columns=sorted(df.columns))

        if linear_regression_fail:
            if test_type == 'transient':
                phdp_log.logwrite('!!! Emissions Cycle %d Validation Linear Regression Failed !!!\n' % ecn)
            else:
                phdp_log.logwrite('!!! Mode Number %d Validation Linear Regression Failed !!!\n' % ecn)

    for idx, validation_data in enumerate(best_validation['validation_data']):
        if validation_data is not None:
            pd.DataFrame(validation_data).to_csv(phdp_globals.options.output_folder + output_prefix +
                                                 '-validation_data_cycle_%s.csv' %
                                                 best_validation['description'][idx], index=False)

    print('CYCLE VALID = %s' % all(cycle_valid))

    return all(cycle_valid), best_validation


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

            if test_name in ('GHGTRNS'):
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
