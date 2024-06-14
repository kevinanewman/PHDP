"""
PHDP report generators

"""

import os

import numpy as np

path = os.path.dirname(os.path.abspath(__file__))

print('path is "%s"' % path)

from phdp import *

from modal_report import *
from transient_report import *
from general_report import *


def calc_emscalresults_drift_check(report_df, emissions_cycle_number, check_phase, limit_pct, driftline, components):
    """
    Calculate pass/fail drift check tests for data in EmsCalResults

    Args:
        report_df (DataFrame): the drift check report data
        emissions_cycle_number (int): emissions cycle number to process
        check_phase (str): 'PRE' or 'POST'
        limit_pct (float): The pass/fail limit tolerance, in percent
        driftline (str): the emissions drift line e.g. 'DILUTE', 'HOT', etc
        components ([str]): list of strings, one or more emissions components, e.g. ['CO2', 'NOX', ...]

    Returns:
        Nothing, updates ``report_df`` with measured and calculated values

    """
    source = phdp_globals.test_data['EmsCalResults']

    for component in components:
        for measurement_type, offset in zip(
                ['DriftSpanValue_ppm', 'DriftSpanRange_ppm', 'DriftZeroMeasured_ppm', 'DriftSpanMeasured_ppm',
                 'DriftZero2Measured_ppm'],
                [5, 1, 2, 6, 9]):

            row_select = ((source['DriftLine'] == driftline) & (source['DriftComponent'] == component)
                          & (source['EmissionsCycleNumber_Integer'] == emissions_cycle_number))

            if any(row_select):
                value = source.loc[row_select, measurement_type].item()

                if measurement_type == 'DriftSpanValue_ppm':
                    nominal_span_ppm = value

                set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), value, col_offset=offset)

                if measurement_type in ['DriftZeroMeasured_ppm', 'DriftZero2Measured_ppm', 'DriftSpanMeasured_ppm']:
                    if measurement_type in ['DriftZeroMeasured_ppm', 'DriftZero2Measured_ppm']:
                        offset_pct_span = value / nominal_span_ppm * 100
                    elif measurement_type == 'DriftSpanMeasured_ppm':
                        offset_pct_span = (value - nominal_span_ppm) / nominal_span_ppm * 100

                    set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), [offset_pct_span],
                                 col_offset=offset + 1)

                    pass_fail = pass_fail_range(offset_pct_span, [-limit_pct, limit_pct])
                    set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), pass_fail,
                                 col_offset=offset + 2)


def calc_driftcheck_drift_check(report_df, emissions_cycle_number, check_phase, limit_pct, driftline, components):
    """
    Calculate pass/fail drift check tests for data in DriftCheck

    Args:
        report_df (DataFrame): the drift check report data
        emissions_cycle_number (int): emissions cycle number to process
        check_phase (str): 'PRE' or 'POST'
        limit_pct (float): The pass/fail limit tolerance, in percent
        driftline (str): the emissions drift line e.g. 'DILUTE', 'HOT', etc
        components ([str]): list of strings, one or more emissions components, e.g. ['CO2', 'NOX', ...]

    Returns:
        Nothing, updates ``report_df`` with measured and calculated values

    """
    source = phdp_globals.test_data['DriftCheck']

    for component in components:
        for measurement_type, offset in zip(
                ['DriftSpanValue_ppm', 'DriftSpanRange_ppm', 'DriftZeroMeasured_ppm', 'DriftSpanMeasured_ppm'],
                [5, 1, 2, 6]):

            row_select = ((source['DriftLine'] == driftline) & (source['DriftComponent'] == component)
                          & (source['EmissionsCycleNumber_Integer'] == emissions_cycle_number))

            if any(row_select):
                value = source.loc[row_select, measurement_type].item()

                if measurement_type == 'DriftSpanValue_ppm':
                    nominal_span_ppm = value

                set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), value, col_offset=offset)

                if measurement_type in ['DriftZeroMeasured_ppm', 'DriftSpanMeasured_ppm']:
                    if measurement_type in ['DriftZeroMeasured_ppm']:
                        offset_pct_span = value / nominal_span_ppm * 100
                    elif measurement_type == 'DriftSpanMeasured_ppm':
                        offset_pct_span = (value - nominal_span_ppm) / nominal_span_ppm * 100

                    set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), [offset_pct_span],
                                 col_offset=offset + 1)

                    pass_fail = pass_fail_range(offset_pct_span, [-limit_pct, limit_pct])
                    set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), pass_fail,
                                 col_offset=offset + 2)


def calc_bagdata_drift_check(report_df, emissions_cycle_number, check_phase, limit_pct, driftline, components):
    """
    Calculate pass/fail drift check tests for data in BagData

    Args:
        report_df (DataFrame): the drift check report data
        emissions_cycle_number (int): emissions cycle number to process
        check_phase (str): 'PRE' or 'POST'
        limit_pct (float): The pass/fail limit tolerance, in percent
        driftline (str): the emissions drift line e.g. 'DILUTE', 'HOT', etc
        components ([str]): list of strings, one or more emissions components, e.g. ['CO2', 'NOX', ...]

    Returns:
        Nothing, updates ``report_df`` with measured and calculated values

    """
    source = phdp_globals.test_data['BagData']

    for component in components:
        for measurement_type, offset in zip(
                ['RbSpanValue_ppm', 'RbRange_ppm', 'RbZeroCalMeasured_ppm', 'RbSpanCalMeasured_ppm',
                 'RbZero2CalMeasured_ppm'],
                [5, 1, 2, 6, 9]):

            row_select = ((source['RbLine'] == driftline) & (source['RbComponent'] == component)
                          & (source['EmissionsCycleNumber_Integer'] == emissions_cycle_number))

            if any(row_select):
                value = source.loc[row_select, measurement_type].item()

                if measurement_type == 'RbSpanValue_ppm':
                    nominal_span_ppm = value

                set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), value, col_offset=offset)

                if measurement_type in ['RbZeroCalMeasured_ppm', 'RbZero2CalMeasured_ppm', 'RbSpanCalMeasured_ppm']:
                    if measurement_type in ['RbZeroCalMeasured_ppm', 'RbZero2CalMeasured_ppm']:
                        offset_pct_span = value / nominal_span_ppm * 100
                    elif measurement_type == 'RbSpanCalMeasured_ppm':
                        offset_pct_span = (value - nominal_span_ppm) / nominal_span_ppm * 100

                    set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), [offset_pct_span],
                                 col_offset=offset + 1)

                    pass_fail = pass_fail_range(offset_pct_span, [-limit_pct, limit_pct])
                    set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), pass_fail,
                                 col_offset=offset + 2)


def calc_bagdriftcheck_drift_check(report_df, emissions_cycle_number, check_phase, limit_pct, driftline, components):
    """
    Calculate pass/fail drift check tests for data in BagDriftCheck

    Args:
        report_df (DataFrame): the drift check report data
        emissions_cycle_number (int): emissions cycle number to process
        check_phase (str): 'PRE' or 'POST'
        limit_pct (float): The pass/fail limit tolerance, in percent
        driftline (str): the emissions drift line e.g. 'DILUTE', 'HOT', etc
        components ([str]): list of strings, one or more emissions components, e.g. ['CO2', 'NOX', ...]

    Returns:
        Nothing, updates ``report_df`` with measured and calculated values

    """
    source = phdp_globals.test_data['BagDriftCheck']

    for component in components:
        for measurement_type, offset in zip(
                ['DriftSpanValue_ppm', 'DriftSpanRange_ppm', 'DriftZeroMeasured_ppm', 'DriftSpanMeasured_ppm'],
                [5, 1, 2, 6]):

            row_select = ((source['DriftLine'] == driftline) & (source['DriftComponent'] == component) &
                          (source['EmissionsCycleNumber_Integer'] == emissions_cycle_number))

            if any(row_select):
                value = source.loc[row_select, measurement_type].item()

                if measurement_type == 'DriftSpanValue_ppm':
                    nominal_span_ppm = value

                set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), value, col_offset=offset)

                if measurement_type in ['DriftZeroMeasured_ppm', 'DriftSpanMeasured_ppm']:
                    if measurement_type in ['DriftZeroMeasured_ppm']:
                        offset_pct_span = value / nominal_span_ppm * 100
                    elif measurement_type == 'DriftSpanMeasured_ppm':
                        offset_pct_span = (value - nominal_span_ppm) / nominal_span_ppm * 100

                    set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), [offset_pct_span],
                                 col_offset=offset + 1)

                    pass_fail = pass_fail_range(offset_pct_span, [-limit_pct, limit_pct])
                    set_value_at(report_df, '%s_%s_%s' % (driftline, component, check_phase), pass_fail,
                                 col_offset=offset + 2)


def generate_driftcheck_report(report_filename, results, test_type, test_name):
    """
    Generate drift check reports, one per emissions cycle

    Args:
        report_filename (str): the Excel report filename, a new sheet will be added
        results (dict): a dictionary containing the results of the emissions calculations.
        test_type (str): test type, i.e. 'transient' or 'modal'
        test_name (str): test type name, e.g. 'FTP', 'RMC', etc

    """
    report_df = pd.read_csv(path + os.sep + 'report_templates' + os.sep + 'drift_check_report_template.csv',
                            encoding='UTF-8', header=None)

    report_df = report_df.fillna('')

    if test_type == 'modal':
        emissions_cycles = [1]
    else:
        emissions_cycles = [results['tadsummary'][i]['EmissionsCycleNumber_Integer'].iloc[0]
                            for i in range(0, len(results['tadsummary']))]

    for emissions_cycle_number in emissions_cycles:
        for check_phase, limit_pct in zip(['PRE', 'POST'], [1, 2]):
            if check_phase == 'PRE':
                calc_method = calc_emscalresults_drift_check
                if test_name == 'RMC':
                    ecn = 0  # for some reason the pre-check is cycle 0 for RMC, not cycle 1 like the rest of the tests
                else:
                    ecn = emissions_cycles[0]
            else:
                calc_method = calc_driftcheck_drift_check
                ecn = emissions_cycles[-1]

            driftline = 'DIRECT'
            components = ['COH', 'CO2', 'CH4', 'THC', 'NOX', 'O2']
            calc_method(report_df, ecn, check_phase, limit_pct, driftline, components)

            driftline = 'HOT'
            components = ['NH3']
            calc_method(report_df, ecn, check_phase, limit_pct, driftline, components)

            driftline = 'DILUTE'
            components = ['COL', 'CO2', 'CH4', 'THC', 'NOX', 'N2O']
            calc_method(report_df, ecn, check_phase, limit_pct, driftline, components)

            driftline = 'EGR'
            components = ['CO2']
            calc_method(report_df, ecn, check_phase, limit_pct, driftline, components)

        for check_phase, limit_pct in zip(['PRE', 'POST'], [1, 2]):
            if check_phase == 'PRE':
                calc_method = calc_bagdata_drift_check
            else:
                calc_method = calc_bagdriftcheck_drift_check

            driftline = 'BAG'
            components = ['COL', 'CO2', 'CH4', 'THC', 'NOX', 'N2O']
            calc_method(report_df, emissions_cycle_number, check_phase, limit_pct, driftline, components)

        with pd.ExcelWriter(
                report_filename,
                mode="a",
                engine="openpyxl",
                if_sheet_exists="replace",
        ) as writer:
            report_df.to_excel(writer, index=False, header=False,
                               sheet_name='Drift Check %d' % emissions_cycle_number)


def CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, value_name):
    """
    Update pre test check template CFR1065EMS-related values and perform pass/fail checks

    Args:
        report_df (DataFrame): the report template dataframe
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format
        row_name (str): name of the row from the report template
        row_select (Series): a Boolean vector that selects a row from ``cfrdata``
        value_name (str): the column name of the target value , e.g. 'Hangup_ppmC'

    Returns:
        Nothing, updates values in ``report_df``

    """
    cfr1065ems = phdp_globals.test_data['CFR1065EMS']

    if any(row_select):
        value = cfr1065ems[value_name].loc[row_select].item()
        set_value_at(report_df, row_name, [value], col_offset=2)
        set_value_at(report_df, row_name, pass_fail_range(value, [-0.5, 0.5]), col_offset=3)

        CFR1065_datetimes(report_df, test_datetime, cfr1065ems, row_name, row_select)


def CFR1065_datetimes(report_df, test_datetime, cfrdata, row_name, row_select):
    """
    Handle date/time and elapsed time values in the pre test check report

    Args:
        report_df (DataFrame): the report template dataframe
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format
        cfrdata (DataFrame): CFR1065XXX data
        row_name (str): name of the row from the report template
        row_select (Series): a Boolean vector that selects a row from ``cfrdata``

    Returns:
        Nothing, updates values in ``report_df``

    """
    from datetime import datetime, timedelta, timezone

    test_timestamp = phdp_globals.test_data['Header']['Time_Date'].iloc[0] * 24 * 3600

    filename_timestamp = datetime.strptime(test_datetime, '%Y%m%d%H%M').timestamp()

    time_offset = filename_timestamp - test_timestamp

    test_time = phdp_globals.test_data['Workstation']['Time_Date'].iloc[0].item() * 24 * 3600 + time_offset
    check_time = cfrdata['TestTime_Date'].loc[row_select].item() * 24 * 3600 + time_offset
    elapsed_time_secs = test_time - check_time

    set_value_at(report_df, row_name,
                 datetime.fromtimestamp(test_time).strftime('%m/%d/%Y %H:%M'), col_offset=5)
    set_value_at(report_df, row_name,
                 datetime.fromtimestamp(check_time).strftime('%m/%d/%Y %H:%M'), col_offset=6)
    set_value_at(report_df, row_name,
                 str(timedelta(seconds=elapsed_time_secs)).rsplit(':', 1)[0], col_offset=7)
    set_value_at(report_df, row_name,
                 pass_fail_range(elapsed_time_secs, [0, 8 * 3600]), col_offset=10)


def generate_pre_test_check_report(report_filename, test_datetime):
    """
        Generate pre test check report

        Args:
            report_filename (str): the Excel report filename, a new sheet will be added
            test_datetime (str): the date and time of the test in YYYMMDDhhmm format

    """
    report_df = pd.read_csv(path + os.sep + 'report_templates' + os.sep + 'pre_test_check_template.csv',
                            encoding='UTF-8', header=None)

    report_df = report_df.fillna('')

    pre_test = phdp_globals.test_data['PreTest'].iloc[0]

    set_value_at(report_df, 'Cell Ambient Temperature', pre_test['CellAmbTempAvg_°C'], col_offset=2)
    set_value_at(report_df, 'Cell Ambient Temperature', pass_fail_range(pre_test['CellAmbTempAvg_°C'].item(), [20, 30]), col_offset=5)

    set_value_at(report_df, 'PowerMap Baro Pressure Diff', pre_test['BarometricPressDiffAvg_kPa'], col_offset=2)
    set_value_at(report_df, 'PowerMap Baro Pressure Diff', pass_fail_range(pre_test['BarometricPressDiffAvg_kPa'].item(), [-5, 5]), col_offset=5)

    set_value_at(report_df, 'CVS Dilution Air', pre_test['CVSDilAirTempAvg_°C'], col_offset=2)
    set_value_at(report_df, 'CVS Dilution Air', pass_fail_range(pre_test['CVSDilAirTempAvg_°C'].item(), [20, 30]), col_offset=5)

    set_value_at(report_df, 'Eng Coolant Out Temperature', pre_test['EngCoolOutTempAvg_°C'], col_offset=2)

    if pd.notna(pre_test['AfterTreatmentTempAvg_°C'].item()):
        set_value_at(report_df, 'After Treatment Temperature', pre_test['AfterTreatmentTempAvg_°C'], col_offset=2)

    set_value_at(report_df, 'Oil Temperature', pre_test['OilGalleryTempAvg_°C'], col_offset=2)

    cfr1065ems = phdp_globals.test_data['CFR1065EMS']

    # LEAK checks
    check_name = 'LeakRate_%'
    row_select = (cfr1065ems['Line'] == 'DIRECT') & (cfr1065ems['Component'] == 'CO2')
    row_name = 'DIRECT CO2 Cold'
    CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DIRECT') & (cfr1065ems['Component'] == 'THC') & (cfr1065ems['ColdHot'] == 'Hot')
    row_name = 'DIRECT THC Hot'
    CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'HOT') & (cfr1065ems['Component'] == 'NH3')
    row_name = 'HOT NH3'
    CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DILUTE') & (cfr1065ems['Component'] == 'CO2')
    row_name = 'DILUTE CO2 Cold'
    CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DILUTE') & (cfr1065ems['Component'] == 'THC') & (cfr1065ems['ColdHot'] == 'Hot')
    row_name = 'DILUTE THC Hot'
    CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DILUTE') & (cfr1065ems['Component'] == 'N2O')
    row_name = 'DILUTE N2O Cold'
    CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, check_name)

    # CFR1065EMS HANGUP checks
    check_name = 'Hangup_ppmC'
    row_select = (cfr1065ems['Line'] == 'DILUTE') & (cfr1065ems['Component'] == 'THC') & (cfr1065ems['Check'] == 'HCHangup')
    row_name = 'DILUTE THC'
    CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DIRECT') & (cfr1065ems['Component'] == 'THC') & (cfr1065ems['Check'] == 'HCHangup')
    row_name = 'DIRECT THC'
    CFR1065EMS_checks(report_df, test_datetime, row_name, row_select, check_name)

    pass_fail_dict = {'Pass': 'pass', 'Fail': 'FAIL'}

    cfr1065cvs = phdp_globals.test_data['CFR1065CVS']

    if cfr1065cvs is not None:
        check_name = 'CheckPassFail'
        row_select = cfr1065cvs['Component'] == 'BSU Fill Line'
        row_name = 'BSU Fill Line Leak Check'
        set_value_at(report_df, row_name, pass_fail_dict[cfr1065cvs[check_name].loc[row_select].item()], col_offset=3)
        CFR1065_datetimes(report_df, test_datetime, cfr1065cvs, row_name, row_select)

        row_select = cfr1065cvs['Component'] == 'BSU Read Line'
        row_name = 'BSU Read Line Leak Check'
        set_value_at(report_df, row_name, pass_fail_dict[cfr1065cvs[check_name].loc[row_select].item()], col_offset=3)
        CFR1065_datetimes(report_df, test_datetime, cfr1065cvs, row_name, row_select)

        for bagpair in [1, 2, 3]:
            check_name = 'AmbBagHangup_ppmC'
            row_select = (cfr1065cvs['Check'] == 'HCHangup') & (cfr1065cvs['Component'] == 'THC') & (cfr1065cvs['HCHangupBagPair_Integer'] == bagpair)
            row_name = 'Amb Bag HC Hangup Test %d' % bagpair
            value = cfr1065cvs[check_name].loc[row_select].item()
            set_value_at(report_df, row_name, [value], col_offset=2)
            set_value_at(report_df, row_name, pass_fail_range(value, [-np.inf, 2]), col_offset=3)
            CFR1065_datetimes(report_df, test_datetime, cfr1065cvs, row_name, row_select)

            check_name = 'SampleBagHangup_ppmC'
            row_select = (cfr1065cvs['Check'] == 'HCHangup') & (cfr1065cvs['Component'] == 'THC') & (cfr1065cvs['HCHangupBagPair_Integer'] == bagpair)
            row_name = 'Sample Bag HC Hangup Test %d' % bagpair
            value = cfr1065cvs[check_name].loc[row_select].item()
            set_value_at(report_df, row_name, [value], col_offset=2)
            set_value_at(report_df, row_name, pass_fail_range(value, [-np.inf, 2]), col_offset=3)
            CFR1065_datetimes(report_df, test_datetime, cfr1065cvs, row_name, row_select)

    cfr1065pm = phdp_globals.test_data['CFR1065PM']

    if cfr1065pm is not None:
        check_name = 'CheckPassFail'
        row_name = 'DLS/PSU Leak Check'
        set_value_at(report_df, row_name, pass_fail_dict[cfr1065pm[check_name].item()], col_offset=3)
        CFR1065_datetimes(report_df, test_datetime, cfr1065pm, row_name, 0)

    with pd.ExcelWriter(
            report_filename,
            mode="w",
            engine="openpyxl",
    ) as writer:
        report_df.to_excel(writer, index=False, header=False,
                           sheet_name='PreTestCheck')


def generate_cycle_validation_report(report_filename, validation_results):
    """
    Generate cycle validation report

    Args:
        report_filename (str): the Excel report filename, a new sheet will be added
        validation_results (dict): dict of results from best cycle validation (min error, min omits, min shift)

    """
    report_df = pd.read_csv(path + os.sep + 'report_templates' + os.sep + 'cycle_validation_report_template.csv',
                            encoding='UTF-8', header=None)

    report_df = report_df.fillna('')

    pass_fail_dict = {True: 'pass', False: 'FAIL'}

    for ecn in validation_results['regression_results']:
        regression_results = validation_results['regression_results'][ecn]

        set_value_at(report_df, 'Cycle Validation', regression_results['descriptor'])
        set_value_at(report_df, 'Shift', regression_results['time_shift'])

        for validation_prefix, validation_name in (
                zip(['speed_rpm', 'torque_Nm', 'power_kW'], ['Speed', 'Torque', 'Power'])):
            for check_name, check_description in zip(
                    ['StdErr', 'Slope', 'Rsq', 'Intercept'],
                    ['Standard Error (SEE)', 'Slope (m)', 'R2', 'Intercept (b)']):
                row_name = '%s %s' % (validation_name, check_description)
                set_value_at(report_df, row_name,
                             [regression_results['%s_%s' % (validation_prefix, check_name)]])
                set_value_at(report_df, row_name,
                             [regression_results['%s_%s_limit_min' % (validation_prefix, check_name)]], col_offset=2)
                set_value_at(report_df, row_name,
                             [regression_results['%s_%s_limit_max' % (validation_prefix, check_name)]], col_offset=3)
                set_value_at(report_df, row_name,
                             pass_fail_dict[regression_results['%s_%sOK' % (validation_prefix, check_name)]], col_offset=5)
            for check_name, check_description in zip(
                    ['Points'], ['Records Used for Analysis (#)']):
                row_name = '%s %s' % (validation_name, check_description)
                set_value_at(report_df, row_name,
                             int(regression_results['%s_%s' % (validation_prefix, check_name)]))

            for check_name, check_description in zip(
                    ['CheckFailCount'], ['Validation']):
                row_name = '%s %s' % (validation_name, check_description)
                set_value_at(report_df, row_name,
                             pass_fail_dict[regression_results['%s_%s' % (validation_prefix, check_name)] == 0])

        with pd.ExcelWriter(
                report_filename,
                mode="a",
                engine="openpyxl",
                if_sheet_exists="replace",
        ) as writer:
            report_df.to_excel(writer, index=False, header=False,
                               sheet_name='Cycle Validation %d' % ecn)
