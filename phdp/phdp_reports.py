"""
PHDP report generators

"""

import os
path = os.path.dirname(os.path.abspath(__file__))

from phdp import *
from common.file_io import file_exists
from constants import constants

from metpy.calc import specific_humidity_from_dewpoint
from metpy.units import units

initial_report = True


def generate_transient_report(output_prefix, calc_mode, results, test_datetime, test_type, test_num, test_site,
                              vehicle_test, pm_mass_mg):
    """
    Generate a transient test report.

    Args:
        output_prefix (str): the prefix to be added to the filename of the generated report file.
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        results (dict): a dictionary containing the results of the emissions calculations.
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format.
        test_type (str): test type, i.e. 'transient' or 'modal'
        test_num (int): the number assigned to the test, e.g. '00139'
        test_site (str): the name of the site where the test was performed, e.g. 'HD02'
        vehicle_test (bool): ``True`` if test has an associated vehicle speed trace
        pm_mass_mg (float): particulate matter mass (mg)

    Returns:
        Report file name for use by subsequent reports

    """
    global initial_report

    for i in range(0, len(results['tadsummary'])):
        emissions_cycle_number = results['tadsummary'][i]['EmissionsCycleNumber_Integer'].iloc[0]

        # create report header
        report_df = pd.read_csv(path + os.sep + 'transient_report_template.csv', encoding='UTF-8', header=None)
        report_df = report_df.fillna('')

        set_value_at(report_df, 'Test Date',
                     '%s/%s/%s' % (test_datetime[4:6], test_datetime[6:8], test_datetime[0:4]))

        set_value_at(report_df, 'Test Cell', test_site)
        set_value_at(report_df, 'Test Number', test_num)
        set_value_at(report_df, 'Test Type', test_type)

        if max([results['tadsummary'][i]['EmissionsCycleNumber_Integer'].iloc[0]
                for i in range(len(results['tadsummary']))]) > 1:
            cycle_name = phdp_globals.test_data['TestDetails'].loc[
                phdp_globals.test_data['TestDetails']['EmissionsCycleNumber_Integer'] == emissions_cycle_number,
                'CycleName'].item()

            set_value_at(report_df, 'Cycle ID', cycle_name)
        else:
            set_value_at(report_df, 'Cycle ID', phdp_globals.test_data['TestDetails']['CycleName'].iloc[0])

        set_value_at(report_df, 'Calculation Mode', calc_mode)

        set_value_at(report_df, 'Original Concentration', results['tadsummary'][i].iloc[0, 0:7].values)
        set_value_at(report_df, 'Corrected Concentration', results['dctadsummary'][i].iloc[0, 0:7].values)
        set_value_at(report_df, 'Original Mass', results['tadsummary'][i].iloc[0, 9:17].values)

        set_value_at(report_df, 'Corrected Mass', results['dctadsummary'][i].iloc[0, 9:17].values)

        set_value_at(report_df, 'Original Background Mass', results['tadsummary'][i].iloc[0, 17:25].values)
        set_value_at(report_df, 'Corrected Background Mass', results['dctadsummary'][i].iloc[0, 17:25].values)

        original_net_mass = results['tadsummary'][i].iloc[0, 25:33].values
        corrected_net_mass = results['dctadsummary'][i].iloc[0, 25:33].values

        set_value_at(report_df, 'Original Net Mass', original_net_mass)
        set_value_at(report_df, 'Corrected Net Mass', corrected_net_mass)

        # if calc_mode == 'dilute':
        #     net_mass_drift_check_pct = ((original_net_mass - corrected_net_mass) /
        #                                 np.maximum(original_net_mass, sys.float_info.epsilon) * 100)
        #     set_value_at(report_df, 'Net Mass Drift Check %', net_mass_drift_check_pct)

        cycle_work_kWh = results['tadsummary'][i]['cycle_work_kWh']

        # TODO: PM calculations
        if pm_mass_mg is not None:
            set_value_at(report_df, 'Total PM Mass', [pm_mass_mg])
            set_value_at(report_df, 'BSPM', [pm_mass_mg / 1000 / cycle_work_kWh.item()])
        else:
            set_value_at(report_df, 'Total PM Mass', 'N/A')
            set_value_at(report_df, 'BSPM', 'N/A')

        set_value_at(report_df, 'Cycle Work', cycle_work_kWh)
        set_value_at(report_df, 'Total Dilution Flow', results['tadsummary'][i]['total_dilute_flow_mol'])

        original_brake_specific_emissions = results['tadsummary'][i].iloc[0, 33:41].values
        corrected_brake_specific_emissions = results['dctadsummary'][i].iloc[0, 33:41].values
        brake_specific_emissions_drift_check_pct = (
                (original_brake_specific_emissions - corrected_brake_specific_emissions) /
                np.maximum(original_brake_specific_emissions, sys.float_info.epsilon) * 100)

        set_value_at(report_df, 'Original Brake Specific Emissions', original_brake_specific_emissions)
        set_value_at(report_df, 'Corrected Brake Specific Emissions', corrected_brake_specific_emissions)
        set_value_at(report_df, 'Brake Specific Emissions Drift Check %', brake_specific_emissions_drift_check_pct)

        # 1065 Drift Check Pass/Fails
        pass_fail = pass_fail_range(brake_specific_emissions_drift_check_pct[0], [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=3)  # CO2

        co_value = min(brake_specific_emissions_drift_check_pct[1],
                       ((original_brake_specific_emissions[1] - corrected_brake_specific_emissions[1]))/3.5 * 100)
        pass_fail = pass_fail_range(co_value, [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=4)  # CO

        pass_fail = pass_fail_range(brake_specific_emissions_drift_check_pct[2], [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=5)  # NOx

        pass_fail = pass_fail_range(brake_specific_emissions_drift_check_pct[3], [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=6)  # HC

        pass_fail = pass_fail_range(brake_specific_emissions_drift_check_pct[6], [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=9)  # NMHC

        pass_fail = pass_fail_range(brake_specific_emissions_drift_check_pct[7], [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=10)  # NMHC+NOx

        set_value_at(report_df, 'ϵrC', results['1036_calculations'][i]['erC_rel_err_%'])
        set_value_at(report_df, 'ϵaC', results['1036_calculations'][i]['eaC_g'])
        set_value_at(report_df, 'ϵaCrate', results['1036_calculations'][i]['eaCrate_g/h'])

        set_value_at(report_df, 'ϵaC', '±%.3f' % results['1036_calculations'][i]['eaC_g_limit'].iloc[0],
                     col_offset=2)
        set_value_at(report_df, 'ϵaCrate', '±%.3f' % results['1036_calculations'][i]['eaCrate_g/h_limit'].iloc[0],
                     col_offset=2)

        set_value_at(report_df, 'ϵrC', results['1036_calculations'][i]['erC_rel_err_%_check'],
                     col_offset=3)
        set_value_at(report_df, 'ϵaC', results['1036_calculations'][i]['eaC_g_check'],
                     col_offset=3)
        set_value_at(report_df, 'ϵaCrate', results['1036_calculations'][i]['eaCrate_g/h_check'],
                     col_offset=3)

        set_value_at(report_df, 'CO2 Energy_corr', results['1036_calculations'][i]['CO2 Energy_Corr g/kWh'],
                     col_offset=2)
        set_value_at(report_df, 'mfuelcor', results['1036_calculations'][i]['mfuelcor_meas'])
        set_value_at(report_df, 'mfuelcor_dil', results['1036_calculations'][i]['mfuelcor_dil'])

        if vehicle_test and calc_mode != 'dilute-bag':
            set_value_at(report_df, 'Simulation Average Vehicle Speed',
                         results['1036_calculations'][i]['simulation_average_vehicle_speed_mps'], col_offset=2)
            set_value_at(report_df, 'CycleAverageEngineWork',
                         results['1036_calculations'][i]['CycleAverageEngineWork_kWh'], col_offset=2)
            set_value_at(report_df, 'CycleAverageIdleSpeed',
                         results['1036_calculations'][i]['CycleAverageIdleSpeed_rpm'], col_offset=2)
            set_value_at(report_df, 'CycleAverageIdleTorque',
                         results['1036_calculations'][i]['CycleAverageTorque_Nm'], col_offset=2)
            set_value_at(report_df, 'EngineToVehicleSpeedRatio',
                         results['1036_calculations'][i]['EngineToVehicleSpeedRatio_rev/mi'], col_offset=2)
        else:
            set_value_at(report_df, 'Simulation Average Vehicle Speed', 'NA', col_offset=2)
            set_value_at(report_df, 'CycleAverageEngineWork', 'NA', col_offset=2)
            set_value_at(report_df, 'CycleAverageIdleSpeed', 'NA', col_offset=2)
            set_value_at(report_df, 'CycleAverageIdleTorque', 'NA', col_offset=2)
            set_value_at(report_df, 'EngineToVehicleSpeedRatio', 'NA', col_offset=2)

        report_filename = (phdp_globals.options.output_folder_base + output_prefix + 'report.xlsx')

        if initial_report:
            mode = 'w'
            initial_report = False
        else:
            mode = 'a'

        with pd.ExcelWriter(
                report_filename,
                mode=mode,
                engine="openpyxl",
        ) as writer:
            report_df.to_excel(writer, index=False, header=False, sheet_name='%s Emissions %d' %
                                                                             (calc_mode, emissions_cycle_number))

    return report_filename


def generate_modal_report(output_prefix, calc_mode, results, test_datetime, test_type, test_num, test_site):
    """
    Generate a modal report and calculated mode-weighted mass and brake-specific emissions.

    Args:
        output_prefix (str): the prefix to be added to the filename of the generated report file.
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        results (dict): a dictionary containing the results of the emissions calculations.
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format.
        test_type (str): test type, i.e. 'transient' or 'modal'
        test_num (int): the number assigned to the test, e.g. '00139'
        test_site (str): the name of the site where the test was performed, e.g. 'HD02'

    Returns:
        Report file name for use by subsequent reports

    """
    # create report header
    report_df = pd.read_csv(path + os.sep + 'modal_report_template.csv', encoding='UTF-8', header=None)
    report_df = report_df.fillna('')

    set_value_at(report_df, 'Test Date', '%s/%s/%s' % (test_datetime[4:6], test_datetime[6:8], test_datetime[0:4]))
    set_value_at(report_df, 'Test Cell', test_site)
    set_value_at(report_df, 'Test Number', test_num)
    set_value_at(report_df, 'Test Type', test_type)
    set_value_at(report_df, 'Cycle ID', '%s' % phdp_globals.test_data['TestDetails']['CycleName'].iloc[0])
    set_value_at(report_df, 'Calculation Mode', calc_mode)

    for i in range(0, len(results['dctadsummary'])):
        mode_number = results['dctadsummary'][i]['ModeNumber_Integer'].iloc[0]
        col_offset = mode_number + 1

        if mode_number > 4:
            # extend report_df
            report_df.insert(mode_number+5, ' ' * mode_number, '')

        set_value_at(report_df, 'Mode', mode_number, col_offset=col_offset)

        set_value_at(report_df, 'Sample CO2', results['dctadsummary'][i]['mCO2_g'],
                     col_offset=col_offset)
        set_value_at(report_df, 'Sample CO', results['dctadsummary'][i]['mCO_g'],
                     col_offset=col_offset)
        set_value_at(report_df, 'Sample NOX', results['dctadsummary'][i]['mNOx_g'],
                     col_offset=col_offset)
        set_value_at(report_df, 'Sample HC', results['dctadsummary'][i]['mTHC_g'],
                     col_offset=col_offset)
        set_value_at(report_df, 'Sample CH4', results['dctadsummary'][i]['mCH4_g'],
                     col_offset=col_offset)
        set_value_at(report_df, 'Sample N2O', results['dctadsummary'][i]['mN2O_g'],
                     col_offset=col_offset)
        set_value_at(report_df, 'Sample NMHC', results['dctadsummary'][i]['mNMHC_g'],
                     col_offset=col_offset)
        set_value_at(report_df, 'Sample NMHC+NOx', results['dctadsummary'][i]['mNMHC_g+mNOx_g'],
                     col_offset=col_offset)

        if calc_mode == 'dilute':
            set_value_at(report_df, 'Background CO2', results['dctadsummary'][i]['mCO2bkgrnd_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Background CO', results['dctadsummary'][i]['mCObkgrnd_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Background NOX', results['dctadsummary'][i]['mNOxbkgrnd_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Background HC', results['dctadsummary'][i]['mTHCbkgrnd_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Background CH4', results['dctadsummary'][i]['mCH4bkgrnd_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Background N2O', results['dctadsummary'][i]['mN2Obkgrnd_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Background NMHC', results['dctadsummary'][i]['mNMHCbkgrnd_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Background NMHC+NOx', results['dctadsummary'][i]['mNMHCbkgrnd_g+mNOxbkgrnd_g'],
                         col_offset=col_offset)

            set_value_at(report_df, 'Net CO2', results['dctadsummary'][i]['mCO2net_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Net CO', results['dctadsummary'][i]['mCOnet_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Net NOX', results['dctadsummary'][i]['mNOxnet_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Net HC', results['dctadsummary'][i]['mTHCnet_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Net CH4', results['dctadsummary'][i]['mCH4net_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Net N2O', results['dctadsummary'][i]['mN2Onet_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Net NMHC', results['dctadsummary'][i]['mNMHCnet_g'],
                         col_offset=col_offset)
            set_value_at(report_df, 'Net NMHC+NOx', results['dctadsummary'][i]['mNMHCnet_g+mNOxnet_g'],
                         col_offset=col_offset)

        set_value_at(report_df, 'Mode Time', results['dctad'][i]['elapsed_time_s'],
                     col_offset=col_offset)
        set_value_at(report_df, 'Power', results['dctad'][i]['Power_kW'],
                     col_offset=col_offset)

        set_value_at(report_df, 'mfuelcor_meas', results['1036_calculations'][i]['mfuelcor_meas'],
                     col_offset=col_offset)
        set_value_at(report_df, 'mfuelcor_dil', results['1036_calculations'][i]['mfuelcor_dil'],
                     col_offset=col_offset)

        set_value_at(report_df, 'ϵrC', results['1036_calculations'][i]['erC_rel_err_%'],
                     col_offset=mode_number)
        set_value_at(report_df, 'ϵaC', results['1036_calculations'][i]['eaC_g'],
                     col_offset=mode_number)
        set_value_at(report_df, 'ϵaCrate', results['1036_calculations'][i]['eaCrate_g/h'],
                     col_offset=mode_number)

        set_value_at(report_df, 'ϵrC Limit', '±2',
                     col_offset=mode_number)
        set_value_at(report_df, 'ϵaC Limit', '±%.3f' % results['1036_calculations'][i]['eaC_g_limit'].iloc[0],
                     col_offset=mode_number)
        set_value_at(report_df, 'ϵaCrate Limit', '±%.3f' % results['1036_calculations'][i]['eaCrate_g/h_limit'].iloc[0],
                     col_offset=mode_number)

        set_value_at(report_df, 'ϵrC Check', results['1036_calculations'][i]['erC_rel_err_%_check'],
                     col_offset=mode_number)
        set_value_at(report_df, 'ϵaC Check', results['1036_calculations'][i]['eaC_g_check'],
                     col_offset=mode_number)
        set_value_at(report_df, 'ϵaCrate Check', results['1036_calculations'][i]['eaCrate_g/h_check'],
                     col_offset=mode_number)

    if 'WeightingFactor_Fraction' in phdp_globals.test_data['CycleDefinition']:
        weighted_power = 0
        for i in range(0, len(results['dctadsummary'])):
            weight = phdp_globals.test_data['CycleDefinition']['WeightingFactor_Fraction'].iloc[i]
            set_value_at(report_df, 'Weighting Factor', weight, col_offset=i+2)
            weighted_power += results['dctad'][i]['Power_kW'] * weight

        set_value_at(report_df, 'Power (kW)', weighted_power)

        # calculate weighted values
        for idx, signal in enumerate(('CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O', 'NMHC')):
            weighted_mass_emissions = 0
            for i in range(0, len(results['dctadsummary'])):
                weighted_mass_emissions += (
                        results['dctadsummary'][i]['m%s_g' % signal] *
                        phdp_globals.test_data['CycleDefinition']['WeightingFactor_Fraction'].iloc[i] /
                        results['dctad'][i]['elapsed_time_s'] * 3600)
            weighted_specific_emissions = weighted_mass_emissions / weighted_power
            set_value_at(report_df, 'Mass Emissions (g/h)', weighted_mass_emissions, col_offset=idx + 1)
            set_value_at(report_df, 'Specific Emsissions (g/kWh)', weighted_specific_emissions, col_offset=idx + 1)

    report_filename = (phdp_globals.options.output_folder_base + output_prefix + 'report.xlsx')

    with pd.ExcelWriter(
            report_filename,
            mode="w",
            engine="openpyxl",
    ) as writer:
        report_df.to_excel(writer, index=False, header=False, sheet_name='Emissions')

    return report_filename


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
    report_df = pd.read_csv(path + os.sep + 'drift_check_report_template.csv', encoding='UTF-8', header=None)
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


def set_average_min_max(report_df, dctad, value_name, signal_name, col_offset, scale=1):
    if signal_name in dctad:
        set_value_at(report_df, value_name, [dctad[signal_name].mean() * scale], col_offset=col_offset)
        set_value_at(report_df, value_name, [dctad[signal_name].min() * scale], col_offset=col_offset+1)
        set_value_at(report_df, value_name, [dctad[signal_name].max() * scale], col_offset=col_offset+2)

        return True
    else:
        return False


def generate_general_report(report_filename, calc_mode, results, test_type, test_datetime, test_site):
    """
    Generate general test report

    Args:
        report_filename (str): the Excel report filename, a new sheet will be added
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        results (dict): a dictionary containing the results of the emissions calculations.
        test_type (str): test type, i.e. 'transient' or 'modal'
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format.
        test_site (str): the name of the site where the test was performed, e.g. 'HD02'

    """
    report_df = pd.read_csv(path + os.sep + 'general_report_template.csv', encoding='UTF-8', header=None)
    report_df = report_df.fillna('')

    if test_type == 'modal':
        emissions_cycles = [1]
    else:
        emissions_cycles = [results['tadsummary'][i]['EmissionsCycleNumber_Integer'].iloc[0]
                            for i in range(0, len(results['tadsummary']))]

    for emissions_cycle_number in emissions_cycles:
        # write header
        set_value_at(report_df, 'Cell', test_site)
        set_value_at(report_df, 'Operator', phdp_globals.test_data['Workstation']['Operator'])
        set_value_at(report_df, 'Test Date', '%s/%s/%s' % (test_datetime[4:6], test_datetime[6:8], test_datetime[0:4]))
        set_value_at(report_df, 'Power Map ID', phdp_globals.test_data['MapResults']['PowerMapID'])
        set_value_at(report_df, 'Engine Model', phdp_globals.test_data['EngineData']['EngModel'])
        set_value_at(report_df, 'Engine Number', phdp_globals.test_data['EngineData']['EngNumber'])
        set_value_at(report_df, 'Engine Family', phdp_globals.test_data['Header']['EngFamily'])

        # drift corrected time-aligned data
        dctad = results['dctad'][emissions_cycle_number-1]
        dctadsummary = results['dctadsummary'][emissions_cycle_number-1]

        set_average_min_max(report_df, dctad, 'Barometric Pressure', 'pCellAmbient_Avg_kPa', col_offset=2)
        pass_fail = (dctad['pCellAmbient_Avg_kPa'].min() >= 92.68 and
                     dctad['pCellAmbient_Avg_kPa'].max() <= 102.68)
        set_value_at(report_df, 'Barometric Pressure', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        set_average_min_max(report_df, dctad, 'Intake Air Temperature', 'tIntakeAir_Avg_°C', col_offset=2)
        pass_fail = (dctad['tIntakeAir_Avg_°C'].min() >= 20 and
                     dctad['tIntakeAir_Avg_°C'].max() <= 30)
        set_value_at(report_df, 'Intake Air Temperature', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        set_average_min_max(report_df, dctad, 'Intake Air Pressure', 'IntakeAirPress_Avg_kPa', col_offset=2)
        set_average_min_max(report_df, dctad, 'Intake Air Dew Point', 'tCellDewPt_Avg_°C', col_offset=2)
        set_average_min_max(report_df, dctad, 'Fuel Flow', 'qmFuel_Avg_g/h', col_offset=2, scale=1/1000)
        set_average_min_max(report_df, dctad, 'Intake Air Flow', 'qmIntakeAir_Avg_kg/h', col_offset=2)
        set_average_min_max(report_df, dctad, 'DEF Flow', 'DEFMassFlowRate_Avg_g/h', col_offset=2, scale=1/1000)
        set_average_min_max(report_df, dctad, 'Exhaust Back-pressure', 'ExhaustBackPressure_kPa', col_offset=2)
        set_average_min_max(report_df, dctad, 'Fuel Pump Inlet Temperature', 'tFuel_Avg_°C', col_offset=2)
        set_average_min_max(report_df, dctad, 'Fuel Pump Supply Pressure', 'pFuelSupply_kPa', col_offset=2)
        set_average_min_max(report_df, dctad, 'Fuel Pump Return Pressure', 'pFuelReturn_kPa', col_offset=2)

        if set_average_min_max(report_df, dctad, 'CA Coolant Temperature', 'tCoolantCA_°C', col_offset=2):
            pass_fail = (dctad['tCoolantCA_°C'].min() >= 20 and
                         dctad['tCoolantCA_°C'].max() <= 30)
            set_value_at(report_df, 'CA Coolant Temperature', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        set_average_min_max(report_df, dctad, 'Coolant Temperature', 'tCoolantIn_°C', col_offset=2)
        # set_average_min_max(report_df, dctad, 'PM Trap Face Temperature', '?', col_offset=2)
        set_average_min_max(report_df, dctad, 'Oil Sump Temperature', 'tOilSump_°C', col_offset=2)
        set_average_min_max(report_df, dctad, 'NOx Humidity Correction (Kh)', 'Kh', col_offset=2)

        set_value_at(report_df, 'FuelType', phdp_globals.test_data['EngineData']['FuelTypeName'], col_offset=2)
        set_value_at(report_df, 'Fuel H/C Ratio',
                     phdp_globals.test_data['EngineData']['FuelHTCRAT_ratio'], col_offset=2)
        set_value_at(report_df, 'Fuel O/C Ratio',
                     phdp_globals.test_data['EngineData']['FuelOTCRAT_ratio'], col_offset=2)
        set_value_at(report_df, 'Fuel Low Heating Value',
                     phdp_globals.test_data['EngineData']['FuelLowHeatingValue_MJ/kg'], col_offset=2)

        # Emissions data

        EmsCalResults = phdp_globals.test_data['EmsCalResults']

        signals = [
            'conRawCO2_Avg_%vol', 'conRawHCO_Avg_ppm', 'conRawNOX_Avg_ppm', 'conRawTHC_Avg_ppmC',
            'conRawCH4cutter_Avg_ppmC', 'conRawO2_Avg_%vol', 'conRawNH3_Avg_ppm',
            'conCO2_Avg_%vol', 'conLCO_Avg_ppm', 'conNOX_Avg_ppm', 'conTHC_Avg_ppmC', 'conCH4cutter_Avg_ppmC',
            'conEGRCO2_Avg_%vol', 'conN2O_Avg_ppm',
        ]

        for signal in signals:
            component, driftline, scale_factor = handle_emscal_driftline(signal)
            emscal_data = (
                EmsCalResults)[
                (EmsCalResults['DriftComponent'] == component) & (EmsCalResults['DriftLine'] == driftline)]

            if not emscal_data.empty:
                signal_prefix = signal.rsplit('_')[0]
                set_value_at(report_df, signal_prefix, emscal_data['DriftZero2Range_ppm'] / scale_factor, col_offset=2)  # Range
                set_value_at(report_df, signal_prefix, emscal_data['DriftSpanValue_ppm'] / scale_factor, col_offset=3)  # Span
                five_pct_over = emscal_data['DriftSpanValue_ppm'].item() / scale_factor * 1.05
                set_value_at(report_df, signal_prefix, [five_pct_over], col_offset=4)  # 5% over

                if set_average_min_max(report_df, dctad, signal_prefix, signal, col_offset=5):
                    pass_fail = dctad[signal].max() < five_pct_over
                    set_value_at(report_df, signal_prefix, pass_fail_range(pass_fail, [True, True]), col_offset=8)

        # CVS tunnel data
        set_average_min_max(report_df, dctad, 'CVS Molar Flow', 'CVSMolarFlow_Avg_mol/s', col_offset=2)
        set_average_min_max(report_df, dctad, 'CVS Volume Flow', 'CVSMolarFlow_Avg_mol/s', col_offset=2, scale=0.024055)

        set_average_min_max(report_df, dctad, 'CVSDilAirTemp', 'CVSDilAirTemp_Avg_°C', col_offset=2)
        pass_fail = (dctad['CVSDilAirTemp_Avg_°C'].min() >= 20 and
                     dctad['CVSDilAirTemp_Avg_°C'].max() <= 30)
        set_value_at(report_df, 'CVSDilAirTemp', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        set_average_min_max(report_df, dctad, 'CVSDilExhTemp', 'CVSDilExhTemp_Avg_°C', col_offset=2)
        set_average_min_max(report_df, dctad, 'CVSDilAirDPTemp', 'CVSDilAirDPTemp_Avg_°C', col_offset=2)
        if 'CVSDilAirRH_Avg_%' in dctad:
            set_average_min_max(report_df, dctad, 'CVSDilAirRH', 'CVSDilAirRH_Avg_%', col_offset=2)

        dctad['cvs_dilution_air_humidity'] = \
            specific_humidity_from_dewpoint(dctad['pCellAmbient_Avg_kPa'].values * units.kPa,
                                            dctad['tCellDewPt_Avg_°C'].values * units.degC).to('g/kg').magnitude
        set_average_min_max(report_df, dctad, 'CVS Dilution Air Humidity', 'cvs_dilution_air_humidity', col_offset=2)

        if set_average_min_max(report_df, dctad, 'CVS Pressure at Exh Entry', 'pTailpipe_Avg_kPa', col_offset=2):
            pass_fail = (dctad['pTailpipe_Avg_kPa'].min() >= -1.2 and
                         dctad['pTailpipe_Avg_kPa'].max() <= 1.2)
            set_value_at(report_df, 'CVS Pressure at Exh Entry', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        if 'BagFillProportionality' in dctadsummary:
            set_value_at(report_df, 'Bag Fill Proportionality (SEE/mean)', dctadsummary['BagFillProportionality'], col_offset=2)
            set_value_at(report_df, 'Bag Fill Proportionality (SEE/mean)',
                         pass_fail_range(dctadsummary['BagFillProportionality'].item(), [0, 3.5]), col_offset=7)

        skip_secs = 5
        if phdp_globals.test_data['CVSDLSSampleResults'] is not None:
            cvs_dlsresults = phdp_globals.test_data['CVSDLSSampleResults'].iloc[
                             int(skip_secs / constants['SamplePeriod_s']):]
            set_average_min_max(report_df, cvs_dlsresults, 'PM Filter Face Velocity', 'FilterFaceVelocity_cm/s',
                                col_offset=2, scale=1/100)
            pass_fail = (cvs_dlsresults['FilterFaceVelocity_cm/s'].min() / 100 >= 0.35 and
                         cvs_dlsresults['FilterFaceVelocity_cm/s'].max() / 100 <= 1)
            set_value_at(report_df, 'PM Filter Face Velocity', pass_fail_range(pass_fail, [True, True]), col_offset=7)

            set_average_min_max(report_df, cvs_dlsresults, 'CVS Tunnel Residence Time', 'CVSResidenceTime_s',
                                col_offset=2)
            set_average_min_max(report_df, cvs_dlsresults, 'PM Probe Residence Time', 'ProbeResidenceTime_s',
                                col_offset=2)

            set_average_min_max(report_df, cvs_dlsresults, 'PM Tunnel Residence Time', 'ResidenceTime_s',
                                col_offset=2)
            pass_fail = cvs_dlsresults['ResidenceTime_s'].min() >= 0.5
            set_value_at(report_df, 'PM Tunnel Residence Time', pass_fail_range(pass_fail, [True, True]), col_offset=7)

            set_average_min_max(report_df, cvs_dlsresults, 'Overall Residence Time', 'TotalResidenceTime_s',
                                col_offset=2)
            pass_fail = (cvs_dlsresults['TotalResidenceTime_s'].min() >= 1 and
                         cvs_dlsresults['TotalResidenceTime_s'].max() <= 5.5)
            set_value_at(report_df, 'Overall Residence Time', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        dctad_trunc = dctad.iloc[int(skip_secs / constants['SamplePeriod_s']):].copy()

        if set_average_min_max(report_df, dctad_trunc, 'CVS Tunnel Dilution Ratio', 'CVS Dilution Ratio', col_offset=2):
            pass_fail = dctad_trunc['CVS Dilution Ratio'].min() >= 2.0
            set_value_at(report_df, 'CVS Tunnel Dilution Ratio', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        set_average_min_max(report_df, dctad_trunc, 'PSU Dilution Ratio', 'PSU Dilution Ratio', col_offset=2)

        if set_average_min_max(report_df, dctad_trunc, 'Overall Dilution Ratio', 'Overall Dilution Ratio', col_offset=2):
            pass_fail = (5.0 <= dctad_trunc['Overall Dilution Ratio'].min() <= 7.0)
            set_value_at(report_df, 'Overall Dilution Ratio', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        if 'Proportionality Check' in dctad_trunc:
            set_value_at(report_df, 'Proportionality Check (SEE/mean flow)', [dctad_trunc['Proportionality Check'].mean()],
                         col_offset=2)
            pass_fail = dctad_trunc['Proportionality Check'].mean() <= 3.5
            set_value_at(report_df, 'Proportionality Check (SEE/mean flow)', pass_fail_range(pass_fail, [True, True]),
                         col_offset=7)

        if 'FilterPressureDrop_kPa' in dctad:
            set_value_at(report_df, 'Filter Pressure Drop Increase', [dctad['FilterPressureDrop_kPa'].iloc[-1] -
                                                                      dctad['FilterPressureDrop_kPa'].iloc[0]],
                         col_offset=2)

        with pd.ExcelWriter(
                report_filename,
                mode="a",
                engine="openpyxl",
                if_sheet_exists="replace",
        ) as writer:
            report_df.to_excel(writer, index=False, header=False,
                               sheet_name='%s General %d' % (calc_mode, emissions_cycle_number))


def generate_pre_test_check_report(report_filename):
    """
        Generate pre test check report

        Args:
            report_filename (str): the Excel report filename, a new sheet will be added

    """
    report_df = pd.read_csv(path + os.sep + 'pre_test_check_template.csv', encoding='UTF-8', header=None)
    report_df = report_df.fillna('')

    pre_test = phdp_globals.test_data['PreTest']

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
    CFR1065EMS_checks(report_df, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DIRECT') & (cfr1065ems['Component'] == 'THC') & (cfr1065ems['ColdHot'] == 'Hot')
    row_name = 'DIRECT THC Hot'
    CFR1065EMS_checks(report_df, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'HOT') & (cfr1065ems['Component'] == 'NH3')
    row_name = 'HOT NH3'
    CFR1065EMS_checks(report_df, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DILUTE') & (cfr1065ems['Component'] == 'CO2')
    row_name = 'DILUTE CO2 Cold'
    CFR1065EMS_checks(report_df, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DILUTE') & (cfr1065ems['Component'] == 'THC') & (cfr1065ems['ColdHot'] == 'Hot')
    row_name = 'DILUTE THC Hot'
    CFR1065EMS_checks(report_df, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DILUTE') & (cfr1065ems['Component'] == 'N2O')
    row_name = 'DILUTE N2O Cold'
    CFR1065EMS_checks(report_df, row_name, row_select, check_name)

    # CFR1065EMS HANGUP checks
    check_name = 'Hangup_ppmC'
    row_select = (cfr1065ems['Line'] == 'DILUTE') & (cfr1065ems['Component'] == 'THC') & (cfr1065ems['Check'] == 'HCHangup')
    row_name = 'DILUTE THC'
    CFR1065EMS_checks(report_df, row_name, row_select, check_name)

    row_select = (cfr1065ems['Line'] == 'DIRECT') & (cfr1065ems['Component'] == 'THC') & (cfr1065ems['Check'] == 'HCHangup')
    row_name = 'DIRECT THC'
    CFR1065EMS_checks(report_df, row_name, row_select, check_name)

    pass_fail_dict = {'Pass': 'pass', 'Fail': 'FAIL'}

    cfr1065cvs = phdp_globals.test_data['CFR1065CVS']

    if cfr1065cvs is not None:
        check_name = 'CheckPassFail'
        row_select = cfr1065cvs['Component'] == 'BSU Fill Line'
        row_name = 'BSU Fill Line Leak Check'
        set_value_at(report_df, row_name, pass_fail_dict[cfr1065cvs[check_name].loc[row_select].item()], col_offset=3)
        CFR1065_datetimes(report_df, cfr1065cvs, row_name, row_select)

        row_select = cfr1065cvs['Component'] == 'BSU Read Line'
        row_name = 'BSU Read Line Leak Check'
        set_value_at(report_df, row_name, pass_fail_dict[cfr1065cvs[check_name].loc[row_select].item()], col_offset=3)
        CFR1065_datetimes(report_df, cfr1065cvs, row_name, row_select)

        for bagpair in [1, 2, 3]:
            check_name = 'AmbBagHangup_ppmC'
            row_select = (cfr1065cvs['Check'] == 'HCHangup') & (cfr1065cvs['Component'] == 'THC') & (cfr1065cvs['HCHangupBagPair_Integer'] == bagpair)
            row_name = 'Amb Bag HC Hangup Test %d' % bagpair
            value = cfr1065cvs[check_name].loc[row_select].item()
            set_value_at(report_df, row_name, [value], col_offset=2)
            set_value_at(report_df, row_name, pass_fail_range(value, [-np.inf, 2]), col_offset=3)
            CFR1065_datetimes(report_df, cfr1065cvs, row_name, row_select)

            check_name = 'SampleBagHangup_ppmC'
            row_select = (cfr1065cvs['Check'] == 'HCHangup') & (cfr1065cvs['Component'] == 'THC') & (cfr1065cvs['HCHangupBagPair_Integer'] == bagpair)
            row_name = 'Sample Bag HC Hangup Test %d' % bagpair
            value = cfr1065cvs[check_name].loc[row_select].item()
            set_value_at(report_df, row_name, [value], col_offset=2)
            set_value_at(report_df, row_name, pass_fail_range(value, [-np.inf, 2]), col_offset=3)
            CFR1065_datetimes(report_df, cfr1065cvs, row_name, row_select)

    cfr1065pm = phdp_globals.test_data['CFR1065PM']

    if cfr1065pm is not None:
        check_name = 'CheckPassFail'
        row_name = 'DLS/PSU Leak Check'
        set_value_at(report_df, row_name, pass_fail_dict[cfr1065pm[check_name].item()], col_offset=3)
        CFR1065_datetimes(report_df, cfr1065pm, row_name, 0)

    with pd.ExcelWriter(
            report_filename,
            mode="a",
            engine="openpyxl",
            if_sheet_exists="replace",
    ) as writer:
        report_df.to_excel(writer, index=False, header=False,
                           sheet_name='PreTestCheck')


def CFR1065_datetimes(report_df, cfrdata, row_name, row_select):
    """
    Handle date/time and elapsed time values in the pre test check report

    Args:
        report_df (DataFrame): the report template dataframe
        cfrdata (DataFrame): CFR1065XXX data
        row_name (str): name of the row from the report template
        row_select (Series): a Boolean vector that selects a row from ``cfrdata``

    Returns:
        Nothing, updates values in ``report_df``

    """
    from datetime import datetime, timedelta

    time_offset = datetime(1899, 12, 29, 22, 26, 57).timestamp()

    test_time = cfrdata['Time_Date'].loc[row_select].item() * 24 * 3600 + time_offset
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


def CFR1065EMS_checks(report_df, row_name, row_select, value_name):
    """
    Update pre test check template CFR1065EMS-related values and perform pass/fail checks

    Args:
        report_df (DataFrame): the report template dataframe
        row_name (str): name of the row from the report template
        row_select (Series): a Boolean vector that selects a row from ``cfrdata``
        value_name (str): the column name of the target value , e.g. 'Hangup_ppmC'

    Returns:
        Nothing, updates values in ``report_df``

    """
    from datetime import datetime, timedelta

    cfr1065ems = phdp_globals.test_data['CFR1065EMS']

    if any(row_select):
        value = cfr1065ems[value_name].loc[row_select].item()
        set_value_at(report_df, row_name, [value], col_offset=2)
        set_value_at(report_df, row_name, pass_fail_range(value, [-0.5, 0.5]), col_offset=3)

        CFR1065_datetimes(report_df, cfr1065ems, row_name, row_select)


def generate_cycle_validation_report(report_filename, validation_results):
    """
    Generate cycle validation report

    Args:
        report_df (DataFrame): the report template dataframe
        validation_results (dict): dict of results from best cycle validation (min error, min omits, min shift)

    """
    report_df = pd.read_csv(path + os.sep + 'cycle_validation_report_template.csv', encoding='UTF-8', header=None)
    report_df = report_df.fillna('')

    pass_fail_dict = {True: 'pass', False: 'FAIL'}

    for i in range(len(validation_results['regression_results'])):
        regression_results = validation_results['regression_results'][i]
        ecn = regression_results['Emissions Cycle Number']

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
