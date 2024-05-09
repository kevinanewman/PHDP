"""
PHDP report generators

"""

import os
path = os.path.dirname(os.path.abspath(__file__))
#
# import pandas as pd


from phdp import *


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

    """
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

        if calc_mode == 'dilute':
            net_mass_drift_check_pct = ((original_net_mass - corrected_net_mass) /
                                        np.maximum(original_net_mass, sys.float_info.epsilon) * 100)
            set_value_at(report_df, 'Net Mass Drift Check %', net_mass_drift_check_pct)

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

        with pd.ExcelWriter(
                report_filename,
                mode="w",
                engine="openpyxl",
        ) as writer:
            report_df.to_excel(writer, index=False, header=False, sheet_name='Emissions %d' % emissions_cycle_number)

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


def calc_emscalresults_drift_check(report_df, check_phase, limit_pct, driftline, components):
    source = phdp_globals.test_data['EmsCalResults']

    for component in components:
        for measurement_type, offset in zip(
                ['DriftSpanValue_ppm', 'DriftSpanRange_ppm', 'DriftZeroMeasured_ppm', 'DriftSpanMeasured_ppm',
                 'DriftZero2Measured_ppm'],
                [5, 1, 2, 6, 9]):
            if any((source['DriftLine'] == driftline) & (source['DriftComponent'] == component)):
                value = source.loc[(source['DriftLine'] == driftline) & (source['DriftComponent'] == component),
                        measurement_type].item()

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


def calc_driftcheck_drift_check(report_df, check_phase, limit_pct, driftline, components):
    source = phdp_globals.test_data['DriftCheck']

    for component in components:
        for measurement_type, offset in zip(
                ['DriftSpanValue_ppm', 'DriftSpanRange_ppm', 'DriftZeroMeasured_ppm', 'DriftSpanMeasured_ppm'],
                [5, 1, 2, 6]):
            if any((source['DriftLine'] == driftline) & (source['DriftComponent'] == component)):
                value = source.loc[(source['DriftLine'] == driftline) & (source['DriftComponent'] == component),
                        measurement_type].item()

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


def calc_bagdata_drift_check(report_df, check_phase, limit_pct, driftline, components):
    source = phdp_globals.test_data['BagData']

    for component in components:
        for measurement_type, offset in zip(
                ['RbSpanValue_ppm', 'RbRange_ppm', 'RbZeroCalMeasured_ppm', 'RbSpanCalMeasured_ppm',
                 'RbZero2CalMeasured_ppm'],
                [5, 1, 2, 6, 9]):
            if any((source['RbLine'] == driftline) & (source['RbComponent'] == component)):
                value = source.loc[(source['RbLine'] == driftline) & (source['RbComponent'] == component),
                        measurement_type].item()

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


def calc_bagdriftcheck_drift_check(report_df, check_phase, limit_pct, driftline, components):
    source = phdp_globals.test_data['BagDriftCheck']

    for component in components:
        for measurement_type, offset in zip(
                ['DriftSpanValue_ppm', 'DriftSpanRange_ppm', 'DriftZeroMeasured_ppm', 'DriftSpanMeasured_ppm'],
                [5, 1, 2, 6]):
            if any((source['DriftLine'] == driftline) & (source['DriftComponent'] == component)):
                value = source.loc[(source['DriftLine'] == driftline) & (source['DriftComponent'] == component),
                        measurement_type].item()

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


def generate_driftcheck_report(report_filename, results, test_datetime):
    """
    Generate a transient test report.

    Args:
        report_filename (str): the Excel report filename, a new sheet will be added
        results (dict): a dictionary containing the results of the emissions calculations.
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format.

    """
    report_df = pd.read_csv(path + os.sep + 'drift_check_report_template.csv', encoding='UTF-8', header=None)
    report_df = report_df.fillna('')

    for check_phase, limit_pct in zip(['PRE', 'POST'], [1, 2]):
        if check_phase == 'PRE':
            calc_method = calc_emscalresults_drift_check
        else:
            calc_method = calc_driftcheck_drift_check

        driftline = 'DIRECT'
        components = ['COH', 'CO2', 'CH4', 'THC', 'NOX', 'O2']
        calc_method(report_df, check_phase, limit_pct, driftline, components)

        driftline = 'HOT'
        components = ['NH3']
        calc_method(report_df, check_phase, limit_pct, driftline, components)

        driftline = 'DILUTE'
        components = ['COL', 'CO2', 'CH4', 'THC', 'NOX', 'N2O']
        calc_method(report_df, check_phase, limit_pct, driftline, components)

        driftline = 'EGR'
        components = ['CO2']
        calc_method(report_df, check_phase, limit_pct, driftline, components)

    for check_phase, limit_pct in zip(['PRE', 'POST'], [1, 2]):
        if check_phase == 'PRE':
            calc_method = calc_bagdata_drift_check
        else:
            calc_method = calc_bagdriftcheck_drift_check

        driftline = 'BAG'
        components = ['COL', 'CO2', 'CH4', 'THC', 'NOX', 'N2O']
        calc_method(report_df, check_phase, limit_pct, driftline, components)

        with pd.ExcelWriter(
                report_filename,
                mode="a",
                engine="openpyxl",
                if_sheet_exists="replace",
        ) as writer:
            report_df.to_excel(writer, index=False, header=False, sheet_name='Drift Check')
