import os
path = os.path.dirname(os.path.abspath(__file__))

from phdp import *
from common.file_io import file_exists

import pandas as pd


def generate_modal_report(report_filename, calc_mode, results, validation_results,
                          test_datetime, test_type, test_num, test_site):
    """
    Generate a modal report and calculated mode-weighted mass and brake-specific emissions.

    Args:
        report_filename (str): the Excel report filename, a new sheet will be added
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        results (dict): a dictionary containing the results of the emissions calculations.
        validation_results (dict): dict of results from best cycle validation (min error, min omits, min shift)
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format.
        test_type (str): test type, i.e. 'transient' or 'modal'
        test_num (int): the number assigned to the test, e.g. '00139'
        test_site (str): the name of the site where the test was performed, e.g. 'HD02'

    Returns:
        Report file name for use by subsequent reports

    """
    # create report header
    report_df = pd.read_csv(path + os.sep + 'report_templates' + os.sep + 'modal_report_template.csv',
                            encoding='UTF-8', header=None)

    report_df = report_df.fillna('')

    cycle_id = phdp_globals.test_data['TestDetails']['CycleName'].iloc[0]

    set_value_at(report_df, 'Test Date', '%s/%s/%s' % (test_datetime[4:6], test_datetime[6:8], test_datetime[0:4]))
    set_value_at(report_df, 'Test Cell', test_site)
    set_value_at(report_df, 'Test Number', test_num)
    set_value_at(report_df, 'Test Type', test_type)
    set_value_at(report_df, 'Cycle ID', '%s' % cycle_id)
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
            if 'mPM_g' in results['dctadsummary'][i]:
                set_value_at(report_df, 'Total PM', results['dctadsummary'][i]['mPM_g'], col_offset=col_offset)
            else:
                set_value_at(report_df, 'Total PM', [0], col_offset=col_offset)

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

        # Mode time validation
        stablization_time_s = results['dctad'][i]['stabilization_time_s']
        set_value_at(report_df, 'Stabilization Time', stablization_time_s, col_offset=col_offset)
        if cycle_id == 'Steady State Fuel Use':
            pass_fail = pass_fail_range(stablization_time_s.item(), [69, 71])
        elif cycle_id == 'Idle Fuel Consumption':
            pass_fail = pass_fail_range(stablization_time_s.item(), [179, 181])
        else:
            pass_fail = pass_fail_range(stablization_time_s.item(), [300, np.inf])
        set_value_at(report_df, 'Stabilization Time Check', pass_fail, col_offset=col_offset)

        mode_time_s = results['dctad'][i]['elapsed_time_s']
        set_value_at(report_df, 'Sample Time', mode_time_s, col_offset=col_offset)
        if cycle_id == 'Steady State Fuel Use':
            pass_fail = pass_fail_range(mode_time_s.item(), [29, 31])
        elif cycle_id == 'Idle Fuel Consumption':
            pass_fail = pass_fail_range(mode_time_s.item(), [599, 601])
        else:
            pass_fail = pass_fail_range(mode_time_s.item(), [60, np.inf])
        set_value_at(report_df, 'Sample Time Check', pass_fail, col_offset=col_offset)

        # Speed / torque validation
        speed_demand_rpm = results['dctad'][i]['DynoSpeedDemand_Avg_rev/min']
        speed_max_rpm = results['dctad'][i]['spDyno_Avg_rev/min_max']
        speed_min_rpm = results['dctad'][i]['spDyno_Avg_rev/min_min']
        set_value_at(report_df, 'Speed Demand', speed_demand_rpm, col_offset=col_offset)
        set_value_at(report_df, 'Average Speed', results['dctad'][i]['spDyno_Avg_rev/min_mean'], col_offset=col_offset)
        set_value_at(report_df, 'Maximum Speed', speed_max_rpm, col_offset=col_offset)
        set_value_at(report_df, 'Minimum Speed', speed_min_rpm, col_offset=col_offset)

        if cycle_id == 'Steady State Fuel Use' or (cycle_id == 'Idle Fuel Consumption' and mode_number > 1):
            speed_tolerance_rpm = 0.01 * phdp_globals.test_data['MapResults']['EngSpeedNhi50_rev/min'].item()
            speed_limit = '%.3f ± %.2f' % (speed_demand_rpm.item(), speed_tolerance_rpm)
        else:
            speed_limit = ''

        if speed_limit:
            set_value_at(report_df, 'Speed Limit', speed_limit, col_offset=col_offset)

            if ((speed_max_rpm.item() <= speed_demand_rpm.item() + speed_tolerance_rpm) and
                    (speed_min_rpm.item() >= speed_demand_rpm.item() - speed_tolerance_rpm)):
                pass_fail = 'pass'
            else:
                pass_fail = 'FAIL'

            set_value_at(report_df, 'Speed Check', pass_fail, col_offset=col_offset)

        torque_demand_Nm = results['dctad'][i]['DynoTorqueDemand_Avg_Nm']
        torque_max_Nm = results['dctad'][i]['tqShaft_Avg_Nm_max']
        torque_min_Nm = results['dctad'][i]['tqShaft_Avg_Nm_min']
        set_value_at(report_df, 'Torque Demand', torque_demand_Nm, col_offset=col_offset)
        set_value_at(report_df, 'Average Torque', results['dctad'][i]['tqShaft_Avg_Nm_mean'], col_offset=col_offset)
        set_value_at(report_df, 'Maximum Torque', torque_max_Nm, col_offset=col_offset)
        set_value_at(report_df, 'Minimum Torque', torque_min_Nm, col_offset=col_offset)

        if cycle_id == 'Steady State Fuel Use' or (cycle_id == 'Idle Fuel Consumption' and mode_number > 1):
            torque_tolerance_Nm = 0.05 * phdp_globals.test_data['MapResults']['EngMaxTorque_Nm'].item()
            torque_limit = '%.3f ± %.2f' % (torque_demand_Nm.item(), torque_tolerance_Nm)
        elif cycle_id == 'Idle Fuel Consumption' and mode_number == 1:
            torque_tolerance_Nm = 25
            torque_limit = '%.3f ± %.2f' % (torque_demand_Nm.item(), torque_tolerance_Nm)
        else:
            torque_limit = ''

        if torque_limit:
            set_value_at(report_df, 'Torque Limit', torque_limit, col_offset=col_offset)

            if ((torque_max_Nm.item() <= torque_demand_Nm.item() + torque_tolerance_Nm) and
                    (torque_min_Nm.item() >= torque_demand_Nm.item() - torque_tolerance_Nm)):
                pass_fail = 'pass'
            else:
                pass_fail = 'FAIL'

            set_value_at(report_df, 'Torque Check', pass_fail, col_offset=col_offset)

    if 'WeightingFactor_Fraction' in phdp_globals.test_data['CycleDefinition']:
        weighted_power = 0
        for i in range(0, len(results['dctadsummary'])):
            weight = phdp_globals.test_data['CycleDefinition']['WeightingFactor_Fraction'].iloc[i]
            set_value_at(report_df, 'Weighting Factor', weight, col_offset=i+2)
            weighted_power += results['dctad'][i]['Power_kW'] * weight

        set_value_at(report_df, 'Power (kW)', weighted_power)

        # calculate weighted values
        for idx, signal in enumerate(('CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O', 'NMHC', 'PM')):
            weighted_mass_emissions = pd.Series(0)
            for i in range(0, len(results['dctadsummary'])):
                if 'm%s_g' % signal in results['dctadsummary'][i]:
                    weighted_mass_emissions += (
                            results['dctadsummary'][i]['m%s_g' % signal] *
                            phdp_globals.test_data['CycleDefinition']['WeightingFactor_Fraction'].iloc[i] /
                            results['dctad'][i]['elapsed_time_s'] * 3600)
            weighted_specific_emissions = weighted_mass_emissions / weighted_power
            set_value_at(report_df, 'Mass Emissions (g/h)', weighted_mass_emissions, col_offset=idx + 1)
            set_value_at(report_df, 'Specific Emsissions (g/kWh)', weighted_specific_emissions, col_offset=idx + 1)

    if file_exists(report_filename):
        mode = 'a'
    else:
        mode = 'w'

    with pd.ExcelWriter(
            report_filename,
            mode=mode,
            engine="openpyxl",
    ) as writer:
        report_df.to_excel(writer, index=False, header=False, sheet_name='Emissions %s' % calc_mode)
