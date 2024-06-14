import os
path = os.path.dirname(os.path.abspath(__file__))

from phdp import *

import pandas as pd


def generate_transient_report(report_filename, calc_mode, results, validation_results,
                              test_datetime, test_type, test_num, test_site, vehicle_test):
    """
    Generate a transient test report.

    Args:
        report_filename (str): the Excel report filename, a new sheet will be added
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        results (dict): a dictionary containing the results of the emissions calculations.
        validation_results (dict): dict of results from best cycle validation (min error, min omits, min shift)
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format.
        test_type (str): test type, i.e. 'transient' or 'modal'
        test_num (int): the number assigned to the test, e.g. '00139'
        test_site (str): the name of the site where the test was performed, e.g. 'HD02'
        vehicle_test (bool): ``True`` if test has an associated vehicle speed trace

    Returns:
        Report file name for use by subsequent reports

    """
    for i in range(0, len(results['tadsummary'])):
        emissions_cycle_number = results['tadsummary'][i]['EmissionsCycleNumber_Integer'].iloc[0]

        # create report header
        report_df = pd.read_csv(path + os.sep + 'report_templates' + os.sep + 'transient_report_template.csv',
                                encoding='UTF-8', header=None)

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

        if validation_results['PM_results']:
            set_value_at(report_df, 'Total PM Mass',
                         [validation_results['PM_results'][emissions_cycle_number]['pm_mass_g'] * 1000])
            set_value_at(report_df, 'BSPM',
                         [validation_results['PM_results'][emissions_cycle_number]['pm_mass_g'] /
                                             cycle_work_kWh.item()])
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
                       (original_brake_specific_emissions[1] - corrected_brake_specific_emissions[1])/3.5 * 100)
        pass_fail = pass_fail_range(co_value, [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=4)  # CO

        nox_value = min(brake_specific_emissions_drift_check_pct[2],
                       (original_brake_specific_emissions[2] - corrected_brake_specific_emissions[2])/0.27 * 100)
        pass_fail = pass_fail_range(nox_value, [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=5)  # NOx

        hc_value = min(brake_specific_emissions_drift_check_pct[3],
                       (original_brake_specific_emissions[3] - corrected_brake_specific_emissions[3])/0.19 * 100)
        pass_fail = pass_fail_range(hc_value, [-4, 4])
        set_value_at(report_df, '1065 Drift Check', pass_fail, col_offset=6)  # HC

        nmhc_value = min(brake_specific_emissions_drift_check_pct[6],
                       (original_brake_specific_emissions[6] - corrected_brake_specific_emissions[6])/0.19 * 100)
        pass_fail = pass_fail_range(nmhc_value, [-4, 4])
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

        with pd.ExcelWriter(
                report_filename,
                mode='a',
                engine="openpyxl",
        ) as writer:
            report_df.to_excel(writer, index=False, header=False, sheet_name='Emissions %s %d' %
                                                                             (calc_mode, emissions_cycle_number))
