from phdp import *

import math
from scipy.stats import linregress


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

    if test_name in ['GHGTRNS', 'GHGCRSE']:
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
            elif test_name in ['GHGTRNS', 'GHGCRSE']:
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
                    if test_name in ['GHGTRNS', 'GHGCRSE']:
                        start_condition = (continuous_data['EmissionsCycleNumber_Integer'] == ecn)
                    else:
                        start_condition = ((continuous_data['ModeNumber_Integer'] == 1) &
                                           (continuous_data['EmissionsCycleNumber_Integer'] == ecn))

                    if test_name == 'RMC':
                        start_condition_index = 9
                    elif test_name in ['GHGTRNS', 'GHGCRSE']:
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
                        .index[end_condition_index])) - 1

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

                    if stp == 'speed_rpm':
                        reference_speed_range = (ref.max() - ref.min()) / ref.mean()

                    if stats['slope'] is None:
                        pass_fail[stp] = 'FAIL'
                        linear_regression_fail = True
                    else:
                        if stp == 'speed_rpm' and reference_speed_range <= 0.1:
                            # only need to check SEE
                            pass_fail[stp] = pass_fail_range(stats['SEE'], limits[stp]['SEE']) == 'pass'
                            fail_count += sum([int(pass_fail_range(stats['SEE'], limits[stp]['SEE']) == 'FAIL')])
                        else:
                            pass_fail[stp] = all([pass_fail_range(stats[k], limits[stp][k]) == 'pass' for k in stats])
                            fail_count += sum([int(pass_fail_range(stats[k], limits[stp][k]) == 'FAIL') for k in stats])

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
