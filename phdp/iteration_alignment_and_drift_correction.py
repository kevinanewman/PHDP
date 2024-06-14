from phdp import *
import test_sites
from constants import constants


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
