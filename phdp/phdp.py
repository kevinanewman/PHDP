"""

PHDP top level code

----

**CODE**

"""

import sys, os

import pandas as pd

path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(path, '..'))  # picks up omega_model sub-packages

from phdp import *
from constants import constants
import test_sites

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
        columns = pd.read_csv(filename, header=None, nrows=1,  encoding=encoding, encoding_errors='strict')
        if units_nrows > 0:
            units = pd.read_csv(filename, header=None, skiprows=1, nrows=units_nrows, encoding=encoding, encoding_errors='strict')
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


def load_data(test_site):
    """
        Load test data into phdp_globals.test_data dict

    Args:
        test_site (str): test site ID, e.g. 'HD02'

    """
    # # read the Horiba file (not used for now)
    # phdp_log.logwrite('reading %s...' % phdp_globals.options.horiba_file)
    # unitized_columns = get_unitized_columns(phdp_globals.options.horiba_file, sheet_name='ContinuousData1')
    # phdp_globals.test_data['horiba_data'] = \
    #     pd.read_excel(phdp_globals.options.horiba_file, names=unitized_columns, header=1, skiprows=0,
    #                   sheet_name='ContinuousData1')

    # read in the rest of the files just in case
    os.chdir(file_io.get_filepath(phdp_globals.options.horiba_file))
    input_files = sorted([f for f in os.listdir() if f.endswith('.csv')])
    for input_file in input_files:
        file_name = input_file.rsplit('.', 2)[-2]
        if file_name != 'Processing':
            phdp_log.logwrite('reading %s...' % input_file)
            if 'tad' not in file_name:
                unitized_columns = get_unitized_columns(input_file, encoding=phdp_globals.options.encoding[test_site])
                phdp_globals.test_data[file_name] = \
                    pd.read_csv(input_file, names=unitized_columns,
                                encoding=phdp_globals.options.encoding[test_site], encoding_errors='strict',
                                header=1, skiprows=0)
            else:
                get_unitized_columns(input_file, units_nrows=0, encoding=phdp_globals.options.output_encoding)  # dump tad columns to logfile for reference
                phdp_globals.test_data[file_name] = pd.read_csv(input_file, encoding=phdp_globals.options.output_encoding)


def time_align_continuous_data(test_site, vehicle_test, sampled_crank, emissions_cycle_number, min_mode_number):
    """
    Time-align continuous data, add Vehicle Moving from cycle definition

    Args:
        test_site (str): e.g. 'HD02'
        vehicle_test (bool): ``True`` if test has an associated vehicle speed trace
        sampled_crank (bool): ``True`` if test has sampled engine cranking and startup
        emissions_cycle_number (int): emissions cycle number to process
        min_mode_number (int): the starting mode number of the data

    Returns:
        Dataframe of time-aligned data

    """
    from test_sites import site_info

    test_sites.init_site_info(test_site)

    SamplePeriod_s = phdp_globals.test_data['TestParameters']['ContinuousLoggerPeriod_s'].item()

    constants['SamplePeriod_s'] = SamplePeriod_s

    time_aligned_data = pd.DataFrame(index=phdp_globals.test_data['ContinuousData'].index)

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

    time_aligned_data = time_aligned_data.iloc[test_start_index:test_end_index]

    time_aligned_data = time_aligned_data.fillna(0)  # TODO: for now, need to confirm...

    return time_aligned_data


def pre_chemical_balance_calculations(time_aligned_data, test_type):
    """
    Calculate values required for chemical balance iteration

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        test_type (str): 'transient' or 'modal'

    Returns:
        Nothing, updates time_aligned_data

    """
    time_aligned_data['BagFillFlow_Avg_m³/s'] = time_aligned_data['BagFillFlow_Avg_l/min'] / 60000

    if test_type == 'transient':
        if 'CVSMolarFlow_Avg_mol/s' not in time_aligned_data:
            time_aligned_data['CVSFlow_mol/s'] = \
                (time_aligned_data['CVSFlow_Avg_m³/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] +
                phdp_globals.test_data['TestParameters']['DiluteSampleVolumeFlow_m³/s'].item()) / 0.024055
        else:
            # TODO: need to verify if this is correct for HD05 transient (FTP) test
            time_aligned_data['CVSFlow_mol/s'] = (
                    time_aligned_data['CVSMolarFlow_Avg_mol/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] / 0.024055 +
                    phdp_globals.test_data['TestParameters']['DiluteSampleMolarFlow_mol/s'].item())
    else:
        if 'CVSMolarFlow_Avg_mol/s' not in time_aligned_data:
            time_aligned_data['qvCVS_Avg_m³/s'] = time_aligned_data['qvCVS_Avg_m³/min'] / 60

            time_aligned_data['CVSFlow_mol/s'] = (
                (time_aligned_data['qvCVS_Avg_m³/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] +
                 phdp_globals.test_data['TestParameters']['DiluteSampleVolumeFlow_l/s'] / 1000) / 0.024055)
        else:
            time_aligned_data['CVSFlow_mol/s'] = (
                    time_aligned_data['CVSMolarFlow_Avg_mol/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] / 0.024055 +
                    phdp_globals.test_data['TestParameters']['DiluteSampleMolarFlow_mol/s'])

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
                time_aligned_data['pH2Odil_kPa'] / time_aligned_data['pCellAmbient_kPa']
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
            time_aligned_data['xH2Odil_mol/mol'] = time_aligned_data['pH2Odilsat_kPa'] / time_aligned_data['pCellAmbient_kPa']
        else:
            # 1065.645-3:
            time_aligned_data['xH2Odil_mol/mol'] = time_aligned_data['pH2Odilsat_kPa'] / time_aligned_data['pCellAmbient_Avg_kPa']

    from constants import constants, update_constants
    update_constants()  # update constants that rely on test fuel properties, etc
    time_aligned_data['Power_kW'] = \
        np.maximum(0, time_aligned_data['tqShaft_Nm'] * time_aligned_data['spDyno_rev/min'] / 9548.8)

    if 'DEFMassFlowRate_Avg_g/h' not in time_aligned_data:
        time_aligned_data['DEFMassFlowRate_Avg_g/h'] = 0

    # 1065.655-20:
    time_aligned_data['alpha'] = \
        CFR1065.alpha(time_aligned_data['qmFuel_Avg_g/h'], time_aligned_data['DEFMassFlowRate_Avg_g/h'])

    # 1065.655-21:
    time_aligned_data['beta'] = \
        CFR1065.beta(time_aligned_data['qmFuel_Avg_g/h'], time_aligned_data['DEFMassFlowRate_Avg_g/h'])

    # 1065.655-23:
    time_aligned_data['delta'] = \
        CFR1065.delta(time_aligned_data['qmFuel_Avg_g/h'], time_aligned_data['DEFMassFlowRate_Avg_g/h'])

    time_aligned_data['gamma'] = 0  # no gamma for now

    time_aligned_data['Tint_K'] = time_aligned_data['tIntakeAir_°C'] + 273.15
    time_aligned_data['Tdewint_°C'] = time_aligned_data['tCellDewPt_°C']

    # 1065.645-1:
    time_aligned_data['Tdewint_K'] = time_aligned_data['Tdewint_°C'] + 273.15

    # 1065.645-1:
    time_aligned_data['pH2Oamb_kPa'] = \
        CFR1065.vapor_pressure_of_water_kPa(time_aligned_data['Tdewint_K'])

    # 1065.645-3:
    time_aligned_data['xH2Oint_mol/mol'] = \
        time_aligned_data['pH2Oamb_kPa'] / time_aligned_data['pCellAmbient_kPa']


def iterate_chemical_balance(time_aligned_data, emissions_cycle_number, drift_corrected=False):
    """
    Iterate the chemical balance equations

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        emissions_cycle_number (int): emissions cycle number to process
        drift_corrected (bool): use drift-correct background if ``True``

    Returns:
        Nothing, updates time_aligned_data

    """
    time_aligned_data['xDil/Exh_mol/mol'] = 0.8
    time_aligned_data['xH2Oexh_mol/mol'] = 2 * time_aligned_data['xH2Oint_mol/mol']
    time_aligned_data['xCcombdry_mol/mol'] = \
        time_aligned_data['conCO2_Avg_%vol'] / 100 + (time_aligned_data['conTHC_Avg_ppmC'] +
                                                      time_aligned_data['conLCO_Avg_ppm']) / 1e6
    time_aligned_data['xH2dry_μmol/mol'] = 0
    time_aligned_data['xint/exhdry_mol/mol'] = 0
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
            phdp_globals.test_data['EmsComponents']['InputName'] == 'conLCO', 'ResidualH2O_%vol'].item()

        time_aligned_data['xCOdry_μmol/mol'] = \
            time_aligned_data['conLCO_Avg_ppm'] / (1 - residual_H2O_pctvol / 100)

        # 1065.655-15
        residual_H2O_pctvol = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['InputName'] == 'conCO2', 'ResidualH2O_%vol'].item()

        time_aligned_data['xCO2dry_%'] = \
            time_aligned_data['conCO2_Avg_%vol'] / (1 - residual_H2O_pctvol / 100)

        # 1065.655-16
        residual_H2O_pctvol = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['InputName'] == 'conNOX', 'ResidualH2O_%vol'].item()
        time_aligned_data['xNOdry_μmol/mol'] = \
            time_aligned_data['conNOX_Avg_ppm'] * 0.75 / (1 - residual_H2O_pctvol / 100)

        # 1065.655-17
        time_aligned_data['xNO2dry_μmol/mol'] = \
            time_aligned_data['conNOX_Avg_ppm'] * 0.25 / (1 - residual_H2O_pctvol / 100)

        if not drift_corrected:
            BagData = phdp_globals.test_data['BagData']
        else:
            BagData = phdp_globals.test_data['drift_corrected_BagData']

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
            time_aligned_data['conTHC_Avg_ppmC'] / (1 - time_aligned_data['xH2Oexh_mol/mol'])

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

        if iteration == 0:
            converged = False
        else:
            converged = ((time_aligned_data - time_aligned_data_prior).abs() <=
                         phdp_globals.options.chemical_balance_convergence_tolerance *
                         time_aligned_data.abs()).all().all()

        time_aligned_data_prior = time_aligned_data.copy()

        if iteration > 10:
            print(emissions_cycle_number, iteration, time_aligned_data['xCcombdry_mol/mol'].iloc[0])

        iteration = iteration + 1

    print('Converged in %d iterations' % iteration)

    if iteration > 10:
        print()


def post_chemical_balance_calculations(time_aligned_data):
    """
    Calculate values after chemical balance iteration

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data

    Returns:
        Nothing, updates time_aligned_data

    """
    from constants import constants

    # CFR 1065.659-1
    xH2OCO2dilmeas = phdp_globals.test_data['EmsComponents'].loc[
        phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilCO2_System']['ResidualH2O_%vol'].item()
    time_aligned_data['xCO2exh_%mol'] = time_aligned_data['conCO2_Avg_%vol'] * \
                                        ((1 - time_aligned_data['xH2Oexh_mol/mol']) / (1 - xH2OCO2dilmeas / 100))
    # CFR 1065.659-1
    time_aligned_data['xTHCexh_μmol/mol'] = \
        time_aligned_data['conTHC_Avg_ppmC'] - \
        phdp_globals.test_data['TestParameters']['InitialDilTHC_ppmC'].item()

    # CFR 1065.660-9
    DiluteRFPFC2H6_Fraction = phdp_globals.test_data['TestParameters']['DiluteRFPFC2H6_Fraction'].item()
    DiluteRFCH4_Fraction = phdp_globals.test_data['TestParameters']['DiluteRFCH4_Fraction'].item()
    time_aligned_data['xCH4exh_μmol/mol'] = \
        (time_aligned_data['conCH4cutter_Avg_ppmC'] -
         time_aligned_data['conTHC_Avg_ppmC'] * DiluteRFPFC2H6_Fraction) / \
        (1 - DiluteRFPFC2H6_Fraction * DiluteRFCH4_Fraction)

    # CFR 1065.660-4
    time_aligned_data['xNMHCexh_μmol/mol'] = \
        (time_aligned_data['xTHCexh_μmol/mol'] -
         time_aligned_data['xCH4exh_μmol/mol'] * DiluteRFCH4_Fraction) / \
        (1 - DiluteRFPFC2H6_Fraction * DiluteRFCH4_Fraction)

    # CFR 1065.659-1
    xH2OCOdilmeas = phdp_globals.test_data['EmsComponents'].loc[
        phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilCO_System']['ResidualH2O_%vol'].item()
    time_aligned_data['xCOexh_μmol/mol'] = \
        time_aligned_data['conLCO_Avg_ppm'] * \
        ((1 - time_aligned_data['xH2Oexh_mol/mol']) / (1 - xH2OCOdilmeas / 100))

    # CFR 1065.659-1
    xH2ONOxdilmeas = phdp_globals.test_data['EmsComponents'].loc[
        phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilNOx_System']['ResidualH2O_%vol'].item()
    time_aligned_data['xNOexh_μmol/mol'] = \
        time_aligned_data['conNOX_Avg_ppm'] * 0.75 * \
        ((1 - time_aligned_data['xH2Oexh_mol/mol']) / (1 - xH2ONOxdilmeas / 100))

    # CFR 1065.659-1
    time_aligned_data['xNO2exh_μmol/mol'] = \
        time_aligned_data['conNOX_Avg_ppm'] * 0.25 * \
        ((1 - time_aligned_data['xH2Oexh_mol/mol']) / (1 - xH2ONOxdilmeas / 100))

    # CFR 1065.670-1
    time_aligned_data['xNOxcorrected_μmol/mol'] = \
        (time_aligned_data['xNOexh_μmol/mol'] + time_aligned_data['xNO2exh_μmol/mol']) * \
        (9.953 * time_aligned_data['xH2Oint_mol/mol'] + 0.832)

    # CFR 1065.650-5
    time_aligned_data['mCO2_g/sec'] = \
        time_aligned_data['xCO2exh_%mol'] / 100 * constants['MCO2_g/mol'] * time_aligned_data['CVSFlow_mol/s']

    # CFR 1065.650-5
    time_aligned_data['mCO_g/sec'] = \
        time_aligned_data['xCOexh_μmol/mol'] / 1e6 * constants['MCO_g/mol'] * time_aligned_data['CVSFlow_mol/s']

    # CFR 1065.650-5
    time_aligned_data['mCH4_g/sec'] = \
        time_aligned_data['xCH4exh_μmol/mol'] / 1e6 * constants['MCH4_g/mol'] * time_aligned_data['CVSFlow_mol/s']

    # CFR 1065.650-5
    time_aligned_data['mTHC_g/sec'] = \
        time_aligned_data['xTHCexh_μmol/mol'] / 1e6 * constants['MTHC_g/mol'] * time_aligned_data['CVSFlow_mol/s']

    # CFR 1065.650-5
    time_aligned_data['mNMHC_g/sec'] = \
        time_aligned_data['xNMHCexh_μmol/mol'] / 1e6 * constants['MNMHC_g/mol'] * time_aligned_data['CVSFlow_mol/s']

    # CFR 1065.650-5
    time_aligned_data['mNOxuncorrected_g/sec'] = \
        (time_aligned_data['xNOexh_μmol/mol'] + time_aligned_data['xNO2exh_μmol/mol']) / 1e6 * \
        constants['MNOx_g/mol'] * time_aligned_data['CVSFlow_mol/s']

    # CFR 1065.650-5
    time_aligned_data['mNOxcorrected_g/sec'] = \
        (time_aligned_data['xNOxcorrected_μmol/mol']) / 1e6 * \
        constants['MNOx_g/mol'] * time_aligned_data['CVSFlow_mol/s']

    # CFR 1065.650-5
    time_aligned_data['mN2O_g/sec'] = \
        time_aligned_data['conN2O_Avg_ppm'] / 1e6 * constants['MN2O_g/mol'] * time_aligned_data['CVSFlow_mol/s']

    # CFR 1065.640-9
    time_aligned_data['Mmix_intake_g/mol'] = \
        constants['Mair_g/mol'] * (1 - time_aligned_data['xH2Oint_mol/mol']) + \
        constants['MH2O_g/mol'] * time_aligned_data['xH2Oint_mol/mol']
    time_aligned_data['nint_mol/sec'] = \
        time_aligned_data['qmIntakeAir_Avg_kg/h'] / 3.6 / time_aligned_data['Mmix_intake_g/mol']

    # CFR 1065.655-26
    time_aligned_data['nexh_mol/sec'] = \
        (time_aligned_data['xraw/exhdry_mol/mol'] - time_aligned_data['xint/exhdry_mol/mol']) * \
        (1 - time_aligned_data['xH2Oexh_mol/mol']) * time_aligned_data['CVSFlow_mol/s'] + \
        time_aligned_data['nint_mol/sec']

    # CFR 1065.667(C)
    time_aligned_data['ndil_mol/sec'] = time_aligned_data['CVSFlow_mol/s'] - time_aligned_data['nexh_mol/sec']


def drift_correct_continuous_data(time_aligned_data, signal_name):
    """
    Perform drift correction for continuous data, per CFR1065.672-1

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        signal_name (str): signal name, e.g. 'conRawCO2_Avg_%vol'

    Returns:
        Nothing, updates dataframe with drift corrected signal.

    """
    rename = {
        'LCO': 'COL',
        'HCO': 'COH',
        'CH4cutter': 'CH4'
    }

    xrefzero = 0

    signal, signal_type, unit = signal_name.split('_')

    if signal.startswith('conRaw'):
        component = signal.replace('conRaw', '')
        driftline = 'DIRECT'
        if component == 'NH3':
            driftline = 'HOT'
    elif signal.startswith('conEGRCO2'):
        component = 'CO2'
        driftline = 'EGR'
    else:
        component = signal.replace('con', '')
        driftline = 'DILUTE'

    if component in rename:
        component = rename[component]

    EmsCalResults = phdp_globals.test_data['EmsCalResults']
    xpre_data = (
        EmsCalResults)[(EmsCalResults['DriftComponent'] == component) & (EmsCalResults['DriftLine'] == driftline)]

    if not xpre_data.empty:
        xrefspan = xpre_data['DriftSpanValue_ppm'].item()
        xprezero = xpre_data['DriftZero2Measured_ppm'].item()
        xprespan = xpre_data['DriftSpanMeasured_ppm'].item()

        DriftCheck = phdp_globals.test_data['DriftCheck']
        xpost_data = (
            DriftCheck)[(DriftCheck['DriftComponent'] == component) & (DriftCheck['DriftLine'] == driftline)]

        xpostzero = xpost_data['DriftZeroMeasured_ppm'].item()
        xpostspan = xpost_data['DriftSpanMeasured_ppm'].item()

        if unit == '%vol':
            scale_factor = 10 ** 4
        else:
            scale_factor = 1

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
                            (DriftCheck['EmissionsCycleNumber_Integer'] == bag_data.loc[idx, 'EmissionsCycleNumber_Integer'])]

    xpostzero = xpost_data['DriftZeroMeasured_ppm'].item()
    xpostspan = xpost_data['DriftSpanMeasured_ppm'].item()

    bag_data.loc[idx, 'RbSmpConc_ppm'] = xrefzero + (xrefspan - xrefzero) * (
            2 * bag_data.loc[idx, 'RbSmpConc_ppm'] - (xprezero + xpostzero)) / (
            (xprespan + xpostspan) - (xprezero + xpostzero))

    bag_data.loc[idx, 'RbAmbConc_ppm'] = xrefzero + (xrefspan - xrefzero) * (
            2 * bag_data.loc[idx, 'RbAmbConc_ppm'] - (xprezero + xpostzero)) / (
            (xprespan + xpostspan) - (xprezero + xpostzero))


def calc_summary_results(time_aligned_data, emissions_cycle_number, test_type, drift_corrected=False):
    """
    Calculate summary results for the given time aligned data

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data
        emissions_cycle_number (int): emissions cycle number to process
        test_type (str): 'transient' or 'modal'
        drift_corrected (bool): use drift-correct background if ``True``

    Returns:
        Summary results in a pandas Series

    """
    from constants import constants

    # calculate summary values
    summary_results = pd.DataFrame(index=[0])
    summary_results['avg_xCO2exh_%mol'] = time_aligned_data['xCO2exh_%mol'].mean()
    summary_results['avg_xCOexh_μmol/mol'] = time_aligned_data['xCOexh_μmol/mol'].mean()
    summary_results['avg_xNOxcorrected_μmol/mol'] = time_aligned_data['xNOxcorrected_μmol/mol'].mean()
    summary_results['avg_xTHCexh_μmol/mol'] = time_aligned_data['xTHCexh_μmol/mol'].mean()
    summary_results['avg_xCH4exh_μmol/mol'] = time_aligned_data['xCH4exh_μmol/mol'].mean()
    summary_results['avg_xN2O_μmol/mol'] = time_aligned_data['conN2O_Avg_ppm'].mean()
    summary_results['avg_xNMHCexh_μmol/mol'] = time_aligned_data['xNMHCexh_μmol/mol'].mean()

    SamplePeriod_s = constants['SamplePeriod_s']

    summary_results['total_dilute_flow_mol'] = time_aligned_data['ndil_mol/sec'].sum() * SamplePeriod_s
    summary_results['cycle_work_kWh'] = (
            ((time_aligned_data['Power_kW'] > 0) * time_aligned_data['Power_kW'] *
             SamplePeriod_s).sum() / 3600)

    # calculate sample mass grams
    summary_results['mCO2_g'] = time_aligned_data['mCO2_g/sec'].sum() * SamplePeriod_s
    summary_results['mCO_g'] = time_aligned_data['mCO_g/sec'].sum() * SamplePeriod_s
    summary_results['mNOx_g'] = time_aligned_data['mNOxcorrected_g/sec'].sum() * SamplePeriod_s
    summary_results['mTHC_g'] = time_aligned_data['mTHC_g/sec'].sum() * SamplePeriod_s
    summary_results['mCH4_g'] = time_aligned_data['mCH4_g/sec'].sum() * SamplePeriod_s
    summary_results['mN2O_g'] = time_aligned_data['mN2O_g/sec'].sum() * SamplePeriod_s
    summary_results['mNMHC_g'] = time_aligned_data['mNMHC_g/sec'].sum() * SamplePeriod_s
    summary_results['mNMHC_g+mNOx_g'] = summary_results['mNMHC_g'] + summary_results['mNOx_g']

    # calculate background mass grams
    from constants import constants

    if not drift_corrected:
        BagData = phdp_globals.test_data['BagData']
    else:
        BagData = phdp_globals.test_data['drift_corrected_BagData']

    for background_component in ['CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O']:
        background_conc = (
            BagData.loc[(BagData['RbComponent'] == str.upper(background_component)) &
                        (BagData['EmissionsCycleNumber_Integer'] == emissions_cycle_number), 'RbAmbConc_ppm'].item())
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
    DiluteRFCH4_Fraction = phdp_globals.test_data['TestParameters']['DiluteRFCH4_Fraction'].item()

    summary_results['mNMHCbkgrnd_g'] = \
        (summary_results['total_dilute_flow_mol'] * constants['MNMHC_g/mol'] *
         (THC_background_conc - CH4_background_conc * DiluteRFCH4_Fraction) / 10 ** 6)

    summary_results['mNMHCbkgrnd_g+mNOxbkgrnd_g'] = summary_results['mNMHCbkgrnd_g'] + summary_results['mNOxbkgrnd_g']

    # calculate net mass grams
    for component in ['CO2', 'CO', 'NOx', 'THC', 'CH4', 'N2O', 'NMHC']:
        summary_results['m%snet_g' % component] = \
            (summary_results['m%s_g' % component] -
             summary_results['m%sbkgrnd_g' % component])

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

    return summary_results


def calc_1036_results(drift_corrected_time_aligned_data, drift_corrected_time_aligned_data_summary_results,
                      emissions_cycle_number, test_type, vehicle_test):
    """
    Calculate 1036 summary and carbon balance error check results

    Args:
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
    mDEF_g = drift_corrected_time_aligned_data[
                 'DEFMassFlowRate_Avg_g/h'].sum() * SamplePeriod_s / 3600
    mdot_avg_CO2DEF = (mDEF_g * constants['MCO2_g/mol'] * constants['wCH4N2O_Mass Fraction of urea in DEF'] /
                       constants['MCH4N2O_g/mol'])

    # CFR 1036.540-5
    mfuel_term = drift_corrected_time_aligned_data['CVSFlow_mol/s'] * drift_corrected_time_aligned_data[
        'xCcombdry_mol/mol'] / (1 + drift_corrected_time_aligned_data['xH2Oexhdry_mol/mol']) * SamplePeriod_s

    mfuel_cycle = constants['MC_g/mol'] / wCmeas * (sum(mfuel_term) - mdot_avg_CO2DEF / constants['MCO2_g/mol'])

    # CFR 1036.535-4
    mfuel_g = drift_corrected_time_aligned_data['qmFuel_Avg_g/h'].sum() * SamplePeriod_s / 3600
    calculations_1036['mfuelcor_meas'] = mfuel_g * EmfuelCmeas / constants['EmfuelCref	MJ/kg'] / constants['wCref']

    # CFR 1036.535-4
    calculations_1036['mfuelcor_dil'] = mfuel_cycle * EmfuelCmeas / constants['EmfuelCref	MJ/kg'] / constants['wCref']

    if test_type == 'modal':
        # units g/s instead of g:
        calculations_1036['mfuelcor_meas'] = calculations_1036['mfuelcor_meas'] / SamplePeriod_s
        calculations_1036['mfuelcor_dil'] = calculations_1036['mfuelcor_dil'] / SamplePeriod_s

    TestDetails = phdp_globals.test_data['TestDetails']

    if test_type == 'transient' and vehicle_test:
        simulation_average_vehicle_speed_mps = \
            TestDetails[TestDetails['EmissionsCycleNumber_Integer'] == emissions_cycle_number][
                'CycleAverageVehicleSpeed_m/s'].item()

        calculations_1036['simulation_average_vehicle_speed_mps'] = simulation_average_vehicle_speed_mps

        calculations_1036['CycleAverageEngineWork_kWh'] = sum(
            drift_corrected_time_aligned_data['Power_kW'] * drift_corrected_time_aligned_data[
                'VehicleMoving_Logical'] * (
                    drift_corrected_time_aligned_data['Power_kW'] > 0)) * SamplePeriod_s / 3600

        calculations_1036['CycleAverageIdleSpeed_rpm'] = (
            (drift_corrected_time_aligned_data['spDyno_rev/min'][
                drift_corrected_time_aligned_data['VehicleMoving_Logical'] == 0]).mean())

        calculations_1036['CycleAverageTorque_Nm'] = (
            (drift_corrected_time_aligned_data['tqShaft_Nm'][
                drift_corrected_time_aligned_data['VehicleMoving_Logical'] == 0]).mean())

        calculations_1036['EngineToVehicleSpeedRatio_rev/mi'] = (
                (drift_corrected_time_aligned_data['spDyno_rev/min'][drift_corrected_time_aligned_data[
                                                                         'VehicleMoving_Logical'] == 1]).mean()
                / 60 / simulation_average_vehicle_speed_mps)

    # Carbon Balance Error Check calculations:
    # CFR 1065.643-6
    calculations_1036['mCexh_g'] = (
            constants['MC_g/mol'] *
            (drift_corrected_time_aligned_data_summary_results['mCO2net_g'] / constants['MCO2_g/mol'] +
             drift_corrected_time_aligned_data_summary_results['mCOnet_g'] / constants['MCO_g/mol'] +
             drift_corrected_time_aligned_data_summary_results['mTHCnet_g'] / constants['MTHC_g/mol']))

    # CFR 1065.643-1
    nint_mol = drift_corrected_time_aligned_data['nint_mol/sec'].sum() * SamplePeriod_s
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
    if test_type == 'transient':
        if 'CycleDuration_s' in TestDetails:
            t_s = TestDetails[TestDetails['EmissionsCycleNumber_Integer'] == emissions_cycle_number][
                'CycleDuration_s'].item()
        else:
            t_s = TestDetails[TestDetails['EmissionsCycleNumber_Integer'] == emissions_cycle_number][
                'PMSampleTime_s'].item()

    else:
        t_s = SamplePeriod_s
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
                    filedialog.askopenfilename(title='Select any test file', filetypes=[('csv', '*.csv')])

            horiba_filename = file_io.get_filename(phdp_globals.options.horiba_file)

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

            results = \
                {'1036_calculations': [], 'tad': [], 'tadsummary': [], 'dctad': [], 'dctadsummary': []}

            phdp_log.logwrite('\nProcessing test %s (%s) from %s...\n' % (test_num, test_type, test_site))

            # load all raw data, even if not all required for now
            load_data(test_site)

            if test_type == 'transient':
                emissions_cycles = \
                    [ecn for ecn in phdp_globals.test_data['ContinuousData']['EmissionsCycleNumber_Integer'].unique()
                     if ecn >= 1]
            else:
                emissions_cycles = phdp_globals.test_data['ModalTestData']['ModeNumber_Integer'].values

            for ecn in emissions_cycles:
                print('\nProcessing %s %d ...' % (test_type, ecn))

                # pull in raw data and time align as necessary
                if test_type == 'transient':
                    emissions_cycle_number = ecn
                    time_aligned_data = time_align_continuous_data(test_site, vehicle_test, sampled_crank,
                                                                   emissions_cycle_number, min_mode_number)
                else:
                    emissions_cycle_number = 1
                    mode_number = ecn
                    time_aligned_data = phdp_globals.test_data['ModalTestData'].loc[
                        phdp_globals.test_data['ModalTestData']['ModeNumber_Integer'] == mode_number].copy()
                    time_aligned_data['SampleTime_s'] = (
                        phdp_globals.test_data)['ModeValidationResults']['SampleTime_s'][mode_number-1]
                    constants['SamplePeriod_s'] = time_aligned_data['SampleTime_s'].item()                    
                    non_numeric_columns = [c for c in time_aligned_data.columns
                                           if pd.api.types.is_object_dtype(time_aligned_data[c])]
                    time_aligned_data = time_aligned_data.drop(non_numeric_columns, axis=1)
                    time_aligned_data['tqShaft_Nm'] = time_aligned_data['tqShaft_Avg_Nm']
                    time_aligned_data['spDyno_rev/min'] = time_aligned_data['spDyno_Avg_rev/min']
                    # time_aligned_data['qmFuel_g/h'] = time_aligned_data['qmFuel_Avg_g/h']
                    time_aligned_data['tIntakeAir_°C'] = time_aligned_data['tIntakeAir_Avg_°C']
                    time_aligned_data['tCellDewPt_°C'] = time_aligned_data['tCellDewPt_Avg_°C']
                    time_aligned_data['pCellAmbient_kPa'] = time_aligned_data['pCellAmbient_Avg_kPa']
                    time_aligned_data = time_aligned_data.dropna(axis=1).reset_index(drop=True)

                # add calculated values
                pre_chemical_balance_calculations(time_aligned_data, test_type)

                # chemical balance iteration to calculate xDil/Exh_mol/mol, xH2Oexh_mol/mol and xCcombdry_mol/mol
                iterate_chemical_balance(time_aligned_data, emissions_cycle_number)

                post_chemical_balance_calculations(time_aligned_data)

                time_aligned_data_summary = (
                    calc_summary_results(time_aligned_data, emissions_cycle_number, test_type))

                drift_corrected_time_aligned_data = time_aligned_data.copy()

                # drift-correct concentrations
                for signal_name in [col for col in drift_corrected_time_aligned_data.columns if col.startswith('con')]:
                    drift_correct_continuous_data(drift_corrected_time_aligned_data, signal_name)

                # drift-correct bag values
                phdp_globals.test_data['drift_corrected_BagData'] = phdp_globals.test_data['BagData'].copy()
                for idx in phdp_globals.test_data['drift_corrected_BagData'].index:
                    if phdp_globals.test_data['drift_corrected_BagData'].loc[idx, 'RbComponent'] != 'NMHC':
                        drift_correct_bag_data(phdp_globals.test_data['drift_corrected_BagData'], idx)

                # add calculated values
                pre_chemical_balance_calculations(drift_corrected_time_aligned_data, test_type)

                # chemical balance iteration to calculate xDil/Exh_mol/mol, xH2Oexh_mol/mol and xCcombdry_mol/mol
                iterate_chemical_balance(drift_corrected_time_aligned_data, emissions_cycle_number, drift_corrected=True)

                post_chemical_balance_calculations(drift_corrected_time_aligned_data)

                drift_corrected_time_aligned_data_summary = (
                    calc_summary_results(drift_corrected_time_aligned_data, emissions_cycle_number, test_type,
                                         drift_corrected=True))

                calculations_1036 = calc_1036_results(drift_corrected_time_aligned_data,
                                                      drift_corrected_time_aligned_data_summary,
                                                      emissions_cycle_number, test_type, vehicle_test)

                if test_type == 'transient':
                    time_aligned_data_summary['EmissionsCycleNumber_Integer'] = emissions_cycle_number
                    drift_corrected_time_aligned_data_summary['EmissionsCycleNumber_Integer'] = emissions_cycle_number
                    calculations_1036['EmissionsCycleNumber_Integer'] = emissions_cycle_number
                else:
                    time_aligned_data_summary['ModeNumber_Integer'] = mode_number
                    drift_corrected_time_aligned_data_summary['ModeNumber_Integer'] = mode_number
                    calculations_1036['ModeNumber_Integer'] = mode_number

                results['tad'].append(time_aligned_data)
                results['dctad'].append(drift_corrected_time_aligned_data)
                results['tadsummary'].append(time_aligned_data_summary)
                results['dctadsummary'].append(drift_corrected_time_aligned_data_summary)
                results['1036_calculations'].append(calculations_1036)

            print('\nSaving results...\n')

            output_prefix = horiba_filename.rsplit('.', 1)[0] + '-'

            if test_type == 'modal':
                index_name = 'ModeNumber_Integer'
            else:
                index_name = 'EmissionsCycleNumber_Integer'
                phdp_globals.test_data['drift_corrected_BagData'].to_csv(
                    phdp_globals.options.output_folder_base + output_prefix + 'dcbagdata.csv', index=False,
                    encoding=phdp_globals.options.output_encoding, errors='replace')

            pd.concat(results['tad']).set_index(index_name).to_csv(
                phdp_globals.options.output_folder_base + output_prefix + 'tad.csv',
                encoding=phdp_globals.options.output_encoding, errors='replace')

            pd.concat(results['dctad']).set_index(index_name).to_csv(
                phdp_globals.options.output_folder_base + output_prefix + 'dctad.csv',
                encoding=phdp_globals.options.output_encoding, errors='replace')

            pd.concat([pd.DataFrame(ts) for ts in results['tadsummary']]).set_index(index_name).to_csv(
                phdp_globals.options.output_folder_base + output_prefix + 'tadsummary.csv',
                encoding=phdp_globals.options.output_encoding, errors='replace')

            pd.concat([pd.DataFrame(ts) for ts in results['dctadsummary']]).set_index(index_name).to_csv(
                phdp_globals.options.output_folder_base + output_prefix + 'dctadsummary.csv',
                encoding=phdp_globals.options.output_encoding, errors='replace')

            pd.concat([pd.DataFrame(ts) for ts in results['1036_calculations']]).set_index(index_name).to_csv(
                phdp_globals.options.output_folder_base + output_prefix + '1036_calculations.csv',
                encoding=phdp_globals.options.output_encoding, errors='replace')

            if test_type == 'transient':
                generate_transient_report(output_prefix, results, test_datetime, test_type, test_num,
                                          test_site, vehicle_test)
            else:
                generate_modal_report(output_prefix, results, test_datetime, test_type, test_num,
                                          test_site, vehicle_test)

            print('done!')

            return results

    except:
        phdp_log.logwrite("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        print("### Check PHDP log for error messages ###")
        phdp_log.end_logfile("\nSession Fail")


def generate_transient_report(output_prefix, results, test_datetime, test_type, test_num, test_site,
                              vehicle_test):
    """

    Args:
        output_prefix:
        results:
        test_datetime:
        test_type:
        test_num:
        test_site:
        vehicle_test:

    Returns:

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
                for i in range(len(results['tadsummary']))] )> 1:
            cycle_name = phdp_globals.test_data['TestDetails'].loc[
                phdp_globals.test_data['TestDetails']['EmissionsCycleNumber_Integer'] == emissions_cycle_number,
                'CycleName'].item()

            set_value_at(report_df, 'Cycle ID', '%s' % cycle_name)
        else:
            set_value_at(report_df, 'Cycle ID', '%s' %
                         phdp_globals.test_data['TestDetails']['CycleName'].iloc[0])

        set_value_at(report_df, 'Original Concentration', results['tadsummary'][i].iloc[0, 0:7].values)
        set_value_at(report_df, 'Corrected Concentration', results['dctadsummary'][i].iloc[0, 0:7].values)
        set_value_at(report_df, 'Original Mass', results['tadsummary'][i].iloc[0, 9:17].values)

        set_value_at(report_df, 'Corrected Mass', results['dctadsummary'][i].iloc[0, 9:17].values)
        set_value_at(report_df, 'Original Background Mass', results['tadsummary'][i].iloc[0, 17:25].values)
        set_value_at(report_df, 'Corrected Background Mass', results['dctadsummary'][i].iloc[0, 17:25].values)

        original_net_mass = results['tadsummary'][i].iloc[0, 25:33].values
        corrected_net_mass = results['dctadsummary'][i].iloc[0, 25:33].values
        net_mass_drift_check_pct = (original_net_mass - corrected_net_mass) / original_net_mass * 100

        set_value_at(report_df, 'Original Net Mass', original_net_mass)
        set_value_at(report_df, 'Corrected Net Mass', corrected_net_mass)
        set_value_at(report_df, 'Net Mass Drift Check %', net_mass_drift_check_pct)

        # TODO: PM calculations
        set_value_at(report_df, 'Total PM Mass', None)
        set_value_at(report_df, 'BSPM', None)

        set_value_at(report_df, 'Cycle Work', results['tadsummary'][i]['cycle_work_kWh'])
        set_value_at(report_df, 'Total Dilution Flow', results['tadsummary'][i]['total_dilute_flow_mol'])

        original_brake_specific_emissions = results['tadsummary'][i].iloc[0, 33:41].values
        corrected_brake_specific_emissions = results['dctadsummary'][i].iloc[0, 33:41].values
        brake_specific_emissions_drift_check_pct = (
                (original_brake_specific_emissions - corrected_brake_specific_emissions) /
                original_brake_specific_emissions * 100)

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

        if vehicle_test:
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

        report_df.to_csv(phdp_globals.options.output_folder_base + output_prefix +
                         '%d-report.csv' % emissions_cycle_number,
                         encoding='UTF-8', index=False, header=False, float_format='%.1g')


def generate_modal_report(output_prefix, results, test_datetime, test_type, test_num, test_site,
                              vehicle_test):
    """

    Args:
        output_prefix:
        results:
        test_datetime:
        test_type:
        test_num:
        test_site:
        vehicle_test:

    Returns:

    """

    # create report header
    report_df = pd.read_csv(path + os.sep + 'modal_report_template.csv', encoding='UTF-8', header=None)
    report_df = report_df.fillna('')

    set_value_at(report_df, 'Test Date', '%s/%s/%s' % (test_datetime[4:6], test_datetime[6:8], test_datetime[0:4]))
    set_value_at(report_df, 'Test Cell', test_site)
    set_value_at(report_df, 'Test Number', test_num)
    set_value_at(report_df, 'Test Type', test_type)
    set_value_at(report_df, 'Cycle ID', '%s' % phdp_globals.test_data['TestDetails']['CycleName'].iloc[0])

    for i in range(0, len(results['dctadsummary'])):
        mode_number = results['dctadsummary'][i]['ModeNumber_Integer'].iloc[0]
        col_offset = mode_number + 1

        if mode_number > 4:
            # extend report_df
            report_df.insert(mode_number+5, ' ' * mode_number, '')

        set_value_at(report_df, 'Mode', mode_number, col_offset=col_offset)

        set_value_at(report_df, 'CVS CO2', results['dctadsummary'][i]['mCO2_g'], col_offset=col_offset)
        set_value_at(report_df, 'CVS CO', results['dctadsummary'][i]['mCO_g'], col_offset=col_offset)
        set_value_at(report_df, 'CVS NOX', results['dctadsummary'][i]['mNOx_g'], col_offset=col_offset)
        set_value_at(report_df, 'CVS HC', results['dctadsummary'][i]['mTHC_g'], col_offset=col_offset)
        set_value_at(report_df, 'CVS CH4', results['dctadsummary'][i]['mCH4_g'], col_offset=col_offset)
        set_value_at(report_df, 'CVS N2O', results['dctadsummary'][i]['mN2O_g'], col_offset=col_offset)
        set_value_at(report_df, 'CVS NMHC', results['dctadsummary'][i]['mNMHC_g'], col_offset=col_offset)
        set_value_at(report_df, 'CVS NMHC+NOx', results['dctadsummary'][i]['mNMHC_g+mNOx_g'], col_offset=col_offset)

        set_value_at(report_df, 'Background CO2', results['dctadsummary'][i]['mCO2bkgrnd_g'], col_offset=col_offset)
        set_value_at(report_df, 'Background CO', results['dctadsummary'][i]['mCObkgrnd_g'], col_offset=col_offset)
        set_value_at(report_df, 'Background NOX', results['dctadsummary'][i]['mNOxbkgrnd_g'], col_offset=col_offset)
        set_value_at(report_df, 'Background HC', results['dctadsummary'][i]['mTHCbkgrnd_g'], col_offset=col_offset)
        set_value_at(report_df, 'Background CH4', results['dctadsummary'][i]['mCH4bkgrnd_g'], col_offset=col_offset)
        set_value_at(report_df, 'Background N2O', results['dctadsummary'][i]['mN2Obkgrnd_g'], col_offset=col_offset)
        set_value_at(report_df, 'Background NMHC', results['dctadsummary'][i]['mNMHCbkgrnd_g'], col_offset=col_offset)
        set_value_at(report_df, 'Background NMHC+NOx', results['dctadsummary'][i]['mNMHCbkgrnd_g+mNOxbkgrnd_g'], col_offset=col_offset)

        set_value_at(report_df, 'Net CO2', results['dctadsummary'][i]['mCO2net_g'], col_offset=col_offset)
        set_value_at(report_df, 'Net CO', results['dctadsummary'][i]['mCOnet_g'], col_offset=col_offset)
        set_value_at(report_df, 'Net NOX', results['dctadsummary'][i]['mNOxnet_g'], col_offset=col_offset)
        set_value_at(report_df, 'Net HC', results['dctadsummary'][i]['mTHCnet_g'], col_offset=col_offset)
        set_value_at(report_df, 'Net CH4', results['dctadsummary'][i]['mCH4net_g'], col_offset=col_offset)
        set_value_at(report_df, 'Net N2O', results['dctadsummary'][i]['mN2Onet_g'], col_offset=col_offset)
        set_value_at(report_df, 'Net NMHC', results['dctadsummary'][i]['mNMHCnet_g'], col_offset=col_offset)
        set_value_at(report_df, 'Net NMHC+NOx', results['dctadsummary'][i]['mNMHCnet_g+mNOxnet_g'], col_offset=col_offset)

        set_value_at(report_df, 'Mode Time', results['dctad'][i]['SampleTime_s'], col_offset=col_offset)
        set_value_at(report_df, 'Power', results['dctad'][i]['Power_kW'], col_offset=col_offset)

        set_value_at(report_df, 'mfuelcor_meas', results['1036_calculations'][i]['mfuelcor_meas'], col_offset=col_offset)
        set_value_at(report_df, 'mfuelcor_dil', results['1036_calculations'][i]['mfuelcor_dil'], col_offset=col_offset)

        set_value_at(report_df, 'ϵrC', results['1036_calculations'][i]['erC_rel_err_%'], col_offset=mode_number)
        set_value_at(report_df, 'ϵaC', results['1036_calculations'][i]['eaC_g'], col_offset=mode_number)
        set_value_at(report_df, 'ϵaCrate', results['1036_calculations'][i]['eaCrate_g/h'], col_offset=mode_number)

        set_value_at(report_df, 'ϵrC Limit', 2, col_offset=mode_number)
        set_value_at(report_df, 'ϵaC Limit', results['1036_calculations'][i]['eaC_g_limit'], col_offset=mode_number)
        set_value_at(report_df, 'ϵaCrate Limit', results['1036_calculations'][i]['eaCrate_g/h_limit'], col_offset=mode_number)

        set_value_at(report_df, 'ϵrC Check', results['1036_calculations'][i]['erC_rel_err_%_check'], col_offset=mode_number)
        set_value_at(report_df, 'ϵaC Check', results['1036_calculations'][i]['eaC_g_check'], col_offset=mode_number)
        set_value_at(report_df, 'ϵaCrate Check', results['1036_calculations'][i]['eaCrate_g/h_check'], col_offset=mode_number)

    report_df.to_csv(phdp_globals.options.output_folder_base + output_prefix + 'report.csv', encoding='UTF-8',
                     index=False, header=False, float_format='%.1g')


if __name__ == "__main__":
    try:
        results = run_phdp(PHDPSettings())
    except:
        print("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        os._exit(-1)
