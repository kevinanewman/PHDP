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
            if file_name != 'tad':
                unitized_columns = get_unitized_columns(input_file, encoding=phdp_globals.options.encoding[test_site])
                phdp_globals.test_data[file_name] = \
                    pd.read_csv(input_file, names=unitized_columns,
                                encoding=phdp_globals.options.encoding[test_site], encoding_errors='strict',
                                header=1, skiprows=0)
            else:
                get_unitized_columns(input_file, units_nrows=0, encoding=phdp_globals.options.output_encoding)  # dump tad columns to logfile for reference
                phdp_globals.test_data[file_name] = pd.read_csv(input_file, encoding=phdp_globals.options.output_encoding)


def time_align_continuous_data(test_site, emissions_cycle_number):
    """
    Time-align continuous data, add Vehicle Moving from cycle definition

    Args:
        test_site (str): e.g. 'HD02'
        emissions_cycle_number (int): emissions cycle number to process

    Returns:
        Dataframe of time-aligned data

    """
    from test_sites import test_sites

    ContinuousLoggerPeriod_s = phdp_globals.test_data['TestParameters']['ContinuousLoggerPeriod_s'].item()
    time_aligned_data = pd.DataFrame(index=phdp_globals.test_data['ContinuousData'].index)
    for source in test_sites[test_site]['signals_and_delays'].keys():
        print(source)
        for signal in test_sites[test_site]['signals_and_delays'][source]:
            print(signal)
            delay_s = test_sites[test_site]['signals_and_delays'][source][signal]
            delay_samples = round(delay_s / ContinuousLoggerPeriod_s)
            time_aligned_data = pd.concat([time_aligned_data,
                                           pd.DataFrame({signal: phdp_globals.test_data[source][signal]
                                                        .iloc[delay_samples:].values})], axis=1)
    time_aligned_data = \
        time_aligned_data[time_aligned_data['EmissionsCycleNumber_Integer'] == emissions_cycle_number].reset_index(drop=True)
    # add vehicle moving flag
    test_cycle_definition = phdp_globals.test_data['CycleDefinition'][phdp_globals.test_data['CycleDefinition']
                                                                      ['EmissionsCycleNumber_Integer'] == emissions_cycle_number]
    vehicle_moving_int = test_cycle_definition['VehicleMoving_Logical']
    time_aligned_data = pd.concat([time_aligned_data,
                                   pd.DataFrame({'VehicleMoving_Logical': vehicle_moving_int.values})], axis=1)
    time_aligned_data['VehicleMoving_Logical'] = 1 * time_aligned_data['VehicleMoving_Logical'].fillna(False)

    return time_aligned_data


def pre_chemical_balance_calculations(time_aligned_data):
    """
    Calculate values required for chemical balance iteration

    Args:
        time_aligned_data (dataframe): time-aligned continuous test data

    Returns:
        Nothing, updates time_aligned_data

    """
    time_aligned_data['BagFillFlow_Avg_m³/s'] = time_aligned_data['BagFillFlow_Avg_l/min'] / 60000

    time_aligned_data['CVSFlow_mol/s'] = \
        (time_aligned_data['CVSFlow_Avg_m³/s'] + time_aligned_data['BagFillFlow_Avg_m³/s'] +
         phdp_globals.test_data['TestParameters']['DiluteSampleVolumeFlow_m³/s'].item()) / 0.024055

    time_aligned_data['Tsat_K'] = time_aligned_data['CVSDilAirTemp_Avg_°C'] + 273.15

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

    # 1065.645 - 3:
    time_aligned_data['xH2Odil_mol/mol'] = \
        time_aligned_data['pH2Odil_kPa'] / time_aligned_data['pCellAmbient_kPa']

    from constants import constants, update_constants
    update_constants()  # update constants that rely on test fuel properties, etc
    time_aligned_data['Power_kW'] = \
        np.maximum(0, time_aligned_data['tqShaft_Nm'] * time_aligned_data['spDyno_rev/min'] / 9548.8)

    # 1065.655-20:
    time_aligned_data['alpha'] = \
        CFR1065.alpha(time_aligned_data['qmFuel_g/h'], time_aligned_data['DEFMassFlowRate_Avg_g/h'])

    # 1065.655-21:
    time_aligned_data['beta'] = \
        CFR1065.beta(time_aligned_data['qmFuel_g/h'], time_aligned_data['DEFMassFlowRate_Avg_g/h'])

    # 1065.655-23:
    time_aligned_data['delta'] = \
        CFR1065.delta(time_aligned_data['qmFuel_g/h'], time_aligned_data['DEFMassFlowRate_Avg_g/h'])

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
            phdp_globals.test_data['EmsComponents']['InputName'] == 'conRawHCO', 'ResidualH2O_%vol'].item()

        time_aligned_data['xCOdry_μmol/mol'] = \
            time_aligned_data['conLCO_Avg_ppm'] / (1 - residual_H2O_pctvol / 100)

        # 1065.655-15
        residual_H2O_pctvol = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['InputName'] == 'conRawCO2', 'ResidualH2O_%vol'].item()

        time_aligned_data['xCO2dry_%'] = \
            time_aligned_data['conCO2_Avg_%vol'] / (1 - residual_H2O_pctvol / 100)

        # 1065.655-16
        residual_H2O_pctvol = phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['InputName'] == 'conRawNOX', 'ResidualH2O_%vol'].item()
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

        print(iteration, time_aligned_data['xCcombdry_mol/mol'][0])
        iteration = iteration + 1


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

    signal, type, unit = signal_name.split('_')

    if signal.startswith('conRaw'):
        component = signal.replace('conRaw', '')
        driftline = 'DIRECT'
        if component == 'NH3':
            driftline = 'HOT'
    else:
        component = signal.replace('con', '')
        driftline = 'DILUTE'

    if component in rename:
        component = rename[component]

    EmsCalResults = phdp_globals.test_data['EmsCalResults']
    xpre_data = (
        EmsCalResults)[(EmsCalResults['DriftComponent'] == component) & (EmsCalResults['DriftLine'] == driftline)]

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
                # filedialog.askopenfilename(title='Select Horiba Tn File', filetypes=[('xlsm', '*.xlsm')])

            horiba_filename = file_io.get_filename(phdp_globals.options.horiba_file)

            test_site, test_datetime, test_num, test_type, _ = horiba_filename.replace('.Tn', '').split('.')

            phdp_log.logwrite('\nProcessing test %s (%s) from %s...\n' % (test_num, test_type, test_site))

            # load all raw data, even if not all required for now
            load_data(test_site)

            emissions_cycle_number = 1  # for now, eventually will loop over cycles, I imagine

            # pull in raw data and time align as necessary
            time_aligned_data = time_align_continuous_data(test_site, emissions_cycle_number)

            # add calculated values
            pre_chemical_balance_calculations(time_aligned_data)

            # chemical balance iteration to calculate xDil/Exh_mol/mol, xH2Oexh_mol/mol and xCcombdry_mol/mol
            iterate_chemical_balance(time_aligned_data, emissions_cycle_number)

            post_chemical_balance_calculations(time_aligned_data)

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
            pre_chemical_balance_calculations(drift_corrected_time_aligned_data)

            # chemical balance iteration to calculate xDil/Exh_mol/mol, xH2Oexh_mol/mol and xCcombdry_mol/mol
            iterate_chemical_balance(drift_corrected_time_aligned_data, emissions_cycle_number, drift_corrected=True)

            post_chemical_balance_calculations(drift_corrected_time_aligned_data)

            # just for development, I think:
            phdp_globals.options.output_folder_base = file_io.get_filepath(phdp_globals.options.horiba_file) + os.sep

            time_aligned_data.to_csv(phdp_globals.options.output_folder_base + 'tad.csv', index=False,
                                     encoding=phdp_globals.options.output_encoding)

            drift_corrected_time_aligned_data.to_csv(phdp_globals.options.output_folder_base + 'dctad.csv', index=False,
                                     encoding=phdp_globals.options.output_encoding)

            phdp_globals.test_data['drift_corrected_BagData'].to_csv(
                phdp_globals.options.output_folder_base + 'dcbagdata.csv', index=False,
                encoding=phdp_globals.options.output_encoding)

            print('done!')

            return time_aligned_data

    except:
        phdp_log.logwrite("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        print("### Check PHDP log for error messages ###")
        phdp_log.end_logfile("\nSession Fail")


if __name__ == "__main__":
    try:
        time_aligned_data = run_phdp(PHDPSettings())
    except:
        print("\n#RUNTIME FAIL\n%s\n" % traceback.format_exc())
        os._exit(-1)
