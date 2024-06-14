import os
path = os.path.dirname(os.path.abspath(__file__))

from phdp import *
from constants import constants

import pandas as pd

from metpy.calc import specific_humidity_from_dewpoint
from metpy.units import units


def set_average_min_max(report_df, dctad, value_name, signal_name, col_offset, scale=1):
    found = False

    if type(dctad) is pd.Series:
        values = dctad
        found = True
    else:
        if signal_name in dctad:
            values = dctad[signal_name]
            found = True
    if found:
        set_value_at(report_df, value_name, [values.mean() * scale], col_offset=col_offset)
        set_value_at(report_df, value_name, [values.min() * scale], col_offset=col_offset+1)
        set_value_at(report_df, value_name, [values.max() * scale], col_offset=col_offset+2)

    return found


def generate_general_report(report_filename, calc_mode, results, validation_results,
                            test_type, test_datetime, test_site):
    """
    Generate general test report

    Args:
        report_filename (str): the Excel report filename, a new sheet will be added
        calc_mode (str): the mode used during the emissions calculations - 'dilute', 'raw', etc.
        results (dict): a dictionary containing the results of the emissions calculations
        validation_results (dict): dict of results from best cycle validation (min error, min omits, min shift)
        test_type (str): test type, i.e. 'transient' or 'modal'
        test_datetime (str): the date and time of the test in YYYMMDDhhmm format.
        test_site (str): the name of the site where the test was performed, e.g. 'HD02'

    """
    report_df = pd.read_csv(path + os.sep + 'report_templates' + os.sep + 'general_report_template.csv',
                            encoding='UTF-8', header=None)

    report_df = report_df.fillna('')

    if test_type == 'modal':
        emissions_cycles = [1]
        modes = [results['tadsummary'][i]['ModeNumber_Integer'].iloc[0]
                            for i in range(0, len(results['tadsummary']))]
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

        if test_type == 'transient':
            # drift corrected time-aligned data
            dctad = results['dctad'][emissions_cycle_number-1]
        else:
            dctad = pd.concat([results['dctad'][mode-1] for mode in modes])

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
        set_average_min_max(report_df, dctad, 'Exhaust Back-pressure', 'ExhaustBackPressure_Avg_kPa', col_offset=2)
        set_average_min_max(report_df, dctad, 'Fuel Pump Inlet Temperature', 'tFuel_Avg_°C', col_offset=2)
        set_average_min_max(report_df, dctad, 'Fuel Pump Supply Pressure', 'pFuelSupply_Avg_kPa', col_offset=2)
        set_average_min_max(report_df, dctad, 'Fuel Pump Return Pressure', 'pFuelReturn_Avg_kPa', col_offset=2)

        if set_average_min_max(report_df, dctad, 'CA Coolant Temperature', 'tCoolantCA_Avg_°C', col_offset=2):
            pass_fail = (dctad['tCoolantCA_Avg_°C'].min() >= 20 and
                         dctad['tCoolantCA_Avg_°C'].max() <= 30)
            set_value_at(report_df, 'CA Coolant Temperature', pass_fail_range(pass_fail, [True, True]), col_offset=7)

        set_average_min_max(report_df, dctad, 'Coolant Temperature', 'tCoolantIn_Avg_°C', col_offset=2)
        set_average_min_max(report_df, dctad, 'PM Trap Face Temperature', 'tEngPmTrapFace_Avg_°C', col_offset=2)
        set_average_min_max(report_df, dctad, 'Oil Sump Temperature', 'tOilSump_Avg_°C', col_offset=2)
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

        EmsCalResults_ecns = EmsCalResults['EmissionsCycleNumber_Integer'].unique()

        if len(EmsCalResults_ecns) == 1:
            emscal_emissions_cycle_number = EmsCalResults_ecns[0]
        else:
            emscal_emissions_cycle_number = emissions_cycle_number

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
                (EmsCalResults['DriftComponent'] == component) & (EmsCalResults['DriftLine'] == driftline) &
                (EmsCalResults['EmissionsCycleNumber_Integer'] == emscal_emissions_cycle_number)]

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

        if validation_results['PM_results']:
            if test_type == 'modal':
                cvs_dilution_ratios = \
                    pd.Series([validation_results['PM_results'][i]['CVS Dilution Ratio'].item()
                               for i in validation_results['PM_results']])
            else:
                cvs_dilution_ratios = validation_results['PM_results'][emissions_cycle_number]['CVS Dilution Ratio']
            set_average_min_max(report_df, cvs_dilution_ratios, 'CVS Tunnel Dilution Ratio',
                                'CVS Dilution Ratio', col_offset=2)
            pass_fail = cvs_dilution_ratios.min() >= 2.0
            set_value_at(report_df, 'CVS Tunnel Dilution Ratio', pass_fail_range(pass_fail, [True, True]), col_offset=7)

            if test_type == 'modal':
                psu_dilution_ratios = \
                    pd.Series([validation_results['PM_results'][i]['PSU Dilution Ratio']
                               for i in validation_results['PM_results']])
            else:
                psu_dilution_ratios = pd.Series(validation_results['PM_results'][
                                                    emissions_cycle_number]['PSU Dilution Ratio'])
            set_average_min_max(report_df, psu_dilution_ratios, 'PSU Dilution Ratio', 'PSU Dilution Ratio',
                                col_offset=2)

            if test_type == 'modal':
                overall_dilution_ratios = \
                    pd.Series([validation_results['PM_results'][i]['Overall Dilution Ratio']
                               for i in validation_results['PM_results']])
            else:
                overall_dilution_ratios = \
                    validation_results['PM_results'][emissions_cycle_number]['Overall Dilution Ratio']

            set_average_min_max(report_df, overall_dilution_ratios, 'Overall Dilution Ratio',
                                   'Overall Dilution Ratio', col_offset=2)
            pass_fail = (5.0 <= overall_dilution_ratios.min() <= 7.0)
            set_value_at(report_df, 'Overall Dilution Ratio', pass_fail_range(pass_fail, [True, True]), col_offset=7)

            set_value_at(report_df, 'Proportionality Check (SEE/mean flow)',
                         [validation_results['PM_results'][1]['Proportionality_pct'].mean() / 100], col_offset=2)
            pass_fail = validation_results['PM_results'][1]['Proportionality_pct'].mean() <= 3.5
            set_value_at(report_df, 'Proportionality Check (SEE/mean flow)', pass_fail_range(pass_fail, [True, True]),
                         col_offset=7)

            set_value_at(report_df, 'Filter Pressure Drop Increase',
                         [validation_results['PM_results'][1]['FilterPressureDrop_kPa'][-1] -
                          validation_results['PM_results'][1]['FilterPressureDrop_kPa'][0]],
                         col_offset=2)

        with pd.ExcelWriter(
                report_filename,
                mode="a",
                engine="openpyxl",
                if_sheet_exists="replace",
        ) as writer:
            report_df.to_excel(writer, index=False, header=False,
                               sheet_name='General %s %d' % (calc_mode, emissions_cycle_number))

