"""
    Test site specific data

"""

from phdp import *

test_sites = dict()

test_sites['HD02'] = {
    'signals_and_delays':
        {
            'ContinuousData':
            {
                'ModeNumber_Integer': 0,
                'EmissionsCycleNumber_Integer': 0,
                'pCellAmbient_kPa': 0,
                'spDyno_rev/min': 0,
                'tqShaft_Nm': 0,

                'qmIntakeAir_Avg_kg/h': phdp_globals.test_data['TestParameters']['AirFlowDelay_s'].item(),
                'qmFuel_g/h': phdp_globals.test_data['TestParameters']['FuelFlowDelay_s'].item(),
                'DEFMassFlowRate_Avg_g/h': 0,
                'tIntakeAir_°C': 0,
                'IntakeAirPress_Avg_kPa': 0,
                'tCellDewPt_°C': 0,

                'CVSDilAirRH_Avg_%': 0,
                'CVSDilAirTemp_Avg_°C': 0,
                'conRawCO2_Avg_%vol': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawCO2_System']['SampleDelay_s'].item(),
                'conRawHCO_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawCO_System']['SampleDelay_s'].item(),
                'conRawNOX_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawNOx_System']['SampleDelay_s'].item(),
                'conRawTHC_Avg_ppmC': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawTHC_System']['SampleDelay_s'].item(),
                'conRawCH4cutter_Avg_ppmC': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawCH4_System']['SampleDelay_s'].item(),
                'conRawO2_Avg_%vol': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawO2_System']['SampleDelay_s'].item(),
                'conRawNH3_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawNH3_System']['SampleDelay_s'].item(),
                'conCO2_Avg_%vol': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilCO2_System']['SampleDelay_s'].item(),
                'conLCO_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilCO_System']['SampleDelay_s'].item(),
                'conNOX_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilNOx_System']['SampleDelay_s'].item(),
                'conN2O_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilN2O_System']['SampleDelay_s'].item(),
                'conTHC_Avg_ppmC': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilTHC_System']['SampleDelay_s'].item(),
                'conCH4cutter_Avg_ppmC': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilCH4_System']['SampleDelay_s'].item(),
                'CVSFlow_Avg_m³/s': 0,
                'BagFillFlow_Avg_l/min': 0,
            },
        },  # signals_and_delays
}  # HD02

test_sites['HD05'] = {
    'signals_and_delays':
        {
            'ContinuousData':
            {
                'ModeNumber_Integer': 0,
                'EmissionsCycleNumber_Integer': 0,
                'pCellAmbient_kPa': 0,
                'spDyno_rev/min': 0,
                'tqShaft_Nm': 0,

                'qmIntakeAir_Avg_kg/h': phdp_globals.test_data['TestParameters']['AirFlowDelay_s'].item(),
                'qmFuel_g/h': phdp_globals.test_data['TestParameters']['FuelFlowDelay_s'].item(),
                'DEFMassFlowRate_Avg_g/h': 0,
                'tIntakeAir_°C': 0,
                'IntakeAirPress_Avg_kPa': 0,
                'tCellDewPt_°C': 0,

                'CVSDilAirRH_Avg_%': 0,
                'CVSDilAirTemp_Avg_°C': 0,
                'conRawCO2_Avg_%vol': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawCO2_System']['SampleDelay_s'].item(),
                'conRawHCO_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawCO_System']['SampleDelay_s'].item(),
                'conRawNOX_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawNOx_System']['SampleDelay_s'].item(),
                'conRawTHC_Avg_ppmC': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawTHC_System']['SampleDelay_s'].item(),
                'conRawCH4cutter_Avg_ppmC': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawCH4_System']['SampleDelay_s'].item(),
                'conRawO2_Avg_%vol': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawO2_System']['SampleDelay_s'].item(),
                'conRawNH3_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'RawNH3_System']['SampleDelay_s'].item(),
                'conCO2_Avg_%vol': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilCO2_System']['SampleDelay_s'].item(),
                'conLCO_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilCO_System']['SampleDelay_s'].item(),
                'conNOX_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilNOx_System']['SampleDelay_s'].item(),
                'conN2O_Avg_ppm': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilN2O_System']['SampleDelay_s'].item(),
                'conTHC_Avg_ppmC': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilTHC_System']['SampleDelay_s'].item(),
                'conCH4cutter_Avg_ppmC': phdp_globals.test_data['EmsComponents'].loc[phdp_globals.test_data['EmsComponents']['ParameterName'] == 'DilCH4_System']['SampleDelay_s'].item(),
                'CVSFlow_Avg_m³/s': 0,
                'BagFillFlow_Avg_l/min': 0,
            },
        },  # signals_and_delays
}  # HD05