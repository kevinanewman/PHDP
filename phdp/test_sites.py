"""
    Test site specific data

"""

from phdp import *

site_info = {'signals_and_delays': dict()}


def get_test_parameters_sample_delay(parameter_name, optional=False):
    """
    Get sample delay from TestParameters data
    
    Args:
        parameter_name (str): signal name, e.g. 'AirFlowDelay_s'
        optional (bool) = ``True`` if signal is optional

    Returns:
        Sample delay, in seconds, or ``None`` if signal is not available

    """
    if not optional or parameter_name in phdp_globals.test_data['TestParameters'].tolist():
        return phdp_globals.test_data['TestParameters'][parameter_name].item()
    else:
        # if optional and not in test data:
        return None


def get_ems_sample_delay(parameter_name, optional=False):
    """
    Get sample delay from EmsComponents data

    Args:
        parameter_name (str): signal name, e.g. 'AirFlowDelay_s'
        optional (bool) = ``True`` if signal is optional

    Returns:
        Sample delay, in seconds, or ``None`` if signal is not available

    """
    if not optional or parameter_name in phdp_globals.test_data['EmsComponents']['ParameterName'].tolist():
        return phdp_globals.test_data['EmsComponents'].loc[
            phdp_globals.test_data['EmsComponents']['ParameterName'] == parameter_name]['SampleDelay_s'].item()
    else:
        # if optional and not in test data:
        return None


def init_site_info(test_site):
    """
    Init ``site_info`` for the given test site

    Args:
        test_site (str): test site name, e.g. 'HD02'

    Returns:
        Nothing, updates ``site_info`` dict

    """
    global site_info

    if test_site in ('HD02', 'HD05'):
        site_info['signals_and_delays']['ContinuousData'] = \
            {
                'ModeNumber_Integer': 0,
                'EmissionsCycleNumber_Integer': 0,
                'pCellAmbient_kPa': 0,
                'spDyno_rev/min': 0,
                'tqShaft_Nm': 0,

                'qmIntakeAir_Avg_kg/h': get_test_parameters_sample_delay('AirFlowDelay_s'),
                'qmFuel_Avg_g/h': get_test_parameters_sample_delay('FuelFlowDelay_s'),
                'DEFMassFlowRate_Avg_g/h': 0,
                'tIntakeAir_°C': 0,
                'IntakeAirPress_Avg_kPa': 0,
                'tCellDewPt_°C': 0,

                'CVSDilAirTemp_Avg_°C': 0,
                'conRawCO2_Avg_%vol': get_ems_sample_delay('RawCO2_System', optional=True),
                'conRawHCO_Avg_ppm': get_ems_sample_delay('RawCO_System', optional=True),
                'conRawNOX_Avg_ppm': get_ems_sample_delay('RawNOx_System', optional=True),
                'conRawTHC_Avg_ppmC': get_ems_sample_delay('RawTHC_System', optional=True),
                'conRawCH4cutter_Avg_ppmC': get_ems_sample_delay('RawCH4_System', optional=True),
                'conRawO2_Avg_%vol': get_ems_sample_delay('RawO2_System', optional=True),
                'conRawNH3_Avg_ppm': get_ems_sample_delay('RawNH3_System', optional=True),
                'conCO2_Avg_%vol': get_ems_sample_delay('DilCO2_System'),
                'conLCO_Avg_ppm': get_ems_sample_delay('DilCO_System'),
                'conNOX_Avg_ppm': get_ems_sample_delay('DilNOx_System'),
                'conN2O_Avg_ppm': get_ems_sample_delay('DilN2O_System'),
                'conTHC_Avg_ppmC': get_ems_sample_delay('DilTHC_System'),
                'conCH4cutter_Avg_ppmC': get_ems_sample_delay('DilCH4_System'),
                'BagFillFlow_Avg_l/min': 0,
                'EngDynoMode': 0,
            }
    if test_site == 'HD02':
        site_info['signals_and_delays']['ContinuousData']['CVSDilAirRH_Avg_%'] = 0
        site_info['signals_and_delays']['ContinuousData']['CVSFlow_Avg_m³/s'] = 0

    if test_site == 'HD05':
        site_info['signals_and_delays']['ContinuousData']['CVSDilAirDPTemp_Avg_°C'] = 0
        site_info['signals_and_delays']['ContinuousData']['CVSMolarFlow_Avg_mol/s'] = 0
