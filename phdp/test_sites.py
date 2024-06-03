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


def get_optional_continuous_signal(parameter_name, default_value, sample_delay):
    """
    Get sample delay optional data

    Args:
        parameter_name (str): signal name, e.g. 'AirFlowDelay_s'
        default_value (numeric): the default value for the optional signal if not present
        sample_delay (numeric): signal delay

    Returns:
        sample delay or ``None`` if signal is not available

    """
    if parameter_name not in phdp_globals.test_data['ContinuousData']:
        phdp_globals.test_data['ContinuousData'][parameter_name] = default_value

    return sample_delay


def init_site_info(test_site, test_type):
    """
    Init ``site_info`` for the given test site

    Args:
        test_site (str): test site name, e.g. 'HD02'
        test_type (str): test type, i.e. 'transient' or 'modal'

    Returns:
        Nothing, updates ``site_info`` dict

    """
    global site_info

    if test_site in ('HD02', 'HD05'):
        if test_type == 'transient':
            site_info['signals_and_delays']['ContinuousData'] = \
                {
                    'Time_Date': 0,
                    'ModeNumber_Integer': 0,
                    'EmissionsCycleNumber_Integer': 0,
                    'pCellAmbient_Avg_kPa': 0,
                    'spDyno_Avg_rev/min': 0,
                    'tqShaft_Avg_Nm': 0,

                    'qmIntakeAir_Avg_kg/h': get_test_parameters_sample_delay('AirFlowDelay_s'),
                    'qmFuel_Avg_g/h': get_test_parameters_sample_delay('FuelFlowDelay_s'),
                    'DEFMassFlowRate_Avg_g/h': get_optional_continuous_signal('DEFMassFlowRate_Avg_g/h',
                                                                              default_value=0, sample_delay=0),
                    'tIntakeAir_Avg_°C': 0,
                    'IntakeAirPress_Avg_kPa': 0,
                    'tCellDewPt_Avg_°C': 0,

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
                    'ExhaustBackPressure_Avg_kPa': 0,
                    'CVSDilExhTemp_Avg_°C': 0,
                    'pTailpipe_Avg_kPa': 0,
                    'tFuel_Avg_°C': 0,

                    # use _Avg_ signals if available, fallback to non-_Avg_
                    'pFuelSupply_Avg_kPa': get_optional_continuous_signal('pFuelSupply_Avg_kPa',
                                                                          default_value=phdp_globals.test_data[
                                                                              'ContinuousData']['pFuelSupply_kPa'],
                                                                          sample_delay=0),

                    'pFuelReturn_Avg_kPa': get_optional_continuous_signal('pFuelReturn_Avg_kPa',
                                                                          default_value=phdp_globals.test_data[
                                                                              'ContinuousData']['pFuelReturn_kPa'],
                                                                          sample_delay=0),

                    'tCoolantCA_Avg_°C': get_optional_continuous_signal('tCoolantCA_Avg_°C',
                                                                          default_value=phdp_globals.test_data[
                                                                              'ContinuousData']['tCoolantCA_°C'],
                                                                          sample_delay=0),

                    'tCoolantIn_Avg_°C': get_optional_continuous_signal('tCoolantIn_Avg_°C',
                                                                          default_value=phdp_globals.test_data[
                                                                              'ContinuousData']['tCoolantIn_°C'],
                                                                          sample_delay=0),

                    'tOilSump_Avg_°C': get_optional_continuous_signal('tOilSump_Avg_°C',
                                                                          default_value=phdp_globals.test_data[
                                                                              'ContinuousData']['tOilSump_°C'],
                                                                          sample_delay=0),

                    'tEngPmTrapFace_Avg_°C': get_optional_continuous_signal('tEngPmTrapFace_Avg_°C',
                                                                          default_value=get_optional_continuous_signal('tEngPmTrapFace_°C', default_value=0, sample_delay=0),
                                                                          sample_delay=0),
                }

            if test_site == 'HD02':
                site_info['signals_and_delays']['ContinuousData']['CVSDilAirRH_Avg_%'] = 0
                site_info['signals_and_delays']['ContinuousData']['CVSFlow_Avg_m³/s'] = 0

            if test_site == 'HD05':
                site_info['signals_and_delays']['ContinuousData']['CVSDilAirDPTemp_Avg_°C'] = 0
                site_info['signals_and_delays']['ContinuousData']['CVSMolarFlow_Avg_mol/s'] = 0

        else:
            site_info['optional_modal_signals'] = {
                'DEFMassFlowRate_Avg_g/h': get_optional_continuous_signal('DEFMassFlowRate_Avg_g/h',
                                                                      default_value=0, sample_delay=0)
        }
