"""
    Test site specific data

"""

test_sites = dict()

test_sites['HD02'] = {
    'signals_and_delays':
        {'horiba_data':
            {
                'Running_Logical': 0,
                'ModeNumber_Integer': 0,
                'CycleNumber_Integer': 0,
                'pCellAmbient_kPa': 0,
                'spDyno_rev/min': 0,
                'tqShaft_Nm': 0,
                'qmIntakeAir_kg/h': 0,
                'qmFuel_g/h': phdp_globals.test_data['TestParameters']['FuelFlowDelay_s'].item(),
                'DEFMassFlowRate_g/h': 0,
                'tIntakeAir_°C': 0,
                'IntakeAirPress_kPa': 0,
                'tCellDewPt_°C': 0,
                'CVSDilExhTemp_°C': 0,
                'CVSDilAirRH_%': 0,
                'CVSDilAirTemp_°C': 0,
                'conRawCO2_%mol': 7.1,
                'conRawHCO_µmol/mol': 7.5,
                'conRawNOX_µmol/mol': 5.5,
                'conRawTHC_µmol/mol': 4.1,
                'conRawCH4cutter_µmol/mol': 4.3,
                'conRawO2_%mol': 0,  # or 7.6?
                'conRawNH3_µmol/mol': 0,  # or 5.0?
                'conCO2_%mol': 9.3,
                'conLCO_µmol/mol': 10.3,
                'conNOX_µmol/mol': 6.5,
                'conN2O_µmol/mol': 9.0,
                'conTHC_µmol/mol': 6.3,
                'conCH4cutter_µmol/mol': 7.5,
                'CVSFlow_m³/s': phdp_globals.test_data['TestParameters']['AirFlowDelay_s'].item(),
            },
            'CycleDefinition': {
                'VehicleMoving_Logical': 0,
            }
        }
}