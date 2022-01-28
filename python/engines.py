
import numpy as np

class Engines():
    """ 
    Turboshaft engines as aircraft components.

    Attributes
    ----------
    number_of_engines : int
        Number of engines.
    power_available : float
        Power available per engine; W
    sfc : float
        Specific fuel consumption; kg/Wh
    a, b : float
        Engine parameters A and B, determining the fuel consumption.
    
    """

    def __init__(self, engine_data:dict):
        """ 
        Load the engine variables.

        Parameters
        ----------
        engine_data : dict
            Engine section of the configuration file.

        """
        self.number_of_engines = engine_data.get('Number of engines', 1)
        self.power_available   = engine_data.get('Power available')
        self.fuel_safety       = engine_data.get('Fuel safety factor', 1)
        self.sfc               = engine_data.get('SFC')
        self.a                 = engine_data.get('A')
        self.b                 = engine_data.get('B')


    def get_sfc(self, temperature_ratio, pressure_ratio, power):
        """ 
        Calculate the power-dependent specific fuel consumption based on the 
        engine parameters A and B. [p.310]

        Parameters
        ----------
        temperature_ratio : float
            Temperature ratio relative to mean sea level.
        pressure_ratio : float
            Pressure ratio relative to mean sea level.
        power : float
            Total power; W
        
        Returns
        -------
        float
            Specific fuel consumption; kg/Wh
        
        """
        # SFC [kg/Wh]
        return (
            self.a * pressure_ratio * np.sqrt(temperature_ratio) / power
            + self.b)
