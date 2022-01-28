
class Fuselage:
    """
    Fuselage as an aircraft component.

    Attributes
    ----------
    drag_area : float
        Drag area of the fuselage.
    download_stationary : float
        Relative increase in thrust due to obstructed downwash (stationary).
    number_of_seats : int
        Number of seats.

    """

    def __init__(self, fuselage_data: dict):
        """
        Load the fuselage variables.

        Parameters
        ----------
        fuselage_data : dict
            Fuselage section of the configuration file.

        """
        self.drag_area           = fuselage_data.get('Drag area')
        self.download_stationary = fuselage_data.get('Download factor', 0)
        self.number_of_seats     = fuselage_data.get('Number of seats')


    def get_fuselage_drag(self, density, flight_speed):
        """ 
        Calculate the fuselage drag based on a constant drag area.

        Parameters
        ----------
        density : float
            Density of the surrounding air; kg/m³
        flight_speed : float
            Flight speed; m/s
        
        Returns
        -------
        float
            Fuselage drag; N        
        
        """
        # Fuselage drag [N]
        return 0.5 * density * flight_speed ** 2 * self.drag_area


    def get_parasite_power(self, density, flight_speed):
        """ 
        Calculate the parasite power due to fuselage drag in forward flight.

        Parameters
        ----------
        density : float
            Density of the surrounding air; kg/m³
        flight_speed : float
            Flight speed; m/s
        
        Returns
        -------
        float
            Parasite power; W

        """
        # Parasite power [W]
        return 0.5 * density * flight_speed ** 3 * self.drag_area


    def get_download_factor_in_flight(self, advance_ratio):
        """ 
        Calculate the download factor in forward flight, assuming a linear 
        decline until advance ratio 0.5. 
        
        Parameters
        ----------
        advance_ratio : float
            Advance ratio.
        
        Returns
        -------
        float
            Download factor in forward flight.
        
        """
        # Fast forward flight
        if abs(advance_ratio) > 0.5:
            return 0

        # Download factor [-]
        return (1 - abs(advance_ratio) / 0.5) * self.download_stationary
