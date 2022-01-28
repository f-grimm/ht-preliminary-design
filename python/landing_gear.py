
class LandingGear:
    """
    Landing gear as an aircraft component.

    Attributes
    ----------
    type_ : str
        Type of landing gear.
    number_of_legs : int
        Number of legs.

    """

    def __init__(self, landing_gear_data: dict):
        """
        Load the landing gear variables.

        Parameters
        ----------
        landing_gear_data : dict
            Landing gear section of the configuration file.

        """
        self.type_          = landing_gear_data.get('Type')
        self.number_of_legs = landing_gear_data.get('Number of legs')
