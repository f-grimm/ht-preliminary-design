
import numpy as np
import yaml

class Aircraft:
    """ 
    Vertical flight aircraft.

    Attributes
    ----------
    name : str
        Name of the aircraft.
    mtow : float
        Maximum take-off weight; kg
    empty_weight_ratio : float
        Ratio of empty weight to maximum take-off weight.
    special_equipment : float
        Special equipment mass; kg
    accessory_power : float
        Accessory power; W
    eta_transmission : float
        Efficiency of the transmission.
    
    """

    def __init__(self, filename: str):
        """ 
        Load the configuration. 

        Parameters
        ----------
        filename : str
            Path and name of the YAML file containing the configuration.

        """
        with open(filename) as file:
            data = yaml.safe_load(file)

        self.name               = data['Name']
        self.mtow               = data['Misc'].get('MTOW')
        self.empty_weight_ratio = data['Misc'].get('Empty weight ratio')
        self.special_equipment  = data['Misc'].get('Special equipment', 0)
        self.accessory_power    = data['Misc'].get('Accessory power', 0)
        self.eta_transmission   = data['Misc'].get('Transm. efficiency', 1.0)
        

    def get_thrust_and_alpha(self, gross_weight, drag, gamma, gravity):
        """ 
        Determine the thrust and angle of attack by means of a force balance 
        with drag and weight vectors (isolated rotor).

        Parameters
        ----------
        gross_weight : float
            Aircraft mass; kg
        drag : float
            Drag force acting on the aircraft; N
        gamma : float
            Climb angle; rad
        gravity : float
            Gravitational acceleration; m/s²
        
        Returns
        ------- 
        float
            Thrust; N
        float
            Angle of attack; rad

        """
        if drag == 0:
            # In hover
            alpha = 0.0
            thrust = gross_weight * gravity
        else:
            # 2D force balance
            weight_vector = np.array([0, -gross_weight * gravity])
            drag_vector = drag * np.array([np.cos(gamma), -np.sin(gamma)])
            thrust_vector = -weight_vector - drag_vector
            thrust = np.linalg.norm(thrust_vector)
            
            # Normalization
            drag_unit_vector = drag_vector / drag
            thrust_unit_vector = thrust_vector / thrust
            
            # Angle between drag and thrust
            theta = np.arccos(np.clip(np.dot(
                drag_unit_vector, thrust_unit_vector), -1.0, 1.0))

            # Angle of attack
            alpha = np.pi / 2 - theta
        
        # Thrust [N], angle of attack [rad]        
        return thrust, alpha
