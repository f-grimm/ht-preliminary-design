
import yaml
import numpy as np
import matplotlib.pyplot as plt

class Mission:
    """
    Mission profile. 
    
    Attributes
    ----------
    name : str
        Name of the mission.
    duration : list[float]
        Duration of each segment; h
    payload : list[float]
        Payload; kg
    crew_mass : list[float]
        Crew mass; kg
    flight_speed : list[float]
        Flight speed; m/s
    height : list[float]
        Height above ground; m, length n + 1 as defined between segments. In 
        case of climb/descent, the half-way height is used to calculate 
        density, temperature, and pressure.
    gravity : list[float]
        Gravitational acceleration; m/s²
    density : list[float]
        Density of the surrounding air; kg/m³
    temperature : list[float]
        Temperature of the surrounding air; K
    pressure : list[float]
        Pressure of the surrounding air; Pa
    climb_angle : list[float]
        Climb angle; rad

    """

    def __init__(self, filename: str):
        """
        Load the mission data.

        Parameters
        ----------
        filename : str
            Path and name of the YAML file containing the mission.

        """
        with open(filename) as file:
            data = yaml.safe_load(file)

        self.name         = data['Name']
        self.duration     = data['Mission']['Duration']
        self.payload      = data['Mission']['Payload']
        self.crew_mass    = data['Mission']['Crew mass']
        self.flight_speed = data['Mission']['Flight speed']
        self.height       = data['Mission']['Height']
        
        # Calculation
        climb = np.diff(data['Mission']['Height'])
        temp_offset = data['Mission']['Temperature offset']
        height_mid = np.array(self.height[:-1]) + climb / 2
        self.gravity = (
            data['Mission'].get('Gravity', [9.81] * len(self.duration)))
        self.density, self.temperature, self.pressure = tuple(
            list(values) for values in zip(*map(
                atmosphere, height_mid, temp_offset)))
        self.climb_angle = list(map(
            get_climb_angle, self.flight_speed, climb, self.duration))
        
        
    def plot_mission(self):
        """
        Plot the height and payload over the duration of the mission.

        """
        # Time, double elements to plot discontinuous lines
        time = [0] + [t for t in np.cumsum(self.duration) for _ in (0, 1)][:-1]

        # Payload, constant in each segment
        payload = [m for m in self.payload for _ in (0, 1)]
        
        # Height, linear in each segment
        height = [h for h in self.height for _ in (0, 1)][1:-1]

        # Plot
        fig, ax = plt.subplots(figsize=(8,5))
        ax2 = ax.twinx()
        ax.set_title('Mission: ' + self.name, fontweight='bold')
        ax.plot(time, height, color='tab:red', label='Height')
        ax.set_xlabel('Time [h]')
        ax.set_ylabel('Height [m]')
        ax.legend(loc='lower left')
        ax.set_ylim(bottom=0)
        ax.grid()
        ax2.plot(time, payload, color='tab:blue', label='Payload')
        ax2.set_ylim(bottom=0, top=max(payload) / 0.6)
        ax2.legend(loc='lower right')
        ax2.set_ylabel('Mass [kg]')
        fig.tight_layout()


def atmosphere(height, temp_offset):
    """
    Calculate the density, temperature, and pressure at a given height 
    based on the international standard atmosphere (ISA). Deviation from 
    the ISA is considered via a temperature offset. [p. 278]
    
    Parameters
    ----------
    height : float
        Height; m
    temp_offset : float
        Temperature offset compared to the ISA; K or °C

    Returns
    -------
    float
        Density; kg/m³
    float
        Temperature; K
    float
        Pressure; Pa
    
    """
    TEMPERATURE_MSL = 288.15
    PRESSURE_MSL = 101325
    RADIUS_EARTH = 6356766
    GAS_CONST = 287.05     
    GAMMA = -0.0065     
    N = 1.235
    
    geopotential_height = RADIUS_EARTH * height / (RADIUS_EARTH + height)
    temperature_isa = TEMPERATURE_MSL + GAMMA * geopotential_height
    temperature = temperature_isa + temp_offset
    pressure = (
        PRESSURE_MSL * (1 + GAMMA / TEMPERATURE_MSL * geopotential_height) 
        ** (N / (N - 1)))
    density = pressure / (GAS_CONST * temperature)

    # Density [kg/m³], temperature [K], pressure [Pa]
    return density, temperature, pressure


def get_climb_angle(flight_speed, climb, duration):
    """
    Calculate the climb angle.
    
    Parameters
    ----------
    flight_speed : float
        Flight speed; m/s
    climb : float
        Vertical distance; m
    duration : float
        Duration; h

    Returns
    -------
    float
        Climb angle; rad

    """
    # In hover
    if flight_speed == 0: 
        return 0.0

    # Climb velocity [m/s]
    climb_velocity = climb / (duration * 60 ** 2)
    
    if abs(climb_velocity) > abs(flight_speed): 
        print('\nWarning: Flight speed insufficient.')

    # Climb angle [rad]
    return np.arcsin(np.clip(climb_velocity / flight_speed, -1, 1))
