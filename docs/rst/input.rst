
Input
=====

Missions and helicopter configurations are provided as ``.yaml``-files in '*python/data*'. Examples are given below.

Missions
--------

.. code-block:: yaml

    Name: name

    Mission:
        
        # Constant within segment
        Duration: [0.1, 0.1] # h
        Payload: [0, 100] # kg
        Crew mass: [0, 0] # kg
        Flight speed: [10, 10] # m/s
        Temperature offset: [0, 0] # °C (deviation from ISA)
        Gravity: [9.76, 9.80] # m/s² (optional, default 9.81)
        
        # Linear with defined points between segments
        Height: [0, 0, 0] # m


Configurations
--------------

.. code-block:: yaml

    Name: name

    Main rotor:
        
        Number of blades: 4
        Kappa: 1.15 # (real/ideal induced power)
        Zero-lift drag coeff.: 0.011 # (avg., for profile power)
        Tip velocity: 213 # m/s
        Installation height: 2.5 # m (for in-ground effect)

    Tail rotor:
        
        # Other rotor attributes are possible, but not used.
        Power fraction: 0.05

    Engines:

        Number of engines: 1
        Power available: 500_000 # W (per engine)
        SFC: 0.00035 # kg/Wh 
        A: 40 # kg/h (A, B are optional, for power-dependent SFC)
        B: 0.00025 # kg/Wh

    Fuselage:

        Download factor: 0.05 # (thrust increase due to fuselage download)
        Drag area: 1.5 # m²
        Number of seats: 3 # (used in empty weight estimation)

    Landing gear:

        Type: Rigid wheels # {Skids, Retractable wheels, Rigid wheels}
        Number of legs: 2 # (used in empty weight estimation)

    Misc:

        Empty weight ratio: 0.5 # (EW / MTOW)
        Accessory power: 48_000 # W [p. 271]
        Transmission efficiency: 0.98 # (default 1.0)
        MTOW: 1500 # kg (not needed for sizing)