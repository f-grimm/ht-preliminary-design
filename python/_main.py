
from mission import Mission
from helicopter import Helicopter
import matplotlib.pyplot as plt
import pandas as pd

def main():

    # Load mission from YAML file
    mission = Mission('data/missions/exercise_9.yaml') 

    # Load helicopter from YAML file
    concept = Helicopter('data/configurations/concept.yaml')

    # Premilinary design
    df = concept.preliminary_design(
        mission, 
        radius_range  = (3.0, 7.0, 7), 
        chord_range   = (0.2, 0.4, 5), 
        use_ew_models = True, 
        status        = True)

    # Options
    pd.set_option(
        'display.precision', 3, 'display.max_rows', 20, 
        'display.max_columns', None, 'display.width', 55)

    # Print
    print(df.sort_values(by='Fuel [kg]'))
    
    # Plots
    concept.plot_results(df, 'MTOW [kg]')
    concept.plot_results(df, 'Fuel [kg]')
    plt.show()


# Run main function
if __name__ == '__main__':
    main()
