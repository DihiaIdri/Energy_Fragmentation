from GradyEnergyCriterion import GradyEnergyCriterion
import pandas as pa


def metal_properties():
    metal_data = pa.read_excel('/Users/dihiaidrici/Desktop/SecondPaper/ShockInMetalProperties.xlsx', sheet_name='Fracture_Properties', header=1, index_col=0)
    strain_rate = 10**5  # s**-1

    # remember that indexing begins at 0
    Al = metal_data.iloc[0]
    Ti = metal_data.iloc[1]
    Zn_compact = metal_data.iloc[2]  # Joe Hooper: Impact fragmentation of brittle metal compact - zinc powder
    Zr = metal_data.iloc[3]

    # Al_results = GradyEnergyCriterion(Al, strain_rate)
    # Zn_results = GradyEnergyCriterion(Zn_compact, strain_rate)
    Ti_results = GradyEnergyCriterion(Ti, strain_rate)


if __name__ == '__main__':
    metal_properties()


