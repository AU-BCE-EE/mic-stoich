source('core.R')
source('constants.R')

# Respiration
micstoich( donor = 'C6H12O6', acceptor = 'O2')
# With biomass
micstoich( donor = 'C6H12O6', acceptor = 'O2', bioform = 'C5H7O2N', fs = 0.1)
# Mehane
micstoich( donor = 'CH4', acceptor = 'O2')

# Methanogenesis
micstoich(donor = 'C6H12O6', acceptor = 'CO2')
micstoich(donor = 'C6H20O6', acceptor = 'CO2')

# No acceptor means fermentation
micstoich(donor = 'C6H12O6', prod = 'CH3COOH')
micstoich(donor = 'C6H12O6', prod = 'CH3CH2OH')
# With biomass
micstoich(donor = 'C6H12O6', prod = 'CH3CH2OH', bioform = 'C5H7O2N', fs = 0.1)
# H2
micstoich(donor = 'C6H12O6', prod = 'H2')
# Mixed products
micstoich(donor = 'C6H12O6', prod = 'CH3COOH CH2CH2COH (H2)3')
micstoich(donor = 'C6H12O6', prod = 'H2', bioform = 'C5H7O2N', fs = 0.1)

# Mixed reactants
micstoich( donor = 'CH3COOH H2', acceptor = 'O2')

# In some cases the same donor can produce different products
# Then the product is needed
micstoich( donor = 'C6H12O6', acceptor = 'NO3-')
micstoich( donor = 'C6H12O6', acceptor = 'NO3- N2')
micstoich( donor = 'C6H12O6', acceptor = 'NO3- NH3')
micstoich( donor = 'C6H12O6', acceptor = 'NO3- NO2-')

# Errors
source('core.R')
micstoich(donor = 'C6H12O6')
micstoich(donor = 'C6H12O6', acceptor = 'Fe')
