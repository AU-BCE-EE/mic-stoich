source('core.R')

micstoich( donor = 'C6H12O6', acceptor = 'O2', bioform = 'C5H7O2N', fs = 0.1, order = 'sort')
micstoich( donor = 'C6H12O6', acceptor = 'CO2', bioform = 'C5H7O2N', fs = 0.1, order = 'sort')
micstoich( donor = 'C6H12O6', prod = 'CH3COOH', bioform = 'C5H7O2N', fs = 0.1, order = 'sort')


micstoich(
  donor = donor, 
  acceptor , 
  prod = NULL,
  bioform = 'C5H7O2N',
  fs, 
  elements = c('C', 'H', 'O', 'N'),
  dropzero, 
  dropsub, 
  order, 
  tol
) 


