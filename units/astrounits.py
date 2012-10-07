from quantities import *


"""=========================
   Adding other units
   ========================="""

# astronomical units
###Distances
kpc = UnitLength('kilo parsec' , 10**3 * parsec , symbol='kpc' , aliases=['kilparsec'])
Mpc = UnitLength('mega parsec' , 10**6 * parsec , symbol='Mpc' , aliases =['megaparsec'])

###Solar units
lsun = LSun = UnitLuminousIntensity('solar luminosity' , 3.839e26 * W    , symbol='LSun'   , aliases = ['lsun' , 'Lsun'])
msun = solMass = Msun =  UnitMass( 'Solar Mass'        , 1.98892e30 * kg , symbol='Msun'   , aliases=['msun'   , 'Msun'] )
RSun = UnitLength('solar radius'                       , 6.955e8 * m     , symbol = 'RSun' , aliases = ['Rsun' , 'rsun'])

###Time
Myr = myr = UnitTime('million year'        , 1000000 * yr    , symbol='Myr'   , aliases = ['myr' , 'Myr'])
Gyr = gyr = UnitTime('giga (billion) year' , 1000000000 * yr , symbol = 'Gyr' , aliases = ['gyr' , 'Gyr'])

percent = UnitQuantity('percent', 0.01 * dimensionless, symbol='%')
metallicity = UnitQuantity('metallicity', dimensionless, symbol = 'metallicity')

