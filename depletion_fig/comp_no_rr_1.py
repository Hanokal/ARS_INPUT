import math
from math import pi, sqrt
import numpy as np


dx = 0.125
dy = 1
dr = 1
dt = 0.121
nx = 6144
ny = 256
nr = 256



l0 = 2.0*math.pi  # wavelength in normalized units
t0 = l0           # optical cycle in normalized units

Lx = nx*dx # viz !
Ly = ny*dy

n_patch = 256


#odhady pro rad_reakci
# quick access databases
chi0    = 0.010
g0      = 1.e3
photon_energy_max = 10.*g0*chi0*0.4

Niterations = 2280*10*8 # 2280*10*5
Lmu  = 1

def my_filter(particles):
    return (particles.px>30.)


delka =128*l0 #hrozny ja vim --- priste uz to snad prepisu ... 

n0 = 0.005 #n0             = 6*0.0005 
rel_dens=8.0125 
xi=13.5*delka
Li= 12.2*delka
xr= 22.758218047 *delka
Lr= 3.9*delka 
n     = 20


def prof_and_lens(x,r):
    
    return n0*(np.exp(-((x-xi)/Li)**n) +rel_dens*np.exp(-((x-xr)/Lr)**n))
        

Main(
    geometry = "AMcylindrical",
    number_of_AM = 2,
    interpolation_order = 2,
    
    cell_length = [dx,dy],
    grid_length  = [Lx,Ly], #256*10
    
    number_of_patches = [ n_patch, 32], #256
    
    timestep = dt,
    simulation_time = Niterations*dt,
    EM_boundary_conditions = [["PML","PML"],["PML","PML"],],
    number_of_pml_cells  = [[40,40],[40,40]],
    reference_angular_frequency_SI = 2.*pi * 3.e8/(Lmu*1.e-6),
    print_every = 100.,
    use_BTIS3_interpolation = True,
)


LoadBalancing(
    every = 20,
)

MovingWindow(
    time_start = Lx*1.15, #1.20
    velocity_x = 0.9997
)


a0=4
tau = 11.25*t0
omega_p = np.sqrt(n0)
lambda_p = 2*np.pi/omega_p
laser_waist = np.sqrt(a0)*lambda_p/np.pi
#print(laser_waist)
#laser_fwhm = 50.
#laser_waist = 180.

LaserGaussianAM(
   box_side        = "xmin",
   omega           = 1.,              # normalized units
   a0              = a0,        # normalized units
   focus           = [Lx/2, 0], # normalized units
   waist           = laser_waist,     # normalized units
   
   time_envelope   = tgaussian(fwhm=tau,center=3*tau, )
)


  


Species(
    name                        = "eon",
    
    position_initialization     = "regular",
    momentum_initialization     = "cold",
    
    particles_per_cell          = 16,
    regular_number = [2,1,8],
    mass                        = 1.0,
    charge                      = -1.0,
    
    charge_density = prof_and_lens,
    mean_velocity = [0., 0., 0.],
    
    pusher = "borisBTIS3",
    time_frozen = 0.0,
    
    boundary_conditions = [
        ["remove", "remove"],
        ["reflective", "remove"],
    ],
    radiation_model = "LL"
)


Species(
    name = 'ion',
    
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    
    particles_per_cell      = 8,
    regular_number          = [1,1,8],
    mass = 1836.,    # normalized units
    charge = 1., # normalized units
    
    number_density = prof_and_lens,#trapezoidal(n0,xvacuum=vacuum_length,xslope1 = slope,xplateau = plateau_length ),
    #gaussian(n0,xvacuum=vacuum_length,xlength=2*plateau_length,yorder=0,xfwhm =2/3 *plateau_length, xcenter = critical_point), #trapezoidal(n0,xvacuum=vacuum_length,xplateau=plateau_length),
    #temperature = [0.001], # normalized units
    boundary_conditions = [
        ["remove", "remove"],
        ["reflective", "remove"],
    ],
    time_frozen = Niterations,
    pusher = "borisBTIS3", 
)


#Species(
#    name = 'nitrogen',
#    ionization_model = 'tunnel',
#    ionization_electrons = 'electron',
#    atomic_number = 7,
#    position_initialization = 'regular',
#    momentum_initialization = 'cold',
#    particles_per_cell = 16,
#    regular_number = [2,1,8],
#    time_frozen = 100000000,
#    mass = 14.*1836.0,
#    charge = 0.0,
#    number_density = nnitr,
#    pusher = "borisBTIS3",
#    boundary_conditions = [
#        ["remove", "remove"],
#        ["reflective", "remove"],
#    ],
#)







####### STANDART ######
globalEvery =400    



DiagScalar(
    every = globalEvery
)

fields_probes               = ['Ex','Ey','Rho',]


DiagProbe(
    every                   = globalEvery,
    origin                  = [0., 1*dy, 1*dy],
    corners                 = [[Lx, 0, 0],],
    number                  = [2500],
    fields                  = ['Ey','Rho',]
 )

"""
DiagProbe(
    every = globalEvery,
    origin = [0., -Ly, 0],
    corners = [ [Lx,-Ly,0.], [0,Ly,0.] ],
    number = [700, 700],
    fields = ['Ex','Ey','Rho_eon',]
)

"""
####### ELECTRON ######
"""
DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["y",    -Ly,    Ly,    800],
             ["py",    -80,    80,    750] ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["z",    -Ly,    Ly,    800],
             ["pz",    -80,    80,    750] ]
)


DiagParticleBinning(
    deposited_quantity = "weight_charge",
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["y",    -Ly,    Ly,    1000],
             ["py",   -80,    80,    1001],
             ["gamma", 0.,    100.,    2,  "edge_inclusive"]]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["z",    -Ly,    Ly,    800],
             ["pz",    -80,    80,    750],
             ["gamma", 0.,    100.,    2,  "edge_inclusive"]]
    )
"""
DiagParticleBinning(
    deposited_quantity = "weight_charge",	
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["moving_x",    0.,    Lx,    1000],
             ["gamma",    -10,    1500,    1001] ]
)


"""
DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["y",    -Ly,    Ly,    800],
             ["py",    -100,    100,    750],
             ["gamma", 200.,    600.,    1, ]])


DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["z",    -Ly,    Ly,    800],
             ["pz",    -100,    100,    750],
             ["gamma", 200.,    600.,    1, ]])
"""

####### RADIATION ######

"""
DiagRadiationSpectrum(
    every = globalEvery,
    time_average = globalEvery ,
    species = ["eon"],
    photon_energy_axis = [photon_energy_max/1.e8,photon_energy_max, 2000,'logscale'],
    axes = []
)


DiagRadiationSpectrum(
    every = Niterations,
    time_average = Niterations,
    species = ["eon"],
    photon_energy_axis = [photon_energy_max/1.e8,photon_energy_max, 2000, 'logscale'],
    axes = []
)

from numpy import arctan2, pi

def angle(p):
    return arctan2(p.py,p.px)

DiagRadiationSpectrum(
    every = globalEvery//2,
    time_average = globalEvery//2 ,
    species = ["eon"],
    photon_energy_axis = [0., 1000., 2000, 'logscale'],
    axes = [
        [angle,-pi,pi,90]
    ]
)

DiagTrackParticles(
    species = "eon",
    every = 250,
    filter = my_filter,
) """
