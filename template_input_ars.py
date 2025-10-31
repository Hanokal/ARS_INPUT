import math
from math import pi, sqrt
import numpy as np

# grid settings
dx = 0.125
dy = 1
dr = 1
dt = 0.121
nx = 6144
ny = 256
nr = 256



l0 = 2.0*math.pi  # wavelength in normalized units
t0 = l0           # optical cycle in normalized units

Lx = nx*dx 
Ly = ny*dy

n_patch = 256



# quick access databases
chi0    = 0.010
g0      = 1.e3
photon_energy_max = 10.*g0*chi0*0.4

Niterations = 2280*10*12 # number of iterations



delka =128*l0   #used for scalling

xr=  18.30249281554*delka #radiator center -- BO 
rel_dens=13.688487208568  # relative density of radiator --BO


n0 = 0.005      # accelerator density
xi=13.5*delka   # accelerator center
Li= 12.2*delka  # accelerator FWHM
Lr= 3.9*delka   # radiator FWHM
m = 10          # SG steepness



def prof_and_lens(x,r):
    """plasma profile"""
    return n0*(np.exp(-((x-xi)/Li)**(2*m) +rel_dens*np.exp(-((x-xr)/Lr)**(2*m)))
        

Main(
    geometry = "AMcylindrical",
    number_of_AM = 2,
    interpolation_order = 2,
    
    cell_length = [dx,dy],
    grid_length  = [Lx,Ly], 
    
    number_of_patches = [ n_patch, 32], 
    
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
    time_start = Lx*1.15, 
    velocity_x = 0.9997
)

RadiationReaction(
   Niel_computation_method = "table",
   # Radiation parameters
   minimum_chi_continuous = 1.5e-4,
   #minimum_chi_discontinuous = 1e-3,
)

a0=4
tau = 11.25*t0
omega_p = np.sqrt(n0)
lambda_p = 2*np.pi/omega_p
laser_waist = np.sqrt(a0)*lambda_p/np.pi


LaserGaussianAM(
   box_side        = "xmin",
   omega           = 1.,                # normalized units
   a0              = a0,                # normalized units
   focus           = [Lx/2, 0],         # normalized units
   waist           = laser_waist,       # normalized units
   
   time_envelope   = tgaussian(fwhm=tau,center=3*tau, )
)


  


Species(                        #electrons
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


Species(                    #ions
    name = 'ion',
    
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    
    particles_per_cell      = 8,
    regular_number          = [1,1,8],
    mass = 1836.,    # normalized units
    charge = 1., # normalized units
    
    number_density = prof_and_lens,
    
    boundary_conditions = [
        ["remove", "remove"],
        ["reflective", "remove"],
    ],
    time_frozen = Niterations,
    pusher = "borisBTIS3", 
)





####### Diagnostics ######
globalEvery =1000    

DiagProbe(
    every = globalEvery,
    origin = [0., -Ly, 0],
    corners = [ [Lx,-Ly,0.], [0,Ly,0.] ],
    number = [700, 700],
    fields = ['Ex','Ey','Rho_eon',]
)


####### ELECTRON ######

DiagParticleBinning(
    deposited_quantity = "weight",
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["y",    -Ly,    Ly,    800],
             ["py",   -80,    80,    750],
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

DiagParticleBinning(
    deposited_quantity = "weight",	
    every = globalEvery,
    #time_average = 1,
    species = ["eon",],
    axes = [ ["moving_x",    0.,    Lx,    500],
             ["px",    -10,    2500,    500] ]
)

DiagPerformances(
    every = 50
)


####### RADIATION ######

DiagRadiationSpectrum(
    every = globalEvery,
    time_average = globalEvery ,
    species = ["eon"],
    photon_energy_axis = [photon_energy_max/1.e8,photon_energy_max, 2000,'logscale'],
    axes = []
)

