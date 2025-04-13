import effectiveness_factor as ef
from typing import List, Optional
from enum import Enum

type Vector = List[float]
type Matrix = List[Vector]


class DiffusivityComputation(Enum):
    Precomputed = "precomputed",
    Ideal = "ideal"
    NonIdeal = "nonideal",

class ReactionType(Enum):   
    TypeI = "TypeI",
    TypeII = "TypeII",
    TypeIII = "TypeIII",
    TypeIV = "TypeIV",
    TypeV = "TypeV",
    TypeVI = "TypeVI",
    TypeVII = "TypeVII",

class CatalystGeometry(Enum):
    Slab = "Slab",
    Cylinder = "Cylinder",
    Sphere = "Sphere",
    Other = "Other"


def effectiveness_factor_wrapper(reaction_type: ReactionType, catalyst_geometry: CatalystGeometry, catalyst_characteristic_length: float, catalyst_porosity: float, catalyst_density: float, kinetic_constant: float, equilibrium_constant: float, equilibrium_Ca: float, concentrations: Vector, diffusivity_computation: DiffusivityComputation, catalyst_tortuosity: Optional[float] = None, diffusivities: Optional[Vector] = None, molecular_weights: Optional[Vector] = None, association_factors: Optional[Vector] = None, molar_volumes_bp: Optional[Vector] = None, viscosities: Optional[Vector] = None, temperature: Optional[float] = None, unifac_R: Optional[Vector] = None, unifac_Q: Optional[Vector] = None, unifac_a: Optional[Matrix] = None, unifac_nu: Optional[Matrix] = None) -> float:
    if diffusivity_computation == DiffusivityComputation.Precomputed: 
        if diffusivities is None:
            raise ValueError("Diffusivities must be provided for precomputed diffusivity calculation.")
        return ef.ef_precomputed_diffusivities(reaction_type.value[0], catalyst_geometry.value[0], catalyst_characteristic_length, catalyst_porosity, catalyst_density, catalyst_tortuosity, kinetic_constant, equilibrium_constant, equilibrium_Ca, concentrations, diffusivities)
    if diffusivity_computation == DiffusivityComputation.Ideal:
        if any(map(lambda x: x is None, [molecular_weights, association_factors, molar_volumes_bp, viscosities, temperature])):
            raise ValueError("All parameters must be provided for ideal diffusivity calculation.")
        return ef.ef_diffusivities_ideal(reaction_type.value[0], catalyst_geometry.value[0], catalyst_characteristic_length, catalyst_porosity, catalyst_density, catalyst_tortuosity, kinetic_constant, equilibrium_constant, equilibrium_Ca, concentrations, molecular_weights, association_factors, molar_volumes_bp, viscosities, temperature)
    if diffusivity_computation == DiffusivityComputation.NonIdeal:
        if any(map(lambda x: x is None, [molecular_weights, association_factors, molar_volumes_bp, viscosities, temperature, unifac_Q, unifac_a, unifac_nu])):
            raise ValueError("All parameters must be provided for non-ideal diffusivity calculation.")
        return ef.ef_diffusivities_new(reaction_type.value[0], catalyst_geometry.value[0], catalyst_characteristic_length, catalyst_porosity, catalyst_density, catalyst_tortuosity, kinetic_constant, equilibrium_constant, equilibrium_Ca, concentrations, molecular_weights, association_factors, molar_volumes_bp, viscosities, temperature, unifac_R, unifac_Q, unifac_a, unifac_nu)



if __name__ == "__main__":
    # Example usage - precomputed diffusivities
    reaction_type = ReactionType.TypeI
    catalyst_geometry = CatalystGeometry.Sphere
    catalyst_characteristic_length = 0.744/2.0 * 1.0e-2 # dm
    catalyst_porosity = 0.4886 
    catalyst_tortuosity = 1.3 
    catalyst_density = 600.0 # g/dm3 
    kinetic_constant = 4.5747e-5 # dm6 mol-1 min-1 gcat-1
    equilibrium_constant = 2.6700 
    equilibrium_Ca = 3.2450 # mol/dm3 
    concentrations = [8.3747, 8.3747, 0.0, 0.0] # mol/dm3
    diffusivities = [0.3172e-4, 0.2203e-4, 0.2009e-4, 0.3684e-4] # dm2/min
    result = effectiveness_factor_wrapper(reaction_type, catalyst_geometry, catalyst_characteristic_length, catalyst_porosity, catalyst_density, kinetic_constant, equilibrium_constant, equilibrium_Ca, concentrations, DiffusivityComputation.Precomputed, catalyst_tortuosity=catalyst_tortuosity, diffusivities=diffusivities)
    print(f"eta = {result:.3f} (pre-computed)")

    # Example usage - ideal diffusivity
    # The diffusivities are in units of cm2/s, so convert the other parameters to the same units. 
    catalyst_characteristic_length = 0.744/2.0 * 1.0e-1 # cm
    catalyst_porosity = 0.4886 
    catalyst_tortuosity = 1.3 
    catalyst_density = 0.6 # g/cm3
    kinetic_constant = 4.5747e-5 * 1e6/60 # cm6 mol-1 s-1 gcat-1
    equilibrium_constant = 2.6700 
    # viscosities in cP
    viscosities = [0.5728194655847797,
                   0.4408845116586904,
                   0.2503508526367409,
                   0.3649544959435543] 
    association_factors = phi = [1, 1.5, 1, 2.6];
    molecular_weights = [60.0520, 46.0684, 88.1051, 18.0153]
    molar_volumes_bp = [66.0, 60.85, 106.3, 18.8]
    temperature = 273.15 + 78.0 # K
    equilibrium_Ca = 3.2450e-3 # mol/cm3 
    concentrations = [3.2482e-3,    3.2482e-3,    5.2844e-3,    5.2844e-3]
    result = effectiveness_factor_wrapper(reaction_type, catalyst_geometry, catalyst_characteristic_length, catalyst_porosity, catalyst_density, kinetic_constant, equilibrium_constant, equilibrium_Ca, concentrations, DiffusivityComputation.Ideal, catalyst_tortuosity=catalyst_tortuosity, molecular_weights=molecular_weights, association_factors=association_factors, molar_volumes_bp=molar_volumes_bp, viscosities=viscosities, temperature=temperature)
    print(f"eta = {result:.3f} (ideal)")

    # Example usage - non-ideal diffusivity
    # Using same inputs as before, but specifying UNIFAC parameters.
    R = [0.9011, 0.6744, 1., 0.92, 1.9031, 1.3013]
    Q = [0.848, 0.540, 1.2, 1.4, 1.728, 1.224]
    a = [
        [0.0, 0., 986.5, 1318., 232.1, 663.5],
        [0., 0., 986.5, 1318., 232.1, 663.5], 
        [156.4, 156.4, 0., 353.5, 101.1, 199.],
        [300., 300., -229.1, 0., 72.87, -14.09],
        [114.8, 114.8, 245.4, 200.8, 0., 660.2],
        [315.3, 315.3, -151., -66.17, -256.3, 0.]
    ]
    nu = [
        [1., 0., 0., 0., 0., 1.],
        [1., 1., 1., 0., 0., 0.],
        [1., 1., 0., 0., 1., 0.],
        [0., 0., 0., 1., 0., 0.]
    ]

    result = effectiveness_factor_wrapper(reaction_type, catalyst_geometry, catalyst_characteristic_length, catalyst_porosity, catalyst_density, kinetic_constant, equilibrium_constant, equilibrium_Ca, concentrations, DiffusivityComputation.NonIdeal, catalyst_tortuosity=catalyst_tortuosity, molecular_weights=molecular_weights, association_factors=association_factors, molar_volumes_bp=molar_volumes_bp, viscosities=viscosities, temperature=temperature, unifac_R=R, unifac_Q=Q, unifac_a=a, unifac_nu=nu)
    print(f"eta = {result:.3f} (non-ideal)")