import effectiveness_factor as ef
from typing import List, Optional
from enum import Enum
from flask import Flask, request, render_template


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

    @staticmethod
    def from_str(label): 
        if label.lower() == "typei" or label.lower() == "type1":
            return ReactionType.TypeI
        if label.lower() == "typeii" or label.lower() == "type2":
            return ReactionType.TypeII
        if label.lower() == "typeiii" or label.lower() == "type3":
            return ReactionType.TypeIII
        if label.lower() == "typeiv" or label.lower() == "type4":
            return ReactionType.TypeIV
        if label.lower() == "typev" or label.lower() == "type5":
            return ReactionType.TypeV
        if label.lower() == "typevi" or label.lower() == "type6":
            return ReactionType.TypeVI
        if label.lower() == "typevii" or label.lower() == "type7":
            return ReactionType.TypeVII
        raise ValueError(f"Unknown reaction type: {label}")
    
#  // type 2, type 4, type 6 => disable B
# // type 3, type 5, type 6 => disable D
    @staticmethod
    def components(reaction_type):
        if reaction_type in [ReactionType.I, ReactionType.VII]: 
            return ['A', 'B', 'C', 'D']
        if reaction_type in [ReactionType.II, ReactionType.IV]:
            return ['A', 'C', 'D']
        if reaction_type in [ReactionType.III, ReactionType.V]:
            return ['A', 'B', 'C']
        return ['A', 'C']


class CatalystGeometry(Enum):
    Slab = "Slab",
    Cylinder = "Cylinder",
    Sphere = "Sphere",
    Other = "Other"

    @staticmethod
    def from_str(label):
        if label.lower() == "slab":
            return CatalystGeometry.Slab
        if label.lower() == "cylinder":
            return CatalystGeometry.Cylinder
        if label.lower() == "sphere":
            return CatalystGeometry.Sphere
        return CatalystGeometry.Other

def extract_float_default(form_data: dict, key: str, default: float | None = None) -> float | None:
    # Helper function to get a value from the form but avoid parsing errors (ValueError) due to the key existing but having an empty string as its value ('')
    value_str = form_data.get(key, str(default))
    try:
        value = float(value_str)
    except ValueError:
        value = default
    return value

def extract_float_error(form_data: dict, key: str) -> float:
    # Helper function to get a value from the form but raise an error if parsing fails
    value_str = form_data.get(key)
    if value_str is None:
        raise ValueError(f"Missing required field: {key}")
    try:
        value = float(value_str)
    except ValueError:
        raise ValueError(f"Invalid numeric value for field: {key}")
    return value

def extract_data_list(form_data: dict, base_key: str, mixture_components: List) -> Vector:
    # Helper function to extract a list of float values from the form data based on a base key and mixture components
    values = []
    for component in mixture_components:
        key = f"{base_key}{component}"
        value = extract_float_error(form_data, key)
        values.append(value)
    return values

def extract_unifac_groups(form_data: dict) -> Vector:
    # Groups have key in the form "group1-A", "group1-B", etc. 
    # Since component "A" always exists, we can use it to get the total number of UNIFAC structural groups entered by the user by counting how many "groupX-A" keys exist  
    groups_ids = []
    for key in form_data.keys():
        if key.startswith("group") and key.endswith("-A"):
            group_id = key.split("-")[0].split("group")[1]
            groups_ids.append(group_id)
    return groups_ids

def init_matrix(nrows: int, ncols: int) -> Matrix:
    return [[0.0 for _ in range(ncols)] for _ in range(nrows)]

def extract_unifac_params(form_data: dict, mixture_components: List):
    groups = extract_unifac_groups(form_data)
    unifac_nu = init_matrix(len(mixture_components), len(groups))
    unifac_a = init_matrix(len(groups), len(groups))
    unifac_R = []
    unifac_Q = []

    # Initialize empty lists for each component
    for i, group_label in enumerate(groups):
        key_R = f"group{group_label}-Rk"
        key_Q = f"group{group_label}-Qk"
        R_value = extract_float_error(form_data, key_R)
        Q_value = extract_float_error(form_data, key_Q)
        unifac_R.append(R_value)
        unifac_Q.append(Q_value)
        for j, component in enumerate(mixture_components):
            key_nu = f"group{group_label}-{component}"
            nu_value = extract_float_error(form_data, key_nu)
            unifac_nu[j][i] = nu_value
        for j, group_label_inner in enumerate(groups):
            key_a = f"energy{group_label}-{group_label_inner}"
            a_value = extract_float_error(form_data, key_a)
            unifac_a[i][j] = a_value
    return unifac_R, unifac_Q, unifac_nu, unifac_a


app = Flask(__name__)
@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")

# TODO see how to handle errors / return errors with information about what caused the failure...

@app.route("/calculate", methods=["POST"])
def calculate():
    # data: a dictionary containing all input parameters
    data = request.get_json()
    
    # NOTE for debugging
    print(data)


    # Main information: diffusivity computation mode, reaction kinetics, catalyst geometry (strings/enums)
    diffusivity_computation = DiffusivityComputation[data.get("diffusivityType")]
    reaction_type = ReactionType.from_str(data.get("rxType"))
    catalyst_geometry = CatalystGeometry.from_str(data.get("geometryType"))
    components = ReactionType.components(reaction_type)

    # Catalyst properties single value data
    catalyst_porosity = extract_float_default(data, "porosity", 0.5)
    catalyst_characteristic_length = extract_float_default(data, "charlength", 0.1)
    catalyst_density = extract_float_default(data, "particledensity", 600.0)
    catalyst_tortuosity = extract_float_default(data, "tortuosity", None)

    # Reaction properties single value data
    kinetic_constant = extract_float_default(data, "kconst", 1.0)
    equilibrium_constant = extract_float_default(data, "keq", 1.0)
    equilibrium_Ca = extract_float_default(data, "concAeq", 1.0)
    temperature = extract_float_default(data, "temperature", None)

    # Numeric list data 
    concentrations = extract_data_list(data, "conc", components)

    # Breakpoint for different diffusivity computation modes
    if diffusivity_computation == DiffusivityComputation.Precomputed:
        diffusivities = extract_data_list(data, "diffusivity", components)
        result = ef.ef_precomputed_diffusivities(
            reaction_type.value[0],
            catalyst_geometry.value[0], 
            catalyst_characteristic_length, 
            catalyst_porosity,
            catalyst_density,
            catalyst_tortuosity,
            kinetic_constant,
            equilibrium_constant,
            equilibrium_Ca,
            concentrations,
            diffusivities
            )
        return {"result": result}

    # ideal or nonideal, needs Wilke-Chang data
    molecular_weights = extract_data_list(data, "mw", components)
    association_factors = extract_data_list(data, "phi", components)
    molar_volumes_bp = extract_data_list(data, "vbp", components)
    viscosities = extract_data_list(data, "viscosity", components)

    if diffusivity_computation == DiffusivityComputation.Ideal:
        result = ef.ef_diffusivities_ideal(
            reaction_type.value[0],
            catalyst_geometry.value[0], 
            catalyst_characteristic_length, 
            catalyst_porosity, 
            catalyst_density, 
            catalyst_tortuosity, 
            kinetic_constant, 
            equilibrium_constant, 
            equilibrium_Ca, 
            concentrations, 
            molecular_weights, 
            association_factors, 
            molar_volumes_bp, 
            viscosities, 
            temperature
        )
        return {"result": result}
    
    # nonideal, needs UNIFAC data
    unifac_R, unifac_Q, unifac_nu, unifac_a = extract_unifac_params(data, components)

    result = ef.ef_diffusivities_new(
            reaction_type.value[0],
            catalyst_geometry.value[0], 
            catalyst_characteristic_length, 
            catalyst_porosity, 
            catalyst_density, 
            catalyst_tortuosity, 
            kinetic_constant, 
            equilibrium_constant, 
            equilibrium_Ca, 
            concentrations, 
            molecular_weights, 
            association_factors, 
            molar_volumes_bp, 
            viscosities, 
            temperature, 
            unifac_R, 
            unifac_Q, 
            unifac_a, 
            unifac_nu
    )
    return {"result": result}

