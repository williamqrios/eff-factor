import effectiveness_factor as ef
from typing import List, Optional, Tuple
from enum import Enum
from flask import Flask, request, render_template, jsonify

class ValidationError(Exception):
    def __init__(self, message):
        self.message = message

def validate_parse_float(float_as_str: Optional[str], key: str) -> float:
    if float_as_str is None or float_as_str.strip() == "": 
        raise ValidationError(f"{key} cannot be empty.")

    try:
        value = float(float_as_str.strip())
    except ValueError:
        raise ValidationError(f"{key} has an invalid numeric value: {float_as_str}")

    return value 


def validate_float(value: float, key: str, zero_allowed: bool = True, negative_allowed: bool = True, limits: Optional[Tuple[float, float]] = None) -> float:
    
    if not zero_allowed and value == 0.0:
        raise ValidationError(f"{key} cannot be zero.")

    if not negative_allowed and value < 0.0:
        raise ValidationError(f"{key} cannot be negative.")

    if limits is not None:
        lower_limit, upper_limit = limits
        if not (lower_limit < value < upper_limit):
            raise ValidationError(f"{key} must be higher than {lower_limit} and lower than {upper_limit}.")

    return value

def validate_result(value):
    if value is None or (isinstance(value, float) and (value != value)):  # Check for NaN
        raise ValidationError("The calculation resulted in an invalid number (NaN). Please check your input parameters.")
    return value


type Vector = List[float]
type Matrix = List[Vector]


class DiffusivityComputation(Enum):
    Precomputed = "precomputed",
    Ideal = "ideal"
    NonIdeal = "nonideal",

    @staticmethod
    def from_str(label: str):
        if label.lower() == "precomputed":
            return DiffusivityComputation.Precomputed
        if label.lower() == "ideal":
            return DiffusivityComputation.Ideal
        if label.lower() == "nonideal":
            return DiffusivityComputation.NonIdeal
        raise ValueError(f"Unknown diffusivity computation type: {label}")

class ReactionType(Enum):   
    TypeI = "TypeI",
    TypeII = "TypeII",
    TypeIII = "TypeIII",
    TypeIV = "TypeIV",
    TypeV = "TypeV",
    TypeVI = "TypeVI",
    TypeVII = "TypeVII",

    @staticmethod
    def from_str(label: str): 
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
    
    @staticmethod
    def components(reaction_type: "ReactionType") -> List[str]:
        if reaction_type in [ReactionType.TypeI, ReactionType.TypeVII]: 
            return ['A', 'B', 'C', 'D']
        if reaction_type in [ReactionType.TypeII, ReactionType.TypeIV]:
            return ['A', 'C', 'D']
        if reaction_type in [ReactionType.TypeIII, ReactionType.TypeV]:
            return ['A', 'B', 'C']
        return ['A', 'C']


class CatalystGeometry(Enum):
    Slab = "Slab",
    Cylinder = "Cylinder",
    Sphere = "Sphere",
    Other = "Other"

    @staticmethod
    def from_str(label: str):
        if label.lower() == "slab":
            return CatalystGeometry.Slab
        if label.lower() == "cylinder":
            return CatalystGeometry.Cylinder
        if label.lower() == "sphere":
            return CatalystGeometry.Sphere
        return CatalystGeometry.Other

def extract_float_default(form_data: dict, key: str, default: Optional[float] = None, zero_allowed: bool = True, negative_allowed: bool = True, limits: Optional[Tuple[float, float]] = None) -> Optional[float]:
    # Helper function to get a value from the form but avoid parsing errors (ValueError) due to the key existing but having an empty string as its value ('')
    value_str = form_data.get(key, str(default))
    try:
        value = float(value_str)
    except ValueError:
        value = default
    else:
        # Runs if try block succeeds
        value = validate_float(value, key, zero_allowed=zero_allowed, negative_allowed=negative_allowed, limits=limits)
    
    return value

def extract_float_error(form_data: dict, key: str, zero_allowed: bool = True, negative_allowed: bool = True, limits: Optional[Tuple[float, float]] = None) -> float:
    # Helper function to get a value from the form but raise an error if parsing fails
    value_str = form_data.get(key)
    value_float = validate_parse_float(value_str, key)
    value = validate_float(value_float, key, zero_allowed=zero_allowed, negative_allowed=negative_allowed, limits=limits)
    return value

def extract_data_list(form_data: dict, base_key: str, mixture_components: List, zero_allowed: bool = True, negative_allowed: bool = True, limits: Optional[Tuple[float, float]] = None) -> Vector:
    # Helper function to extract a list of float values from the form data based on a base key and mixture components
    values = []
    for component in mixture_components:
        key = f"{base_key}{component}"
        value = extract_float_error(form_data, key, zero_allowed=zero_allowed, negative_allowed=negative_allowed, limits=limits)
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

def extract_unifac_params(form_data: dict, mixture_components: List[str]) -> tuple:
    groups = extract_unifac_groups(form_data)
    unifac_nu = init_matrix(len(mixture_components), len(groups))
    unifac_a = init_matrix(len(groups), len(groups))
    unifac_R = []
    unifac_Q = []

    # Initialize empty lists for each component
    for i, group_label in enumerate(groups):
        key_R = f"group{group_label}-Rk"
        key_Q = f"group{group_label}-Qk"
        R_value = extract_float_error(form_data, key_R, zero_allowed=True, negative_allowed=False, limits=None)
        Q_value = extract_float_error(form_data, key_Q, zero_allowed=True, negative_allowed=False, limits=None)
        unifac_R.append(R_value)
        unifac_Q.append(Q_value)
        for j, component in enumerate(mixture_components):
            key_nu = f"group{group_label}-{component}"
            nu_value = extract_float_error(form_data, key_nu, zero_allowed=True, negative_allowed=False, limits=None)
            unifac_nu[j][i] = nu_value
        for j, group_label_inner in enumerate(groups):
            key_a = f"energy{group_label}-{group_label_inner}"
            a_value = extract_float_default(form_data, key_a, 0.0, zero_allowed=True, negative_allowed=True, limits=None)
            unifac_a[i][j] = a_value
    return unifac_R, unifac_Q, unifac_nu, unifac_a


app = Flask(__name__)
@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")

@app.route("/calculate", methods=["POST"])
def calculate():
    try: 
        # data: a dictionary containing all input parameters
        data = request.get_json()
        
        # NOTE for debugging
        # print(data)

        # Main information: diffusivity computation mode, reaction kinetics, catalyst geometry (strings/enums)

        diffusivity_computation = DiffusivityComputation.from_str(data.get("diffusivityType"))
        reaction_type = ReactionType.from_str(data.get("rxType"))
        catalyst_geometry = CatalystGeometry.from_str(data.get("geometryType"))
        components = ReactionType.components(reaction_type)


        # Catalyst properties single value data
        catalyst_porosity = extract_float_error(data, "porosity", zero_allowed=False, negative_allowed=False, limits=(0.0, 1.0))
        catalyst_characteristic_length = extract_float_error(data, "charlength", zero_allowed=False, negative_allowed=False, limits=None)
        catalyst_density = extract_float_error(data, "particledensity", zero_allowed=False, negative_allowed=False, limits=None)
        catalyst_tortuosity = extract_float_default(data, "tortuosity", None, zero_allowed=False, negative_allowed=False)

        # Reaction properties single value data
        kinetic_constant = extract_float_error(data, "kconst", zero_allowed=False, negative_allowed=False, limits=None)
        equilibrium_constant = extract_float_error(data, "keq", zero_allowed=False, negative_allowed=False, limits=None)
        equilibrium_Ca = extract_float_error(data, "concAeq", zero_allowed=False, negative_allowed=False, limits=None)
        temperature = extract_float_error(data, "temperature", zero_allowed=False, negative_allowed=False, limits=None)

        # Numeric list data 
        concentrations = extract_data_list(data, "conc", components, zero_allowed=True, negative_allowed=False, limits=None)

        # Breakpoint for different diffusivity computation modes
        if diffusivity_computation == DiffusivityComputation.Precomputed:
            diffusivities = extract_data_list(data, "diffusivity", components, zero_allowed=False, negative_allowed=False, limits=None)
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
            result = validate_result(result)
            return jsonify({"result": result}), 200

        # ideal or nonideal, needs Wilke-Chang data
        molecular_weights = extract_data_list(data, "mw", components, zero_allowed=False, negative_allowed=False, limits=None)
        association_factors = extract_data_list(data, "phi", components, zero_allowed=False, negative_allowed=False, limits=None)
        molar_volumes_bp = extract_data_list(data, "vbp", components, zero_allowed=False, negative_allowed=False, limits=None)
        viscosities = extract_data_list(data, "viscosity", components, zero_allowed=False, negative_allowed=False, limits=None)

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
            result = validate_result(result)
            return jsonify({"result": result}), 200
        
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
        result = validate_result(result)
        return jsonify({"result": result}), 200
    # Errors raised by validation functions in this module
    except ValidationError as e:
        return jsonify({"error": str(e.message)}), 400
    # Errors raised by calling the Rust functions with invalid parameters not caught by the above validation 
    except ValueError as e: 
        return jsonify({"error": str(e)}), 400
    # Any other error due to server issues
    except Exception as e:
        return jsonify({"error": f"An unexpected error occurred: {str(e)}"}), 500

