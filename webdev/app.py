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

class CatalystGeometry(Enum):
    Slab = "Slab",
    Cylinder = "Cylinder",
    Sphere = "Sphere",
    Other = "Other"

app = Flask(__name__)
@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")
