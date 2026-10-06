# matrix_models

Codes for various matrix models: Langevin simulations of bosonic gauged
matrix quantum mechanics.

## Contents

| Folder | Model |
|---|---|---|
| [`D4MM_v06.10.26`](D4MM_v06.10.26) | bosonic D=4 matrix model (3 matrices + gauge holonomy) | 

| What it does |
generates configurations by Langevin dynamics and measures 
the Polyakov loop, 
the extent of space and 
the energy on every saved configuration; parallel run script included |

Each folder has its own README with the model, how to build and run, the output format and the known limits.

## Quick start (D4MM)

    cd D4MM_v06.10.26
    python3 run.py list
    python3 run.py run seeds --tier smoke

Needs a C++17 compiler and Python 3.

## Created by

Arpith Kumar, Postdoctoral Researcher, Bielefeld University (10.2026)
