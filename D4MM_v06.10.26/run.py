#!/usr/bin/env python3
"""Run the D=4 matrix model for a whole set of parameters (temperatures, lattice sizes, seeds), in parallel.

    python3 run.py list                                   show the parameter sets
    python3 run.py run scan_L20 --tier smoke              compile and run every point of a set (tier = how long)
    python3 run.py run scan_L20 --tier standard --jobs 4  four runs at the same time
    python3 run.py status scan_L20 --tier standard        which runs are finished

The program has N, L, the temperature and lambda as compile-time constants (src/model.h), so for every run this script copies src/
into the run folder, sets those constants, compiles, writes params.txt and starts the program there. Nothing in src/ is changed.

Result of a run: runs/<set>/<tier>/<name>/  with  meas/extent.txt  meas/polyakov.txt  meas/energy.txt  (one line per measured
configuration), run.log, params.txt, run.json (all settings, timing, status) and src/ (the exact code that was compiled).
Only the standard library of Python 3 and a C++17 compiler (g++ or clang++) are needed.

How the file is organised: settings (BASE_STEP, TIERS, PARAM_SETS) -> small helpers (step_for, build_runs, set_constants)
-> run_one (compile and run ONE point) -> main (command line, then run all points). The physics lives in src/, not here.
"""

import argparse                      # reads the command line
import concurrent.futures            # runs several simulations at the same time
import json                          # run.json: a machine-readable record of every run
import os                            # number of cpu cores
import re                            # search-and-replace of the constants in model.h
import shutil                        # copy / delete folders, find the compiler
import subprocess                    # start the compiler and the simulation
import sys                           # exit with a message
import time                          # timing of a run
from pathlib import Path             # file paths

HERE = Path(__file__).resolve().parent            # the folder that contains this script (and src/)

# ---------------------------------------------------------------- settings
BASE_STEP = 2e-4                                  # FIXED step of a run: one value per run, same at every step. step_for() lowers it for high T or fine lattices
ADAPT_FACTOR = 3.0                                # ADAPTIVE part, acts inside a run: a step is shortened if its largest force exceeds 3 x the running mean (0 = off)
STEPS_PER_SECOND_N4_L20 = 22000                   # speed of the program for N=4, L=20 on a laptop; only used for the time estimate

# How long a run is. therm = steps before the first measurement, sweeps = steps with measurements, gap = measure every `gap` steps.
TIERS = {
    "smoke": dict(therm=2000, sweeps=20000, gap=20),            # a few seconds: only to see that everything works
    "standard": dict(therm=50000, sweeps=1600000, gap=100),     # about 100 s for N=4, L=20
    "long": dict(therm=100000, sweeps=10000000, gap=100),       # about 6 times longer than standard
}

# The parameter sets. A point is a dict with N, L, T and optionally seed (default 1) or step (default step_for).
# Add your own sets here; `what` is the line printed by `python3 run.py list`.
PARAM_SETS = {
    "scan_L20": dict(
        what="N=4, L=20: 17 temperatures T = 0.2 ... 2.5, denser around T = 1.1 (the transition region)",
        points=[dict(N=4, L=20, T=t) for t in (0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5)],
    ),
    "scan_L20_L40": dict(
        what="N=4, L=20 and L=40: nine temperatures T = 0.2 ... 4 each (compares two lattice spacings)",
        points=[dict(N=4, L=length, T=t) for length in (20, 40) for t in (0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 4.0)],
    ),
    "seeds": dict(
        what="N=4, L=20, T=2.5 with four random seeds: the spread between seeds shows the true statistical error",
        points=[dict(N=4, L=20, T=2.5, seed=s) for s in (1, 2, 3, 4)],
    ),
}


# ---------------------------------------------------------------- small helpers
def step_for(temperature, length):
    """Langevin step: BASE_STEP, but at most 2 percent of the lattice spacing a = 1/(T L), because the force grows like 1/a."""
    spacing = 1.0 / (temperature * length)        # lattice spacing a = beta / L = 1 / (T L)
    return min(BASE_STEP, 0.02 * spacing)         # the smaller of the two


def build_runs(set_name, tier, outdir):
    """Full description of every run of a set: parameters, the statistics of the tier, and the folder of the run."""
    runs = []                                                          # one dict per run
    for point in PARAM_SETS[set_name]["points"]:                       # every point of the set
        run = {"seed": 1, **TIERS[tier], **point}                      # start with seed 1 and the tier; the point's own entries win
        run.setdefault("step", step_for(run["T"], run["L"]))           # no step given in the point: choose it by the rule above
        name = f"N{run['N']}_L{run['L']}_T{run['T']:g}"                # folder name from N, L, T ...
        if "seed" in point:                                            # ... plus the seed, if this point set one
            name += f"_seed{run['seed']}"
        run["name"] = name
        run["folder"] = str(Path(outdir) / set_name / tier / name)     # runs/<set>/<tier>/<name>
        runs.append(run)
    return runs


def set_constants(header_text, run):
    """Put N, L and T of this run into a copy of model.h. Each of the three lines must be found exactly once (lambda stays 1)."""
    replacements = (                                                   # (what to look for, what to put instead)
        (r"const int L = \d+;", f"const int L = {run['L']};"),
        (r"const int NCOLOR = \d+;", f"const int NCOLOR = {run['N']};"),
        (r"const double TEMP = [0-9.eE+-]+;", f"const double TEMP = {run['T']!r};"),
    )
    for pattern, replacement in replacements:
        header_text, count = re.subn(pattern, replacement, header_text)    # replace and count the replacements
        if count != 1:                                                     # the line was missing or appears twice: model.h was edited
            raise RuntimeError(f"model.h: expected exactly one line matching {pattern!r}, found {count}")
    return header_text


def count_lines(path):
    """Number of lines of a text file, 0 if it does not exist (used to check that a run produced all its measurements)."""
    if not path.exists():
        return 0
    with path.open() as handle:
        return sum(1 for _ in handle)


# ---------------------------------------------------------------- one run
def run_one(run, force):
    """Compile and run ONE point of a set; returns a one-line message for the screen."""
    folder = Path(run["folder"])
    record = folder / "run.json"                                       # the record of this run

    # --- skip a finished run (unless --force) -------------------------------------------------------------
    if record.exists() and json.loads(record.read_text()).get("status") == "done" and not force:
        return f"{run['name']}: already done, skipped"
    if folder.exists():
        shutil.rmtree(folder)                                          # start clean: an earlier attempt was partial, or --force

    # --- copy the sources, with N, L, T of this run written into model.h ------------------------------------
    (folder / "src").mkdir(parents=True)
    for path in (HERE / "src").iterdir():                              # every file of src/
        text = path.read_text()
        if path.name == "model.h":
            text = set_constants(text, run)                            # only model.h gets the run's constants
        (folder / "src" / path.name).write_text(text)

    # --- params.txt: step, thermalisation, production, measure_every, adaptive factor, seed (see README) --------
    parameters = [repr(run["step"]), run["therm"], run["sweeps"], run["gap"], ADAPT_FACTOR, run["seed"]]
    (folder / "params.txt").write_text("\n".join(str(value) for value in parameters) + "\n")

    # --- record that the run has started (bookkeeping, not important for the physics) -------------------------
    started = time.time()
    status = dict(run, status="running", started=time.strftime("%Y-%m-%d %H:%M:%S"))
    record.write_text(json.dumps(status, indent=2) + "\n")

    # --- compile ---------------------------------------------------------------------------------------------
    compiler = shutil.which("g++") or shutil.which("clang++")          # first compiler found
    sources = [str(p) for p in sorted((folder / "src").glob("*.cpp"))] # all .cpp files of the copy
    compile_command = [compiler, "-O2", "-std=c++17", *sources, "-o", str(folder / "d4mm")]    # -O2: optimise; -std=c++17: needed for <filesystem>
    if compiler is None or subprocess.run(compile_command, capture_output=True).returncode != 0:
        status["status"] = "compile failed"
        record.write_text(json.dumps(status, indent=2) + "\n")
        return f"{run['name']}: COMPILE FAILED"

    # --- run the simulation in the run folder (it reads params.txt there and writes meas/ there) ----------------
    with (folder / "run.log").open("w") as log:                        # everything the program prints goes into run.log
        code = subprocess.run([str(folder / "d4mm")], cwd=folder, stdout=log, stderr=subprocess.STDOUT).returncode

    # --- check that the run produced all its measurements -----------------------------------------------------
    expected = run["sweeps"] // run["gap"]                             # number of measurements a complete run has
    rows = {name: count_lines(folder / "meas" / f"{name}.txt") for name in ("extent", "polyakov", "energy")}
    complete = code == 0 and all(count == expected for count in rows.values())    # exit code 0 and the right number of lines in all files

    # --- how often the adaptive step shortened a step (printed by the program at the end of run.log) ---------------
    report = re.search(r"steps shortened (\d+) of (\d+)", (folder / "run.log").read_text())
    shortened = ""
    if report:
        shortened = f", {100 * int(report.group(1)) / int(report.group(2)):.2f}% of steps shortened"

    # --- final record and message (bookkeeping) -----------------------------------------------------------------
    status.update(status="done" if complete else "failed", seconds=round(time.time() - started, 1), rows=rows, expected_rows=expected)
    record.write_text(json.dumps(status, indent=2) + "\n")
    outcome = "done" if complete else "FAILED"
    return f"{run['name']}: {outcome} in {status['seconds']} s ({rows['polyakov']} measurements{shortened})"


# ---------------------------------------------------------------- command line
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=("list", "run", "status"))                       # what to do
    parser.add_argument("set", nargs="?", choices=sorted(PARAM_SETS))                      # which parameter set (not needed for `list`)
    parser.add_argument("--tier", choices=sorted(TIERS), default="standard")               # how long each run is
    parser.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1), help="runs at the same time (one core each)")
    parser.add_argument("--outdir", default=str(HERE / "runs"), help="where the run folders are created")
    parser.add_argument("--force", action="store_true", help="redo runs that are already finished")
    parser.add_argument("--dry-run", action="store_true", help="show the runs and a time estimate, start nothing")
    args = parser.parse_args()

    # `list`: print the parameter sets and stop
    if args.action == "list":
        for name, info in PARAM_SETS.items():
            print(f"{name:13s} {len(info['points']):2d} runs  {info['what']}")
        return
    if args.set is None:
        parser.error("name a parameter set (see `list`)")

    # the runs of the chosen set, and a time estimate (bookkeeping: the program speed scales roughly like 1/L)
    runs = build_runs(args.set, args.tier, args.outdir)
    tier = TIERS[args.tier]
    print(f"{args.set}, tier {args.tier} (therm {tier['therm']}, sweeps {tier['sweeps']}, gap {tier['gap']})")
    seconds = [(r["therm"] + r["sweeps"]) / (STEPS_PER_SECOND_N4_L20 * 20 / r["L"]) for r in runs]     # steps / (steps per second)
    hours = sum(seconds) / 3600
    print(f"  {len(runs)} runs, about {hours:.2f} h on one core, about {hours / args.jobs:.2f} h with {args.jobs} jobs")

    # `status`: show which runs are finished and stop
    if args.action == "status":
        for run in runs:
            record = Path(run["folder"]) / "run.json"
            state = json.loads(record.read_text()).get("status") if record.exists() else "not started"
            print(f"  {run['name']:30s} {state}")
        return

    # `--dry-run`: show every run with its step and seed, start nothing
    if args.dry_run:
        for run in runs:
            print(f"  {run['name']:30s} T={run['T']:<5g} step={run['step']:.2e} seed={run['seed']}")
        return

    # `run`: start all runs, `jobs` of them at a time
    if shutil.which("g++") is None and shutil.which("clang++") is None:
        sys.exit("no C++ compiler found: install g++ or clang++ first")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:      # each thread waits for one simulation process
        for message in pool.map(lambda run: run_one(run, args.force), runs):        # results are printed in the order of the runs
            print(message, flush=True)


if __name__ == "__main__":
    main()
