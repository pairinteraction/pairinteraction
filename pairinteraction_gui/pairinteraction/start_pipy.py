#!/usr/bin/env python3
"""The main function will be called by the pairinteraction GUI to start the calculation.
"""
import argparse
import concurrent.futures
import itertools
import json
import os
import pickle
import sys
from functools import partial

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import pipy  # noqa


def main(args):
    # Load params from conf.json file
    with open(args.path_config) as f:
        conf = json.load(f)

    # convert conf to settings and save as settings.json
    settings = conf_to_settings(conf, args.path_cache)
    with open(args.path_config.replace("conf.json", "settings.json"), "w") as f:
        json.dump(settings, f, indent=4)

    # start the calculation
    config = settings["config"]
    scriptoptions = settings["scriptoptions"]
    output(f"{'>>TYP':5}{settings['button_id']:7}")
    output(f"{'>>TOT':5}{scriptoptions['numBlocks']*scriptoptions['listoptions']['steps']:7}")
    if config["nAtoms"] == 2:
        for bn, syms in enumerate(scriptoptions["symmetries_list"]):
            config.update(syms)
            settings["blocknumber"] = bn
            do_simulations(settings)
    else:
        settings["blocknumber"] = 1
        do_simulations(settings)
    output(f"{'>>END':5}")


def do_simulations(settings, pass_atom="direct"):
    # TODO decide, wether pass_atom direct or path is better (faster, memory efficient,...)
    param_list = get_param_list(settings)
    ip_list = list(range(len(param_list)))

    atom = pipy.atom_from_config(settings["config"])
    atom.system.buildHamiltonian()

    if atom.nAtoms == 2:
        allQunumbers = [[*qns[:, 0], *qns[:, 1]] for qns in np.array(atom.allQunumbers)]
    else:
        allQunumbers = atom.allQunumbers
    _name = f"{'one' if atom.nAtoms == 1 else 'two'}_{atom.config.toHash()}"
    path_basis = f"/tmp/basis_{_name}_blocknumber_{settings['blocknumber']}.csv"
    basis = np.insert(allQunumbers, 0, np.arange(len(allQunumbers)), axis=1)
    np.savetxt(path_basis, basis, delimiter="\t", fmt=["%d"] + ["%d", "%d", "%.1f", "%.1f"] * atom.nAtoms)
    output(f"{'>>BAS':5}{len(basis):7}")
    output(f"{'>>STA':5} {path_basis} ")

    if pass_atom == "direct":
        config = {"atom": atom}
    elif pass_atom == "path":
        atom_path = path_basis.replace("basis_", "atom_").replace(".csv", ".pkl")
        with open(atom_path, "wb") as file:
            pickle.dump(atom, file)
        config = {"atom_path": atom_path}
    else:
        config = settings["config"]

    p_one_run = partial(one_run, config, param_list)

    num_pr = settings.get("runtimeoptions", {}).get("NUM_PROCESSES", 1)
    num_pr = os.cpu_count() if num_pr in [0, -1] else num_pr

    if num_pr > 1 and not atom.config.isReal():
        raise NotImplementedError(
            "Parallelization for complex calculations is not implemented yet, since there is some bug occuring."
        )

    if num_pr > 1:
        with concurrent.futures.ProcessPoolExecutor(num_pr) as executor:
            futures = [executor.submit(p_one_run, i) for i in ip_list]
            for future in concurrent.futures.as_completed(futures):
                print_completed(settings, future.result())
    else:
        for i in ip_list:
            result = p_one_run(i)
            print_completed(settings, result)

    if "atom_path" in config:
        os.remove(config["atom_path"])


def one_run(config, param_list, ip):
    # get default Atom from config (either given, load it or newly construct it)
    if "atom" in config:
        atom = config["atom"]
    elif "atom_path" in config:
        with open(config["atom_path"], "rb") as f:
            atom = pickle.load(f)
    else:
        atom = pipy.atom_from_config(config)
    config = atom.config

    atom.updateFromParams(param_list[ip])
    dimension = len(atom.basisQunumbers)

    real_complex = "real" if config.isReal() else "complex"
    pathCacheMatrix = config.pathCache() + f"cache_matrix_{real_complex}_new/"
    os.makedirs(pathCacheMatrix, exist_ok=True)
    _name = f"{'one' if atom.nAtoms == 1 else 'two'}_{config.toHash()}"
    filename = pathCacheMatrix + _name + ".pkl"

    if not os.path.exists(filename):
        output(f"Calculating {ip}, {dimension}x{dimension}")
        atom.calcEnergies()
        energies0 = config.getEnergiesPair() if atom.nAtoms == 2 else config.getEnergiesSingle()
        data = {
            "energies": atom.energies - np.mean(energies0),
            "basis": atom.vectors,
            "params": config.toOutDict(),
        }
        with open(filename, "wb") as f:
            pickle.dump(data, f)
        with open(filename.replace(".pkl", ".json"), "w") as f:
            json.dump(data["params"], f, indent=4)

    result = {
        "ip": ip,
        "dimension": dimension,
        "filename": filename,
    }
    return result


def get_param_list(settings):
    scriptoptions = settings.setdefault("scriptoptions", {})
    if "param_list" in scriptoptions:
        return scriptoptions["param_list"]
    elif "listoptions" in scriptoptions and "steps" in scriptoptions["listoptions"]:
        listoptions = scriptoptions["listoptions"]
        steps = listoptions["steps"]
        k_lists = {}
        for k in ["Bx", "By", "Bz", "Ex", "Ey", "Ez", "distance"]:
            if listoptions.get("min" + k) is not None and listoptions.get("max" + k, None) is not None:
                k_lists[k] = np.linspace(listoptions["min" + k], listoptions["max" + k], steps)
        param_list = [{k: v[i] for k, v in k_lists.items()} for i in range(steps)]
        scriptoptions["param_list"] = param_list
    else:
        raise NotImplementedError("TODO: create param list from individual lists, including symmetries")
    return param_list


def print_completed(settings, result):
    output(f"{'>>DIM':5}{result['dimension']:7}")
    numBlocks = settings["scriptoptions"]["numBlocks"]
    totalstep = 1 + result["ip"] * numBlocks + settings["blocknumber"]
    output(
        f"{'>>OUT':5}{totalstep:7}{result['ip']:7}{numBlocks:7}" f"{settings['blocknumber']:7} {result['filename']} "
    )


def conf_to_settings(conf, pathCache):
    config = {}
    for k, v in conf.items():
        if (
            k.startswith("delta")
            or any(k == q + i for q in ["species", "n", "l", "j", "m"] for i in ["", "1", "2"])
            or k in ["diamagnetism", "samebasis", "exponent"]
        ):
            config[k] = v

    config["method"] = ["WHITTAKER", "NUMEROV"][conf["missingCalc"]]
    config["pathCache"] = pathCache
    if conf["button_id"] == 2:
        config["pair.conserveMomenta"] = conf.get("conserveM", False)
        config["pair.useQunumbersSinglePair"] = True
        config["nAtoms"] = 2
    elif conf["button_id"] in [0, 1]:
        config["nAtoms"] = 1
        for k in ["species", "n", "l", "j", "m"]:
            config[k] = config.pop(k + "1", config.pop(k + "2", None))
    else:  # 3
        config["nAtoms"] = 1
    config.update({k: conf.get("min" + k, 0) for k in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]})

    listoptions = {
        "onlySameSector": conf.get("sametrafo", None),
        "steps": conf.get("steps", None),
        **{mm + "distance": conf.get(mm + "R", None) for mm in ["min", "max"]},
        **{mm + key: conf.get(mm + key, None) for key in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"] for mm in ["min", "max"]},
    }
    no_steps = all(
        listoptions["min" + key] == listoptions["max" + key] for key in ["Ex", "Ey", "Ez", "Bx", "By", "Bz", "distance"]
    )
    if no_steps:
        listoptions["steps"] = 1
    config["isReal"] = all([listoptions.get(k, 0) == 0 for k in ["minBy", "maxBy", "minEy", "maxEy"]])

    if config["nAtoms"] == 2:
        symmetries = {}
        for s in ["inv", "per", "ref"]:
            v = [eo for eo in ["O", "E"] if conf.get(s + eo, False)]
            symmetries[s] = v if len(v) > 0 else [None]
        symmetries_list = []
        for inv, per, ref in itertools.product(*[symmetries[k] for k in ["inv", "per", "ref"]]):
            if listoptions["onlySameSector"]:
                # In case of inversion and permutation symmetry: the orbital parity is conserved and we
                # just use the blocks with the same parity as the inital state
                L = config["l1"] + config["l2"]
                if (f"{inv}{per}" in ["EO", "OE"] and L % 2 == 0) or (f"{inv}{per}" in ["EE", "OO"] and L % 2 == 1):
                    continue
                # TODO add more symmetry conditions here, see HamiltonianTwo.cpp
            symmetries_list.append({"inversion": inv, "permutation": per, "reflection": ref})
    else:
        symmetries_list = [{k: False for k in ["inversion", "permutation", "reflection"]}]

    runtimeoptions = {"NUM_PROCESSES": int(os.environ.get("NUM_PROCESSES", 0))}

    # touch unimportant options
    for k in ["dd", "dq", "qq", "missingWhittaker", "precision", "zerotheta"]:
        conf.get(k, None)

    settings = {
        "config": config,
        "scriptoptions": {
            "listoptions": listoptions,
            "symmetries_list": symmetries_list,
            "numBlocks": len(symmetries_list),
        },
        "runtimeoptions": runtimeoptions,
        "button_id": conf["button_id"],
    }
    return settings


def output(x):
    print(x, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scripts")
    parser.add_argument("--run_gui", action="store_true", help="Use this only to run a calculation from the GUI!")
    parser.add_argument("-c", "--path_config", type=str, help="Path to the config.")
    parser.add_argument("-o", "--path_cache", type=str, help="Path to the cache.")
    args = parser.parse_args()

    if args.run_gui:
        main(args)
    else:
        raise ValueError(
            "For now this script can only be called from the GUI and there is no other use calling this script!"
        )
