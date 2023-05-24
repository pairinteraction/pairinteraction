"""The main function will be called by the pairinteraction GUI to start the calculation.
"""
import itertools
import json
import multiprocessing
import os
import pickle
import sys
import tempfile

import numpy as np


# FIXME we should use from pairinteraction_gui import pipy here,
# however, this will lead to bugs in the pyinstall on mac and windows
# probably you have to fix this by adapting cmake and/or the .spec file
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import pipy  # noqa


def main(paths, kwargs):
    kwargs.setdefault("printFunction", lambda msg: print(msg, flush=True, end=""))

    path_conf, path_cache = paths["path_conf"], paths["path_cache"]
    # Load params from conf.json file
    with open(path_conf) as f:
        conf = json.load(f)

    # convert conf to settings and save as settings.json
    settings = conf_to_settings(conf, path_cache)
    with open(path_conf.replace("conf.json", "settings.json"), "w") as f:
        json.dump(settings, f, indent=4)

    # start the calculation
    config = settings["config"]
    scriptoptions = settings["scriptoptions"]
    output(f"{'>>TYP':5}{settings['button_id']:7}", kwargs)
    output(f"{'>>TOT':5}{scriptoptions['numBlocks']*scriptoptions['listoptions']['steps']:7}", kwargs)
    if config["nAtoms"] == 2:
        for bn, syms in enumerate(scriptoptions["symmetries_list"]):
            config.update(syms)
            settings["blocknumber"] = bn
            do_simulations(settings, kwargs)
    else:
        settings["blocknumber"] = 1
        do_simulations(settings, kwargs)
    info("all Hamiltonian processed", kwargs)
    output(f"{'>>END':5}", kwargs)


def do_simulations(settings, kwargs, pass_atom="direct&delete", context="default"):
    assert context in ["default", "fork", "spawn", "forkserver"]
    assert pass_atom in ["direct", "path", "direct&delete", "path&delete", "config"]

    kwargs["atom"] = atom = pipy.atom_from_config(settings["config"])
    info(f"construct {atom}", kwargs)

    # create and save atom basis
    if atom.nAtoms == 2:
        allQunumbers = [[*qns[:, 0], *qns[:, 1]] for qns in np.array(atom.allQunumbers)]
    else:
        allQunumbers = atom.allQunumbers
    info(f"basis size with restrictions: {len(allQunumbers)}", kwargs)
    _name = f"{'one' if atom.nAtoms == 1 else 'two'}_{atom.config.toHash()}"
    path_tmp = tempfile.gettempdir()
    path_basis = os.path.join(path_tmp, f"basis_{_name}_blocknumber_{settings['blocknumber']}.csv")
    basis = np.insert(allQunumbers, 0, np.arange(len(allQunumbers)), axis=1)
    np.savetxt(path_basis, basis, delimiter="\t", fmt=["%d"] + ["%d", "%d", "%.1f", "%.1f"] * atom.nAtoms)
    info(f"save {type(atom).__name__} basis", kwargs)
    output(f"{'>>BAS':5}{len(basis):7}", kwargs)
    output(f"{'>>STA':5} {path_basis} ", kwargs)

    # precalculate matrix elements
    param_list = get_param_list(settings)
    info("precalculate matrix elements", kwargs)
    # make sure also interactions are precalculated (distance must be not None) and all fields, that might occur
    atom.updateFromParams(param_list[0])
    atom.system.buildInteraction()
    atom.updateFromParams(param_list[-1])
    atom.system.buildInteraction()

    # Decide how to pass the atom to the subprocesses
    if "delete" in pass_atom:
        atom.delete()
        # delete the cpp object, so it is not pickled and let the subprocess recreate them
        # Ideally, each subprocess recreates it just once. This seems to work for forked processes,
        # however for spawned processes, the recreation happens for each run - even if the same SpawnPoolWorker is used
        # Why is the state of the "Atom" object not stored in the spawned processes? # FIXME
    if "direct" in pass_atom:
        config = {"atom": atom}
    elif "path" in pass_atom:
        atom_path = path_basis.replace("basis_", "atom_").replace(".csv", ".pkl")
        with open(atom_path, "wb") as file:
            pickle.dump(atom, file)
        config = {"atom_path": atom_path}
    elif pass_atom == "config":
        config = atom.config.toOutDict()

    # Some definitions for the multiprocessing
    ip_list = list(range(len(param_list)))
    run_kwargs = {"config": config, "param_list": param_list}
    dict_list = [(ip, run_kwargs) for ip in ip_list]

    global P_ONE_RUN

    def P_ONE_RUN(ip):
        return one_run(ip, config, param_list)

    num_pr = settings.get("runtimeoptions", {}).get("NUM_PROCESSES", 1)
    num_pr = os.cpu_count() if num_pr in [0, -1] else num_pr

    # The actual calculation
    if num_pr > 1:
        if context == "default":
            kwargs["pool"] = pool = multiprocessing.Pool(num_pr)
        elif context in ["fork", "spawn", "forkserver"]:
            kwargs["pool"] = pool = multiprocessing.get_context(context).Pool(
                num_pr
            )  # mimic windows behaviour with spawn
        if pool._ctx.get_start_method() == "fork":
            results = pool.imap_unordered(P_ONE_RUN, ip_list)
        else:  # this is slower on linux and sometimes had weird bugs when aborting calculation
            results = pool.imap_unordered(one_run_dict, dict_list)
        for result in results:
            print_completed(settings, result, kwargs)
        pool.terminate()  # equivalent to with pool
    else:
        for i in ip_list:
            result = P_ONE_RUN(i)
            print_completed(settings, result, kwargs)

    # Clean up
    atom.delete()
    del kwargs["atom"]
    if "atom_path" in config:
        os.remove(config["atom_path"])


def one_run_dict(args):
    ip, run_kwargs = args
    return one_run(ip, **run_kwargs)


def one_run(ip, config, param_list):
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
    pathCacheMatrix = os.path.join(config.pathCache(), f"cache_matrix_{real_complex}_new")
    os.makedirs(pathCacheMatrix, exist_ok=True)
    _name = f"{'one' if atom.nAtoms == 1 else 'two'}_{config.toHash()}"
    filename = os.path.join(pathCacheMatrix, _name + ".pkl")
    filename_json = os.path.join(pathCacheMatrix, _name + ".json")

    if not os.path.exists(filename) or not os.path.exists(filename_json):
        atom.calcEnergies()
        energies0 = config.getEnergiesPair() if atom.nAtoms == 2 else config.getEnergiesSingle()
        data = {
            "energies": atom.energies - np.mean(energies0),
            "basis": atom.vectors,
            "params": config.toOutDict(),
        }

        print(f"{ip+1}. Hamiltonian diagonalized ({dimension}x{dimension}) ({multiprocessing.current_process().name})")
        with open(filename, "wb") as f:
            pickle.dump(data, f)
        with open(filename_json, "w") as f:
            json.dump(data["params"], f, indent=4)
    else:
        print(f"{ip+1}. Hamiltonian loaded")

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


def print_completed(settings, result, kwargs):
    output(f"{'>>DIM':5}{result['dimension']:7}", kwargs)
    numBlocks = settings["scriptoptions"]["numBlocks"]
    totalstep = 1 + result["ip"] * numBlocks + settings["blocknumber"]
    output(
        f"{'>>OUT':5}{totalstep:7}{result['ip']:7}{numBlocks:7}" f"{settings['blocknumber']:7} {result['filename']} ",
        kwargs,
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

    runtimeoptions = {"NUM_PROCESSES": int(conf.get("NUM_PROCESSES", 0))}

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


def output(msg, kwargs):
    printFunction = kwargs["printFunction"]
    printFunction(msg + "\n")


def info(msg, kwargs):
    printFunction = kwargs["printFunction"]
    atom = kwargs.get("atom", None)
    if atom is not None:
        pi = "pireal" if atom.config.isReal() else "picomplex"
        no = "One-atom" if atom.config.nAtoms() == 1 else "Two-atom"
        msg = f"{pi}: {no} Hamiltonian, {msg}"
    printFunction(msg + "\n")
