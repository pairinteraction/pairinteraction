"""The main function will be called by the pairinteraction GUI to start the calculation.
"""
import itertools
import json
import multiprocessing
import os
import pickle
import sys
import tempfile
import time

import numpy as np
from PyQt5.QtCore import QThread


# FIXME we should use from pairinteraction_gui import pipy here,
# however, this will lead to bugs in the pyinstall on mac and windows
# probably you have to fix this by adapting cmake and/or the .spec file
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import pipy  # noqa


P_ONE_RUN = None


class PipyThread(QThread):
    NO_POOL = "NO_POOL"

    def __init__(self, all_queues, context="default", pass_atom="direct&delete", parent=None):
        super().__init__(parent)
        self.all_queues = all_queues

        assert context in ["default", "fork", "spawn", "forkserver"]
        assert pass_atom in ["direct", "path", "direct&delete", "path&delete", "config"]
        self.context = context
        self.pass_atom = pass_atom

        self.params = None
        self.num_pr = None
        self._pool = None
        self._kwargs = None

    @property
    def pool(self):
        if self._pool is None:
            return self.createPool()
        return self._pool

    def createPool(self):
        print("Creating new pool")
        if self.num_pr == 1:
            self._pool = self.NO_POOL
        elif self.context == "default":
            self._pool = multiprocessing.Pool(self.num_pr)
        elif self.context in ["fork", "spawn", "forkserver"]:
            # mimic windows behaviour with spawn
            self._pool = multiprocessing.get_context(self.context).Pool(self.num_pr)
        return self._pool

    def killPool(self):
        if self._pool == self.NO_POOL:
            pass
            # self.exiting = True
        elif self._pool is not None:
            self._pool.terminate()
            self._pool.join()
            self._pool = None

    def setNumProcesses(self, num_pr):
        num_pr = None if num_pr == 0 else num_pr
        if self.num_pr != num_pr:
            self.killPool()
            self.num_pr = num_pr

    def setParams(self, params):
        self.params = params

    @property
    def kwargs(self):
        self._kwargs = {"pool": self.pool, "printFunction": self.all_queues.processOneLine, "params": self.params}
        return self._kwargs

    def run(self):
        self.start_pipy(self.kwargs)

    def terminate(self):
        # FIXME this is not the cleanest way, maybe using multiprocessing.manager/event/value would be better
        # but also might be slower and introduces bugs, where the gui does not close, althoug the terminal terminated
        current_time = time.time()
        self.killPool()
        super().terminate()
        self.wait()

        # delete files, that where changed during pool.terminate was called
        for real_complex in ["real", "complex"]:
            pathCacheMatrix = os.path.join(self.params["pathCache"], f"cache_matrix_{real_complex}_new")
            if not os.path.isdir(pathCacheMatrix):
                continue
            for fn in os.listdir(pathCacheMatrix):
                f = os.path.join(pathCacheMatrix, fn)
                if not os.path.isfile(f):
                    continue
                if current_time < os.path.getmtime(f):
                    os.remove(f)

        atom = self._kwargs.get("atom", None)
        if atom is not None:
            atom.delete()

    def start_pipy(self, kwargs):
        # convert params to settings and save as settings.json
        params = kwargs["params"]
        settings = params_to_settings(params)
        with open(params["pathConf"].replace("conf.json", "settings.json"), "w") as f:
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
                self.run_simulations(settings, kwargs)
        else:
            settings["blocknumber"] = 1
            self.run_simulations(settings, kwargs)
        info("all Hamiltonian processed", kwargs)
        output(f"{'>>END':5}", kwargs)

    def run_simulations(self, settings, kwargs):
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
        param_list = self.get_param_list(settings)
        info("precalculate matrix elements", kwargs)
        # make sure also interactions are precalculated (distance must be not None) and all fields, that might occur
        atom.updateFromParams(param_list[0])
        atom.system.buildInteraction()
        atom.updateFromParams(param_list[-1])
        atom.system.buildInteraction()

        # Decide how to pass the atom to the subprocesses
        if "delete" in self.pass_atom:
            atom.delete()
            # delete the cpp object, so it is not pickled and let the subprocess recreate them
            # Ideally, each subprocess recreates it just once. This seems to work for forked processes,
            # however for spawned processes, you should check it again! # FIXME
        if "direct" in self.pass_atom:
            config = {"atom": atom}
        elif "path" in self.pass_atom:
            atom_path = path_basis.replace("basis_", "atom_").replace(".csv", ".pkl")
            with open(atom_path, "wb") as file:
                pickle.dump(atom, file)
            config = {"atom_path": atom_path}
        elif self.pass_atom == "config":
            config = atom.config.toOutDict()

        # Some definitions for the multiprocessing
        ip_list = list(range(len(param_list)))
        starmap_list = [(ip, config, param_list) for ip in ip_list]

        global P_ONE_RUN

        def P_ONE_RUN(ip):
            return self.one_run(ip, config, param_list)

        # The actual calculation
        pool = kwargs.get("pool", None)
        if pool == PipyThread.NO_POOL or pool is None:
            for i in ip_list:
                result = P_ONE_RUN(i)
                self.print_completed(settings, result, kwargs)
        else:
            if (
                pool._ctx.get_start_method() != "fork"
            ):  # this is slower on linux and sometimes had weird bugs when aborting calculation
                results = pool.starmap_async(self.one_run, starmap_list).get()
            else:  # fork
                # restart pool, because the pool must be created after defining the function
                # this somehow is still faster than using starmap_async
                self.killPool()
                pool = self.pool
                results = pool.imap_unordered(P_ONE_RUN, ip_list)
                # this is slower, but equivalent to spawn
                # results = pool.starmap_async(self.one_run, starmap_list).get()
            for result in results:
                self.print_completed(settings, result, kwargs)
            # don't terminate pool / use with pool, such that we can reuse it

        # Clean up
        atom.delete()
        del kwargs["atom"]
        if "atom_path" in config:
            os.remove(config["atom_path"])

    @staticmethod
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

            print(
                f"{ip+1}. Hamiltonian diagonalized ({dimension}x{dimension}) ({multiprocessing.current_process().name})"
            )
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

    @staticmethod
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

    @staticmethod
    def print_completed(settings, result, kwargs):
        output(f"{'>>DIM':5}{result['dimension']:7}", kwargs)
        numBlocks = settings["scriptoptions"]["numBlocks"]
        totalstep = 1 + result["ip"] * numBlocks + settings["blocknumber"]
        output(
            f"{'>>OUT':5}{totalstep:7}{result['ip']:7}{numBlocks:7}"
            f"{settings['blocknumber']:7} {result['filename']} ",
            kwargs,
        )


def params_to_settings(params):
    config = {}
    for k, v in params.items():
        if (
            k.startswith("delta")
            or any(k == q + i for q in ["species", "n", "l", "j", "m"] for i in ["", "1", "2"])
            or k in ["diamagnetism", "samebasis", "exponent"]
        ):
            config[k] = v

    config["method"] = ["WHITTAKER", "NUMEROV"][params["missingCalc"]]
    config["pathCache"] = params["pathCache"]
    if params["button_id"] == 2:
        config["pair.conserveMomenta"] = params.get("conserveM", False)
        config["pair.useQunumbersSinglePair"] = True
        config["nAtoms"] = 2
    elif params["button_id"] in [0, 1]:
        config["nAtoms"] = 1
        for k in ["species", "n", "l", "j", "m"]:
            config[k] = config.pop(k + "1", config.pop(k + "2", None))
    else:  # 3
        config["nAtoms"] = 1
    config.update({k: params.get("min" + k, 0) for k in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]})

    listoptions = {
        "onlySameSector": params.get("sametrafo", None),
        "steps": params.get("steps", None),
        **{mm + "distance": params.get(mm + "R", None) for mm in ["min", "max"]},
        **{
            mm + key: params.get(mm + key, None)
            for key in ["Ex", "Ey", "Ez", "Bx", "By", "Bz"]
            for mm in ["min", "max"]
        },
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
            v = [eo for eo in ["O", "E"] if params.get(s + eo, False)]
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

    runtimeoptions = {"NUM_PROCESSES": int(params.get("NUM_PROCESSES", 0))}

    # touch unimportant options
    for k in ["dd", "dq", "qq", "missingWhittaker", "precision", "zerotheta"]:
        params.get(k, None)

    settings = {
        "config": config,
        "scriptoptions": {
            "listoptions": listoptions,
            "symmetries_list": symmetries_list,
            "numBlocks": len(symmetries_list),
        },
        "runtimeoptions": runtimeoptions,
        "button_id": params["button_id"],
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
