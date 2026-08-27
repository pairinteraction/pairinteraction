# AGENTS.md

## Building

Always build via `uv sync`, do not use `cmake` directly.
If `uv` options are unavailable, update `uv`.
If `ccache` is not found, omit the `CMAKE_CXX_COMPILER_LAUNCHER` option.

Build, and rebuild after a change to C++ sources, with:

```bash
uv venv  # only if .venv did not exist before
uv pip install -r .build_requirements.txt  # only if .venv did not exist before
uv sync --inexact --no-build-isolation --reinstall-package pairinteraction -v \
  -Cbuild-dir=.venv/build_pip \
  -Ccmake.define.CMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -Cbuild.tool-args=-j4
```

Add `-v -Cbuild.verbose=true` when a build fails.

## Tips for Testing

If you changed C++ code, rebuild and run:

```bash
uv run --no-project pytest -k run_unit_tests
```

If you changed Python code, run the full test suite:

```bash
uv run --no-project pytest
```

Note that `--no-project` is always needed.

## Tips for Formatting and Linting

```bash
uvx pre-commit run --files $(git diff --name-only HEAD) $(git ls-files --others --exclude-standard)
```

## Container limits

If you are running inside a container, the host's CPU and memory figures overstate what you may actually use.

Read the CPU and memory the container may actually use from the cgroup.
Each value falls back to the host figure when the cgroup file is absent
(cgroup v1, non-Linux) or prints `max`:

```bash
cpus=$(awk '{print ($1=="max") ? 0 : ($1 < $2    ? 1 : int($1/$2))}'   /sys/fs/cgroup/cpu.max    2>/dev/null)
mem=$( awk '{print ($1=="max") ? 0 : ($1 < 2^30 ? 1 : int($1/2^30))}' /sys/fs/cgroup/memory.max 2>/dev/null)
[ "${cpus:-0}" -gt 0 ] || cpus=$(nproc)
[ "${mem:-0}"  -gt 0 ] || mem=$(free -g | awk 'NR==2 {print $2}')
echo "CPUs=$cpus  memory=${mem} GiB"
```

For the `-Cbuild.tool-args=-j4` argument above, use at most one job per CPU, and at most one per
2 GiB, because the template-heavy translation units (Eigen, duckdb, nanobind) peak around 2 GiB each.
