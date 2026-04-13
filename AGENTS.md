# AGENTS.md

## Building

Building from scratch can take a long time. Use the following commands to take advantage of caching. If `uv` options are unavailable, update `uv`. If `ccache` is not found, omit the `CMAKE_CXX_COMPILER_LAUNCHER` option.

```bash
uv venv
uv pip install -r .build_requirements.txt
uv sync --reinstall --inexact --no-build-isolation -v \
  -Cbuild.verbose=true \
  -Cbuild-dir=build_pip \
  -Ccmake.define.CMAKE_CXX_COMPILER_LAUNCHER=ccache
```

Always build via `uv sync`, do not use `cmake` directly.

## Tips for Testing

If you changed C++ code, rebuild and run:

```bash
uv run --no-project pytest -k run_unit_tests --durations=0
```

If you changed Python code, run the test suite:

```bash
uv run --no-project pytest --durations=0
```

Note that `--no-project` is always needed.

## Tips for Formatting and Linting

```bash
uvx pre-commit run --files $(git diff --name-only HEAD)
```
