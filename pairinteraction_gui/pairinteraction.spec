# -*- mode: python -*-
import fnmatch
import os

block_cipher = None

a = Analysis(
    ["pairinteraction_gui/pairinteraction_app/app.py"],
    pathex=[".", "pairinteraction"],
    binaries=[],
    datas=[("pairinteraction_gui/pairinteraction_app/icon.png", "."), ("pairinteraction_gui/conf", "conf")],
    hiddenimports=["scipy.integrate", "scipy._lib.messagestream"],
    hookspath=[],
    runtime_hooks=[],
    excludes=["matplotlib", "OpenGL", "PyQt5.QtOpenGL"],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
)

forbidden_binaries = [
    "api-ms-win*",
    "mkl_*",
    "mfc140u.dll",
    "MSVCP140.dll",
    "VCRUNTIME140.dll",
]
exempted_binaries = [
    "mkl_rt*.dll",
    "mkl_core*.dll",
    "mkl_def*.dll",
    "mkl_mc*.dll",
    "mkl_intel_thread*.dll",
    "mkl_intel_thread*.dll",
    "mkl_vml_avx*.dll",
    "mkl_vml_def*.dll",
    "mkl_vml_mc*.dll",
    "mkl_vml_cmpt*.dll",
    "libimalloc.dll",
]
a.binaries = [
    x
    for x in a.binaries
    if not any(fnmatch.fnmatch(os.path.basename(x[1]), y) for y in forbidden_binaries)
    or any(fnmatch.fnmatch(os.path.basename(x[1]), y) for y in exempted_binaries)
]

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(
    pyz, a.scripts, exclude_binaries=True, name="pairinteraction_gui", debug=False, strip=False, upx=True, console=True
)
coll = COLLECT(exe, a.binaries, a.zipfiles, a.datas, strip=False, upx=True, name="pairinteraction_gui")
