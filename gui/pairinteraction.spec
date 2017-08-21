# -*- mode: python -*-

block_cipher = None

a = Analysis(['pairinteraction'],
             pathex=[],
             binaries=[],
             datas=[],
             hiddenimports=['six', 'scipy.integrate'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['matplotlib', 'OpenGL', 'PyQt5.QtOpenGL'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
a.binaries = [x for x in a.binaries if not not os.path.basename(x[1]).startswith("api-ms-win") and not os.path.basename(x[1]).startswith("mkl_") and not os.path.basename(x[1]) in ["mfc140u.dll", "MSVCP140.dll", "VCRUNTIME140.dll"] or os.path.basename(x[1]) in ["mkl_core.dll", "mkl_def.dll", "mkl_intel_thread.dll"]]
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='pairinteraction',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='pairinteraction')
