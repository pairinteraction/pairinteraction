# -*- mode: python -*-

block_cipher = None

mkl_dlls = [
    ('C:\\Program Files\\Miniconda3\\Library\\bin\\mkl_tbb_thread.dll', '')
]

a = Analysis(['startgui'],
             pathex=['..\\build\\gui'],
             binaries=mkl_dlls,
             datas=None,
             hiddenimports=['six', 'scipy.integrate'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['matplotlib', 'OpenGL', 'PyQt5.QtOpenGL'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='startgui',
          debug=False,
          strip=False,
          upx=True,
          console=True )
