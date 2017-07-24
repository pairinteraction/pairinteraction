# -*- mode: python -*-

block_cipher = None

a = Analysis(['pairinteraction'],
             pathex=[],
             binaries=[],
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
