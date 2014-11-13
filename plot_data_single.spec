a = Analysis(['plot_data_single.py'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)

import mpl_toolkits.basemap
import os

src_basedata = os.path.join(mpl_toolkits.basemap.__path__[0], "data")
tgt_basedata = os.path.join('mpl_toolkits', 'basemap', 'data')

exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas + Tree(src_basedata, prefix=tgt_basedata),
          name='plot_data_single',
          debug=False,
          strip=None,
          upx=True,
          console=True )
