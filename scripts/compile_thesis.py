"""Compile the two assembled Typst documents using an installed Typst SDK."""
import argparse
import json
from pathlib import Path
import time
import typst

parser=argparse.ArgumentParser()
parser.add_argument('bundle',type=Path)
parser.add_argument('--font',default='Liberation Serif')
parser.add_argument('--package-path',default='/tmp/tcc-packages')
parser.add_argument('--package-cache-path',default='/tmp/tcc-typst-cache')
args=parser.parse_args()
root=args.bundle.resolve()
records=[]
for stem in ('resultados-conclusao','tcc-com-resultados'):
    started=time.monotonic()
    typst.compile(str(root/f'{stem}.typ'),output=str(root/f'{stem}.pdf'),root=str(root),
                  package_path=args.package_path,package_cache_path=args.package_cache_path,
                  sys_inputs={'fonte':args.font})
    records.append({'source':f'{stem}.typ','output':f'{stem}.pdf','duration_s':time.monotonic()-started})
    print(f'Compiled {stem}.pdf',flush=True)
(root/'compilacao.json').write_text(json.dumps({'font':args.font,'typst_version':typst.__version__,
    'package':'@preview/abntyp:0.1.5','documents':records},indent=2)+'\n')
