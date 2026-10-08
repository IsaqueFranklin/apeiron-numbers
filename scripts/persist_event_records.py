"""Replace delivery symlinks with verified persistent CSV copies after space is freed."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

parser=argparse.ArgumentParser()
parser.add_argument('results',type=Path)
args=parser.parse_args()
root=args.results.resolve()
links=[p for p in root.glob('*/events.csv') if p.is_symlink()]
required=sum(p.stat().st_size for p in links)+64*1024*1024
if shutil.disk_usage(root).free<required:
    raise SystemExit(f'Need {required/1e9:.2f} GB free to preserve complete events; existing data unchanged.')
verification=json.loads((root/'comparacao'/'verificacao-integridade.json').read_text())
for link in links:
    source=link.resolve(strict=True)
    temporary=link.with_name('events.csv.persisting')
    if temporary.exists(): raise FileExistsError(temporary)
    digest=hashlib.sha256()
    with source.open('rb') as reader, temporary.open('xb') as writer:
        while block:=reader.read(8*1024*1024):
            writer.write(block); digest.update(block)
    expected=verification['scenarios'][link.parent.name]['events_sha256']
    if digest.hexdigest()!=expected:
        raise SystemExit(f'Hash mismatch in {temporary}; original link preserved.')
    temporary.replace(link)
    print(f'Preserved {link} ({source.stat().st_size/1e6:.1f} MB)',flush=True)
manifest_path=root/'manifest.json'
manifest=json.loads(manifest_path.read_text())
manifest.setdefault('delivery',{})['all_event_records_persistent']=True
manifest_path.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
