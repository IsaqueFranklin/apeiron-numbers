"""Recalculate the existing bits with the corrected block protocol, preserving history."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from apeiron.analysis import nist_subset
from apeiron.cli import provenance, write_json
from apeiron.comparison import compare_results
from apeiron.plotting import plot_results


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            result.update(chunk)
    return result.hexdigest()


def reanalyse(root, plot=True):
    root = Path(root).resolve()
    manifest = json.loads((root/'manifest.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('cannot reanalyse incomplete simulations')
    names = [s['name'] for s in manifest['scenarios']]
    snapshots = {name: json.loads((root/name/'summary.json').read_text()) for name in names}
    if any(s['status'] != 'complete' for s in snapshots.values()):
        raise ValueError('cannot reanalyse incomplete scenarios')
    history = root/'historico'/'avaliacao-m128'
    if history.exists():
        raise FileExistsError(f'history already exists: {history}')
    history.mkdir(parents=True)
    # Retain the original summaries, report, figures, audit and source references.
    shutil.copy2(root/'manifest.json', history/'manifest.json')
    shutil.copytree(root/'comparacao', history/'comparacao', ignore=shutil.ignore_patterns('figuras'))
    for name in names:
        folder = history/name
        folder.mkdir()
        shutil.copy2(root/name/'summary.json', folder/'summary.json')
        for path in (root/name).glob('diagnostics.*'):
            shutil.copy2(path, folder/path.name)
    started = time.monotonic()
    source = root/'statistical-analysis-source'
    source.mkdir(exist_ok=True)
    project = Path(__file__).resolve().parents[1]
    metadata = provenance()
    hashes = {}
    for relative in (*metadata['source_sha256'], 'scripts/reanalyse_nist.py'):
        target = source/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(project/relative, target)
        hashes[relative] = digest(target)
    record = {'status': 'running', 'utc_started': datetime.now(timezone.utc).isoformat(),
              'command': list(sys.argv), 'history': 'historico/avaliacao-m128',
              'protocol': {'block_size_rule': 'max(20, ceil(n/50))', 'alpha': .01,
                           'suite': 'SP 800-22 Rev. 1a monobit, block frequency, runs'},
              'provenance': metadata, 'source_sha256': hashes,
              'raw_simulation_data_changed': False, 'scenarios': {}}
    record_path = root/'comparacao'/'reavaliacao-nist.json'
    write_json(record_path, record)
    try:
        for name, summary in snapshots.items():
            folder = root/name
            old_block_size = summary['nist_sp800_22_subset']['block_size']
            before = {file: digest(folder/file) for file in ('bits.txt', 'waiting_times.csv')}
            raw = (folder/'bits.txt').read_bytes().strip()
            if set(raw)-{48,49}:
                raise ValueError(f'{name}: invalid bits')
            bits = np.frombuffer(raw, dtype=np.uint8)-48
            if len(bits) != summary['digitization']['bits']:
                raise ValueError(f'{name}: changed bit count')
            tests = nist_subset(bits)
            if not all(tests['block_recommendations'].values()):
                raise ValueError(f'{name}: block recommendation failed')
            summary['nist_sp800_22_subset'] = tests
            summary['statistical_reanalysis'] = {
                'record': '../comparacao/reavaliacao-nist.json', 'old_block_size': old_block_size,
                'block_size_rule': 'max(20, ceil(n/50))', 'same_bits_and_waiting_times': True}
            # Original simulation provenance and parameters remain recorded separately.
            write_json(folder/'summary.json', summary)
            if plot:
                plot_results(folder)
            after = {file: digest(folder/file) for file in before}
            if before != after:
                raise ValueError(f'{name}: raw arrays changed during reanalysis')
            record['scenarios'][name] = {'bits': len(bits), 'before_sha256': before, 'after_sha256': after,
                                         'block_size': tests['block_size'], 'complete_blocks': tests['complete_blocks'],
                                         'residual_bits': tests['residual_bits'], 'tests': tests['tests']}
            print(f'{name}: {len(bits):,} original bits, M={tests["block_size"]}, {tests["complete_blocks"]} blocks', flush=True)
        if plot:
            # Retain the original files instead of deleting generated or user-added figures.
            (root/'comparacao'/'figuras').rename(history/'comparacao'/'figuras')
        compare_results(root, plot=plot)
        record.update({'status': 'complete', 'duration_s': time.monotonic()-started,
                       'utc_finished': datetime.now(timezone.utc).isoformat()})
        write_json(record_path, record)
        manifest['statistical_reanalysis'] = {
            'status': 'complete', 'record': 'comparacao/reavaliacao-nist.json',
            'block_size_rule': record['protocol']['block_size_rule'], 'simulation_rerun': False}
        write_json(root/'manifest.json', manifest)
    except BaseException as exc:
        record.update({'status': 'incomplete', 'failure': f'{type(exc).__name__}: {exc}'})
        write_json(record_path, record)
        raise
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('results', type=Path)
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args()
    reanalyse(args.results, plot=not args.no_plot)
