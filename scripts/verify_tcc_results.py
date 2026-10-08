"""Audit completed TCC runs, reading all events/times/bits in their original order."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from apeiron.analysis import bits_from_waiting_times, nist_subset
from apeiron.model import K_B


def verify(root, expected=1_000_000):
    started=time.monotonic()
    root=Path(root)
    manifest=json.loads((root/'manifest.json').read_text())
    if manifest['status'] != 'complete': raise ValueError('manifest is not complete')
    source_checks={}
    for name,digest in manifest['provenance']['source_sha256'].items():
        actual=hashlib.sha256((root/'source'/name).read_bytes()).hexdigest()
        if actual != digest: raise ValueError(f'source hash mismatch: {name}')
        source_checks[name]=actual
    result={'status':'complete','expected_events_per_scenario':expected,'source_checks':source_checks,'scenarios':{}}
    initial=[]
    for entry in manifest['scenarios']:
        name=entry['name']; directory=root/name
        summary=json.loads((directory/'summary.json').read_text())
        if summary['status'] != 'complete': raise ValueError(f'{name}: incomplete')
        initial.append(summary['lattice']['initial_occupancy_sha256'])
        times=[]
        with (directory/'waiting_times.csv').open(newline='') as handle:
            for i,row in enumerate(csv.DictReader(handle),1):
                if int(row['event']) != i: raise ValueError(f'{name}: waiting-time order')
                times.append(float(row['dt_s']))
        times=np.asarray(times,dtype=float)
        raw=(directory/'bits.txt').read_bytes().strip()
        if len(times) != expected or len(raw) != expected: raise ValueError(f'{name}: counts')
        bits,threshold=bits_from_waiting_times(times)
        if raw != (bits+48).tobytes(): raise ValueError(f'{name}: bit extraction')
        if threshold != summary['digitization']['median_dt_s']: raise ValueError(f'{name}: threshold')
        recorded_tests = summary['nist_sp800_22_subset']
        if nist_subset(bits, block_size=recorded_tests['block_size']) != recorded_tests:
            raise ValueError(f'{name}: test summaries')
        cumulative=0.; sampled=0; offered=0
        with (directory/'events.csv').open(newline='') as handle:
            for count,event in enumerate(csv.DictReader(handle),1):
                if int(event['event']) != count: raise ValueError(f'{name}: event order')
                if float(event['dt_s']) != times[count-1]: raise ValueError(f'{name}: event/time mismatch')
                cumulative+=float(event['dt_s'])
                if cumulative != float(event['time_s']): raise ValueError(f'{name}: clock mismatch')
                barriers=json.loads(event['barriers_eV'])
                if len(barriers) != 12 or min(barriers)<0: raise ValueError(f'{name}: barriers')
                offered+=len(barriers)
                if count == 1 or count % 10000 == 0 or count == expected:
                    sampled+=1
                    saddles=json.loads(event['saddles_eV']); rates=json.loads(event['rates_per_s'])
                    w=float(event['source_well_eV']); total=float(event['total_rate_per_s'])
                    if any(e != s-w for e,s in zip(barriers,saddles)): raise ValueError(f'{name}: barrier formula')
                    expected_rates=[summary['parameters']['attempt_frequency']*math.exp(-e/(K_B*summary['parameters']['temperature'])) for e in barriers]
                    if expected_rates != rates: raise ValueError(f'{name}: rate formula')
                    if math.fsum(rates) != total: raise ValueError(f'{name}: total rate')
                    if -math.log(float(event['waiting_u']))/total != float(event['dt_s']): raise ValueError(f'{name}: residence time')
        if count != expected or offered != 12*expected: raise ValueError(f'{name}: event/offered counts')
        if cumulative != summary['total_time_s']: raise ValueError(f'{name}: total time')
        for key in ('mu_s','mu_w','sigma_s','sigma_w'):
            if summary['parameters'][key] != summary['preflight']['parameters'][key]: raise ValueError(f'{name}: preflight parameters')
        preflight=np.load(directory/'preflight_barriers.npy',mmap_mode='r')
        if len(preflight)!=summary['preflight']['count'] or len(preflight)!=expected: raise ValueError(f'{name}: preflight count')
        if not math.isclose(float(preflight.mean()),summary['preflight']['mean_eV'],rel_tol=1e-12,abs_tol=1e-14): raise ValueError(f'{name}: preflight mean')
        if not math.isclose(float(preflight.std()),summary['preflight']['std_eV'],rel_tol=1e-12,abs_tol=1e-14): raise ValueError(f'{name}: preflight deviation')
        if summary['lattice']['initial_composition'] != summary['lattice']['final_composition']: raise ValueError(f'{name}: composition')
        digest = hashlib.sha256()
        with (directory/'events.csv').open('rb') as binary:
            while block := binary.read(8*1024*1024):
                digest.update(block)
        result['scenarios'][name]={'events':count,'waiting_times':len(times),'bits':len(raw),'offered_barriers':offered,
            'preflight_barriers':len(preflight),'selected_barriers':summary['selected_barriers']['count'],
            'events_sha256':digest.hexdigest(),
            'all_events_order_clock_barrier_counts_verified':True,'all_waiting_time_bits_verified':True,
            'sampled_full_rate_formula_events':sampled,'nist_summary_reproduced':True}
        print(f'{name}: {count:,} events, {len(raw):,} bits; all checks passed',flush=True)
    if len(set(initial)) != 1: raise ValueError('different initial chemical configurations')
    result['same_initial_chemical_configuration']=True
    result['duration_s']=time.monotonic()-started
    path=root/'comparacao'/'verificacao-integridade.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    parser.add_argument('--expected-events',type=int,default=1_000_000)
    args=parser.parse_args()
    verify(args.output,args.expected_events)
