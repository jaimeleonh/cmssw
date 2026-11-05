#!/usr/bin/env python3
"""Summarize recorded cmsRun cases and pair the separately inspected scores."""
import argparse
import gzip
import json
from pathlib import Path
import math


def scores(path):
    if path.exists():return [json.loads(x) for x in path.read_text().splitlines()]
    with gzip.open(str(path)+'.gz','rt') as f:return [json.loads(x) for x in f]


def compare(directory,left,right):
    a=scores(directory/(left+'.scores.jsonl'));b=scores(directory/(right+'.scores.jsonl'))
    error=[]
    for x,y in zip(a,b):
        assert x['event']==y['event'] and x['jets']==y['jets']
        assert x['mask_population_histogram']==y['mask_population_histogram']
        assert x['nonfinite']==0 and y['nonfinite']==0
        assert len(x['scores'])==len(y['scores'])==x['jets']
        error.extend(abs(i-j) for i,j in zip(x['scores'],y['scores']))
    assert len(a)==len(b)
    return {'left':left,'right':right,'events':len(a),'jets':len(error),
            'max_abs_error':max(error),'mean_abs_error':sum(error)/len(error),
            'rmse':math.sqrt(sum(x*x for x in error)/len(error))}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path,nargs='?',default=Path(__file__).resolve().parents[1]/'performance_results/cmssw')
    a=p.parse_args(); rows=[]
    for path in sorted(a.directory.glob('*.json')):
        x=json.loads(path.read_text())
        if 'settings' not in x or x.get('returncode')!=0:continue
        sync=x.get('steady',{}).get('SoftJetIdSC4',{}).get('mean_ms')
        completion=x.get('completion_throughput',{})
        rows.append({'case':path.stem,'settings':x['settings'],'measurement_note':x.get('measurement_note'),
            'tagger_diagnostic_ms':sync,'pipeline_ms':completion.get('ms_per_event'),
            'events_per_second':completion.get('events_per_second'),
            'gpu_peak_MiB':x['gpu_peak_sampled_bytes']/2**20,
            'gpu_above_idle_MiB':(x['gpu_peak_sampled_bytes']-x['gpu_idle_bytes'])/2**20,
            'rss_peak_MiB':x['rss_peak_sampled_bytes']/2**20})
    sample=scores(a.directory/'validate_deepsets_graph.scores.jsonl')
    h=sample[0]['mask_population_histogram'];n=sum(h)
    checks=[]
    for left,right in [('validate_deepsets_graph','validate_deepsets_torch'),
        ('validate_part_graph','validate_part_torch'),('validate_part_int8','validate_part_torch'),
        ('validate_qat_graph','validate_qat_compiled')]:
        try:checks.append(compare(a.directory,left,right))
        except FileNotFoundError:pass
    output={'cases':rows,'score_comparisons':checks,
        'input':{'jets_per_inspected_event':n,'inspected_events':len(sample),'mask_population_histogram':h,
            'mean_valid_particles':sum(i*k for i,k in enumerate(h))/n,
            'one_particle_fraction':h[1]/n,'full_16_particle_fraction':h[16]/n,
            'distinct_score_vectors_in_ten_inspected_events':len(set(tuple(x['scores']) for x in sample))}}
    (a.directory/'summary.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output['score_comparisons'],indent=2))
    for x in rows:print(x['case'],x['pipeline_ms'],x['tagger_diagnostic_ms'],round(x['gpu_peak_MiB'],1))
if __name__=='__main__':main()
