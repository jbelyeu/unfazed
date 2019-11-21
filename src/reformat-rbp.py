import sys
import toolshed as ts
from collections import defaultdict

for i,l in enumerate(ts.reader(sys.argv[1], header='ordered')):
    seen_variants = defaultdict(list)

    if i == 0: 
        keys = l.keys()
        keys.append('phase_rbp')
        print '\t'.join(keys)
    if l['mom_evidence'] == 'none' or l['dad_evidence'] == 'none':
        l['phase_rbp'] = 'na'
        print '\t'.join(l.values())
        continue
    mom_ev, dad_ev = int(l['mom_evidence']), int(l['dad_evidence'])
    if all([ev > 3 for ev in (mom_ev, dad_ev)]): 
        l['phase_rbp'] = 'na'
        print '\t'.join(l.values())
        continue
    if mom_ev + dad_ev < 2:
        l['phase_rbp'] = 'na'
        print '\t'.join(l.values())
        continue
    if mom_ev > dad_ev:
        l['phase_rbp'] = 'maternal'
    elif dad_ev > mom_ev:
        l['phase_rbp'] = 'paternal'
    else:
        l['phase_rbp'] = 'na'

    print '\t'.join(l.values())

