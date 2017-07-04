#!/usr/bin/python
# Timing test about how long it takes to generate a diff list

import qnetu
import sampling
import mytime

net_yaml = """
states:
- name: INITIAL
  queues: [INITIAL]
  successors: [TIER1]
  initial: TRUE
- name: TIER1
  queues: [ TIER1_0 , TIER1_1 ]
  successors: [TIER2]
- name: TIER2
  queues: [ TIER2_0 , TIER2_1 ]
queues:
- { name: INITIAL, service: [G, 2, 0.5 ] }
- { name: TIER1_0, processors: 3, service: [LN, 2.25, 0.2] }
- { name: TIER1_1, processors: 3, service: [LN, 2.25, 0.2] }
- { name: TIER2_0, processors: 2, service: [LN, 1.75, 0.2] }
- { name: TIER2_1, processors: 2, service: [LN, 1.75, 0.2] }
"""

def main ():
    sampling.set_seed (13372)

    net = qnetu.qnet_from_text (net_yaml)
    arrv = net.sample (1000)
    subset = arrv.subset_by_task (0.0)
    ne = arrv.num_events()

    tmr = mytime.timeit()
    N = 1000
    evts = arrv.all_events()

    tot_size = 0
    nl = 0

    for rep in range(N):
        evt = sampling.random_elt (evts)
        q = evt.queue()
        if evt.a != 0:
            dl = q.pyDiffListForArrival (evt, evt.a + 0.01)
            tot_size += len(dl)
            nl += 1

        dl = q.pyDiffListForDeparture (evt, evt.d + 0.01)
        tot_size += len(dl)
        nl += 1

    elapsed = tmr.total ("Time for %d diff lists" % N)

    print "Diff lists per second = %.4f" % (N / elapsed)
    print "Average size = %.4f" % (tot_size / float(nl))

if __name__ == "__main__":
    main()

