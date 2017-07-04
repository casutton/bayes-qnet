#
#  List of useful statistics that could be caluclated from an Arrivals object.
#

import numpy
import arrivals

## utilitise

def averaging (stat):
    """Converts a statistic for a single arrival to one that
       averages over a list of arrivals."""
    return lambda arrvl: numpy.mean ([ stat(a) for a in arrvl ])

def aggregate_by_q (stat, arrv, aggregator=numpy.mean, filter=lambda evt: True):
    """Converts a statitstic that returns a value for a single event
       to a function that maps an Arrivals to a list of the mean statistic
       value for each queue."""
    X = [ list() for i in xrange(arrv.num_queues()) ]
    for evt in arrv:
        if filter(evt):
            X[evt.qid].append (stat(evt))
    return [ aggregator(xs) for xs in X ]

def aggregate_by_task (stat, arrv, aggregator=numpy.mean, filter=lambda tid: True):
    X = []
    for tid in arrv.all_tasks():
        if filter(arrv.events_of_task (tid)):
            X.append (stat(arrv.events_of_task(tid)))
    if len(X) == 0: return 0
    return aggregator (X)


def response_time (evt): return evt.d - evt.a

def is_obs (e):
    return e.obs_a and e.obs_d

## statistics

def utilization (arrv):
    """Returns observed utilization of all queues in arrv."""
    max_d = 0
    tot_time = [0] * arrv.num_queues()
    for evt in arrv:
        max_d = max (evt.d, max_d)
        tot_time[evt.qid] += evt.s
    return [ x / max_d for x in tot_time ]

def mean_wait (arrv):
    """Returns mean waiting time for each queue in arrv."""
    return aggregate_by_q (arrivals.Event.wait, arrv)

def mean_response_time (arrv):
    return aggregate_by_q (response_time, arrv)

def std_response_time (arrv):
    return aggregate_by_q (response_time, arrv, aggregator=numpy.std)

def mean_service (arrv):
    return aggregate_by_q (lambda evt: evt.s, arrv)

def std_service (arrv):
    return aggregate_by_q (lambda evt: evt.s, arrv, aggregator=numpy.std)

def max_service (arrv):
    return aggregate_by_q (lambda evt: evt.s, arrv, aggregator=max)

def mean_obs_response (arrv):
    return aggregate_by_q (response_time, arrv, filter=is_obs)

def std_obs_response (arrv):
    return aggregate_by_q (response_time, arrv, aggregator=numpy.std, filter=is_obs)

def mean_obs_service (arrv):
    return aggregate_by_q (lambda evt: evt.s, arrv, filter=is_obs)

def mean_latent_service (arrv):
    return aggregate_by_q (lambda evt: evt.s, arrv, filter=lambda evt: not is_obs(evt))

def nobs_by_queue (arrv):
    return aggregate_by_q (lambda evt: 1, arrv, filter=is_obs, aggregator=numpy.sum)

def std_obs_service (arrv):
    return aggregate_by_q (lambda evt: evt.s, arrv, aggregator=numpy.std, filter=is_obs)

def mean_wait_percentage (arrv):
    """Returns list, of for each Q, average over evenst of wait time / response time."""
    def wait_pct (evt):
        r = evt.d - evt.a
        if r < 1e-50:
            return 0.0
        else:
            return 1.0 - (evt.s / (evt.d - evt.a))
    return aggregate_by_q (wait_pct, arrv)


def _inner_qsize_for_qid (arrv, qi):
        evts = arrv.events_of_qid (qi)
        times = [ [e.a, 1, e] for e in evts ]
        times.extend ( [e.d, -1, e] for e in evts )
        times.sort()
        
        cum = 0
        for tup in times:
            tup.append (cum)
            cum += tup[1]

        return times

def all_qsize_for_qid (arrv, qi):
    return [ (t,N,evt) for t,foo,evt,N in _inner_qsize_for_qid(arrv,qi)  ]

def mean_qsize_for_qid (arrv, qi):
    times = _inner_qsize_for_qid (arrv, qi)
    return numpy.mean([ tup[3] for tup in times if tup[1] == 1 ])

def mean_qsize_arriving (arrv):
    """Mean queue length as seen by an arriving customer.  Return list with one value for each q in network."""
    return [ mean_qsize_for_qid(arrv,qi) for qi in xrange(arrv.num_queues()) ]
        
def mean_task_time (arrv):
    def time_of_evtl (evtl):
        A = min(e.d for e in evtl)
        D = max(e.d for e in evtl)
        return D-A
    return aggregate_by_task (time_of_evtl, arrv)

def mean_task_time_in_range (arrv, l, r):
    def time_of_evtl (evtl):
        A = min(e.d for e in evtl)
        D = max(e.d for e in evtl)
        return D-A
    def task_arrival_in_range (evtl):
        A = min(e.d for e in evtl)
        return ((l <= A) and (A < r))
    return aggregate_by_task (time_of_evtl, arrv, filter=task_arrival_in_range)

def num_tasks_in_arrival_range (arrv, l, r):
    def task_arrival_in_range (evtl):
        A = min(e.d for e in evtl)
        return ((l <= A) and (A < r))
    return aggregate_by_task (lambda x: 1, arrv, aggregator=len, filter=task_arrival_in_range)

def mean_task_service (arrv):
    def service_of_task (evtl):
        return numpy.sum (e.s for e in evtl)
    return aggregate_by_task (service_of_task, arrv)

def mean_task_length (arrv):
    return aggregate_by_task (len, arrv)
    
def total_service (arrv):
    return numpy.sum (evt.s for evt in arrv)

## output

def write_arrv (f, arrv): 
    evt2n = dict()
    for qi in range(arrv.num_queues()):
        times = _inner_qsize_for_qid (arrv, qi)
        for tup in times:
            if tup[1] == 1:
                evt2n[tup[2]] = tup[3]
    for evt in arrv:
        f.write (evt.as_csv())
        f.write (" ")
        f.write (str(evt2n[evt]))
        f.write ("\n")
    
    
## "exported" stats

STATS = [ \
    utilization, mean_wait, mean_response_time, std_response_time, \
    mean_service, std_service,
    mean_obs_response, std_obs_response,
    mean_obs_service, std_obs_service,
    mean_latent_service, nobs_by_queue,
    mean_wait_percentage, mean_qsize_arriving, mean_task_time \
  ] 

def all_stats_string ():
    return ' '.join([ fn.__name__ for fn in STATS ])

def stringify (val):
    if isinstance(val, list):
        return " ".join (map (str, val))
    else:
        return str(val)


## main method: Run given stats on a given arrivals object

import sys
import re
from optparse import OptionParser

import qnetu

def main():

    parser = OptionParser(description="Compute usage statistics from tasks in a queueing network.")
    parser.add_option("--network", dest="netf",
                      help="YAML file containing network description", metavar="FILE")
    parser.add_option("--arrivals", dest="arrvf",
                      help="YAML file containing job arrivals", metavar="FILE")
    parser.add_option("--use-multif", dest="useMultif", type="int",
                      help="If true, use multif format.")
    parser.add_option("--stats", "--statistics", dest="statsl",
                      help="""List of statistics to report after every iteration.
                          Example: --stat "mean_wait, mean_qsize_arriving"
         Available statistics: %s""" % all_stats_string(), metavar="LIST")
    parser.add_option ("--prefix", dest="prefix",
                       help="Prefix for constructing output files.")

    (opts, args) = parser.parse_args()
    if len(args) > 0:
        parser.error ("Invalid option string %s" % args[0])
        sys.exit(1)

    f = open(opts.netf)
    net = qnetu.qnet_from_text (f)
    f.close()
    
    if opts.useMultif:
        arrv = qnetu.read_multif_of_prefix (opts.arrvf, net)
    else:
        f = open(opts.arrvf)
        arrv = qnetu.load_arrivals (net, f)
        f.close()
 
    if opts.statsl == "ALL":
        stat_fns = STATS[:]
    else:
        comma = re.compile ("[, ]+")
        stat_fns = [ globals()[fn_name] for fn_name in comma.split(opts.statsl) ]
    
    for fn in stat_fns:
        f = open ("%s_%s.txt" % (opts.prefix, fn.__name__), "w")
        f.write (stringify (fn (arrv)))
        f.write ("\n")
        f.close()

if __name__ == "__main__": main()


