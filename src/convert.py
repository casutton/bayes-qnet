#!/usr/bin/python
#  Front end to predict 


from optparse import OptionParser

import sampling
import qnetu
import qstats

import numpy

MULTI = 0
SINGLE = 1
TASK = 2
NQ = 3

def main():
    parser = OptionParser(description="Front-end for greedy bottleneck finding.")
    parser.add_option("--network", "-n", dest="netf",
                      help="YAML file describing network governing the data", metavar="FILE")
    parser.add_option("--input", dest="inf", help="Input arrivals")
    parser.add_option("--output", dest="outf", help="Output arrivals")
    parser.add_option("--multi2single", dest="m2s", action="store_true", default=False)
    parser.add_option("--multi2tasks", dest="m2t", action="store_true", default=False)
    parser.add_option("--multi2nq", dest="m2nq", type="int", default=-1)
    parser.add_option("--single2multi", dest="s2m", action="store_true", default=False)
    parser.add_option("--single2tasks", dest="s2t", action="store_true", default=False)

    (options, args) = parser.parse_args()    

    if options.m2s:
        in_format = MULTI
        out_format = SINGLE
    elif options.s2m:
        in_format = SINGLE
        out_format = MULTI
    elif options.m2t:
        in_format = MULTI
        out_format = TASK
    elif options.s2t:
        in_format = SINGLE
        out_format = TASK
    elif options.m2nq >= 0:
        in_format = MULTI
        out_format = NQ
    else:
        print "Specify --multif2single or --single2multif"
        sys.exit (1)


    f = open(options.netf)
    net = qnetu.qnet_from_text (f)
    f.close()

    arrv = ((in_format == MULTI and qnetu.read_multif_of_prefix (options.inf, net)) or
            (in_format == SINGLE and qnetu.read_from_table (net, options.inf)))
        
    if out_format == MULTI: qnetu.write_multif_to_prefix (options.outf, arrv)
    elif out_format == SINGLE: qnetu.write_to_table (options.outf, arrv)
    elif out_format == TASK: qnetu.write_tasks (options.outf, arrv)
    elif out_format == NQ: write_nq (options.outf, options.m2nq, arrv)

def write_nq (outf, qi, arrv):
    f = open (outf,"w")
    f.write ("T EVT N\n")
    for t,N,evt in qstats.all_qsize_for_qid (arrv, qi):
        f.write ("%.4f %d %d\n" % (t,evt.eid,N))
    f.close()

if __name__ == "__main__":
    main()
