
import yaml

import numpy
import numpy.random
import hmm
import distributions
import qnet
import arrivals
import queues
import estimation

import sampling 
import utils

def qnet_from_fname (fname):
    f = open (fname)
    net = qnet_from_text (f)
    f.close()
    return net

def qnet_from_text (yaml_text):
    ylist = yaml.load (yaml_text)
    
    qname2id = {}
    sname2id = {}
    snames = []
    qnames = []
    s_succ = {}
    s_q = {}
    
    def add_state (name):
        if not name in sname2id:
            snames.append (name)
            sname2id[name] = len(sname2id)
    
    
    def add_queue (name):
        if not name in qname2id:
            qnames.append (name)
            qname2id[name] = len(qname2id)

    # error checking
    if not isinstance (ylist, dict):
        raise Exception ("Bad YAML file.  Needed a dictionary:\n%s" % ylist)

    # read state structure
    slist = ylist['states']
    for s in slist:
        add_state (s['name'])
        sid = sname2id.get (s['name'])

        if 'successors' in s: 
            map (add_state, s['successors'])
            s_succ[sid] = map (sname2id.get, s['successors'])
            
        map (add_queue, s['queues'])
        s_q[sid] = map (qname2id.get, s['queues'])
        
    # create HMM
    ns = len(snames)
    nq = len(qnames)
    T = numpy.zeros (( ns, ns ))  # transition matrix
    O = numpy.zeros (( ns, nq ))  # observation matrix
    
    # populate transition matrix
    for sid,succl in s_succ.iteritems():
        for succ in succl:
            T[sid,succ] = 1./len(succl)  # assume uniform

    # populate emission matrix
    for sid,ql in s_q.iteritems():
        for qid in ql:
            O[sid,qid] = 1./len(ql)  # assume uniform
            
    fsm = hmm.HMM (T, O)
    
    # HMM done. read queue information
    qs = []
    factors = []
    universe = {}
    
    for qhash in ylist['queues']:
        q = queue_from_yhash (qhash, universe)
        qs.append (q)
        factors.append (q.service)
        
    # construct qnet object
    net = qnet.Qnet(qs, factors, universe, qname2id, sname2id, fsm)
    
    return net

def queue_from_text (text, universe):
    return queue_from_yhash (yaml.load (text), universe)

def queue_from_yhash (qhash, universe):
    qname = qhash['name']
    sdist = _tmpl_from_yaml (qhash['service'], qname, universe)
    num_procs = _intuit_procs (qhash)
    
    if 'type' in qhash:
        qtype = qhash['type']
    elif num_procs == 1:
        qtype = "GG1"
    else:
        qtype = "GGk"
    qtype = qtype.upper()    
        
    if qtype == "GG1":
        q = queues.QueueGGk (qname, sdist, 1)
    elif qtype == "GGK":
        q = queues.QueueGGk (qname, sdist, num_procs)
    elif qtype == "DELAY":
        q = queues.DelayStation (qname, sdist)
    elif qtype == "GG1R":
        q = queues.QueueGG1R (qname, sdist)
    elif qtype == "PS":
        q = queues.QueuePS (qname, sdist)
    else:
        raise Exception ("Could not understand queue %s type %s" % (qname, qtype))

    return q

def _intuit_procs (tag):
    if 'processors' in tag.keys():
        return tag['processors']
    else:
        return 1
        
def _parse_args (tag):
    if isinstance (tag,list):
        dist_type = tag[0]
        args = tag[1:]
    else:
        dist_type = tag
        args = None
    return dist_type, args
    
def _tmpl_from_yaml (tag, qname, universe):
    dist_type, args = _parse_args (tag)
    
    if dist_type == "MIX":
        i = 0
        ws = []
        dists = []
        num_components = 0
        while i < len(args):
            ws.append (args[i])
            subtype, subargs = _parse_args (args[i+1])
            dists.append(_dist_from_yaml(subtype, subargs))
            num_components += 1
            i += 2
        var_name = "COMPONENT_%s" % qname
        universe[var_name] = range (num_components)
        return queues.MixtureTemplate (var_name, ws, dists)        
    else:
        return queues.DistributionTemplate (_dist_from_yaml (dist_type, args))


def _dist_from_yaml (dist_type, args):
    if dist_type == 'M':
        scale = args[0] if args else 1.0
        return distributions.Exponential(scale)
    if dist_type == "G":
        if isinstance (args, list):
            shape, scale = args
        elif isinstance (args, dict):
            shape = args['shape']
            scale = args['scale']
        else:
            shape = scale = 1.0
        return distributions.Gamma (shape=shape, scale=scale)
    if dist_type == "LN":
        if isinstance (args, list):
            mu, sig = args
        elif isinstance (args, dict):
            mu = args['mu']
            sig = args['sigma']
        else:
            mu = 0
            sig = 1
        return distributions.LogNormal (mu, sig)
    if dist_type == 'D':
        return distributions.Dirac (args[0])
    else:
        raise Exception ("Could not understand tag " % tag)


def arrv_from_yaml (net, yaml_text):
   ylist = yaml.load (yaml_text)    
   arrv = arrivals.Arrivals(net)
   eid = 0
   
   for l in ylist:
     try: 
        tid = int (l['task'])
        q_name = l['queue']
        q = net.queue_by_name (q_name)
        if q is None:
            raise Exception ("Error in YAML file: Could not find queue $%s$" % q_name)
        qid = net.qid_of_queue (q)
        
        sid = net.sid_by_name (l['state']) if ('state' in l) else -1
            
        a = float (l['arrival'])
        d = float (l['departure'])
        s = 0
        
        arrv.append (arrivals.Event (eid, tid, qid, q, a, s, d, sid, obs_a=1, obs_d=1))
        
        eid = eid + 1
        
     except TypeError, e:
       print "Error reading ", l
       raise e

   arrv.recompute_all_service()

   return arrv
    
def _render_event (evt):
    return { 'task': evt.tid, 'arrival': evt.a, 'departure': evt.d, 'queue': evt.queue().name, 'state': evt.state, 'obs_a': evt.obs_a, 'obs_d': evt.obs_d }
    
    
# interface with YAML package

class QnetLoader (yaml.Loader):
    def __init__ (self, stream, net): 
        yaml.Loader.__init__ (self, stream)
        self.net = net
    def qnet(self): return self.net
    
def event_representer (dumper, evt):
    return dumper.represent_mapping (u'!Event', _render_event(evt))
    
def arrivals_representer (dumper, arrv):
    return dumper.represent_mapping (u'!Arrivals', { 'events': arrv.all_events() })

def event_constructor (loader, node):
    # Create the actual event later.  Just return hash
    mapping = loader.construct_mapping (node)
    return mapping
    
def arrivals_constructor (loader, node):
    m = loader.construct_mapping (node, deep=True)
    evt_dict_lst = m['events']
    
    net = loader.qnet()
    arrv = arrivals.Arrivals (net)
    eid = 0
    
    evt_out = []
    
    for l in evt_dict_lst:
      try:
        tid = int (l['task'])
        q_name = l['queue']
        q = net.queue_by_name (q_name)
        if q is None:
            raise Exception ("Error in YAML file: Could not find queue $%s$" % q_name)
        qid = net.qid_of_queue (q)
        
        if isinstance (l['state'], int):
            sid = l['state']
        else:
            sid = net.sid_by_name (l['state']) if ('state' in l) else -1
            
        a = float (l['arrival'])
        d = float (l['departure'])
        s = 0
        
        obs_a = l['obs_a'] if ('obs_a' in l) else 1
        obs_d = l['obs_d'] if ('obs_d' in l) else 1

        evt = arrivals.Event (eid, tid, qid, q, a, s, d, sid, obs_a=obs_a, obs_d=obs_d)
        evt_out.append (evt)
        
        eid = eid + 1
      except TypeError, e:  
        print "Error reading arrival from ", l
	raise e

    # Make to sure to add events in the order of their arrival
    evt_out.sort (key=lambda e: (e.a, e.d))
    for e in evt_out: 
        arrv.append (e)
    
    arrv.recompute_all_service()
    arrv.validate()

    # Remove non-FIFO tasks
    tasks_deleted = dict()
    for e in evt_out:
        if not (e.tid in tasks_deleted):
            if e.s < 0:                
                prev = e.previous_by_queue()
                print "WARNING: Deleting non-fifo event\n  %s\n  %s\n  %s" % (e.previous_by_queue(), e, e.next_by_queue())
                arrv.delete_task (prev.tid) 
                tasks_deleted[prev.tid] = True
                if e.s < 0:  # if deleting previous task didn't help, delete this guy.
                    print "WARNING: Also deleting e:\n  %s\n  %s" % (e.previous_by_queue(), e)
                    arrv.delete_task (e.tid) 
                    tasks_deleted[e.tid] = True
                    
    arrv.validate()
    
    return arrv
    
yaml.add_representer (arrivals.Event, event_representer)
yaml.add_constructor (u"!Event", event_constructor)

yaml.add_representer (arrivals.Arrivals, arrivals_representer)
yaml.add_constructor (u"!Arrivals", arrivals_constructor)



def load_arrivals (net, stream):
    loader = QnetLoader (stream, net)
    if loader.check_data():
        return loader.get_data()

#  Alternate arrivals format

def write_multif_to_prefix (prefix, arrv, obs_only=True):
    statef = open ("%sstate.txt" % prefix, "w")
    af = open ("%sa.txt" % prefix, "w")
    df = open ("%sd.txt" % prefix, "w")
    arrv = write_multifile_arrv (arrv, statef, af, df, obs_only=obs_only)
    statef.close()
    af.close()
    df.close()

def read_multif_of_prefix (prefix, net):
    statef = open ("%sstate.txt" % prefix)
    af = open ("%sa.txt" % prefix)
    df = open ("%sd.txt" % prefix)
    arrv = read_multifile_arrv (net, statef, af, df)
    statef.close()
    af.close()
    df.close()
    return arrv

def write_multifile_arrv (arrv, statef, af, df, obs_only=True):
    for evt in arrv:
        if (not obs_only) or (evt.obs_a or evt.obs_d):
            next_id = evt.next_by_task().eid if evt.next_by_task() else -1
            statef.write ("%d %d %s %d %d\n" % (evt.tid, evt.eid, evt.queue().name, evt.state, next_id))

    for evt in arrv:
        if (not obs_only) or evt.obs_a:
            af.write ("%d %.15f %d\n" % (evt.eid, evt.a, evt.obs_a))

    for evt in arrv:
        if (not obs_only) or evt.obs_d:
            df.write ("%d %.15f %d\n" % (evt.eid, evt.d, evt.obs_d))
                     

def read_multifile_arrv (net, statef, af, df):
#    statef = fileify (statef)
#    af = fileify (af)
#    df = fileify (df)

    e2a = dict()
    e2d = dict()
    e2nxt = dict()

    evts = []

    for line in af.readlines():
        if line != "\n":
            fields = line.split()
            if len(fields) == 2: fields.append(1) # default OBS value
            eid,a,obs = fields
            eid = int (eid)
            e2a[eid] = (float(a), int(obs))

    for line in df.readlines():
        if line != "\n":
            fields = line.split()
            if len(fields) == 2: fields.append(1) # default OBS value
            eid,d,obs = fields
            eid = int (eid)
            e2d[eid] = (float(d),int(obs))

    for line in statef.readlines():
        if line != "\n":
            tid, eid, qname, sid, next = line.split()
            eid = int(eid)
            tid = int(tid)
            sid = int(sid)
            if int(next) >= 0:
                e2nxt[eid] = int(next)
            else:
                e2nxt[eid] = None
                
            q = net.queue_by_name (qname)
	    if q is None:
              print net.qname2id.keys()
	      raise Exception ("Cannot find queue named %s in\n%s" % (qname, net.as_yaml()))
            qid = net.qid_of_queue (q)

            if not (eid in e2a): raise Exception ("No arrival file entry for event %d" % eid)
            if not (eid in e2d): raise Exception ("No arrival file entry for event %d" % eid)
            
            a, obs_a = e2a.get(eid)
            d, obs_d = e2d.get(eid)

            evts.append ((d, a, qid, tid, eid, q, sid, obs_a, obs_d))

    evts.sort()

    arrv = arrivals.Arrivals(net, cap=2*len(evts))
    for evt in evts:
        d, a, qid, tid, eid, q, sid, obs_a, obs_d = evt
        if a is None:
            a = 0
            obs_a = 0
        if d is None:
            d = 0
            obs_d = 0
        e = arrivals.Event (eid, tid, qid, q, a, 0, d, sid, obs_a=obs_a, obs_d=obs_d)
        arrv.insert (e)

    arrv.recompute_all_service()

    for e,n in e2nxt.iteritems():
        e0 = arrv.event (e)
        e1 = arrv.event (n) if n else None

        e0.set_next_by_task (e1)

        eid0 = e0.eid if e0 else -1
        eid1 = e1.eid if e1 else -1
#        print "DBG: Adding task link %d --> %d" % (eid0, eid1)

        assert e0.next_by_task() == e1
        if e1:
            assert e1.previous_by_task() == e0, "Huh? Couldn't link by task\n  %s\n  %s" % (e0, e1)

    for evt in arrv:
        if evt.qid == 0:
            evt.set_previous_by_task (None)

#     for evt in arrv:
#         prev_byt = evt.previous_by_task()
#         print "++++++\n  %s\n  %s" % (prev_byt, evt)
#         if prev_byt:
#             assert prev_byt.next_by_task() == evt

    arrv.regenerate_all_caches()

    return arrv


def fileify (fname):
    if isinstance (fname, file):
        return fname
    elif isinstance (fname, str):
        return open(str)
    else:
        raise Exception ("What is %s" % fname)

#  Read an arrivals written using arrivals.write_csv

def read_from_table (net, fname):
    f = open(fname)
    arrv = read_csv_from_lines (net, f.readlines())
    f.close()
    return arrv

def read_table_from_string (net, string):
    return read_csv_from_lines (net, string.split("\n"))

def read_csv_from_lines (net, lines):
    evts = []
    for line in lines[1:]:
        fields = line.split()
        if len(fields) == 0: continue
        fi = 0
        eid = int(fields[fi]); fi += 1
        tid = int(fields[fi]); fi += 1
        qid = int(fields[fi]); fi += 1
        a = float(fields[fi]); fi += 1
        d = float(fields[fi]); fi += 1
        s = float(fields[fi]); fi += 1
        fi += 1 # skip w
        state = int(fields[fi]); fi += 1
        proc = int(fields[fi]); fi += 1
        obs_a = int(fields[fi]); fi += 1
        obs_d = int(fields[fi]); fi += 1

        q = net.queues[qid]
        evt = arrivals.Event (eid, tid, qid, q, a, s, d, state, obs_a, obs_d, proc)
        evts.append (evt)

    evts.sort (key=arrivals.Event.arrival)

    arrv = arrivals.Arrivals(net, cap=2*len(evts))
    for evt in evts:
        arrv.insert (evt)
 
    # this handles piddly rounding errors that can occur
    #  in writing and reading
    arrv.recompute_all_service()

    # hack: auxiliaries not initialized.  I'm not sure that they're readable in the file

    return arrv

def write_to_table (fname, arrv):
    f = open(fname, "w")
    arrv.write_csv (f)
    f.close()

def write_tasks (fname, arrv):
    f = open(fname, "w")
    for ti in arrv.all_tasks():
        evts = arrv.events_of_task (ti)
        f.write (utils.as_str (e.d for e in evts))
        f.write ("\n")
    f.close()

#  Configuration

import misc

def read_configuration (stream):
    yml = yaml.load (stream)
    if not isinstance (yml, dict):
        raise Exception ("Could not understand YAML in config file.")
        
    for k in yml.keys():
        if k == "one_d_sampler":
            qnet.set_sampler (yml[k])
        elif k == "debug_integrate":
            misc.DEBUG = yml[k]        
        elif k == "sampling_type":
            estimation.set_sampler (yml[k])
        elif k == "use_rtp":
            qnet.set_use_rtp (int(yml[k]))
        elif k == "proposal":
            prp_names = queues.proposal_names()
            v = yml[k].lower()
            if not v in prp_names:
                raise Exception("Could not understand proposal: %s" % v)
            else:
                queues.set_proposal (prp_names.index(v))        
        elif k == "initializer":
            if yml[k].find('.') < 0:
                value = yml[k]
            else:
                names = yml[k].split('.')
                if len(names) != 2:
                    raise Exception ("Could not understand initializer %s" % yml[k])
                module = __import__(names[0])
                value = getattr(module, names[1])
            qnet.set_initialization_type (value)
        elif k == "slice_sorted":
            v = int(yml[k])
            qnet.set_slice_sorted (v)
        else:
            raise Exception("Could not understand config key %s" % k)
    

# Seeds

def set_seed (seed):
    numpy.random.seed (seed)
    sampling.set_seed (seed)
