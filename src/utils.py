import numpy

def as_str (l):
    return ' '.join (map(str, l))

def randelt (l):
    return l[numpy.random.randint(len(l))]

def append (*l):
    result = []
    for x in l: result.append (x)
    return result

def roll_die (p):
    if sum(p) > 1 + 1e-12:
        raise "Huh? p: %s" % (p,) 
    # Weird numpy thing
    if sum(p) > 1 - 1e-10:
        for i in range(len(p)):
            if p[i] > 1e-10:
                p[i] -= 1e-10
                break
    a = numpy.random.multinomial (1, p)
    return int(numpy.where(a==1)[0])

def delete_all (l, elt):
    ndel = 0
    for x in l:
        if x == elt: ndel+=1
    nrem = len(l) - ndel
    newl = [None] * nrem
    i = 0
    for x in l:
        if x != elt:
            newl[i] = x
            i += 1
    return newl
