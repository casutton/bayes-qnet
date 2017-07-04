#!/usr/bin/python

import sys
import re

def main():
    dotf,start,to = sys.argv[1:]

    graph = {}

    exp = re.compile ("(\w+) -> (\w+)")
    f = open (dotf)
    for line in f:
        match = exp.search(line)
        if match:
            add_arc (graph, match.group(1), match.group(2))
    f.close()

    path = find_path (graph, start, to, [], {})
    if path:
        path.reverse()
        path.insert (0, start)
        path.insert (-1, to)
    print path

def add_arc (graph, start, to):
    if not start in graph:
        graph[start] = []
    if not to in graph:
        graph[to] = []
    graph[start].append (to)
    
def find_path (graph, start, to, the_path, traversed):
    for child in graph[start]:
        if child == to:
            return the_path
        if not child in traversed:
            new_path = [child]
            new_path.extend (the_path)
            new_traversed = dict(traversed)
            new_traversed[start] = True
            path = find_path (graph, child, to, new_path, new_traversed)
            if path: return path
    return None
    
main()
