#!/usr/bin/python -u

import sys
s = 0
n = 0
i = int(sys.argv[1])
minv = 1e+30
maxv = 0

offset = 0
if len(sys.argv) == 4:
    offset = int(sys.argv[3])
for line in sys.stdin:
    v = float(line.split()[i][offset:])
    s += v
    if minv > v: minv = v
    if maxv < v: maxv = v
    n += 1

if n > 0: 
    if s > 0: perc = 100.0 * (maxv - minv) / (s / n)
    else: perc = 0
    if len(sys.argv) == 3:
        if sys.argv[2] == 'e':
            print "%.3e %.3e %.3e %d %d" % (s / n, minv, maxv, n, perc)
        elif sys.argv[2] == '2':
            print "%.2f %.2f %.2f %d" % (s / n, minv, maxv, n)
        elif sys.argv[2] == '3':
            print "%.3f %.3f %.3f %d" % (s / n, minv, maxv, n)
        elif sys.argv[2] == '4':
            print "%.4f %.4f %.4f %d" % (s / n, minv, maxv, n)
    else:    
        print s / n, s, minv, maxv, n
