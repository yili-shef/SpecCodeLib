#! /usr/bin/env python
import os, sys

# This Python script removes a list of files which have same prefix and suffix and
# are labeled by a sequence of numbers with a constant stride.

if len(sys.argv) != 6: 
    print 'Usage: rmfiles.py <start-#> <end-#> <stride> <prefix-of-file-name> <suffix-of-file-name>' 
    sys.exit(1) 

lbnd = int(sys.argv[1])
ubnd = int(sys.argv[2])
stride = int(sys.argv[3])
prfx = sys.argv[4]
sffx = sys.argv[5]
flst = range(lbnd, ubnd+stride, stride)
for x in flst:
    fnm = 'rm -f ' + prfx + str(x) + sffx
    os.system(fnm)

sys.exit(0)
