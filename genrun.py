#! /usr/bin/env python

# This Python script generates a shell script which runs the 
# executables for the DeformationField128 project. 
# Some of the parameters are specific to the project. It can be altered to meet the needs of 
# other projects. 

import sys, string

if len(sys.argv) != 4: 
    print 'Usage: genrun.py <start-file-num> <end-file-num> <output-file-name>' 
    sys.exit(1) 

lbnd = int(sys.argv[1])
ubnd = int(sys.argv[2])
ofname = sys.argv[3]

strubnd = ' ' + str(ubnd) + ' '
strlbnd = ' ' + str(lbnd) + ' '
stride = 5

npnt = ' 2097152 '
dt = ' 0.01225 '
nres = ' 128 '

runlst = range(lbnd+stride, ubnd+stride, stride)

fid = open(ofname, 'w')
for x in runlst:
    prfx = 'xyzbw' + str(x) + 'to' + str(lbnd) + '-'
    nstep = ' ' + str(x - lbnd) + ' '

    str1 = './gridxyz-dp.x ' + nres + str(x) + ' ' + prfx
    str2 = './particletracking-dp.x' + nres + str(x) + nstep + '0' + npnt + dt + '0' + ' ' + prfx
    str3 = './bijongrid-dp.x' + nres + dt + strlbnd + str(x) + ' ' + prfx
    str4 = './rmfiles.py ' + str(lbnd+1) + ' ' + str(x-1) + ' 1 ' \
                           + './out/' + prfx + ' ' \
                           + '.dat'

    fid.write(str1 + '\n')
    fid.write(str2 + '\n')
    fid.write(str3 + '\n')
    fid.write(str4 + '\n')
    fid.write('\n')

fid.close()

