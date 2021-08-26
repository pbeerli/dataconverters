#!/usr/bin/env python
# read a modern migrate file
#
# this is a module for reuse
# (c) Peter Beerli 2016,2021
#
#
import sys

def reader(filename):
    with open(filename, 'r') as infile:
        tmp = infile.read() 
    data = tmp.splitlines()
    #print(data)
    return data

def readpop(data,pop):
    #print("@data",data)
    numind, poptitle = data[0].split(None,1)
    print (numind,poptitle)
    p = []
    for i in range(1,int(numind)+1):
        #print i,data[i][:10]
        p.append([data[i][:10],data[i][10:].strip()])
    #p = [[data[i][:10],data[i][10:].strip()] for i in range(1,int(numind))]
    return int(numind), p, poptitle

def get_header(data):
    newdata = [d for d in data if d[0]!='#'] 
    header = newdata[0]
    tmp = header.split()
    numpop, loci = int(tmp[0]),int(tmp[1])
    return numpop,loci

def split_migrate1(data):
    newdata = [d for d in data if d[0]!='#'] 
    header = newdata[0]
    tmp = header.split()
    numpop, loci = int(tmp[0]),int(tmp[1])
    if len(tmp)>2:
        title = "".join(tmp[2:])
    else:
        title = " "
    sitesline = newdata[1]
    #print (numpop,'@',loci,'@',title)
    #print (sitesline)
    if sitesline[0] != '(':
        print ('failed')
        exit(-1)
    populations = []
    seqloci = []*loci
    start = 2
    for pop in range(numpop):
        numind, pp, poptitle = readpop(newdata[start:],pop)
        populations.append(pp)
        start += numind+1
    return populations,title

def split_migrate(data):
    newdata = [d for d in data if d[0]!='#'] 
    header = newdata[0]
    tmp = header.split()
    numpop, loci = int(tmp[0]),int(tmp[1])
    if len(tmp)>2:
        title = "".join(tmp[2:])
    else:
        title = " "
    sitesline = newdata[1]
    #print numpop,'@',loci,'@',title
    #print sitesline
    if sitesline[0] != '(':
        print ('failed')
        exit(-1)
    populations = []
    seqloci = []*loci
    start = 2
    locations=[]
    for pop in range(numpop):
        numind, pp, poptitle = readpop(newdata[start:],pop)
        populations.append(pp)
        locations.append(poptitle)
        start += numind+1
    return populations,title,locations


if __name__ == '__main__':
    data = reader(sys.argv[1])
    #print data
    populations,title,locations = split_migrate(data)
    print (populations)
    print(locations)

