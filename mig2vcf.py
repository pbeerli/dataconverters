#!/usr/bin/env python3
#
import os
import sys
import re
import datetime

def readline(f):
    x = f.readline()
    while x[0]=='#':
        x = f.readline()
    return x

def readmigline(f):
    x = readline(f).strip()
    while x[0]=='#':
        x = readline(f).strip()
    #print(x)
    return x

def oldstyle(line):
    #instructions=[loci,siteslist]
    l = line.split()
    loci = len(l)
    siteslist = [int(li) for li in l]
    return [loci,siteslist]

def newstyle(line):
    # separate loci that are in different blocks
    line = line.replace('s','')
    commas = line.split(",")
    num_comma = len(commas)
    # split within block
    loci = [b.count("(") for b in commas]
    block = [b.split('(') for b in commas]
    block = [[b1.replace(')','') for b1 in b]  for b in block if b]
    block = [list(filter(None,b)) for b in block]
    block= [[sum([int(ciii) for ciii in cii.split()]) for cii in ci if cii!='' and cii!=' '] for ci in block] 
    return [loci,block]

def read_migratedata(file):
    f = open(file,'r')
    # read the first line pop loci title
    titleline = readmigline(f)
    titleline = titleline.split()
    pop, loci = [int(i) for i in titleline[:2]]
    # read second line: we only do sequences here the could come in two
    # different flavors
    #oldstyle: 340 450 ....  each locus with the #sites, all loci are in
    # different blocks
    #newstyle: (s100 s200) (s300), (s5000) linked (in same () and unlinked ()()
    # and in different blocks ',' with ',' all loci are on the same line!
    line = readmigline(f)
    if "(" in line:
        instructions = newstyle(line)
    else:
        instructions = oldstyle(line)
    #print("instructions",instructions)
    popsequences = []
    names=[]
    for p in range(pop):
        # read the population line this could be many numbers and title
        titleline = readmigline(f).split()
        #print(titleline)
        numind = int(titleline[0])
        #print(numind)
        #sys.exit(-1)
        sequences=[]
        blocks = len(instructions[0])
        #print("blocks",blocks)
        #print("instructions:",instructions)
        #print("numind:", numind)
        for locus in range(blocks):
            nn=[]
            #ss = [[] for i in range(len(instructions[1][locus]))]
            #print("@", instructions[0][locus],len(instructions[1][locus]),instructions[1][locus])
            ss=[]
            for i in range(numind):
                #print ind,locus
                a = readmigline(f)
                nn.append(a[:10])    # individual name
                seq = a[10:].strip() # total sequence in block, could > 1 loci
                #print (locus, seq)
                start = 0
                stop = 0
                #print("reading seq:", start,stop)
                subloci = instructions[1][locus]
                ss2=[]
                for z,j in enumerate(subloci):
                    start = stop
                    stop += j
                    ss2.append(seq[start:stop]) 
                    #print("startstop",start,stop)
                ss.append(ss2)
            xss = list(zip(*ss))
            #for j in range(len(instructions[1][locus])):
            for xssi in xss:
                sequences.append(xssi)
        popsequences.append(sequences)
        names.append(nn)
    #print("@NAMES",names)
    #print(len(popsequences),end=' | ')
    #for s in popsequences:
    #    print(len(s),end=';')
    #    for si in s:
    #        print(len(si),end=' ')
    #    print()
    #print()
    return pop, instructions[0], instructions[1], popsequences, names

def chunks(n,l):
    size = len(l)
    partsize = int(size/n)
    if size > partsize*n:
        lastpartsize = partsize + size-(partsize*n)
    else:
        lastpartsize = partsize
    starters = list(range(0,size,partsize))
    stoppers = list(range(partsize,size,partsize))
    #stoppers.extend(size)
    #print(starters)
    #print(stoppers)
    return [l[start:stop] for start,stop in zip(starters,stoppers)]
        
def to_matrix(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]    


def snps(numpop, sequences):
    mysnps=[] 
    #sequences = [[1,2,3],[1,2,3]]
    #to
    #sequences = [1,1,2,2,3,3]
    #and then [1,2,3]
    x=list(zip(*sequences)) #  [[1,2,3],[1,2,3]] --> [[1,1],[2,2],[3,3]]
    #x = sequences
    y=[]
    for xi in x: #  [[1,1],[2,2],[3,3]] -> [[1],[2],[3]]
        yt=[]
        for xii in xi:
            #print("#",xii)
            yt.extend(xii)
        y.append(yt)
    #print("y:",len(y),len(y[0]),len(y[1]))
    count=-1
    refseq=[]
    
    for yi in y:
        count+=1
        refseq.append(yi[0])
        #print("in snps",len(yi),len(yi[0]))
        ###print(len(yi),f"{[yii[0] for yii in yi]}")
        news = list(map(list, zip(*yi)))
        #print(news)
        for i,ni in enumerate(news):
            nis = sorted(list(set(ni)))
            #print("@", i,nis,"".join(ni))
            if len(nis)>1:
                if nis[0] == refseq[-1][i]:
                    if len(nis)>=3:
                        u = (count, i, nis[0], nis[1]+","+nis[2])
                    else:
                        u = (count, i, nis[0], nis[1])
                elif nis[1] == refseq[-1][i]:
                    if len(nis)>=3:
                        u = (count, i, nis[1], nis[0]+","+nis[2])
                    else:
                        u = (count, i, nis[1], nis[0])
                else:
                    if len(nis)>=3:
                        if nis[2]==refseq[-1][i]:
                            u = (count, i, nis[2], nis[0]+","+nis[1])
                #print(f"snps: {u}")
                mysnps.append(u)
    return mysnps,refseq,y

def vcfwriter(vcfdata,vcfname,reffile,sampledata,names, ploidy):
#    ##fileformat=VCFv4.4
#    ##fileDate=20090805
#    ##source=myImputationProgramV3.1
#    ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
#    ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
#    ##phasing=partial
#    ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
#    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#    ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
#    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
#    ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
#    ##FILTER=<ID=q10,Description="Quality below 10">
#    ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
#    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#    ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#    #CHROM POS     ID        REF    ALT     QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
#    20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
#    20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
#    20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
#    20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
#    20     1234567 microsat1 GTC    G,GTCT  50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3
    if ".vcf" not in vcfname:
        vcfname = vcfname+".vcf"
    f = open(vcfname,'w')
    software="mig2vcf.py"
    today = datetime.date.today()
    totalsize=0
    samplesize = []
    contigs = ""
    for sis,s in enumerate(sampledata):
        totalsize += len(s[0])
        contigs +=f"##contig=<ID={1+sis},length={len(s[0])},assembly=1,md5=none,species=\"x\",taxonomy=x>\n"

        samplesize.append(len(s))
        continue
    #print("@@@@@", samplesize)
    header = f'##fileformat=VCFv4.2\n\
##fileDate={today}\n\
##source={software}\n\
##reference=file://{reffile}\n\
{contigs}\
##phasing=all\n\
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
    f.write(header)
    #print("@",names)
    newnames=[]
    newnames = [ni.strip() for nn in names for ni in nn][:samplesize[0]]
    #print("------",newnames)
    ploidy = int(ploidy)
    stride = ploidy
    if ploidy>1:
        if ":" in newnames[0]:
            newnames=[ni.split(":")[0] for ni in newnames]
        if samplesize[0] % 2 == 0:
            ss = samplesize[0]-1
        else:
            ss = samplesize[0]-1
    else:
        ss = samplesize[0]
    #print("samplesize",ss)
    for i in range(0,ss,stride):
        if i<ss-1:
            tab='\t'
        else:
            tab='\n'
        z = f"{newnames[i]}{tab}"
        #print(z,newnames[i])
        f.write(z)
        #print(z,end="")
    #z = f"{newnames[i]}\n"
    #f.write(z)
    #print(z,end="")
    oldchrom = 0
    count=0
    for line in vcfdata:
        #print("++++",line)
        tab = '\t'
        chrom, pos, ref, alt = line
        if oldchrom!=chrom:
            count = 0
            oldchrom = chrom
        count += 1
        #f.write(str(chrom)+"\t"+str(pos+1)+f"\t{count}\t"+ref+"\t"+alt+"\tPASS\t.\tNS="+str(samplesize[chrom-1])+"\tGT\t")
        f.write(f"chr{chrom}\t{pos}\t{count}\t{ref}\t{alt}\t50\tPASS\tNS={samplesize[chrom-1]}\tGT\t")
        # deal with the individuals
        #print(chrom,len(sampledata[chrom-1]))
#        print("@@",list(range(0,samplesize[chrom-1],ploidy)))
        for sa in range(0,samplesize[chrom],ploidy):
            #print("@",sa)
            if samplesize[chrom]-ploidy == sa:
                tab=f"\n"
            else:
                tab=f"\t"
            pa = sa
            #print(sa,ploidy)
            for p in range(ploidy):
                if p==ploidy-1:
                    delim=''
                else:
                    delim='|'
                #print(f"xxxxx {chrom} {pa} {pos}    {ref} {delim}")
                if sampledata[chrom][pa][pos]==ref+delim:
                    f.write(str(0))
                else:
                    if "," in alt:
                        if sampledata[chrom][pa][pos]==alt[0]:
                            f.write(str(1)+delim)
                        else:
                            f.write(str(2)+delim)
                    else:
                        f.write(str(1)+delim)
                pa += 1
            f.write(tab)
        #f.write("@\n")

def flatten(L):
    for item in L:
        try:
            yield from flatten(item)
        except TypeError:
            yield item
            
def refseqwriter(refname, refseq, title, loci):
    if type(loci) is list:
        #print(loci)
        ll=sum(flatten(loci))
    else:
        ll=int(loci)
    if ll>1:
        # the outcommented section delivers n loci fasta files
        #for l in range(ll):
        #    f=open(refname+f"{l+1:04d}.fasta",'w')
        #    f.write("> Locus "+str(l+1)+" "+title+"\n")
        #    f.write(refseq[l])
        #    f.close()
        f=open(refname+".fasta",'w')
        for l in range(ll):
            f.write("> Locus "+str(l+1)+" "+title+"\n")
            f.write(refseq[l])
            f.write('\n')
        f.close()
    else:
        f=open(refname,'w')
        f.write(">"+title+"\n")
        f.write("".join(refseq))
        f.close()


def parseargs(args):
    if "--help" in args or "-h" in args:
        print("Syntax: mig2vcf -i migratedatafile -o vcfoutfile -d ploidy")
        print("        [this creates also fasta reference files for each locus for the VCF file")
        sys.exit(-1)

    da = {'-p': 1, '-o':'temp.vcf', '-i':'infile'}
    args.pop(0)
    for i,a in enumerate(args):
        if a[0]=='-':
            da[a] = args[i+1]
    try:
        return [da['-i'], da['-o'], da['-o'].replace(".vcf",".fasta"),da['-p']]
    except:
        print("Syntax: mig2vcf -i migratedatafile -o vcfoutfile -d ploidy")
        print("        [this creates also fasta reference files for each locus for the VCF file")
        sys.exit(-1)

if __name__ == '__main__':
    
    datafile, vcffile, reffile, diploid = parseargs(sys.argv)
    try:
        print("Input:    ", datafile, file=sys.stderr)
        print("Output:   ", vcffile, file=sys.stderr)
        print("Reference:", reffile, file=sys.stderr)
        print("Ploidy:   ", diploid, file=sys.stderr)
    except:
        print("Syntax: mig2vcf -i migratedatafile -o vcfoutfile -p ploidy")
        print("        [this creates also fasta reference files for each locus for the VCF file")
        sys.exit(-1)

    numpop, loci, sites, sequences, names = read_migratedata(datafile)
    #print(f"{numpop} {loci} {sites} {len(sequences)}")
    vcfdata, refseq, sampledata = snps(numpop,sequences)
    vcfname=vcffile
    #print(vcfdata)
    vcfwriter(vcfdata,vcfname,reffile,sampledata,names, diploid)
    os.system(f'bgzip -f {vcfname}; tabix {vcfname}.gz')
    #print('done')
    title=f"Reference [{names[0][0].strip()}] for VCF file: "+vcfname
    refseqwriter(reffile, refseq,title,loci)
    #sys.exit(0)
