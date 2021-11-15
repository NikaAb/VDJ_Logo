
from scipy import stats
from IO import Reader
import numpy as np
import itertools
import math
from CM import online_CM


class Stattest():
    def __init__(self):
        pass
    def mutationData(self,alignment1,alignment2,maxHeight=5.0):

        #maxHeight=math.log2(self.alphabet)
        badchar=['-','.','~']
        nrow=len(alignment1)
        ncol=len(alignment1[0])

        native=alignment2[0]

        seq_length=ncol
        lnat=len(native)
        #print(nrow,lnat,ncol)
        #print(native)
        #print(alignment1[0])
        if ncol!=lnat:
            print("Germline seqquence and Alignment length are not same \n")
            seq_length=min([ncol,lnat])
            #return None
        #TODO check if the tow sequences are of equal length or not

        position={}
        germline={}

        for c in range(seq_length):
            gaps=0
            t={}
            total_mut=0
            for r in range(nrow):
                #print(r,c)
                aa=alignment1[r][c]
                #print(r,c,aa)
                #print(aa)
                ga=native[c]

                if aa in badchar:
                    gaps+=1
                    #print (c,gaps)
                    continue
                if aa==ga:
                    continue

                total_mut+=1
                if t.__contains__(aa):
                    f=t[aa]
                    t[aa]=f+1.0

                else:
                    t.__setitem__(aa,1.0)
            #print(c,gaps,nrow,t)
            if gaps==nrow:
                #print ("Eliminated")
                continue
            if t=={}:
                germline[c]={'-':maxHeight}
            else:
                germline[c]={native[c]:total_mut}

            #for key in t:
            #    t[key]=(t[key],

            #t=OrderedDict(sorted(t.items(), key=lambda x:x[1]))

            #print (t)
            position.__setitem__(c,(t,nrow-gaps))
        return position



def findIntersection(d1,d2):
    k1=d1.keys()
    k2=d2.keys()
    return list(set.intersection(set(k1),set(k2)))
def findUnion(d1,d2):
    k1 = d1.keys()
    k2 = d2.keys()
    return list(set.union(set(k1), set(k2)))
def computepValue(r1,r2,alpha=0.05):

    print(r1,r2)
    s=0
    _,p,_,_=stats.chi2_contingency(np.array([r1,r2]))
    if p<=alpha:
        s=1
    return p,s

def combinations(tmp):
    pairs=list(itertools.combinations(tmp,2))
    AApairs=[]
    for pair in pairs:
        p1,n1,r1=pair[0]
        p2,n2,r2=pair[1]
        if p1=={} or p2=={}:
            continue
        #commonAA=findIntersection(p1,p2)
        commonAA=findUnion(p1,p2)
        if commonAA==[]:
            continue
        for AA in commonAA:
            if p1.__contains__(AA):
               pr1=int(p1[AA])
            else:
                pr1=0
            if p2.__contains__(AA):
                pr2 = int(p2[AA])
            else:
                pr2 = 0
            #pr1=int(p2[AA])
            #pr2=int(p2[AA])

            if r1-pr1>0 or r2-pr2>0:
                p_value,s=computepValue([pr1,r1-pr1],[pr2,r2-pr2])
                AApairs.append((AA,(n1,pr1,r1),(n2,pr2,r2),p_value,s))
            else:
                continue
    return AApairs
def getPvalue(alignments,pos):
    tmp=[]
    for all in alignments:
        p,n=all
        d,nr=p[pos] #d={A:10,D:20},nr=100

        tmp.append((d,n,nr))
    return combinations(tmp)

def significanceTest(alignments):
    #Position    Mutant  PIOL DLBCL P-value Significance(0/1)

    #1   D   PIOL-10/35 DLBCL-22/100   0.0001  1
    #1  D   PIOL-10/35  PCNSL-34/35 0.0003  1



    no_alignments=len(alignments)

    p,_=alignments[0]
    positions=p.keys()
    result={}
    for pos in positions:
        pairValue=getPvalue(alignments,pos)
        if pairValue==[]:
            continue
        result.__setitem__(pos,pairValue)
    return result


def read(filename,start_index,length):
    lib=Reader(start_index,length)
    return lib.readCSV(filename)

def writeProportion(alignFiles,germFiles,names,dstfile):

    lib=Stattest()

    raw_proportion=[]
    for f,g,n in zip(alignFiles,germFiles,names):
        mat=read(f,0,-1)
        nat=read(g,0,-1)
        p=lib.mutationData(mat,nat)
        #print(p)
        raw_proportion.append((p,n))
    #for x in raw_proportion:
    #    print(x)
    result=significanceTest(raw_proportion)
    with open(dstfile,'w') as F:
        for pos in result:
            ls=result[pos]
            for l in ls:
                F.write(str(pos)+"\t"+
                        l[0]+"\t"+
                        l[1][0]+'-'+str(l[1][1])+'/'+str(l[1][2])+"\t"+
                        l[2][0]+'-'+str(l[2][1])+'/'+str (l[2][2])+"\t"+
                        str(l[3])+"\t"+
                        str(l[4])+"\n")

            #print (pos,":",result[pos])
    print (result[10])




if __name__ == '__main__':
    #alignFiles=["Dataset/Piol_rev/piol.txt","Dataset/09_03_rev/DLBCL_Seq.txt","Dataset/09_03_rev/PCNS_Seq_V34.txt"]
    #germFiles=["Dataset/Piol_rev/IGVH4_34.txt","Dataset/09_03_rev/DLBCL_Ref.txt","Dataset/09_03_rev/PCNS_Ref_V34.txt"]
    #names=["PIOL","DLBCL","PCNSL"]

    alignFiles=["Dataset/Logo_v34/PIOL_34.txt","Dataset/Logo_v34/PCNSL_34.txt","Dataset/Logo_v34/DLBCL_34.txt","Dataset/Logo_v34/CLL_34.txt"]
    germFiles=["Dataset/Logo_v34/IGVH4_34.txt"]*len(alignFiles)
    names=["PIOL","PCNSL","DLBCL","PIOL"]

    writeProportion(alignFiles,germFiles,names,"Pvalue_PIOL_PCNSL_DLBCL_PIOL_v34.txt")
