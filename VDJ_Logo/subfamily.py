from CM import online_CM,online_CM_fasta,AA2NUM,NUM2AA
import numpy as np
import scipy as sp
import math
from collections import OrderedDict





def I(mat):

    row, col=mat.shape
    const=math.log2(20)
    P = [0 for i in range(row)]
    for i in range(row):
        total=float(sum(mat[i][:20]))
        #print(mat[i][:20])
        P[i]=const+sum([j/total*math.log2(j/total) for j in mat[i][:20] if j!=0 ])
    return  np.array(P)


def normalize_frequency(mat):

    row,col=mat.shape

    for i in range(row):
        total=float(sum(mat[i][:20]))
        if total<=0:
            total=1

        for j in range(col):
            mat[i][j]/=total
    return mat


def correction_factor(srow, rrow):
    total=srow+rrow
    log20=math.log2(20)
    return log20/(log20+(srow/total)*math.log2(srow/total)+(rrow/total)*math.log2(rrow/total))


def compute_height(S_R,fsr,Pi):
    row, col=S_R.shape
    T=np.transpose(S_R)
    T=np.multiply(T,fsr)
    T=np.multiply(T,Pi)
    #print(T.shape)
    return np.transpose(T)




def collect_residue(HeightMat):
    posdic={}
    negdic={}
    row, col=HeightMat.shape

    for i in range(row):
        pTD = {}
        nTD = {}
        for j in range(col):
            if HeightMat[i][j]==0.0:
                continue
            elif HeightMat[i][j]<0.0:
                nTD.update({NUM2AA[j]:abs(HeightMat[i][j])})

            elif HeightMat[i][j]>0.0:
                pTD.update({NUM2AA[j]: HeightMat[i][j]})


        if pTD=={}:
            pTD.__setitem__("-",5.0)
        if nTD=={}:
            nTD.__setitem__("-",5.0)
        pTD = OrderedDict(sorted(pTD.items(), key=lambda x: x[1]))
        nTD = OrderedDict(sorted(nTD.items(), key=lambda x: x[1]))
        posdic.__setitem__(i,pTD)
        negdic.__setitem__(i,nTD)
    return (posdic,negdic)





def subfamily(familyData, subfamilyData, aln_format, germ_format):
    print ("DEGUG::subfamily");
    if aln_format=="fasta" or aln_format=="aln" :
        familyMat = online_CM_fasta(familyData)
    else:
        familyMat=online_CM(familyData)
    if germ_format=="fasta" or germ_format=="aln":
        subfamilyMat = online_CM_fasta(subfamilyData)
    else:
        subfamilyMat=online_CM(subfamilyData)

    fullMat=np.add(familyMat,subfamilyMat)

    Pi=I(fullMat)
    #print(Pi, Pi.shape)
    normFamily=normalize_frequency(familyMat)
    normSubFamily=normalize_frequency(subfamilyMat)

    S_R=np.subtract(normSubFamily,normFamily)

    rrow,_=familyMat.shape
    srow,_=subfamilyMat.shape

    fsr=correction_factor(srow,rrow)
    H=compute_height(S_R,fsr,Pi)
    posdic,negdic=collect_residue(H)
    #print ("DEGUG:: pos neg",  posdic, negdic)
    return (posdic,negdic)














if __name__ == '__main__':
    subfamily("Dataset/test1.txt","Dataset/test2.txt")
    #print (AA2NUM)
    #print(NUM2AA)