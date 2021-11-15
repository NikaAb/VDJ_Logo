"""
Author: Bishnu Sarker, Research Intern at Laboratory of Computational and Quantitative Biology, UPMC, Paris, France.
Email: bishnukuet@gmail.com

"""


import math
import os
import sys
import pickle
from collections import OrderedDict

import re

AminoAcidName=list()
AminoAcidCode3=list()
AminoAcidCode1=list()








def readAlignment(fileName):
    '''
    This method returns the matrix build out of the alignment provided by the file fileName.
    :param fileName: file containing the multiple sequence alignment of fixed length
    :return: return the sequence in matrix format.
    '''
    mat=list()
    for line in open(fileName):
        #print line
        sequence=line.split()[1]
        mat.append(list(sequence))
    return mat

def readAminoAcids(fileName):
    with open(fileName) as F:
        for line in F:
            name,code3,code1=line.split(',')
            AminoAcidName.append(name.strip())
            AminoAcidCode3.append(code3.strip())
            AminoAcidCode1.append(code1.strip())

def computeHeight(alignments):

    Height=dict()
    alignment=readAlignment(alignments)
    fre_dic=computeFrequency(alignment)

    information=I(fre_dic)
    l=len(information)
    for i in range(l):
        fi=fre_dic[i]
        pi=information[i]
        hi=heightAt(fi,pi)
        Height.__setitem__(i,hi)
    return Height

def heightAt(acid_freq,I_pi):
    heightDic=dict()

    for key in acid_freq.keys():
        h=acid_freq[key]*I_pi
        heightDic.__setitem__(key,h)
    return heightDic

def I(SeqFreq,alphabet=20):

    information=dict()
    l=len(SeqFreq)

    #total=float (l)
    def Iat(pos):
        letters=SeqFreq[pos].keys()
        #print letters
        sum=0
        for ch in letters:

            freq=SeqFreq[pos][ch]
            #print freq
            sum+=(freq*math.log(freq,2))
        return sum
    for i in range(l):
        inf=math.log(alphabet,2)+Iat(i)
        information.__setitem__(i,inf)

    return information
def diff_S_R(S,R):


    diff=S
    S_keys=set(S.keys())
    R_keys=set(R.keys())
    #SUR=S_keys.union(R_keys)

    #SAR=S_keys.intersection(R_keys)

    for key in R_keys:
        if diff.__contains__(key):
            #print (diff[key])
            d=diff[key]-R[key]
            diff[key]=d
        else:
            diff.__setitem__(key,-R[key])
    return diff


def correctionFactor(S,R,alphabet=20):

    totalSeq=float(S+R)
    fs=S/totalSeq
    fr=R/totalSeq
    fsr=(math.log(alphabet,2))/((fs*math.log(fs,2)+(fr*math.log(fr,2))))
    return fsr

def heightWithErrorAt(SR,fSR,I_pi):

    height=dict()
    for key in SR.keys():
        h=SR[key]*fSR*I_pi
        height.__setitem__(key,h)
    return height


def computeHeightWithError(trainSeq, testSeq):


    Height=dict()
    S=computeFrequency(readAlignment(trainSeq))
    R=computeFrequency(readAlignment(testSeq))

    #print type(S), type(R)
    #print S
    #print R
    fSR=correctionFactor(len(S),len(R))
    information=I(S)
    #print information
    l=len(R)
    #print l
    for i in range(l):

        S_R=diff_S_R(S[i],R[i])
        Height.__setitem__(i, heightWithErrorAt(S_R,fSR,information[i]))

    return Height



def maxKey(D):

    max=-100000
    maxkey=None
    for k in D:
        #print k
        if D[k]>max:
            max=D[k]
            maxkey=k
    #print maxkey
    return maxkey


def computeConsensus(filename):
    mat=readAlignment(filename)
    freq=computeFrequency(mat)
    seq=">consensus\n"
    for pos in freq:
        #print pos
        temp=freq[pos]
        #print temp
        #newFreq.__setitem__(pos,maxKey(temp))
        seq+=maxKey(temp)
        #print newFreq

    return seq

def computeFrequencyAt(seqmat,j,ncol):

    pos_set=list()
    temp=dict()
    total=float(ncol)
    i=0
    while i<len(seqmat):
        l=seqmat[i][j]
        if temp.__contains__(l):
            temp[l]+=1/total
        else:
            temp.__setitem__(l,1/total)
        i+=1
    return temp

def computeFrequency(mat):

    freq_dic=dict()

    nrow=len(mat)
    ncol=len(mat[0])

    for pos in range(ncol):
        temp=computeFrequencyAt(mat,pos,ncol)
        freq_dic.__setitem__(pos,temp)
    return freq_dic


def saveDataStructure(dataStructure,filename):
    with open(filename,'wb') as F:
        pickle.dump(dataStructure,F)

def run():
    Seq_SRC="sampledata/alignment.txt"

    #print math.log(10,2)
    readAminoAcids("sampledata/amino_acid_code.txt")
    '''
    mat=readAlignment(Seq_SRC)
    SeqFreq=computeFrequency(mat)
    total_count=len(SeqFreq)
    information=I(SeqFreq,total_count)
    '''
    #computeHeight(Seq_SRC)
    print ( computeConsensus("sampledata/alignment.txt"))
    print (computeConsensus("sampledata/alignmentTest.txt"))
    #saveDataStructure(computeHeightWithError("Data/alignment.txt","Data/alignmentTest.txt"),"positions.pkl")

class sequencelib():
    def __init__(self):

        self.alphabet=20
        #print("sequence type is:", sequence_type)


    def setAlphabet(self,sequence_type="protein"):
        if sequence_type.lower()=="protein":
            self.alphabet=20
        elif sequence_type.lower()=="dna":
            self.alphabet=4




        pass

    def __str__(self):
        pass
    def __repr__(self):
        pass




    def readAlinment(self,filename):

        mat=list()
        i=0

        with open(filename) as F:
            for line in F:
                if line[0]==">":
                    if i>0:
                        mat.append(seq)
                    i=i+1
                    seq=list()
                    continue
                else:
                    seq.extend(line.strip('\n'))
        return mat
    def computeEntropy(self,t):
        if not isinstance(t,dict):
            return t
        sum=math.log2(self.alphabet)
        for ch in t:
            sum=sum+t[ch]*math.log2(t[ch])
        return sum
    def sequence2Dict(self,seq):


        position={}
        for c in range(len(seq)):
            position[c]={seq[c]:1}

        return position
    def computeRelativeFrequency(self,alignMat,maxHeight=5.0):
        badchar=['-','.','~']
        nrow=len(alignMat)
        ncol=len(alignMat[0])
        position=dict()
        for c in range(ncol):
            t={}
            toal_mut=0.0
            for r in range(nrow):
                aa=alignMat[r][c]
                if aa in badchar:
                    continue
                toal_mut+=1
                if t.__contains__(aa):
                    f=t[aa]
                    t[aa]=f+1.0/nrow
                else:
                    t.__setitem__(aa,1.0/nrow)

            if t!={}:
                Pi=self.computeEntropy(t)
                for key in t:
                    t[key]=t[key]*Pi#/nrow*maxHeight
            #Pi=self.computeEntropy(t)

                t=OrderedDict(sorted(t.items(),key=lambda x:x[1]))

                position.__setitem__(c,t)
        return position


    def mutationFrequency(self,alignment1,alignment2):
        #print (alignment1)
        #print(alignment2)
        #printmat(alignment1)
        #printmat(alignment2)
        maxHeight=math.log2(self.alphabet)
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
            total_mut=0.0
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
                germline[c]={native[c]:total_mut/nrow*maxHeight}

            for key in t:
                t[key]=t[key]/nrow*maxHeight

            t=OrderedDict(sorted(t.items(), key=lambda x:x[1]))

            #print (t)
            position.__setitem__(c,t)


        return position,germline

    def frequency2Height(self,F):
        pass




def printmat(mat):
    for row in mat:
        print (row)
if __name__ == '__main__':

   lib=sequencelib()