"""
Author: Bishnu Sarker, Research Intern at Laboratory of Computational and Quantitative Biology, UPMC, Paris, France.
Email: bishnukuet@gmail.com

"""


import pickle
import os,sys
from sequencelib import sequencelib, computeHeightWithError
from IO import Reader
from drawinglib import Sequencelogo
from argsparse import parse_arguments
from hightmeasure import read_data
from subfamily import subfamily


from random import randint

SEQ_HEIGHT=dict()
parameterSet={"sequence_type":"Protein",
              "stack_width":0.8,
              "stack_height":5.0,
              "alignment1_index":0,
              "alignment2_index":0,
              "sequence_length":-1,
              "show_X_axis":0,
              "show_Y_axis":0,
              "block_size":0,
              "transparency":0.6,
              "show_position":0,
              "tic_size":1.

              }

Usage="""
        python3 biologo.py --logotype [-p|-m|-s|-mv|-mg] --alignments [-f|-a|-t] filename --germline [-f|-a|-t] filename --output [-e|-s|-p] filename
    """
def restoreDataStructure(filename):

    with open(filename,'rb') as F:
        data=pickle._load(F)
    return data


def setConfiguration(configFile="Configuration.config"):
    defaultparam={}
    with open(configFile) as F:
        for line in F:
            line=line.strip()

            if line.isspace():
                #black lines, just consume
                continue
            if line=='':
                continue
            elif line.startswith("#"):
                #comment lines, Just consume
                continue
            else:
                values=line.split("=")
                k=values[0].strip() #parameterSet[values[0].strip()]
                v=parameterSet[k]
                parameterSet[k]=type(v)(values[1].strip())
                #defaultparam.__setitem__(values[0].strip(),values[1].strip())
    return parameterSet

def reader(alignment,file_format,start_index=0,length=-1):
    lib=Reader(start_index,length)
    mat=None
    if file_format=="aln":
        mat=lib.readAln(alignment)
    elif file_format=="fasta":
        mat=lib.readFasta(alignment)
    elif file_format=="csv":
        mat=lib.readCSV(alignment)
    return mat

def sequenceType(mat):


    r=len(mat)
    c=len(mat[0])
    protein=False
    DNA=False

    for j in range(100):


        ri=randint(0,r)
        rj=randint(0,c)
        l=mat[ri][ri]
        if  l  in ['A','a','T','t','C','c','G','g']:
            protein=False
            DNA=True

        elif l in ['-','.']:
            continue

        else:
            protein=True
            DNA=False
            break
    print (protein,DNA)
    if protein==True:
        return "protein"
    else:
        return "DNA"

#{'in': {'names': ['Name1', 'Name2', 'Name3', 'NameN'], 'formats': ['fasta', 'aln'], 'germline': ['filename'], 'alignments': ['filename1', 'filename2', 'filename3', 'filenameN']},
# 'out': {'output': ['filename'], 'format': ['eps']}}
def build_logo(commands):
    lib=sequencelib()
    parameterset=setConfiguration()
    b=commands["set"]["blocksize"]
    parameterset["block_size"]=b
    l = commands["set"]["length"]
    parameterset["sequence_length"] = l

    s=commands['set']["start"]

    parameterset["sequence_type"]=commands["set"]["seqtype"]
    parameterset["alignment1_index"]=s
    parameterset["alignment2_index"]=s


    logo=Sequencelogo(commands["out"]["output"][0],commands["out"]["format"][0], parameters=parameterset)
    logotype=commands["in"]["Type"]

    if logotype=="PlainLogo":
        filename=commands["in"]["alignments"][0]
        aln_format=commands["in"]["formats"][0]

        mat=reader(filename,aln_format,s,l)

        seqlength=len(mat[0])

        if b>=seqlength:
            logo.set_blocksize(seqlength+1)

        position = lib.computeRelativeFrequency(mat)
        #logo.sequenceLogo(position)
        logo.plainLogo(position)

    elif logotype=="MutationLogo":
        filename = commands["in"]["alignments"][0]
        aln_format = commands["in"]["formats"][0]
        germline=commands["in"]["germline"][0]
        germ_format = commands["in"]["formats"][1]

        seq_mat = reader(filename, aln_format,s,l)
        germ_mat=reader(germline, germ_format,s,l)
        alignment, germline = lib.mutationFrequency(seq_mat, germ_mat)

        seqlength = len(seq_mat[0])

        if b >= seqlength:
            logo.set_blocksize(seqlength + 1)

        logo.mutationLogo(alignment, nativeSeq=germline)

    elif logotype=="SubfamilyLogo":
        filename = commands["in"]["alignments"][0]
        aln_format = commands["in"]["formats"][0]
        germline = commands["in"]["germline"][0]
        germ_format = commands["in"]["formats"][0]

        posdic, negdic = subfamily(filename, germline, aln_format, germ_format)
        seqlength = len(posdic)
        if b >= seqlength:
            logo.set_blocksize(seqlength + 1)

        logo.subfamilyLogo(posdic, negdic)

        '''
        fam_mat = reader(filename, aln_format,s,l)
        sub_mat = reader(germline, germ_format,s,l)

        seqlength = len(fam_mat[0])

        if b >= seqlength:
            logo.set_blocksize(seqlength + 1)

        p1 = lib.computeRelativeFrequency(fam_mat)
        p2 = lib.computeRelativeFrequency(sub_mat)
        logo.mutationLogo(p1, nativeSeq=p2)
        '''
    elif logotype=="MultiviewLogo":
        pos_list=[]
        germ_list=[]
        names=commands['in']['names']

        filenames = commands["in"]["alignments"]
        aln_format = commands["in"]["formats"][0]

        germline = commands["in"]["germline"][0]
        germ_format = commands["in"]["formats"][1]

        germ_mat=reader(germline,germ_format,s,l)
        #print (type(germ_mat))
        seqlength = len(germ_mat[0])

        #if b >= seqlength:
        #    logo.set_blocksize(seqlength + 1)
        logo.set_blocksize(seqlength + 1)

        for filename in filenames:
            aln_mat=reader(filename,aln_format,s,l)
            mr,gr=lib.mutationFrequency(aln_mat,germ_mat)
            pos_list.append(mr)
            germ_list.append(gr)
        logo.compareLogo(pos_list,germ_list,names,x=2.0)

    elif logotype=="MutagenesisLogo":
        filename = commands["in"]["alignments"][0]
        #s=Sequencelogo("DCA",parameters=parameterSet)
        #th=math.log2(0.4)


        data,wildtype=read_data(filename,posth=commands["set"]["posth"],negth=commands["set"]["negth"], start_index=s+1, length=l)
        #writefile("data.txt",data)
        logo.DCALogo(data,wildtype)

def main():
    arguments = sys.argv[1:]
    commands = parse_arguments(arguments)
    print(commands)
    print("Drawing Logo....")
    build_logo(commands)


if __name__ == '__main__':

    main()
    #height=computeHeightWithError("Dataset/subfamilylogo/family.txt","Dataset/subfamilylogo/family.txt")
    #print(height)

#--logotype -s --alignments -t Dataset/test2.txt   --germline -t  Dataset/test1.txt --output -p FL --settings -seqtype p -start 1 -length -1 -blocksize 200
