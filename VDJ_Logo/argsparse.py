


def get_input_format(case):
    format="fasta"
    if case=="-f":
        format="fasta"
    elif case=="-a":
        format="aln"
    elif case=="-t":
        format="csv"
    return format

def get_output_format(case):
    format="svg"
    if case=="-e":
        format="eps"
    elif case=="-p":
        format="pdf"
    elif case=="-s":
        format="svg"
    return format
def get_sequence_type(case):
    seqtype="protein"
    case=case.upper()
    if case=="P" or case=="PROTEIN":
        seqtype="protein"
    elif case=="D" or case=="DNA":
        seqtype="dna"
    elif case=="A" or case=="AUTO":
        seqtype="auto"
    return seqtype
def collect_settings(kk,arguments):
    kk=kk+1
    # parameters: -seqtype, -logorange -blocksize

    pk=get_option(arguments,"-seqtype",kk)
    seqtype=get_sequence_type(arguments[pk+1])

    pk=get_option(arguments,"-start",kk)

    start=int(arguments[pk+1])
    if start>0:
        start=start-1
    pk = get_option(arguments, "-length", kk)
    length = int(arguments[pk + 1])

    pk = get_option(arguments, "-blocksize", kk)
    block=int(arguments[pk+1])
    blocksize=block
    if length>0:
        if block>length:
            blocksize=length
        else:
            blocksize=block
    pk=get_option(arguments, "-posth", kk)
    print(kk, pk)
    if pk==kk or pk==len(arguments):
        posth=-6
    else:
        posth=float(arguments[pk+1])

    pk=get_option(arguments,"-negth", kk)
    if pk == kk or pk==len(arguments):
        negth=-0.6
    else:
        negth=float(arguments[pk+1])

    return {"seqtype":seqtype,"length":length,"blocksize":blocksize,"start":start,"negth":negth,"posth":posth}

def collect_output(kk,arguments):
    outputs=[]
    formats=[]

    if arguments[kk + 1] in ['-e', '-p', '-s']:
        formats.append(get_output_format(arguments[kk + 1]))
        kk += 1
    outputs.append(arguments[kk + 1])
    return {"output":outputs,"format":formats}
def collect_inputs(k, arguments):
    alignments=[]
    germline=[]
    names=[]
    formats=[]
    type=""
    posth = -0.5
    negth = -6
    if arguments[k+1].upper()=="-plainlogo".upper() or arguments[k+1].upper()=="-p".upper() :
        kk=get_option(arguments,"--alignments")
        if arguments[kk + 1] in ['-f', '-a', '-t']:
            formats.append(get_input_format(arguments[kk+1]))
            kk += 1
        alignments.append(arguments[kk+1])
        type="PlainLogo"
    elif arguments[k + 1].upper() == "-mutation".upper() or arguments[k+1].upper()=="-m".upper():
        kk = get_option(arguments, "--alignments")

        if arguments[kk + 1] in ['-f', '-a', '-t']:
            formats.append(get_input_format(arguments[kk+1]))
            kk += 1
        alignments.append(arguments[kk+1])

        kk = get_option(arguments, "--germline")
        if arguments[kk + 1] in ['-f', '-a', '-t']:
            formats.append(get_input_format(arguments[kk+1]))
            kk += 1
        germline.append(arguments[kk + 1])
        type = "MutationLogo"
    elif arguments[k + 1].upper() == "-subfamily".upper() or arguments[k+1].upper()=="-s".upper():
        kk = get_option(arguments, "--alignments")
        if arguments[kk + 1] in ['-f', '-a', '-t']:
            formats.append(get_input_format(arguments[kk+1]))
            kk += 1
        alignments.append(arguments[kk + 1])

        kk = get_option(arguments, "--germline")
        if arguments[kk + 1] in ['-f', '-a', '-t']:
            kk += 1
        germline.append(arguments[kk + 1])
        type = "SubfamilyLogo"
    elif arguments[k + 1].upper() == "-multiview".upper() or arguments[k+1].upper()=="-mv".upper():
        kk = get_option(arguments, "--alignments")

        while True:

            if arguments[kk+1].startswith("--"):
                break
            elif arguments[kk+1] in ['-f','-a','-t']:


                formats.append(get_input_format(arguments[kk+1]))
                kk += 1
                continue
            else:
                kk = kk + 1
                alignments.append(arguments[kk])

        kk = get_option(arguments, "--germline")
        if arguments[kk + 1] in ['-f', '-a', '-t']:
            formats.append(get_input_format(arguments[kk+1]))
            kk += 1
        germline.append(arguments[kk + 1])

        kk=get_option(arguments,"--labels")
        while True:
            if arguments[kk+1].startswith("--"):
                break
            else:
                kk+=1
                names.append(arguments[kk])
        type = "MultiviewLogo"
    elif arguments[k + 1].upper() == "-mutagenesis".upper() or arguments[k+1].upper()=="-mg".upper():
        kk = get_option(arguments, "--input")
        alignments.append(arguments[kk+1])
        #kk==get_option(arguments,"-pth")
        #posth=float(arguments[kk+1])

        #kk == get_option(arguments, "-nth")
        #negth = float(arguments[kk + 1])



        type = "MutagenesisLogo"
    D={"alignments":alignments,"germline":germline,"names":names,"formats":formats,"Type":type,"posth":posth,"negth":negth}
    return D
def get_option(arguments,option,pos=0):

    l = len(arguments)
    i = pos
    while i<l:
        if arguments[i].upper() == option.upper():
            break
        else:
            i+=1
            continue
    return i
def parse_arguments(arguments):

    #print(arguments)
    k=get_option(arguments,'--logotype')

    D_in=collect_inputs(k, arguments)
    kk = get_option(arguments, "--output")
    D_out=collect_output(kk,arguments)
    kk=get_option(arguments,"--settings")
    D_set=collect_settings(kk,arguments)

    return {"in":D_in,"out":D_out,"set":D_set}




if __name__ == '__main__':

    Test1="biologo --logotype -m --alignments -t Dataset/Piol_rev/piol.txt  --germline -t Dataset/Piol_rev/IGVH4_34.txt  --output -p mutation"
    Test2="biologo --logotype -p --alignments -t Dataset/Piol_rev/piol.txt --output -p plain --settings -seqtype 0 -start 2 -length -1 -blocksize 50"
    Test3="biologo --logotype -s --alignments -f Dataset/Piol_rev/piol.txt  --germline -a Dataset/Piol_rev/IGVH4_34.txt  --output -p mutation"
    Test4="biologo --logotype -mv --alignments -t Dataset/Piol_rev/piol.txt Dataset/09_03_rev/DLBCL_Seq.txt Dataset/09_03_rev/PCNS_Seq_V34.txt --labels PIOL DLBCL PCNSL --germline -t Dataset/Piol_rev/IGVH4_34.txt  --output -p MV"
    Test5="biologo --logotype -mg --input filename1 --output -p filename  --settings -posth N -negth N"
    print(parse_arguments(Test2.split()[1:]))