
"""
def executeCommand(command):


    if not command["configuration"]["filename"] is None:
        setConfiguration(configFile=command["configuration"]["filename"])
    #print(parameterSet)
    lib=sequencelib()
    logo=Sequencelogo(command["output"]["filename"],command["output"]["format"],parameters=parameterSet)

    if command["output"]["logotype"]=="-normal":
        mat=reader(command,"alignment1",parameterSet["alignment1_index"],parameterSet["sequence_length"])
        if parameterSet["sequence_type"].lower()=="auto":
            lib.setAlphabet(sequence_type=sequenceType(mat))
        else:
            lib.setAlphabet(sequence_type=parameterSet["sequence_type"])
        position=lib.computeRelativeFrequency(mat)
        logo.sequenceLogo(position)
    elif command["output"]["logotype"]=="-mutation":

        mat=reader(command,"alignment1",parameterSet["alignment1_index"],parameterSet["sequence_length"])
        native=reader(command,"alignment2",parameterSet["alignment2_index"],parameterSet["sequence_length"])

        alignment,germline=lib.mutationFrequency(mat,native)
        #alignment=lib.mutationData(mat,native)

        logo.mutationLogo(alignment,nativeSeq=germline)
    elif command["output"]["logotype"]=="-subfamily":
        mat1=reader(command,"alignment1",parameterSet["alignment1_index"],parameterSet["sequence_length"])
        mat2=reader(command,"alignment2",parameterSet["alignment2_index"],parameterSet["sequence_length"])
        p1=lib.computeRelativeFrequency(mat1)
        p2=lib.computeRelativeFrequency(mat2)
        logo.mutationLogo(p1,nativeSeq=p2)

def readCommandline():

    #'''
    #This function will parse commandline parameters and check the consistancy with predefined format
    #:return: a list of parameters
    #'''
    #lib=sequencelib()
    command=\
        {
        "alignment1":
            {
                "filename":None,
                "format":None,
                "reader":None

            },
        "alignment2":
            {
                "filename":None,
                "format":None,
                "reader":None

            },
        "output":
            {
                "filename":None,
                "format":None,
                "logotype":"normal"
            },
        "configuration":
            {
                "filename":None,
                "format":None
            },

        }
    arguments=sys.argv[1:]
    #print (sys.argv)
    #print (arguments)

    i=0
    l=len(arguments)
    while i<l :
        if arguments[i]=="--input":
            if arguments[i+1]=="-f":

                command["alignment1"]["filename"]=arguments[i+2]
                command["alignment1"]["format"]="fasta"
                if arguments[i+3]=="-f":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="fasta"

                elif arguments[i+3]=="-a":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="aln"
                elif arguments[i+3]=="-t":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="csv"
                else:
                    print(Usage)

            elif arguments[i+1]=="-a":
                command["alignment1"]["filename"]=arguments[i+2]
                command["alignment1"]["format"]="aln"
                if arguments[i+3]=="-f":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="fasta"
                elif arguments[i+3]=="-a":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="aln"
                elif arguments[i+3]=="-t":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="csv"
                else:
                    print(Usage)
            elif arguments[i+1]=="-t":
                command["alignment1"]["filename"]=arguments[i+2]
                command["alignment1"]["format"]="csv"
                if arguments[i+3]=="-f":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="fasta"
                elif arguments[i+3]=="-a":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="aln"
                elif arguments[i+3]=="-t":
                    command["alignment2"]["filename"]=arguments[i+4]
                    command["alignment2"]["format"]="csv"
                else:
                    print(Usage)
            else:
                print (Usage)
            i=i+5
        elif arguments[i]=="--output":
            if arguments[i+1]=="-e":
                command["output"]["filename"]=arguments[i+2]
                command["output"]["format"]="EPS"
            elif arguments[i+1]=="-p":
                command["output"]["filename"]=arguments[i+2]
                command["output"]["format"]="PDF"
            elif arguments[i+1]=="-s":
                command["output"]["filename"]=arguments[i+2]
                command["output"]["format"]="SVG"
            else:
                print (Usage)
            i=i+3
        elif arguments[i]=="--configuration":
            command["configuration"]["filename"]=arguments[i+1]
            i=i+2
        elif arguments[i]=="--logotype":
            command["output"]["logotype"]=arguments[i+1]
            i=i+2
        else:
            print (Usage)
    return command

def run_comparison():
    lib=sequencelib()
    R=Reader(0,-1)
    #logo=Sequencelogo(command["output"]["filename"],command["output"]["format"],parameters=parameterSet)
    '''
    alignFiles=["Dataset/Seqs_am_rev/CLONE_1_am.txt",
                "Dataset/Seqs_am_rev/CLONE_2_am.txt",
                "Dataset/Seqs_am_rev/CLONE_3_am.txt",
                "Dataset/Seqs_am_rev/CLONE_5_am.txt",
                "Dataset/Seqs_am_rev/CLONE_8_am.txt",

                "Dataset/Seqs_am_rev/CLONE_22_am.txt",
                "Dataset/Seqs_am_rev/CLONE_31_am.txt",
                "Dataset/Seqs_am_rev/CLONE_55_am.txt",
                "Dataset/Seqs_am_rev/CLONE_56_am.txt"]
    germFiles=["Dataset/Refs_rev/CLONE_1_ref_am.txt",
               "Dataset/Refs_rev/CLONE_2_ref_am.txt",
               "Dataset/Refs_rev/CLONE_3_ref_am.txt",
               "Dataset/Refs_rev/CLONE_5_ref_am.txt",
               "Dataset/Refs_rev/CLONE_8_ref_am.txt",

               "Dataset/Refs_rev/CLONE_22_ref_am.txt",
               "Dataset/Refs_rev/CLONE_31_ref_am.txt",
               "Dataset/Refs_rev/CLONE_55_ref_am.txt",
               "Dataset/Refs_rev/CLONE_56_ref_am.txt"]
    names=["CLONE1",
           "CLONE2",
           "CLONE3",
           "CLONE5",
           "CLONE8",

           "CLONE22",
           "CLONE31",
           "CLONE55",
           "CLONE56"]
    '''
    alignFiles=["Dataset/Piol_rev/piol.txt","Dataset/09_03_rev/DLBCL_Seq.txt","Dataset/09_03_rev/PCNS_Seq_V34.txt"]
    germFiles=["Dataset/Piol_rev/IGVH4_34.txt","Dataset/09_03_rev/DLBCL_Ref.txt","Dataset/09_03_rev/PCNS_Ref_V34.txt"]
    names=["PIOL","DLBCL","PCNSL"]
    setConfiguration()
    logo=Sequencelogo("CLL","all",parameters=parameterSet)
    pos_list=[]
    germ_list=[]
    for a,g,n in zip(alignFiles,germFiles,names):
        mat=R.readCSV(a)
        nat=R.readCSV(g)
        mr,gr=lib.mutationFrequency(mat,nat)
        #print(mat)
        #print(nat)
        pos_list.append(mr)
        germ_list.append(gr)
    #print(len(pos_list))
    #print (len(germ_list))
    logo.compareLogo(pos_list,germ_list,names,x=2.0)

"""