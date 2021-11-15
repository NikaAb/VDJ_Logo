import numpy as np
from itertools import combinations


AA2NUM={"A":0,  "D":2,   "C":1  , "E":3 ,  "F":4,   "G":5,   "H":6,   "I":7,   "K":8,"L":9,
        "M":10,   "N":11,   "P":12 ,  "Q":13,   "R":14,   "S":15,   "T":16,   "V":17,   "W":18,   "Y":19,".":20,"-":21}
NUM2AA={0:"A",1:"C",2:"D",3:"E",4:"F",5:"G",6:"H",7:"I",8:"K",9:"L",
        10:"M",11:"N",12:"P",13:"Q",14:"R",15:"S",16:"T",17:"V",18:"W",19:"Y",20:".",21:"-"}



def seq2vec(sequence):
    seq_vec = np.zeros((len(sequence), len(AA2NUM)))
    i=0
    for ch in sequence:

        seq_vec[i,AA2NUM[ch]]=1
        i+=1
    return seq_vec
def update_concensus(current, seed):
    #print current.shape, seed.shape
    return np.add(current,seed)
def vec2seq(vec):
    maxs = np.argmax(vec, axis=1)
    con_seq = ""
    for x in maxs:
        con_seq = con_seq + NUM2AA[x]
    return con_seq

def concensus_machine_CM(sequences):
    concensus=np.zeros((len(sequences[0]),len(AA2NUM)))
    for sequence in sequences:
        concensus=update_concensus(concensus,seq2vec(sequence))

    return vec2seq(concensus)

def online_CM_fasta(filename):
    with open(filename) as F:
        seq_length = 0
        for line in F:
            # print(line)
            # print (line)
            if line == "":
                continue
            if line.startswith(">"):
                continue
            lin = line.strip()
            # print (lin)
            seq_length = len(lin)
            break

        concensus = np.zeros((seq_length, len(AA2NUM)))
    with open(filename) as F:
        for line in F:
            # print(line)
            if line == "":
                continue
            if line.startswith(">"):
                continue
            # print(line)
            line = line.strip()
            # print(line)
            concensus = update_concensus(concensus, seq2vec(line))
    return concensus
    pass
def online_CM(filename):

    with open(filename) as F:
        seq_length=0
        for line in F:
            #print(line)
            #print (line)
            if line=="":
                continue
            lin = line.strip().split()[1]
            #print (lin)
            seq_length=len(lin)
            break


        concensus = np.zeros((seq_length, len(AA2NUM)))
    with open(filename) as F:
        for line in F:
            #print(line)
            if line=="":
                continue
            #print(line)
            line=line.strip().split()[1]
            #print(line)
            concensus=update_concensus(concensus,seq2vec(line))
    return concensus
def filter_germline(consensus,germline):
    CM_germ=seq2vec(germline)

    x,y=consensus.shape
    total_seq=sum(consensus[0][:])
    #print (total_seq)

    #print ((x,y),CM_germ.shape)
    for i in range(x):
        for j in range(y):
            if consensus[i][j]==0:
                continue
            else:
                if CM_germ[i][j]==1:
                    consensus[i][j]=0
                else:
                    continue
    return (consensus,total_seq)



def get_con_seq(sequences):

    if len(sequences)==1:
        return sequences[0]['CDR3']
    concensus = np.zeros((len(sequences[0]['CDR3']), len(AA2NUM)))

    for sequence in sequences:
        #print len(sequence['CDR3']), sequence['CDR3']
        concensus = update_concensus(concensus, seq2vec(sequence["CDR3"]))
    return vec2seq(concensus)

def get_distinct_mutant(data_vector):
    y=len(data_vector)
    count=0
    for j in range(y):
        if j==20 or j==21:
            continue
        else:
            if data_vector[j]!=0:
                count+=1
    return count

def reporting_stat(alignments, germline, names=["set1","Set2"],output="reporting_stat.csv"):
    W=open(output,"w")

    AA=[NUM2AA[i] for i in range(22)]
    tite=","+names[0]+",,,,,,,,,,,,,,,,,,,,,,,,,,,"+names[1]+",,,,,,,,,,,,,,,,,,,,,,,,,,\n"
    header="Position,"+ \
           names[2] + ","+str(AA).strip("[]")+","+"No Mutants"+","+"Total Mutations"+","+"Total Sequences"+","+"Percent Mutation"+","+\
           names[2]+","+str(AA).strip("[]")+","+"No Mutants"+","+"Total Mutations"+","+"Total Sequences"+","+"Percent Mutation"+"\n"
    W.write(tite)
    W.write(header)
    CM_Data1,total_seq1=filter_germline(online_CM(alignments[0]),germline)
    CM_Data2,total_seq2 = filter_germline(online_CM(alignments[1]),germline)
    print (names,total_seq1,total_seq2)
    #print(CM_Data1)
    #CM_germline=seq2vec("QVQLQQWGA-GLLKPSETLSLTCAVYGGSF----SGYYWSWIRQPPGKGLEWIGEINHS---GSTNYNPSLK-SRVTISVDTSKNQFSLKLSSVTAADTAVYYC")
    #print (CM_germline)
    x,y=CM_Data2.shape
    #print(x,y)
    for i in range(x):

        total_mutation1=sum(CM_Data1[i][:])-CM_Data1[i][20]-CM_Data1[i][21]
        total_mutation2 = sum(CM_Data2[i][:]) - CM_Data2[i][20] - CM_Data2[i][21]

        no_mutant1=get_distinct_mutant(CM_Data1[i][:])
        no_mutant2 = get_distinct_mutant(CM_Data2[i][:])

        active_seq_count1=(total_seq1-CM_Data1[i][20]-CM_Data1[i][21])
        active_seq_count2 = (total_seq2 - CM_Data2[i][20] - CM_Data2[i][21])

        if active_seq_count1==0:
            active_seq_count1=1.0
        if active_seq_count2==0:
            active_seq_count2=1.0

        mut_percet1=round(100.0*total_mutation1/active_seq_count1,2)
        mut_percet2=round(100.0*total_mutation2/active_seq_count2,2)



        #total_seq1=sum(CM_Data1[i][:])
        #total_seq2=sum(CM_Data2[i][:])
        #print (total_mutation1,total_mutation2)
        #print (total_seq1,total_seq2)

        line = str(i) + ','+germline[i]+','
        for j in range(y):
            #print(NUM2AA[j], CM_germline[i][j])
            if CM_Data1[i][j]==0 :
                line += "\t" + ","
            else:
                line+=str(int(CM_Data1[i][j]))+","
        line += str(no_mutant1)+','+str(int(total_mutation1))+","+str(int(total_seq1))+','+str(mut_percet1)+','+germline[i]+','
        for j in range(y):
            if CM_Data2[i][j] == 0:
                line += "\t" + ","
            else:
                line += str(int(CM_Data2[i][j])) + ","
        line +=str(no_mutant2)+','+ str(int(total_mutation2))+","+str(int(total_seq2))+','+str(mut_percet2)

        W.write(line.strip(',')+"\n")
    W.close()


def find_combinations(files):
    pairs=combinations(files,2)
    return list(pairs)

def massive_reports(alignments, germline,names, directory="reports",suffix="34"):
    pairs=find_combinations(range(len(alignments)))
    for pair in pairs:

        alns=[alignments[pair[0]],alignments[pair[1]]]
        name=[names[pair[0]],names[pair[1]],"IGHV"+"_"+suffix]
        report_name=directory+"/"+name[0]+"_"+name[1]+"_v"+suffix+".csv"

        print (alns,name,report_name)
        reporting_stat(alns,germline,names=name,output=report_name )





if __name__ == '__main__':



    aln_v34=["Dataset/Logo_v34/DLBCL_34.txt", "Dataset/Logo_v34/PIOL_34.txt",
             "Dataset/Logo_v34/PCNSL_34.txt", "Dataset/Logo_v34/CLL_34.txt"]
    germline34 = "QVQLQQWGA-GLLKPSETLSLTCAVYGGSF----SGYYWSWIRQPPGKGLEWIGEINHS---GSTNYNPSLK-SRVTISVDTSKNQFSLKLSSVTAADTAVYYC"
    names34 = ["DLBCL", "PIOL", "PCNSL", "CLL"]

    aln_v7 = ["Dataset/Logo_v7/DLBCL_7.txt", "Dataset/Logo_v7/PIOL_7.txt",
               "Dataset/Logo_v7/PCNSL_7.txt", "Dataset/Logo_v7/CLL_7.txt"]
    names7 = ["DLBCL", "PIOL", "PCNSL", "CLL"]
    germline7="EVQLVESGG-GLVQPGGSLRLSCAASGFTF----SSYWMSWVRQAPGKGLEWVANIKQD--GSEKYYVDSVK-GRFTISRDNAKNSLYLQMNSLRAEDTAVYYC"
    #print (concensus_machine_CM(["AAGDDFWSGYSV",  "ARGYDFWSGYCY",  "ARGYDFWSGYSY","ARGYDFWSGYQN", "AAGYDFWSGYYF", "AGDDFWSGYFGF"]))
    #find_combinations(index_34)
    massive_reports(aln_v34,germline34,names=names34,directory="reports/V4_34_update",suffix="34")
    #massive_reports(aln_v7,germline7,names=names7,directory="reports/V3_7_update",suffix="7",)
    #reporting_stat(["Dataset/Logo_v7/DLBCL_7.txt", "Dataset/Logo_v7/PIOL_7.txt"], germline7,
     #              names=["DLBCL", "PIOL", "IGHV3_7"],
     #              output="reports/DLBCL_PIOL_v8.csv")

    #Cm1=online_CM("Dataset/Logo_v34/CLL_34.txt")
    #print(sum(Cm1[0][:]))
    '''
    reporting_stat(["Dataset/Logo_v7/DLBCL_7.txt","Dataset/Logo_v7/PIOL_7.txt"],germline7,
                   names=["DLBCL","PIOL","IGHV3_7"],
                   output="reports/DLBCL_PIOL_v7.csv")
    reporting_stat(["Dataset/Logo_v34/DLBCL_34.txt", "Dataset/Logo_v34/PIOL_34.txt"], germline7,
                   names=["DLBCL", "PIOL", "IGHV3_7"],
                   output="reports/DLBCL_PIOL_v7.csv")



    [[ 4.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
   0.  0.]
 [ 1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  3.  0.  0.
   0.  0.]
 [ 0.  1.  2.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
   0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  4.  0.  0.  0.  0.  0.  0.  0.
   0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  4.  0.  0.  0.  0.  0.  0.
   0.  0.]
 [ 0.  0.  0.  0.  0.  4.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
   0.  0.]]
    '''

