import os

def readVgene(src_file, dst_file):
        mat=list()
        #nativeSeq=None
        #YYC
        F2=open(dst_file,'w')
        with open(src_file) as F:
            for line in F:
                if line.isspace(): continue
                #print (line)
                t=line.strip("\n ").split()
                if t[1][103]=='C':
                    F2.write(t[0]+'\t'+t[1][:104]+'\n')
                elif t[1][104]=="C":
                    F2.write(t[0]+'\t'+t[1][:104]+'\n')



                #print(t)
                #print(t)
                #if t[2]=="germline":
                #    nativeSeq=t[4]
                #else:
        F2.close()
                #TODO check the length and if given
def convertFiles(src_dir,dst_dir):
    filenames=getFiles(src_dir)
    for filename in filenames:
        readVgene(os.path.join(src_dir,filename),os.path.join(dst_dir,filename))

def getFiles(directory):
    return os.listdir(directory)
if __name__ == '__main__':
    pass
    #convertFiles("Dataset/10_03","Dataset/10_03_rev")
    #readVgene("Dataset/09_03_rev/DLBCL_7.txt","Dataset/09_03_rev/DLBCL_Seq_V7_rev.txt")