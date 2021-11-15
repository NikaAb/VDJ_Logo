import re


class Reader():
    def __init__(self,start_index=0,length=-1):
        self.start_index=start_index
        self.end_index=start_index+length
        if self.start_index>self.end_index:
            self.all=True
        else:
            self.all=False

    def readVgene(self,filename):
        mat=list()
        #nativeSeq=None
        #YYC

        with open(filename) as F:
            for line in F:
                if line.isspace(): continue
                print (line)
                t=line.strip("\n ").split()[1][:104]
                #print(t)
                #print(t)
                #if t[2]=="germline":
                #    nativeSeq=t[4]
                #else:
                #TODO check the length and if given length is greater than the sequence length, what to do?
                if self.all:
                    mat.append(list(t)[self.start_index:])
                else:

                    mat.append(list(t)[self.start_index:self.end_index])
        return mat

    def dic2Mat(self,seqs):
        mat=list()
        for seq in seqs:
            line=list(seqs[seq])
            if self.all:
                mat.append(line[self.start_index:])
            else:
                if len(line)>self.end_index:
                    mat.append(line[self.start_index:self.end_index])
                else:
                    mat.append(line[self.start_index:])

        return mat
    def readCSV(self,filename,sep="\t"):
        mat=list()
        #nativeSeq=None

        with open(filename) as F:
            for line in F:
                if line.isspace(): continue
                t=line.strip("\n ").split()
                if self.all:
                    mat.append(list(t[1])[self.start_index:])
                else:
                    if len(t[1])> self.end_index:
                        mat.append(list(t[1])[self.start_index:self.end_index])
                    else:
                        mat.append(list(t[1])[self.start_index:])

        return mat
    def readFasta(self,filename):
        mat=list()
        with open(filename) as F:

            for line in F:
                line=line.strip().upper()
                if line.isspace():continue
                elif line.startswith(">"):
                    continue
                elif line.startswith(";"):
                    continue
                elif line.startswith("*"):
                    continue
                else:
                    if self.all:
                        mat.append(list(line)[self.start_index:])
                    else:
                        if len(line)>self.end_index:
                            mat.append(list(line)[self.start_index:self.end_index])
                        else:
                            mat.append(list(line)[self.start_index:])

        return mat


    def readAln(self,filename):
        '''
        Rule to process a Aln file:
        1. Line starts with CLUSTAL.*
        2. One or more spaces
        3. A block of lines containing sequences
        4.
        :param filename:
        :return:
        '''
        header_line = re.compile(r'(CLUSTAL.*)$')

    # (sequence_id) (Sequence) (Optional sequence number)
        seq_line   = re.compile(r'(\s*\S+\s+)(\S+)\s*(\d*)\s*$')
        match_line = re.compile(r'([\s:\.\*]*)$')

        seqs=dict()

        with open(filename) as F:
            for line in F:
                line=line.strip()
                #print(line)
                if line.isspace():continue

                else:

                    p1=header_line.match(line)
                    p2=seq_line.match(line)
                    p3=match_line.match(line)
                    if not p1 is None:continue
                    elif not p3 is None: continue
                    elif not p2 is None:
                        #print(p3)
                        seq=p2.group()
                        if seq=='':
                            continue

                        seq=seq.split()
                        #print (seq)
                        id=seq[0]

                        if seqs.__contains__(id):
                            seqs[id]+=seq[1]
                        else:
                            seqs.__setitem__(id,seq[1])

        return self.dic2Mat(seqs)


    def readGeneBank(self,filename):
        pass
    def readFa(self,filename):
        return self.readFasta(filename)
    def readBlast(self,filename):
        pass
    def readMsf(self,filename):
        seqs={}
        number_line=re.compile(r'(\d+\s+\d+)$')
        end_header=re.compile(r'(//)(\s*)$')
        seq_line=re.compile(r'\s*(\S+)\s+([\S\s.?]+)$')
        seq_flag=False
        with open(filename) as F:
            for line in F:
                line=line.strip()
                if line.isspace():continue

                if seq_flag is False:
                    if end_header.match(line) is None:
                        seq_flag=False
                        continue
                    else:
                        seq_flag=True
                else:

                    p_num=number_line.match(line)

                    if not p_num is None: continue
                    p_seq=seq_line.match(line)
                    if not p_seq is None:
                        #print(p)
                        seq=p_seq.group(1,2)

                        id=seq[0]
                        if seqs.__contains__(id):
                            seqs[id]+=seq[1].replace(' ','')
                        else:
                            seqs.__setitem__(id,seq[1].replace(' ',''))
                        seq_flag=True

        return self.dic2Mat(seqs)

    def readGcg(self,filename):
        pass
    def readNbrf(self,filename):
        pass
    def readPhy(self,filename):
        pass
    def readPir(self,filename):
        pass
    def readMat(self,filename):
        pass
    def readPhylip(self,filename):
        pass
    def readPdb(self,filename):
        pass

if __name__ == '__main__':
    reader=Reader()
    filename="Dataset/09_03/DLBCL_Seq.txt"
    mat=reader.readVgene(filename)
    print(mat[0])