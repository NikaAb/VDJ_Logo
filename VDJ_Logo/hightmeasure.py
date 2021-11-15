import math
#from sequencelib import sequencelib


from drawinglib import Sequencelogo
def exp_diff(fitA,th):
    return fitA-th


def read_data(filename,posth=-1,negth=-4, start_index=0, length=-1):
    pos_data={}
    neg_data={}
    wildtype={}
    data={}

    end_index=-1
    if length!=-1:
        end_index=start_index+length
    print(start_index,end_index)

    with open(filename) as F:
        for line in F:
            line=line.strip().split()[0:4]
            pos,wt,mt,fit=int(line[0]),line[1],line[2],float(line[3])

            if data.__contains__(pos):
                t=data[pos]
                t.__setitem__(mt,math.log2(fit))
                data.__setitem__(pos,t)
            else:
                data.__setitem__(pos,{mt:math.log2(fit)})
                wildtype.__setitem__(pos,wt)
    heights={}
    for pos in data:
        if end_index==-1:
            if pos>=start_index:
                heights.__setitem__(pos, compute_height(data[pos], posth, negth))
        else:
            if pos>=start_index and pos<end_index:
                heights.__setitem__(pos,compute_height(data[pos],posth,negth))
    #print(heights)
    #print(heights.keys())
    return heights,wildtype
def computeEntropy(self,t):
    if not isinstance(t,dict):
        return t
    sum=math.log2(self.alphabet)
    for ch in t:
        sum=sum+t[ch]*math.log2(t[ch])
    return sum
def compute_height(t,posth,negth):
    pos={}
    neg={}
    sumpos=0.0
    sumneg=0.0
    entropypos=10.0
    entropyneg=10.0#math.log2(20)
    for k in t:
        if t[k]<=negth:

            step=math.fabs(negth-t[k])
            sumneg=sumneg+step
        elif t[k]>=posth:
            step=math.fabs(posth-t[k])
            sumpos=sumpos+step

    for ki in t:
        #if t[k]<negth and t[k]>posth:

        if t[ki]>=posth:
            step=math.fabs((t[ki]-posth)/sumpos)

            entropypos-=step*math.log2(step)
            pos.__setitem__(ki,step)
        if t[ki]<=negth:
            step=math.fabs((t[ki]-negth)/sumneg)

            entropyneg-=step*math.log2(step)
            neg.__setitem__(ki,step )
    for l in pos:
        pos[l]=(pos[l]*entropypos)*0.15
    for l in neg:
        neg[l]=(neg[l]*entropyneg)*0.15
    return {"pos":pos, "neg":neg}


    #print (data)
    #print(wildtype)
            #print(pos,wt,mt,fit)




def writefile(filename,data):
    with open(filename,"w") as F:
        for pos in data:
            posv=data[pos]["pos"]
            negv=data[pos ]["neg"]
            F.write(str(posv)+"\t"+str(negv)+"\n")



if __name__ == '__main__':
    parameterSet={"sequence_type":"Protein",
              "stack_width":0.8,
              "stack_height":5.0,
              "alignment1_index":0,
              "alignment2_index":0,
              "sequence_length":-1,
              "show_X_axis":1,
              "show_Y_axis":0,
              "block_size":60,
              "transparency":0.6,
              "show_position":0,
              "tic_size":1.

              }
    s=Sequencelogo("DCA",parameters=parameterSet)
    #th=math.log2(0.4)
    data,wildtype=read_data("mutagenesis.txt",posth=-0.5,negth=-6)
    #writefile("data.txt",data)
    s.DCALogo(data,wildtype)

