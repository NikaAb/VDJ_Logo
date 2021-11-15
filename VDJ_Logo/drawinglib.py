"""
Author: Bishnu Sarker, Research Intern at Laboratory of Computational and Quantitative Biology, UPMC, Paris, France.
Email: bishnukuet@gmail.com

"""



from pyx import canvas,unit,path,color,text,trafo,style,graph
from ColorScheme import ColorScheme
from pyx.graph.axis import linear,texter
from subfamily import subfamily
import sys 

import math
from collections import OrderedDict


class Sequencelogo():
    def __init__(self,outfile="Biologo",ext="pdf",xs=3.0,parameters=None):


        unit.set(xscale=xs)
        #text.LatexRunner()
        #text.set(mode="latex")
        text.set(text.LatexRunner)
        #text.preamble(r"\usepackage{courier}")

        self.drawingCanvas=canvas.canvas()

        #setting the parameters
        self.outfile=outfile
        self.format=ext

        if not parameters is None:

            self.block_size=parameters["block_size"]
            self.stack_height=parameters["stack_height"]
            self.stack_width=parameters["stack_width"]
            self.alignment1_index=parameters["alignment1_index"]
            self.alignment2_index=parameters["alignment2_index"]
            self.sequence_type=parameters["sequence_type"]
            self.transparency=parameters["transparency"]
            self.show_X_axis=parameters["show_X_axis"]
            self.show_Y_axis=parameters["show_Y_axis"]
            self.show_position=parameters["show_position"]
            self.tic_size=parameters["tic_size"]
            self.sequence_length=parameters["sequence_length"]


    def set_blocksize(self,blocksize):
        self.block_size = blocksize
    def get_blocksize(self):
        return self.block_size
    def set_seqlength(self,length):
        self.sequence_length=length
    def get_seqlength(self):
        return self.sequence_length
    def set_configuration(self):
        pass

    def read_configuration(self):
        pass

    def show(self,sequenceHeights,x_run=0,y_run=0,maxHeight=5.0):

        x=x_run
        y=y_run
        #for position in sequenceHeights:


        for position in sequenceHeights:
            positiondict=sequenceHeights[position]
            #print (positiondict)
            i,j=self.draw_stack(positiondict,x_run=x,y_run=y)
            x+=i

        x_axis=path.line(-1.0,0.0,x,0.0)
        y_axis=path.line(-1.0,0.0,-1.0,maxHeight)
        self.drawingCanvas.stroke(x_axis,[style.linewidth.thin,color.cmyk.Red])
        self.drawingCanvas.stroke(y_axis,[style.linewidth.thin,color.cmyk.Red,trafo.scale(sx=1,sy=1)])

        self.write(fileFormat="pdf")

    def draw_stack(self,position,x_run=0.,y_run=0.):

        y_neg=y_run
        y_pos=y_run
        c_pos=canvas.canvas()
        c_neg=canvas.canvas()
        cs=ColorScheme()
        for letter in position:
            if position[letter]>=0:
                #letter=r"\texttt{H}"
                yscale=math.fabs(position[letter])
                #print("DEBUG:: draw_stack:: yscale ", yscale)
                dimensions=c_pos.text(x_run,y_pos,r"\textsf{"+letter+"}",
                                      [ cs.getColor(letter),
                                        text.halign.center,
                                        text.size.Huge,
                                        trafo.scale(sx=0.8,sy=yscale)]
                                      )
                y_pos+=dimensions.height*yscale
            else:
                yscale=math.fabs(position[letter])
                #print("DEBUG:: draw_stack:: yscale ", yscale)
                dimensions=c_neg.text(x_run,y_neg,letter,
                                      [cs.getColor(letter),
                                       text.halign.center,
                                       text.size.Huge,
                                       color.transparency(0.7),
                                       trafo.scale(sx=0.8,sy=yscale)]
                                      )
                y_neg+=dimensions.height*yscale

        width=dimensions.width
        cc=canvas.canvas()
        cc.insert(c_pos)
        cc.insert(c_neg,[trafo.mirror()])
        self.drawingCanvas.insert(cc)
        #self.drawingCanvas.insert(c_neg,[trafo.rotate(180)])


        return width,0.0

    def drawStack(self,pos,lettermap,x_pos=0.0,y_pos=0.0,transparency=0.0,maxWidth=0.8):

        if not (isinstance(lettermap,OrderedDict) or isinstance(lettermap,dict)):
            print ("Please provide with valid data.\n "
                   "data should be a  dictionary of a stack of letter as key  containig relative height as valuefollowing format:\n "
                   "{'A':0.50,'T':0.50} \n")
            return None

        #maxWidth=0.0
        pos_canvas=canvas.canvas()
        cs=ColorScheme()
        x_run=x_pos
        y_run=y_pos

        for letter in lettermap:

            yscale=math.fabs(lettermap[letter])

            loc=pos_canvas.text(x_run,y_run,r"\textsf{"+letter+"}",
                                [ cs.getColor(letter),
                                  color.transparency(transparency),
                                  text.halign.center,
                                  text.size.large,
                                  trafo.scale(sx=maxWidth,sy=yscale)])
            #if maxWidth<=loc.width:
            #    maxWidth=loc.width

            y_run+=loc.height*yscale
            #print ("DEBUG:: drawStack y_run", y_run)
        #pos_canvas.text(x_run,1.5,pos,[text.size.small,color.cmyk.Black,trafo.rotate(45)])



        return pos_canvas,maxWidth

    #TODO: user specific Ruler
    def drawRuler(self,x=0.0,y=0.0,height=1.0,width=5.0,tic_size=0.5):

        no_tic=int(height/tic_size)


        #width=width-2
        if self.show_X_axis:
            x_axis=path.line(x-1,y,x+width,y)
            self.drawingCanvas.stroke(x_axis,[style.linewidth.Thick,color.cmyk.Black])

        if self.show_Y_axis:

            self.drawingCanvas.text(x-2.2,y+tic_size*3, "Bits",
                                    [text.size.tiny,
                                     text.valign.top,
                                     text.halign.center,
                                     color.cmyk.Black,
                                     trafo.rotate(90)])
            y_axis=path.line(x-1,y,x-1,y+height)
            self.drawingCanvas.stroke(y_axis,[style.linewidth.Thick,color.cmyk.Black])


            for i in range(0,no_tic):
                y=y+tic_size
                y_tic=path.line(x-1.2,y,x-0.8,y)
        #y_tic2=path.line(x-1.2,y+2,x-0.8,y+2)
                self.drawingCanvas.stroke(y_tic,[style.linewidth.Thick,color.cmyk.Black])
                self.drawingCanvas.text(x-1.5,y+0.2,i+1,
                                        [text.size.tiny,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         color.transparency(0.7),
                                         ])

            y_tic1=path.line(x-1.2,height,x-0.8,height)


            self.drawingCanvas.stroke(y_tic1,[style.linewidth.Thick,color.cmyk.Black])
        #self.drawingCanvas.stroke(box,[style.linewidth.Thick,color.cmyk.Black, trafo.mirror()])

    def mutationLogo(self,alignment,nativeSeq=None,x=0.0,y=0.0):
        blockSize=self.block_size
        maxHeight=self.stack_height
        maxWidth=self.stack_width
        transparency=self.transparency
        aln_len=len(alignment)
        aln_remain=aln_len
        counter=1

        if not isinstance(alignment,dict):
            print ("Please provide with valid data.\n "
                   "data should be a  dictionary following format:\n "
                   "{0:{'A':0.50,'T':0.50}} \n")
            return None
        if not nativeSeq is None:
            if not isinstance(nativeSeq,dict):
                print ("Please provide with valid data.\n "
                   "data should be a  dictionary following format:\n "
                   "{0:{'A':0.50,'T':0.50}} \n")
                return None

        block=0
        x_run=x
        y_run=y
        self.drawRuler(x_run,y_run,height=maxHeight,width=maxWidth*blockSize,tic_size=self.tic_size)
        
       
        
        if not nativeSeq is None:
            #Subfamily logo                      
            #for position in alignment:
            for position in alignment:
                nucleotideNumber = position+self.alignment1_index+1
                #print ("DEBUG:: ", (position + 1),alignment[position],x_run,y_run)
                #writing upper side                
                canvasAl,width1=self.drawStack(position,alignment[position],x_run,y_run)
                #writing down side   
                canvasNat,width2=self.drawStack(position,nativeSeq[position],x_run,math.fabs(y_run),transparency=transparency)

                self.drawingCanvas.insert(canvasAl)
                self.drawingCanvas.insert(canvasNat,[trafo.mirror()])

                #self.drawingCanvas.text(x_run,1.5,position,[text.size.small,color.cmyk.Black,color.transparency(0.8),trafo.scale(sx=1),trafo.rotate(90)])
                if self.show_position:

                    '''self.drawingCanvas.text(x_run,y_run-(maxHeight+1),position+self.alignment2_index+1,
                                        [text.size.tiny,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         color.transparency(transparency),
                                         trafo.rotate(-90)])'''

                    if(nucleotideNumber%5 == 0):

                        self.drawingCanvas.text(x_run,y_run+(maxHeight+1),nucleotideNumber,
                                            [text.size.tiny,
                                             text.valign.top,
                                             text.halign.center,
                                             color.cmyk.Black,
                                             color.transparency(transparency),
                                             trafo.rotate(-90)])

                #block+=1
                if counter%blockSize==0:

                    x_run=x
                    #block=position/blockSize
                    y_run=y_run-(maxHeight+2)-maxHeight-1.5
                    self.drawRuler(x_run,y_run,height=maxHeight,width=maxWidth*blockSize,tic_size=self.tic_size)

                else:
                    x_run+=maxWidth
                counter+=1

        else:
            #Mutation Logo germiline        
            for position in alignment:
                nucleotideNumber=position+self.alignment1_index
                canvasAl,width1=self.drawStack(alignment[position],x_run,y_run)

                #canvasNat,width2=self.drawStack(nativeSeq[position],x_run,y_run,opacity=0.5)
                self.drawingCanvas.insert(canvasAl)
                #self.drawingCanvas.insert(canvasNat,[trafo.mirror()])
                if(nucleotideNumber%5 == 0):
                    self.drawingCanvas.text(x_run,y_run-(maxHeight+1),position+self.alignment1_index,
                                            [text.size.tiny,
                                             text.valign.top,
                                             text.halign.center,
                                             color.cmyk.Black,
                                             color.transparency(transparency),
                                             trafo.rotate(-90)])

                if counter%blockSize==0:


                    x_run=x
                    #block=position/blockSize
                    y_run=y_run-(maxHeight+2)-maxHeight-1
                else:
                    x_run+=maxWidth
                counter+=1


        self.write(fileFormat=self.format)

    def drawSequenceLine(self,alignment,nativeSeq=None,name=None,x=2.0,y=0.0):
        blockSize=self.block_size
        maxHeight=self.stack_height
        maxWidth=self.stack_width
        transparency=self.transparency

        if not isinstance(alignment,dict):
            print ("Please provide with valid data.\n "
                   "data should be a  dictionary following format:\n "
                   "{0:{'A':0.50,'T':0.50}} \n")
            return None
        if not nativeSeq is None:
            if not isinstance(nativeSeq,dict):
                print ("Please provide with valid data.\n "
                   "data should be a  dictionary following format:\n "
                   "{0:{'A':0.50,'T':0.50}} \n")
                return None

        block=0
        x_run=x
        y_run=y
        self.drawRuler(x_run,y_run,height=maxHeight,width=maxWidth*blockSize,tic_size=self.tic_size)

        if not nativeSeq is None:
            for position in alignment:
                canvasAl,width1=self.drawStack(position,alignment[position],x_run,y_run)

                canvasNat,width2=self.drawStack(position,nativeSeq[position],x_run,math.fabs(y_run),transparency=transparency)

                self.drawingCanvas.insert(canvasAl)
                self.drawingCanvas.insert(canvasNat,[trafo.mirror()])

                #self.drawingCanvas.text(x_run,1.5,position,[text.size.small,color.cmyk.Black,color.transparency(0.8),trafo.scale(sx=1),trafo.rotate(90)])
                if self.show_position:

                    self.drawingCanvas.text(x_run,y_run-(maxHeight+1),position+self.alignment2_index+1,
                                        [text.size.tiny,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         color.transparency(transparency),
                                         trafo.rotate(-90)])

                    self.drawingCanvas.text(x_run,y_run+(maxHeight+1),position+self.alignment1_index+1,
                                        [text.size.tiny,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         color.transparency(transparency),
                                         trafo.rotate(-90)])

                #block+=1
                '''
                if position!=0 and (position+1)%blockSize==0:
                    x_run=x
                    #block=position/blockSize
                    y_run=y_run-(maxHeight+2)-maxHeight-1.5
                    self.drawRuler(x_run,y_run,height=maxHeight,width=maxWidth*blockSize,tic_size=self.tic_size)
                '''
                x_run+=maxWidth
        self.drawingCanvas.text(x_run+1,y_run+maxHeight/2,name,[text.size.tiny, text.valign.top, text.halign.center,color.cmyk.Red])
        y_run=y_run-(maxHeight+2)-maxHeight-1.5

        return y_run

    def compareLogo(self,aligns, germlines=None, names=None,x=2.0,y=0.0):
        if not isinstance(aligns,list):
            print ("Invalide data type")
            return -1
        if germlines==None:
            print("Invalid Data Type")
            return -1
        y_run=y
        x_run=x
        for a,g,n in zip(aligns,germlines,names):
            y_run=self.drawSequenceLine(a,g,n,x=x_run,y=y_run)

        self.write(fileFormat='pdf')

    def DCALogo(self,data,wildtype,x=0.0,y=0.0):
        blockSize=self.block_size
        maxHeight=self.stack_height
        maxWidth=self.stack_width
        transparency=self.transparency
        start=self.alignment1_index
        if blockSize>len(data):
            blockSize=len(data)


        #print (start)
        count=1
        x_run=x
        y_run=y
        #self.drawRuler(x_run,y_run,height=maxHeight,width=maxWidth*len(data))
        self.drawingCanvas.text(x_run,y_run-(maxHeight),start+1,
                                        [text.size.tiny,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         trafo.rotate(-90)])

        x_axis=path.line(x_run-1,y_run,maxWidth*blockSize,y_run)
        self.drawingCanvas.stroke(x_axis,[style.linewidth.Thick,color.cmyk.Black])

        x_axis=path.line(x_run-1,y_run-1.2,maxWidth*blockSize,y_run-1.2)
        self.drawingCanvas.stroke(x_axis,[style.linewidth.Thick,color.cmyk.Black])

        for pos in data:
            positive=data[pos]["pos"]
            negative=data[pos]["neg"]
            canvasPos,_=self.drawStack(pos,positive,x_pos=x_run,y_pos=y_run)

            canvasNeg,_=self.drawStack(pos,negative,x_pos=x_run,y_pos=math.fabs(y_run-1.2))

            self.drawingCanvas.insert(canvasPos)

            self.drawingCanvas.insert(canvasNeg,[trafo.mirror()])


            self.drawingCanvas.text(x_run,y_run-0.2,wildtype[pos],
                                        [text.size.small,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         color.transparency(transparency)]
                                    )
            #print (pos)
            if count%10==0:
                self.drawingCanvas.text(x_run,y_run-(maxHeight),pos,
                                        [text.size.tiny,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         trafo.rotate(-90)])

            if count%blockSize==0:
                x_run=x
                y_run=y_run-2*maxHeight+1

                x_axis=path.line(x_run-1,y_run,maxWidth*blockSize,y_run)
                self.drawingCanvas.stroke(x_axis,[style.linewidth.Thick,color.cmyk.Black])

                x_axis=path.line(x_run-1,y_run-1.2,maxWidth*blockSize,y_run-1.2)
                self.drawingCanvas.stroke(x_axis,[style.linewidth.Thick,color.cmyk.Black])

            else:
                x_run+=maxWidth
            count+=1
        self.write(fileFormat="pdf")

    def sequenceLogo(self,alignment,x=0.0,y=0.0):
        #print(self.block_size)
        blockSize=self.block_size
        maxHeight=self.stack_height
        maxWidth=self.stack_width
        transparency=self.transparency

        if not isinstance(alignment,dict):
            print ("Please provide with valid data.\n "
                   "data should be a  dictionary following format:\n "
                   "{0:{'A':0.50,'T':0.50}} \n")
            return None

        x_run=x
        y_run=y

        self.drawRuler(x_run,y_run,height=maxHeight,width=maxWidth*blockSize-2)

        for position in alignment:
                canvasAl,width1=self.drawStack(position,alignment[position],x_run,y_run)

                #canvasNat,width2=self.drawStack(nativeSeq[position],x_run,y_run,opacity=0.5)
                self.drawingCanvas.insert(canvasAl)
                #self.drawingCanvas.insert(canvasNat,[trafo.mirror()])
                if self.show_position:
                    self.drawingCanvas.text(x_run,y_run-1,position+self.alignment1_index,
                                        [text.size.tiny,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         color.transparency(transparency),
                                         trafo.rotate(-90)])

                if position!=0 and (position+1)%blockSize==0:
                    x_run=x
                    #block=position/blockSize
                    y_run=y_run-(maxHeight+2)
                    self.drawRuler(x_run,y_run,height=maxHeight,width=maxWidth*blockSize)
                else:
                    x_run+=maxWidth

        self.write(fileFormat=self.format)

    def plainLogo(self,alignment,x=0.0,y=0.0):
        # print(self.block_size)
        blockSize = self.block_size
        maxHeight = self.stack_height
        maxWidth = self.stack_width
        transparency = self.transparency
        aln_len=len(alignment)
        aln_remain=aln_len
        if not isinstance(alignment, dict):
            print("Please provide with valid data.\n "
                  "data should be a  dictionary following format:\n "
                  "{0:{'A':0.50,'T':0.50}} \n")
            return None

        x_run = x
        y_run = y

        self.drawRuler(x_run, y_run, height=maxHeight, width=maxWidth * blockSize)
        counter=1; 
        for position in alignment:

            canvasAl, width1 = self.drawStack(position, alignment[position], x_run, y_run)

            # canvasNat,width2=self.drawStack(nativeSeq[position],x_run,y_run,opacity=0.5)
            self.drawingCanvas.insert(canvasAl)
            # self.drawingCanvas.insert(canvasNat,[trafo.mirror()])
            if self.show_position:
                self.drawingCanvas.text(x_run, y_run - 1, position + self.alignment1_index+1,
                                        [text.size.tiny,
                                         text.valign.top,
                                         text.halign.center,
                                         color.cmyk.Black,
                                         color.transparency(transparency),
                                         trafo.rotate(-90)])

            #print("DEBUG:counter=", counter, "aln_remain=", aln_remain, "blockSize=", blockSize, "posotion", position)
            if counter%blockSize==0 and counter < len(alignment):
            	#print ("DEBUG:Entrou")
            	if aln_remain < blockSize:
            		x_run = x
            		# block=position/blockSize
            		y_run = y_run - (maxHeight + 2)
            		self.drawRuler(x_run, y_run, height=maxHeight, width=maxWidth * blockSize)
            	else:
            		x_run = x
            		# block=position/blockSize
            		y_run = y_run - (maxHeight + 2)
            		self.drawRuler(x_run, y_run, height=maxHeight, width=maxWidth * blockSize)
            else:
            	x_run += maxWidth
            aln_remain = aln_remain - counter
            counter+=1

            """


            if position != 0 and (position + 1) % blockSize == 0:
                x_run = x
                # block=position/blockSize
                y_run = y_run - (maxHeight + 2)
                self.drawRuler(x_run, y_run, height=maxHeight, width=maxWidth * blockSize)
            else:
                x_run += maxWidth
            """

        self.write(fileFormat=self.format)

    def subfamilyLogo(self,alignment1, alignment2, N=0):
        print ("DEGUG:: subfamilyLogo")    
        y=self.mutationLogo(alignment1,alignment2)
        self.write(fileFormat='pdf')


    def write(self,fileFormat="all"):
        #self.insert(svgfile.svgfile(0, 0, self.outfile)
        self.drawingCanvas.writeSVGfile(self.outfile)
        if fileFormat.lower()=="eps":
            self.drawingCanvas.writeEPSfile(self.outfile)
        elif fileFormat.lower()=="pdf":
            self.drawingCanvas.writePDFfile(self.outfile)
        elif fileFormat.lower()=="svg":
            pass
        else:
            self.drawingCanvas.writeEPSfile(self.outfile)
            self.drawingCanvas.writePDFfile(self.outfile)


if __name__ == '__main__':
    N=Sequencelogo()
    p, n=subfamily("Dataset/test1.txt", "Dataset/test2.txt")
    N.subfamilyLogo(p,n)
