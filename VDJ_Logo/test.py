from pyx import canvas,unit,path,color,text,trafo,style,graph
from ColorScheme import ColorScheme
from pyx.graph.axis import linear,texter

c = canvas.canvas()
c.text(0, 0, "Hello, world!")
c.stroke(path.line(0, 0, 2, 0))
#c.writeEPSfile()
#c.writePDFfile()
c.insert(0,0, "test.svg")

c.writeSVGfile("output.svg")

