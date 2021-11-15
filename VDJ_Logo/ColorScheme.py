"""
Author: Bishnu Sarker, Research Intern at Laboratory of Computational and Quantitative Biology, UPMC, Paris, France.
Email: bishnukuet@gmail.com

"""

from pyx import color

Default_Color={"A":color.cmyk.Periwinkle,
               "R":color.cmyk.BurntOrange,
               "N":color.cmyk.Maroon,
               "D":color.cmyk.Red,
               "B":color.cmyk.Orange,
               "C":color.cmyk.RedViolet,
               "E":color.cmyk.Lavender,
               "Q":color.cmyk.Emerald,
               "Z":color.cmyk.Violet,
               "G":color.cmyk.MidnightBlue,
               "H":color.cmyk.Blue,
               "I":color.cmyk.Cyan,
               "L":color.cmyk.PineGreen,
               "K":color.cmyk.Sepia,
               "M":color.cmyk.Green,
               "F":color.cmyk.Tan,
               "P":color.cmyk.Orchid,
               "S":color.cmyk.Black,
               "T":color.cmyk.Salmon,
               "W":color.cmyk.TealBlue,
               "Y":color.cmyk.OliveGreen,
               "V":color.cmyk.Fuchsia,
               "-":color.cmyk.White,
               ".":color.cmyk.White,
               "X":color.cmyk.Magenta
               }
ClustalScheme={i:color.cmyk.Orange for i in ["G","P","S","T"] }
ClustalScheme.update({i:color.cmyk.Red for i in ["H","K","R"]})
ClustalScheme.update({i:color.cmyk.Blue for i in ["F","W","Y"]})
ClustalScheme.update({i:color.cmyk.Green for i in ["I","L","M","V"]})
ClustalScheme.update({i:color.cmyk.White for i in [".","-"]})

class ColorScheme():
    def __init__(self,setColors=Default_Color):
        self.Colors=setColors
    def getColor(self,AminoAcidCode):
        return self.Colors[AminoAcidCode]

