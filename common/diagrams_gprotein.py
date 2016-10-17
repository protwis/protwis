from common.diagrams import Diagram

from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

from django.utils.safestring import mark_safe

from math import cos, sin, pi, floor,sqrt
from datetime import datetime


class DrawGproteinPlot(Diagram):

    def __init__(self, residue_list, protein_class, protein_name, nobuttons = None):
        self.nobuttons = 'gprotein'
        self.type = 'snakeplot'
        plot_data = {}
        plot_data['direction'] = [0,0, 1, 0, 1, 0, 1, 0]; # 0: EC->IC, 1: IC->EC
        plot_data['helixRadius'] = 70

        self.receptorId = protein_name
        self.family = protein_class
        self.output = ''
        residueType = 'sp'

        # FIXME DO PUREIMAGE
        pureImage = False
        #$pureImage = isset($_GET['pureimage']) && $_GET['pureimage'] == 'TRUE' ? TRUE : FALSE;

        # get sequence, baldwin, and bw information of this receptor

        self.sequence = residue_list
        self.segments = {}
        i = 0
        for r in self.sequence:
            if r.protein_segment:
                segment = str(r.protein_segment.slug)
            elif r.segment_slug: #from family aligment
                segment = str(r.segment_slug)

            if segment not in self.segments: self.segments[segment] = []
            label = ''
            displaylabel = ''
            if r.generic_number:
                label = r.generic_number.label
            elif hasattr(r, 'family_generic_number'):
                label = r.family_generic_number
            if r.display_generic_number: displaylabel = r.display_generic_number.label
            displaylabel = r.amino_acid + str(r.sequence_number) + " \n " + displaylabel
            if hasattr(r, 'frequency'):
                displaylabel = displaylabel + "\n" + r.frequency
            self.segments[segment].append([r.sequence_number,r.amino_acid,label,displaylabel])
            i += 1

        # for helix_num in range(1,2): #FIX for missing generic numbers
        #     rs = self.segments['H5']
        #     for i in range(0,len(rs)):
        #         if not rs[i][2]:
        #             if i+1<len(rs): #if there is a next one
        #                 if rs[i+1][2]: #if it has generic number
        #                     number = str(int(rs[i+1][2].split('x')[1])-1)
        #                     rs[i][2] = str(helix_num) + "x" + number
        #                     print(rs[i][2])

        self.helixWidth = 85           # Width of helix
        self.resNumPerRow = 4          # Residue number per row in helix
        self.angleDeg = 22.0           # Angle size of each helix turn
        self.residue_radius = 12     # Radius of the residue circle

        # svg image padding offset
        self.offsetX = -200
        self.offsetY = -50

        # margin between two helixes
        self.margin = 0

        # highest and lowest bound of this svg
        self.high =0
        self.low = 0

        # keep track of max Y positions of intra/extra loops
        self.maxY = {'intra':0,'extra':0}
        self.maxX = {'left':0,'right':0}

        # helices length
        # helicesLength = Svg::getSnakePlotHelicesLength($baldwin, $helixWidth, $angleDeg) #FIXME

        # top and bottom residue coords in each helix
        self.TBCoords = {}

        self.output = ""
        self.traceoutput = ""
        self.helixoutput = ""

        self.helixoutput += self.drawSnakePlotHelix('H5')


    def __str__(self):

        self.output = "<g id=snake transform='translate(0, " + str(-self.low+ self.offsetY) + ")'>" + self.traceoutput+self.output+self.helixoutput+self.drawToolTip() + "</g>"; #for resizing height
        return mark_safe(self.create(self.output,self.maxX['right']+30,self.high-self.low+self.offsetY*2,"snakeplot", self.nobuttons))

    def drawSnakePlotHelix(self, segment):
        rs = self.segments[segment]
        helix_num = 1
        self.TBCoords[helix_num] = {}

        if helix_num%2!=0: rs.reverse() # reverse direction for even helix because they go from inside to outside

        output_residues = []

        res_num = len(self.segments[segment])
        output_residue_in = ''
        output_residue_out = ''
        output_trace = ''

        startX = self.helixWidth+self.offsetX+(self.margin+self.helixWidth)*(helix_num-1)
        startY = self.offsetY


        row_length = 3
        row_pos = 0
        row = 0
        prevGeneric = '0.0.0'
        bulgeX = 0
        bulgeY = 0
        bulge = 0
        skip = 0
        indentX = -self.residue_radius+3
        indentY = 3
        for i in range(0,res_num):
            prevGeneric_number = prevGeneric.split('.')[2]
            currGeneric_number = rs[i][2].split('.')[2]
            if (helix_num%2==0 and prevGeneric_number+'1'==currGeneric_number) or (helix_num%2!=0 and str(int(prevGeneric_number)-1)+'1'==currGeneric_number):
                bulge = 1
                if row_pos==0:  # if first in row, use space for bulge
                    bulgeY = 5
                    bulgeX = 7
                else:
                    bulgeY = 5
                    bulgeX = 5
                row_length+=1
            elif i!=0 and ((helix_num%2!=0 and int(prevGeneric_number)-1!= int(currGeneric_number)) or (helix_num%2==0 and int(prevGeneric_number)+1!= int(currGeneric_number))):
                skip = 1
                if row_pos!=0 and row_pos+1<row_length:
                    nextX =round(startX-(row_pos+1)*self.residue_radius*1.5+indentX+bulgeX)
                    nextY = round(startY+row*self.residue_radius*2.4+(row_pos+1)*self.residue_radius*0.5+indentY+bulgeY)
                    output_trace += "<line x1="+str(prevX)+" y1="+str(prevY)+" x2="+str(nextX)+" y2="+str(nextY)+" stroke='grey' fill='none' stroke-width='1' stroke-dasharray='1,1' />"
                    row_pos +=1
                elif row_pos+1==row_length:
                    row+=1
                    row_pos=0
                    row_length = 3 if row_length == 4 else 4
                else:
                    row_pos +=1

            # move left as you go down a row
            x = round(startX-row_pos*self.residue_radius*1.6+indentX+bulgeX)

            # Move down with right amount
            y = round(startY+row*self.residue_radius*2.4+row_pos*self.residue_radius*0.5+indentY+bulgeY)
            output_residue = self.DrawResidue(x,y,rs[i][1], rs[i][0], rs[i][3], self.residue_radius)


            if x<self.maxX['left']: self.maxX['left'] = x
            if x>self.maxX['right']: self.maxX['right'] = x

            row_pos += 1
            if bulge==1:
                if row_pos==1:  # if first in row, use space for bulge
                    bulgeY = -3
                    bulgeX = 10
                else:
                    bulgeY = -3
                    bulgeX = 7
                rs[i][2] = prevGeneric # make it the prev one, to catch missing ones correctly
                bulge = 0

            if row_length==3:
                output_residue_in += output_residue
            else:
                output_residue_out += output_residue

            output_residues.append(output_residue)

            if i==0: self.TBCoords[helix_num]['extra'] = [x,y]
            if i==res_num-1: self.TBCoords[helix_num]['intra'] = [x,y]


            if (row_pos==1 and row!=0) or (skip==1 and row_pos==2): # if need for trace
                if row_length==3: points = "M "+str(prevX)+" "+str(prevY)+" Q"+str(prevX-40)+" "+str(prevY+30)+", "+str(x-21)+" "+str(y-8)+" T"+str(x)+" "+str(y)
                if row_length>=4: points = "M "+str(prevX)+" "+str(prevY)+" Q"+str(prevX-40)+" "+str(prevY+30)+", "+str(x-24)+" "+str(y-7)+" T"+str(x)+" "+str(y)
                output_trace += "<path d='" + points + "' stroke='grey' fill='none' stroke-width='2'  />"

            # alternate between 4 and 3 res per row
            if row_length>3 and row_pos>=row_length:
                row_length=3
                row_pos = 0
                row += 1
                bulgeX = 0
                bulgeY = 0
                indentX = -self.residue_radius+3
                indentY = 3
            elif row_length==3 and row_pos>=3:
                row_length=4
                row_pos = 0
                row += 1
                bulgeX = 0
                bulgeY = 0
                indentX = 0
                indentY = 0

            skip = 0
            prevX = x
            prevY = y
            prevGeneric = rs[i][2]

        temp = ''
        if helix_num%2!=0: output_residues.reverse()
        for res in output_residues:
            temp += res

        return output_trace+temp
