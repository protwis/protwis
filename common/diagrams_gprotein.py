from common.diagrams import Diagram
from common.definitions import G_PROTEIN_SEGMENTS

from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

from django.utils.safestring import mark_safe

from math import cos, sin, pi, floor,sqrt
from datetime import datetime
from collections import OrderedDict

class DrawGproteinPlot(Diagram):

    def __init__(self, residue_list, protein_class, protein_name, nobuttons = None):
        self.nobuttons = 'gprotein'
        self.type = 'snakeplot'
        plot_data = {}
        plot_data['direction'] = [0, 0, 1, 0, 1, 0, 1, 0]  # 0: EC->IC, 1: IC->EC
        plot_data['helixRadius'] = 70

        self.receptorId = protein_name
        self.family = protein_class
        self.output = ''

        # FIXME DO PUREIMAGE
        # $pureImage = isset($_GET['pureimage']) && $_GET['pureimage'] == 'TRUE' ? TRUE : FALSE;

        # get sequence, baldwin, and bw information of this receptor

        self.sequence = residue_list
        self.segments = {}
        self.segments_full = OrderedDict()
        i = 0
        for r in self.sequence:
            if r.protein_segment:
                segment = str(r.protein_segment.slug)
            elif r.segment_slug:  #from family aligment
                segment = str(r.segment_slug)

            if segment not in self.segments:
                self.segments[segment] = []
                self.segments_full[segment] = r.protein_segment
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
            self.segments[segment].append([r.sequence_number, r.amino_acid,label,displaylabel])
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
        self.offsetX = 0 #-200
        self.offsetY = 0 #-50

        # margin between two helixes
        self.margin = 10

        # highest and lowest bound of this svg
        self.high =0
        self.low = 0

        # keep track of max Y positions of intra/extra loops
        self.maxY = {'bottom':0,'top':0}
        self.maxX = {'left':0,'right':0}

        # helices length
        # helicesLength = Svg::getSnakePlotHelicesLength($baldwin, $helixWidth, $angleDeg) #FIXME

        # top and bottom residue coords in each helix
        self.TBCoords = {}

        self.output = ""
        self.traceoutput = ""
        self.helixoutput = ""

        # Draw sheets and helices
        self.count = 1
        self.count_sheet = 0
        for s in G_PROTEIN_SEGMENTS['Full']:
            if s in self.segments_full:
                if self.segments_full[s].category=='helix':
                    self.helixoutput += self.drawSnakePlotHelix(s)
                    self.count += 1
                elif self.segments_full[s].category=='sheet':
                    self.helixoutput += self.drawSnakePlotSheet(s)
                    self.count += 1
                    self.count_sheet += 1
        # Draw loops
        self.count = 0
        for s in G_PROTEIN_SEGMENTS['Full']:
            if s in self.segments_full and self.segments_full[s].category=='loop':
                #pass
                self.drawSnakePlotLoop(s)
            else:
                self.count += 1



    def __str__(self):

        self.output = "<g id=snake transform='translate(0, " + str(-self.low+ self.offsetY) + ")'>" + self.traceoutput+self.output+self.helixoutput+self.drawToolTip() + "</g>"; #for resizing height
        return mark_safe(self.create(self.output,self.maxX['right']+30,self.high-self.low+self.offsetY*2,"snakeplot", self.nobuttons))

    def drawSnakePlotHelix(self, segment):
        rs = self.segments[segment]
        helix_num = self.count
        self.TBCoords[helix_num] = {}

        if helix_num%2!=0: rs.reverse() # reverse direction for even helix because they go from inside to outside

        output_residues = []

        res_num = len(self.segments[segment])
        output_residue_in = ''
        output_residue_out = ''
        output_trace = ''

        startX = self.helixWidth+self.offsetX+(self.margin+self.helixWidth)*(helix_num-1)-(self.count_sheet*20)
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
            if ((helix_num%2==0 and prevGeneric_number+'1'==currGeneric_number) or (helix_num%2!=0 and str(int(prevGeneric_number)-1)+'1'==currGeneric_number)) and i!=0:
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

            if i==0: self.TBCoords[helix_num]['top'] = [x,y]
            if i==res_num-1: self.TBCoords[helix_num]['bottom'] = [x,y]


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

    def drawSnakePlotSheet(self, segment):
        rs = self.segments[segment]
        helix_num = self.count
        self.TBCoords[helix_num] = {}

        if helix_num%2!=0: rs.reverse() # reverse direction for even helix because they go from inside to outside

        output_residues = []

        res_num = len(self.segments[segment])
        output_residue_in = ''
        output_residue_out = ''
        output_trace = ''

        startX = 50+self.offsetX+(self.margin+self.helixWidth)*(helix_num-1)-(self.count_sheet*20)
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
                    #output_trace += "<line x1="+str(prevX)+" y1="+str(prevY)+" x2="+str(nextX)+" y2="+str(nextY)+" stroke='grey' fill='none' stroke-width='1' stroke-dasharray='1,1' />"
                    row_pos +=1
                elif row_pos+1==row_length:
                    row+=1
                    row_pos=0
                    row_length = 3 if row_length == 4 else 4
                else:
                    row_pos +=1

            # move left as you go down a row
            x = round(startX) #+indentX+bulgeX

            # Move down with right amount
            y = round(startY+i*self.residue_radius*1.5)
            output_residue = self.DrawResidueSquare(x,y,rs[i][1], rs[i][0], rs[i][3], self.residue_radius)


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

            if i==0: self.TBCoords[helix_num]['top'] = [x,y]
            if i==res_num-1: self.TBCoords[helix_num]['bottom'] = [x,y]


            if (row_pos==1 and row!=0) or (skip==1 and row_pos==2): # if need for trace
                if row_length==3: points = "M "+str(prevX)+" "+str(prevY)+" Q"+str(prevX-40)+" "+str(prevY+30)+", "+str(x-21)+" "+str(y-8)+" T"+str(x)+" "+str(y)
                if row_length>=4: points = "M "+str(prevX)+" "+str(prevY)+" Q"+str(prevX-40)+" "+str(prevY+30)+", "+str(x-24)+" "+str(y-7)+" T"+str(x)+" "+str(y)
                # output_trace += "<path d='" + points + "' stroke='grey' fill='none' stroke-width='2'  />"

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

    def drawSnakePlotLoop(self, segment):

        y_offset = 50
        font_size = 12
        font_family = 'courier'
        bezier_pull = 80
        name = segment
        x_at_max_y = 0

        rs = self.segments[segment] # get residues

        start = 1
        res_before = []
        res_helix = []
        res_after = []

        if self.count % 2 == 0:
            position = 'bottom'
            orientation = 1
        else:
            position = 'top'
            orientation = -1

        if self.count not in self.TBCoords:
            return 0

        if self.count+1 not in self.TBCoords:
            return 0

        # Get positions of two  linking residues from each helix
        x1 = self.TBCoords[self.count][position][0]
        y1 = self.TBCoords[self.count][position][1]
        x2 = self.TBCoords[self.count+1][position][0]
        y2 = self.TBCoords[self.count+1][position][1]

        boxX = (x1+x2)/2 # midway between
        if position=='top':
            boxY = min(y1,y2)-y_offset # over helix
            y_indent = -1*bezier_pull
        if position=='bottom':
            boxY = max(y1,y2)+y_offset # over helix
            y_indent = bezier_pull

        points = str(x1)+","+str(y1)+" "+str(boxX)+","+str(boxY)+" "+str(x2)+","+str(y2)
        points2 = "M "+str(x1)+" "+str(y1)+" Q"+str(boxX)+" "+str(boxY+y_indent)+" "+str(x2)+" "+str(y2)

        # Getting midpoint of Bezier curve http://www.svgbasics.com/curves.html
        Dx = ((x1+boxX)/2)
        Ex = ((x2+boxX)/2)
        Fx = (Dx+Ex)/2

        Dy = ((y1+boxY+y_indent)/2)
        Ey = ((y2+boxY+y_indent)/2)
        Fy = (Dy+Ey)/2

        #JUST SIMPLE
        #self.output += "<path class='"+name+" short' d='" + points2 + "' stroke='black' fill='none' stroke-width='2' />"
        # self.output += "<rect onclick='toggleLoop(\"."+name+"\",\"short\");' class='"+name+" short' x="+str(Fx-18)+" y="+str(Fy-13)+" rx=5 ry=5 width='35' height='20' stroke='black' fill='white' stroke-width='1' style2='fill:red;stroke:black;stroke-width:5;opacity:0.5'/>"
        # self.output += str("<text  onclick='toggleLoop(\"."+name+"\",\"short\");' class='"+name+" short' x="+str(Fx)+" y="+str(Fy)+" text-anchor='middle' font-size="+str(font_size)+" font-family='"+font_family+"'>"+name+"</text>")


        y_indent = y_indent*len(rs)/5 # get an approx need for y_indent for size of loop

        loop_long_length = 0
        super_loop_long_length = 40
        between_residues = 18

        length_of_residues_in_loop = len(rs)*between_residues-self.residue_radius
        length = self.lengthbezier([x1,y1],[boxX,boxY+y_indent],[x2,y2],0.001)

        if len(rs)<super_loop_long_length:
            tries = 0 # adjust size
            while abs(length-length_of_residues_in_loop-70)>5:
                # print(abs(length-length_of_residues_in_loop+100),length,length_of_residues_in_loop,tries)
                if length-length_of_residues_in_loop-70>5:
                    y_indent *=0.9
                else:
                    y_indent *=1.1
                length = self.lengthbezier([x1,y1],[boxX,boxY+y_indent],[x2,y2],0.001)

                tries += 1
                if tries>100:
                    break

        pos = (length-length_of_residues_in_loop)/2 # get start pos

        indentX = 0
        indentY2 = 0
        prev_where = [x1,y1]

        # make rounded arc
        points2 = "M "+str(x1)+" "+str(y1)+" Q"+str(boxX)+" "+str(boxY+y_indent)+" "+str(x2)+" "+str(y2)
        labelbox = self.wherebezier([x1,y1],[boxX,boxY+y_indent],[x2,y2],0.001,length/2)

        labelbox[1][1] += orientation*40

        self.output += "<path class='"+name+"' d='" + points2 + "' stroke='black' fill='none' stroke-width='2' />"

        max_y = y1
        for i in range(0,len(rs)):
            r = rs[i]
            where = self.wherebezier([x1,y1],[boxX,boxY+y_indent],[x2,y2],0.001,pos)

            self.output += self.DrawResidue(where[1][0],where[1][1],r[1], r[0], r[3], self.residue_radius-1,name)
            pos += between_residues

            if where[1][1]>self.high: self.high = where[1][1]
            if where[1][1]<self.low: self.low = where[1][1]
            prev_where = where[1][0],where[1][1]

            if orientation==-1:
                if where[1][1]<self.maxY[position]: self.maxY[position] = where[1][1]
            else:
                if where[1][1]>self.maxY[position]: self.maxY[position] = where[1][1]

            if orientation==-1:
                if where[1][1]<max_y:
                    max_y = where[1][1]
                    x_at_max_y = where[1][0]
            else:
                if where[1][1]>max_y:
                    max_y = where[1][1]
                    x_at_max_y = where[1][0]
            x_at_max_y = where[1][0]

        if orientation==1:
            max_y = max_y+25
        else:
            max_y = max_y-20
        self.output += "<rect onclick='toggleLoop(\"."+name+"\",\"long\");' class='"+name+"' x="+str(x_at_max_y-25)+" y="+str(max_y-13)+" rx=5 ry=5 width='50' height='20' stroke='black' fill='white' stroke-width='1' style2='fill:red;stroke:black;stroke-width:5;opacity:0.5'/>"
        self.output += str("<text  onclick='toggleLoop(\"."+name+"\",\"long\");' class='"+name+"' x="+str(x_at_max_y)+" y="+str(max_y)+" text-anchor='middle' font-size="+str(font_size)+" font-family='"+font_family+"'>"+name+"</text>")
