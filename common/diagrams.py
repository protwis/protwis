from math import cos, sin, tan, pi, sqrt, pow
import string, time, math, random

def uniqid(prefix='', more_entropy=False):
    m = time.time()
    uniqid = '%8x%05x' %(int(math.floor(m)),int((m-math.floor(m))*1000000))
    if more_entropy:
        valid_chars = list(set(string.hexdigits.lower()))
        entropy_string = ''
        for i in range(0,10,1):
            entropy_string += random.choice(valid_chars)
        uniqid = uniqid + entropy_string
    uniqid = prefix + uniqid
    return uniqid

class Diagram:
    def create(self, content,sizex,sizey,name, nobuttons):
        #diagram_js = self.diagramJS()
        if nobuttons=='gprotein' or nobuttons=='arrestin':
            return ("<svg id=\""+name+"\" " +
            "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\""+str(sizex)+"\" height=\""+str(sizey)+"\" " +
            "style='stroke-width: 0px; background-color: white;'>\n"+content+"</svg>" +
            self.drawColorPanel(nobuttons)) #width=\"595\" height=\"430\"
        elif nobuttons:
            return ("<svg id=\""+name+"\" " +
            "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\""+str(sizex)+"\" height=\""+str(sizey)+"\" " +
            "style='stroke-width: 0px; background-color: white;'>\n"+content+"</svg>") #width=\"595\" height=\"430\"
        else:
            return ("<svg id=\""+name+"\" " +
            "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\""+str(sizex)+"\" height=\""+str(sizey)+"\" " +
            "style='stroke-width: 0px; background-color: white;'>\n"+content+"</svg>" +
            self.drawColorPanel()) #width=\"595\" height=\"430\"

    def drawToolTip(self):
        output2 = """<g id='tool-tip-{}' transform='translate(0,0)' visibility='hidden'>
            <rect x='0' y='-40' width='1' height='25' stroke='black' fill='white' stroke-width='1' />
            <text x='0' y='-23' text-anchor='middle' font-family='Arial' font-size='12' fill='black'></text>
            </g>""".format(self.type)
        output= ""

        return output

    def drawColorPanel(self, nobuttons=None):

        boxstyle = """<style>
        .pick-color  {
          display:inline-block;
          width: 35px;
          height: 20px;
          margin: 1px;
          border-radius: 5px;
          border: 2px solid #000;
        }

        .long {
          display: none
        }
        .tooltip-inner {
            white-space:pre-wrap;
            max-width:none;
        }

        </style>
        """

        presetColors = {'D': ['#E60A0A', '#FDFF7B'],'E': ['#E60A0A', '#FDFF7B'],
                                    'K': ['#145AFF', '#FDFF7B'],'R': ['#145AFF', '#FDFF7B'],
                                    'S': ['#A70CC6', '#FDFF7B'],'T': ['#A70CC6', '#FDFF7B'],
                                    'N': ['#A70CC6', '#FDFF7B'],'Q': ['#A70CC6', '#FDFF7B'],
                                    'V': ['#E6E600', '#000000'],'L': ['#E6E600', '#000000'],
                                    'I': ['#E6E600', '#000000'],'A': ['#E6E600', '#000000'],
                                    'M': ['#E6E600', '#000000'],'F': ['#18FF0B', '#000000'],
                                    'Y': ['#18FF0B', '#000000'],'W': ['#0BCF00', '#000000'],
                                    'H': ['#0093DD', '#000000'],'P': ['#CC0099', '#FDFF7B'],
                                    'C': ['#B2B548', '#000000'],'G': ['#FF00F2', '#000000'],
                                    '-': ['#FFFFFF', '#000000']
                                    }
        fillcolors = [['#CCCCCC', '#000000']]
        for key,value in presetColors.items():
            if value not in fillcolors:
                fillcolors.append(value)

        colors = "<div id=\"cp2_"+self.type+"\" class=\"\">"
        for color in fillcolors:
            colors += "<div class='pick-color "+self.type+" selected' id='pick-"+color[0]+"-"+color[1]+"' style='background-color: "+color[0]+";'>&nbsp;</div>"

        colors += ""
        colors += "<div class=\"input-group-addon\" style='width:0;padding:0;border:0;background-color:0;display:inline;position:relative;top:-6px'><i class='pick-color "+self.type+" selected'></i></div>"
        colors += "<input type=\"text\" id='custom_color_"+self.type+"' value=\"#00AABB\" class=\"\" size=8 />"
        colors += ""

        # colors += "<span id=\"cp2_"+self.type+"\" class=\"\">"
        # colors += "<span class=\"pick-color\" "+self.type+" id='pick-"+fillcolors[0][0]+"-"+fillcolors[0][1]+"' style='background-color: "+fillcolors[0][0]+";'><i></i></span>"
        # colors += "<input type=\"text\" value=\"#00AABB\" class=\"\" size=8 />"
        # colors += "</span>"


        output = ("<br>Pick color:" +
            colors +
            "</div>")

        output += '<br><button style="width:120px;" onclick="applyPresentColors(\''+self.type+'\')">Properties</button> <button style="width:120px;" onclick="resetColors(\''+self.type+'\')">Clear</button>'
        if str(self.receptorId)=='family_diagram_preloaded_data':
            output += '<br><button style="width:220px;" onclick="ajaxMutantsPos(\''+self.type+'\');">Show Invitro Mutants</button>'
            output += ' <button style="width:220px;" onclick="ajaxInteractionsPos(\''+self.type+'\')">Show Interactions from structures</button>'
            output += '<br><button style="width:220px;" onclick="ajaxNaturalMutationPos(\''+self.type+'\')">Show Natural Genetic Variations</button>'
            output += ' <button style="width:120px;" onclick="ajaxPTMPos(\''+self.type+'\')">Show PTM sites</button>'
            # output += ' <button style="width:220px;" onclick="ajaxCancerMutationPos(\''+self.type+'\')">Show Cancer Mutations</button>'
            # output += ' <button style="width:220px;" onclick="ajaxDiseaseMutationPos(\''+self.type+'\')">Show Disease Mutations</button>'
        else:
            if nobuttons == 'gprotein':
                output += ' <button style="width:220px;" onclick="ajaxInterface(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Receptor Interface</button>'
                output += '<br><button style="width:220px;" onclick="ajaxNaturalMutation(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Natural Genetic Variations</button>'
                output += ' <button style="width:220px;" onclick="ajaxPTMs(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show PTM sites</button>'
                # output += ' <button style="width:220px;" onclick="ajaxCancerMutation(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Cancer Mutations</button>'
                # output += ' <button style="width:220px;" onclick="ajaxDiseaseMutation(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Disease Mutations</button>'
                if 'human' in self.receptorId: # only show barcode button for human gproteins
                    output += '<br><button style="width:120px;" onclick="ajaxBarcode(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Barcode</button>'
            elif nobuttons == 'arrestin':
                output += ' <button style="width:220px;" onclick="ajaxInterface(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Receptor Interface</button>'
                output += '<br><button style="width:220px;" onclick="ajaxNaturalMutation(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Natural Genetic Variations</button>'
                output += ' <button style="width:220px;" onclick="ajaxPTMs(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show PTM sites</button>'

            else:
                output += '<br><button style="width:220px;" onclick="ajaxMutants(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Invitro Mutants</button>'
                output += ' <button style="width:220px;" onclick="ajaxInteractions(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Interactions from Structures</button>'
                output += '<br><button style="width:220px;" onclick="ajaxNaturalMutation(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Natural Genetic Variations</button>'
                output += ' <button style="width:120px;" onclick="ajaxPTMs(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show PTM sites</button>'
                # output += ' <button style="width:220px;" onclick="ajaxCancerMutation(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Cancer Mutations</button>'
                # output += ' <button style="width:220px;" onclick="ajaxDiseaseMutation(\''+self.type+'\',\''+str(self.receptorId)+'\')">Show Disease Mutations</button>'

        if nobuttons != 'gprotein' and nobuttons != 'arrestin':
            output += '<br><small>Invitro Mutant Data: Increased binding/potency: <font style="color: #000; background-color: #87E88F" color="#87E88F">>5-fold</font>, <font style="color: #000; background-color: #66B36C" color="#66B36C">>10-fold</font>; Reduced binding/potency: <font style="color: #FFF; background-color: #FF7373" color="#FF7373">>5-fold</font>, <font style="color: #FDFF7B; background-color: #FA1111" color="#FA1111">>10-fold</font>; <font style="color: #000; background-color: #F7DA00" color="#F7DA00">No/low effect (<5-fold)</font>; and <font style="color: #000; background-color: #D9D7CE" color="#D9D7CE">N/A</font> </small>'

        return boxstyle+ output

    #Draws a ring of a helical wheel
    def DrawResidue(self, x,y,aa,residue_number,label,radius, resclass = '',cfill="white", precolor = False):
        id = residue_number
        idtext = str(id) + 't'
        tfill = 'black'
        x = round(x)
        y = round(y)
        #if (isset(_GET['precolor']) && _GET['precolor'] == 'TRUE') precolor = TRUE
        # if (precolor) {
        #     iid = str_replace('.', '_', id)
        #     iidtext = str_replace('.', '_', idtext)
        #     cfill = isset(_SESSION['color_pattern'][iid]) ? _SESSION['color_pattern'][iid] : 'white'
        #     tfill = isset(_SESSION['color_pattern'][iidtext]) ? _SESSION['color_pattern'][iidtext] : 'black'
        # }
        output =  """
            <circle class='{} rcircle' cx='{}' cy='{}' r='{}' stroke='black' stroke-width='2' fill='{}'
            fill-opacity='1' id='{}' title='{}' original_title='{}' original_cx='{}' original_cy='{}'/>
            <text x='{}' y='{}' text-anchor='middle' dominant-baseline='middle' font-family='helvetica' font-size='16' fill=''
            id='{}' class='rtext {}' title='{}' original_title='{}' original_x='{}' original_y='{}'> {} </text>
            """.format(resclass,x,y,radius,cfill,id,label,label,x,y,x,y+2,idtext,resclass,label,label,x,y+2,aa) #aa
        return output

    def DrawResidueSquare(self, x,y,aa,residue_number,label,radius, resclass = '',cfill="white", precolor = False):
        id = residue_number
        idtext = str(id) + 't'
        tfill = 'black'
        #if (isset(_GET['precolor']) && _GET['precolor'] == 'TRUE') precolor = TRUE
        # if (precolor) {
        #     iid = str_replace('.', '_', id)
        #     iidtext = str_replace('.', '_', idtext)
        #     cfill = isset(_SESSION['color_pattern'][iid]) ? _SESSION['color_pattern'][iid] : 'white'
        #     tfill = isset(_SESSION['color_pattern'][iidtext]) ? _SESSION['color_pattern'][iidtext] : 'black'
        # }
        output =  """
            <rect class='{} rcircle' x='{}' y='{}' height='{}' width='{}' stroke='black' stroke-width='1' fill='{}'
            fill-opacity='1' id='{}' title='{}' original_title='{}'/>
            <text x='{}' y='{}' text-anchor='middle' font-family='helvetica' font-size='16' fill=''
            id='{}' class='rtext {}' title='{}' original_title='{}'> {} </text>
            """.format(resclass,x-radius*(1.5/2),y-radius*(1.5/2),radius*1.5,radius*1.5,cfill,id,label,label,x,y+6,idtext,resclass,label,label,aa) #aa
        return output
    def deg2rad(self,degrees):
        radians = pi * degrees / 180
        return radians

    def bezier(self,p0,p1,p2,t):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve

        v1x = p1[0]-p0[0]
        v1y = p1[1]-p0[1]

        i1 = [p0[0]+(p1[0]-p0[0])*t,p0[1]+(p1[1]-p0[1])*t]
        i2 = [p1[0]+(p2[0]-p1[0])*t,p1[1]+(p2[1]-p1[1])*t]

        return [i1[0]+(i2[0]-i1[0])*t,i1[1]+(i2[1]-i1[1])*t]

    def bezier_high(self,p0,p1,p2,p3,t):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve

        i1 = self.bezier(p0,p1,p2,t)
        i2 = self.bezier(p1,p2,p3,t)

        return [i1[0]+(i2[0]-i1[0])*t,i1[1]+(i2[1]-i1[1])*t]

    def bezier_high2(self,p0,p1,p2,p3,p4,t):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve
        #https://www.jasondavies.com/animated-bezier/
        i1 = self.bezier_high(p0,p1,p2,p3,t)
        i2 = self.bezier_high(p1,p2,p3,p4,t)

        return [i1[0]+(i2[0]-i1[0])*t,i1[1]+(i2[1]-i1[1])*t]

    def lengthbezier(self,p0,p1,p2,step,p3=False,p4=False):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve

        pos = 0
        length = 0
        p = p0
        while pos <= 1:


            if p3==False:
                xy = self.bezier(p0,p1,p2,pos)
            elif p4==False:
                xy = self.bezier_high(p0,p1,p2,p3,pos)
            elif p4!=False:
                xy = self.bezier_high2(p0,p1,p2,p3,p4,pos)

            length += math.sqrt( (xy[0]-p[0])**2 + (xy[1]-p[1])**2 )
            p = xy
            pos += step

        return round(length)

    def wherebezier(self,p0,p1,p2,step,stop,p3=False,p4=False):
        #https://en.wikipedia.org/wiki/B%C3%A9zier_curve

        pos = 0
        length = 0
        p = p0
        xy = [0,0]

        if stop<0:
            if p3==False:
                stop = self.lengthbezier(p0,p1,p2,step)+stop
            elif p4==False:
                stop = self.lengthbezier(p0,p1,p2,step,p3)+stop
            else:
                stop = self.lengthbezier(p0,p1,p2,step,p3)+stop

        while pos <= 1:

            if length>stop: #stop if it reached the length along the line
                break

            if p3==False:
                xy = self.bezier(p0,p1,p2,pos)
            elif p4==False:
                xy = self.bezier_high(p0,p1,p2,p3,pos)
            else:
                xy = self.bezier_high2(p0,p1,p2,p3,p4,pos)

            length += math.sqrt( (xy[0]-p[0])**2 + (xy[1]-p[1])**2 )
            p = xy
            pos += step

        return pos,xy


    #find slope and y-intercept of a line through two points
    def LineEquation(self,p1, p2):
        # slope
        m = (p2['y']-p1['y'])/(p2['x']-p1['x'])

        # y-intercept
        b = p1['y'] - m*p1['x']

        # direction
        if p1['y'] >= p2['y'] and p1["x"] > p2["x"]:
            x = -1
            y = -1
        elif p1['y'] > p2['y'] and p1["x"] < p2["x"]:
            x = 1
            y = 1
        elif p1['y'] < p2['y'] and p1["x"] > p2["x"]:
            x = -1
            y = -1
        elif p1['y'] <= p2['y'] and p1["x"] < p2["x"]:
            x = 1
            y = 1

        return {'m':m, 'b':b, 'x':x, 'y':y}

    def MoveAlongLine(self,pixels_to_move,m,perpendicular, xDir=1, yDir=1):
        # perpendicular line
        if perpendicular == True:
            if m != 0:
                m = -1/m
            else:
                return {"x":0, "y":pixels_to_move}

        # a^2+b^2=c^2 where a=x and b=mx
        # x^2+mx^2=c^2 (the b of y=mx+b can be omitted as the y-intercept is zero)
        c = pixels_to_move
        x = sqrt( pow(c,2)/(1+pow(m,2)) )
        y = m*x

        return {"x":x*xDir, "y":y*yDir}

    #Draws the full backbone representation for one (box) helix

    def DrawBackbone(self,coordinates):

        # box properties
        numSides = 4

        # start at residue position 1
        cur_res = 1

        # loop through residues and draw backbone
        output = ""
        #for (i=1;i<=count(coordinates);i++) {
        for i in range(1,len(coordinates)+1):
            cur_res = i

            # find next residue
            next_res = cur_res-1
            if next_res<1:
                next_res = 19

            # find prev residue
            prev_res = cur_res+1;
            if prev_res>len(coordinates):
                prev_res = 2

            # next-in-wheel residue orientation
            wheel_next_res = cur_res-numSides
            if wheel_next_res <= (1-numSides):  # if current res is the last res, next is the first res
                wheel_next_res = len(coordinates)
            elif wheel_next_res < 1:
                wheel_next_res += len(coordinates)-1

            # line thickness
            thickness = 6*((i/len(coordinates)/1.2))
            thicknessNext = 6*(((i+1)/len(coordinates)/1.2))
            thicknessPrev = 6*(((i-1)/len(coordinates)/1.2))

            # find residue base points
            cur_points = self.ResiduePoints(cur_res, thickness, coordinates)
            next_points = self.ResiduePoints(next_res, thicknessNext, coordinates)
            prev_points = self.ResiduePoints(prev_res, thicknessPrev, coordinates)

            # lines
            cur_in = self.LineEquation(cur_points[1], cur_points[2])
            cur_out = self.LineEquation(cur_points[4], cur_points[3])
            next_in = self.LineEquation(next_points[1], next_points[2])
            next_out = self.LineEquation(next_points[4], next_points[3])
            prev_in = self.LineEquation(prev_points[1], prev_points[2])

            # line points
            p1 = cur_points[1];
            p2 = self.LineIntercept(cur_in['m'], cur_in['b'], next_out['m'], next_out['b'])
            p3 = self.LineIntercept(cur_in['m'], cur_in['b'], next_in['m'], next_in['b'])
            p4 = self.LineIntercept(cur_in['m'], cur_in['b'], prev_in['m'], prev_in['b'])
            p5 = self.LineIntercept(cur_out['m'], cur_out['b'], prev_in['m'], prev_in['b'])
            p6 = cur_points[4]

            # shorter line for the last residue
            if i == len(coordinates):
                move_in = self.MoveAlongLine(40,cur_in['m'],False,cur_in['x'],cur_in['y']);
                p4 = {'x':p1['x']+move_in['x'], 'y':p1['y']+move_in['y']}
                p5 = {'x':p6['x']+move_in['x'], 'y':p6['y']+move_in['y']}

            # draw line
            points = [p2,p1,p6,p5,p4,p3]
            points_txt = ""
            for coord in points:
                points_txt += str(coord['x'])+","+str(coord['y'])+" "

            # gradient
            gradientId = uniqid()
            angle = 1/tan(self.deg2rad(cur_in['m']))
            output +=  """
                <defs>
                    <linearGradient id='{}' x1='100%' y1='0%' x2='0%' y2='0%' gradientTransform='rotate({})'>
                    <stop offset='0%' stop-color='#00cc00' stop-opacity='1'/>
                    <stop offset='100%' stop-color='#006600' stop-opacity='1'/>
                    </linearGradient>
                </defs>
                """.format(gradientId,angle)

            lineFill = "white"

            # add SVG to output
            output += "<polyline stroke='black' stroke-width='0.5' points='"+points_txt+"' fill='"+lineFill+"'/>"

        return output;

    def LineIntercept(self,m1, b1, m2, b2):
        # line intercept
        intercept = {}
        intercept['x'] = (b2-b1)/(m1-m2)
        intercept['y'] = m1*intercept['x']+b1

        return intercept

    #Finds points on both sides of the residue center and the equation of a line to a reference residue
    def ResiduePoints(self, cur_res, thickness, coordinates):
        # FIXME write a general formula
        ref_residues = {
            20 : 2,
            16 : 2,
            12 : 2, # exception
            8 : 17,
            4 : 13,
            19 : 1,
            15 : 1,
            11 : 20,
            7 : 20,
            3 : 16,
            18 : 4,
            14 : 4,
            10 : 19,
            6 : 19,
            2 : 15,
            17 : 3,
            13 : 3,
            9 : 18,
            5 : 18,
            1 : 14}

        ori = {}

        ref_res = ref_residues[cur_res]

        # next-in-wheel residue orientation
        numSides = 4
        wheel_next_res = cur_res-numSides
        if wheel_next_res <= (1-numSides): # if current res is the last res, next is the first res
            wheel_next_res = len(coordinates)
        elif wheel_next_res < 1:
            wheel_next_res += len(coordinates)-1

        # point orientation
        if (coordinates[cur_res]['y'] >= coordinates[wheel_next_res]['y'] and coordinates[cur_res]['x'] >
            coordinates[wheel_next_res]['x']):
            # SW
            ori['x'] = -1
            ori['y'] = -1

        elif (coordinates[cur_res]['y'] > coordinates[wheel_next_res]['y'] and coordinates[cur_res]['x'] <
            coordinates[wheel_next_res]['x']):
            # NW
            ori['x'] = 1
            ori['y'] = 1

        elif (coordinates[cur_res]['y'] < coordinates[wheel_next_res]['y'] and coordinates[cur_res]['x'] >
            coordinates[wheel_next_res]['x']):
            # SE
            ori['x'] = -1
            ori['y'] = -1

        elif (coordinates[cur_res]['y'] <= coordinates[wheel_next_res]['y'] and coordinates[cur_res]['x'] <
            coordinates[wheel_next_res]['x']):
            # NE
            ori['x'] = 1
            ori['y'] = 1


        # slope of line from current residue to reference residue
        m = (coordinates[ref_res]['y']-coordinates[cur_res]['y'])/(coordinates[ref_res]['x']-
            coordinates[cur_res]['x'])

        # calculate coordinates of perpendicular line
        per_move = self.MoveAlongLine(thickness/2,m,True); # move thickness/2 pixels along a perpendicular line

        # define points
        points = [0];
        points.append({'x':coordinates[cur_res]['x']+per_move['x']*ori['x'], 'y':coordinates[cur_res]['y']+
            per_move['y']*ori['y']})
        points.append({'x':coordinates[ref_res]['x']+per_move['x']*ori['x'], 'y':coordinates[ref_res]['y']+
            per_move['y']*ori['y']})
        points.append({'x':points[2]['x']+per_move['x']*ori['x']*-2, 'y':points[2]['y']+per_move['y']*
            ori['y']*-2})
        points.append({'x':points[1]['x']+per_move['x']*ori['x']*-2, 'y':points[1]['y']+per_move['y']*
            ori['y']*-2})

        return points
