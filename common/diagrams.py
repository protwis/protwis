

from math import cos, sin, pi, sqrt, pow

class Diagram:
    def create(self, content):
        return "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"595\" height=\"430\">\n"+content+"</svg>"

    #Draws a ring of a helical wheel  
    def DrawResidue(self, x,y,aa,residue_number,label,radius,cfill="white", precolor = False):
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
            <circle cx='{}' cy='{}' r='{}' stroke='black' stroke-width='2' fill='{}' 
            fill-opacity='1' id='id' class='rcircle' onclick='residueColor.setColor(evt);'
            onmouseover='showToolTip(x,y,"label",id);' onmouseout='hideToolTip();'/>
            <text x='{}' y='{}' text-anchor='middle' font-family='helvetica' font-size='16' fill='tfill'
            id='idtext' onclick='residueColor.setColor(evt);' class='rtext'
            onmouseover='showToolTip(x,y,"label",id);' onmouseout='hideToolTip();'>{}</text>
            """.format(x,y,radius,cfill,x,y+6,aa) #aa
        return output

    def deg2rad(self,degrees):
        radians = pi * degrees / 180
        return radians

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

# class drawHelix:
#     return 1

# class DrawBackbone(coordinates):
#     return 1

# class DrawBackground(helix_coordinates,$fill_color):
#     return 1

