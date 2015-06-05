from common.diagrams import Diagram

from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

from django.utils.safestring import mark_safe

from math import cos, sin, pi

class DrawHelixBox(Diagram):

    plot_data = {}
    plot_data['direction'] = [0,0, 1, 0, 1, 0, 1, 0]; # 0: EC->IC, 1: IC->EC
    plot_data['helixRadius'] = 70;
    # Class A
    plot_data['Class A'] = {}
    plot_data['Class A']['coordinates'] = [0,[455,240],[385,108],[245,118],[105,95],[75,235],[185,320],[324,303]]
    plot_data['Class A']['helixTopResidues'] = [0,'1x32', '2x63', '3x26', '4x61', '5x39', '6x57', '7x32']
    plot_data['Class A']['helixTopResidues'] = [0,32, 63, 26, 61, 39, 57, 32]
    plot_data['Class A']['rotation'] = [0,340, 320, 290, 125, 40, 180, 130] # in degrees

    def __init__(self, protein):
        self.receptorId = protein
        self.family = protein.get_protein_class()
        self.output = ''

        for i in range(1,len(self.plot_data[self.family]['coordinates'])):
            self.residuelist = Residue.objects.filter(
                protein_segment__slug='TM'+str(i), protein_conformation__protein__entry_name=self.receptorId)

            self.output += self.DrawHelix(self.plot_data[self.family]['coordinates'][i][0],
                self.plot_data[self.family]['coordinates'][i][1],
                self.residuelist,
                self.plot_data['helixRadius'],
                self.plot_data['direction'][i],
                i,
                self.plot_data[self.family]['helixTopResidues'][i],
                self.plot_data[self.family]['rotation'][i]
                )


    def __str__(self):  
        return mark_safe(self.create(self.output))

    def DrawHelix(self, startX,startY,residuelist,radius,direction,helixNum,helixTopResidue,rotation):

        #fetch residues | add sort by sequence_number FIXME
        
        #.prefetch_related(
         #       'protein_conformation__protein', 'protein_conformation__state', 'protein_segment',
          #      'generic_number__scheme', 'display_generic_number__scheme', 'alternative_generic_numbers__scheme')

        sequence = {}
        for r in residuelist:
            #print(r)
            #print(r.amino_acid)
            #print(r.sequence_number)
            sequence[int(r.generic_number.label[2:])] = {'residueType':r.amino_acid,'residueNumber':r.sequence_number}
        print("HelixNum",helixNum,"residueNumbers",sequence.keys())
        # box size
        numResPerSide = 5
        numResInBox = numResPerSide*4
        residueRadius = 12
        
        # start positions of the four sides
        output_residue = ""

        helixSideCoord = {}

        for i in range(1,5):
            helixSideCoord[i] = {}
            helixSideCoord[i]['x'] = startX + radius * cos( self.deg2rad((i-1)*90+rotation) )
            helixSideCoord[i]['y'] = startY + radius * sin( self.deg2rad((i-1)*90+rotation) )

        # starting side
        helixSide = 4 # start on the last side of the helix

        # residue position incrementer
        resPosIncr = 1

        # helix direction indicator
        helixDirection = 1
        if helixNum % 2 != 0:
            helixDirection = -1

        # init residue number (draw lowest residue first)
        currentResidue = helixTopResidue - helixDirection * (numResInBox - 1)

        coordinates = {}

        for i in range(1,numResInBox+1):
            
            # next helix side
            if helixSide == 1:
                # go back to side 1
                nextHelixSide = 4
            else:
                # next side
                nextHelixSide = helixSide - 1

            # find number for next residue (for bulge assignment)
            nextResidue = currentResidue + helixDirection

            # line from the current side to the next
            #print(helixSide)
            print("i:"+str(i),"helixSide",helixSide,"currentResidue",currentResidue,"nextResidue",nextResidue)
            lineEquation = self.LineEquation(
                {
                    "x":helixSideCoord[helixSide]["x"],
                    "y":helixSideCoord[helixSide]["y"]},
                {
                    "x":helixSideCoord[nextHelixSide]["x"],
                    "y":helixSideCoord[nextHelixSide]["y"]})

            # direction of perpendicular move (for bulges), depends on line slope
            perpDir = -1
            if lineEquation["m"] < 0:
                perpDir = 1

            # length of perpendicular move (for bulges), depends on line slope
            perpLen = 5
            if lineEquation["m"] > 1 and lineEquation["x"] > 0:
                perpLen = 15

            x = {}
            y = {}
            coordinates[i] = {}
            # calculate coordinates of the residue (and save them for backbone drawing)
            move = self.MoveAlongLine(residueRadius*1.4, lineEquation['m'], False, lineEquation['x'],
                lineEquation['y'])
            x['0'] = coordinates[i]["x"] = helixSideCoord[helixSide]["x"] + resPosIncr*move["x"]
            y['0'] = coordinates[i]["y"] = helixSideCoord[helixSide]["y"] + resPosIncr*move["y"]

            # init array for residues to be drawn in this position (can be more than on in a bulge)
            coordinateIndices = []
            if currentResidue in sequence.keys():
                coordinateIndices = ['0']

            # check for bulges (residue number suffixed with 1) or constrictions (missing
            # residue number)
            if int(str(nextResidue)+"1") in sequence:
                print("FOUND BULGE")
                coordinateIndices.append('1')
                perpMove = self.MoveAlongLine(perpLen, lineEquation["m"], True, lineEquation['x'],
                    lineEquation['y']);
                x['1'] = (helixSideCoord[helixSide]["x"] + (resPosIncr+0.6)*move["x"] +
                     resPosIncr*perpDir*perpMove["x"])
                y['1'] = (helixSideCoord[helixSide]["y"] + (resPosIncr+0.6)*move["y"] +
                     resPosIncr*perpDir*perpMove["y"])
            
            if coordinateIndices: #not empty
                for coordinateIndex in coordinateIndices:
                    tempCurrentResidue = currentResidue;
                    if coordinateIndex=='1': #not empty
                        tempCurrentResidue = int(str(nextResidue)+str(coordinateIndex))

                    # Get label information of each residue.
                    residue_number = sequence[tempCurrentResidue]['residueNumber'];
                    label = ''
                    #label = self::getResidueLabel(receptorId, helixNum, residue_number); #FIXME IMPLEMENT
                    #label = "m: " . lineEquation["m"] . " x: " . lineEquation["x"];
                    
                    # draw the residue
                    output_residue += self.DrawResidue(x[coordinateIndex],y[coordinateIndex],
                        sequence[tempCurrentResidue]['residueType'],residue_number, label,residueRadius);
                
            

            # residue position incrementer
            if i % 4 == 0:
                resPosIncr += 1

            # increment side
            helixSide = nextHelixSide

            # next residue
            currentResidue = nextResidue
        




        helix_number_svg = "<text x='"+str(startX)+"' y='"+str(startY+7)+"' text-anchor='middle' font-family='helvetica' font-size='20'>"+str(helixNum)+"</text>\n"

        return output_residue+helix_number_svg
