#!/usr/bin/env python
import urllib
from Bio import Phylo
import sys,os


class PrepareTree:
    
    def __init__(self):
        self.proteins = {}
        self.colors = {}
        self.ligand_names={}
        self.crystalized_ligands = []
        self.family_names = {}
        self.sub_names={}
        self.styles = ''
        self.defaultColours={}

##### Convert HSV color to HEX color #####
    def HSV_2_RGB(self,HSV):
        ''' Converts an integer HSV tuple (value range from 0 to 255) to an RGB tuple '''
        # Unpack the HSV tuple for readability
        H, S, V = HSV
        # Check if the color is Grayscale
        if S == 0:
            R = V
            G = V
            B = V
        else:
            # Make hue 0-5
            region = H // 43
            # Find remainder part, make it from 0-255
            remainder = (H - (region * 43)) * 6
            # Calculate temp vars, doing integer multiplication
            P = (V * (255 - S)) >> 8
            Q = (V * (255 - ((S * remainder) >> 8))) >> 8
            T = (V * (255 - ((S * (255 - remainder)) >> 8))) >> 8
            # Assign temp vars based on color cone region
            if region == 0:
                R = V
                G = T
                B = P
            elif region == 1:
                R = Q 
                G = V 
                B = P
            elif region == 2:
                R = P 
                G = V 
                B = T
            elif region == 3:
                R = P 
                G = Q 
                B = V
            elif region == 4:
                R = T 
                G = P 
                B = V
            else: 
                R = V 
                G = P 
                B = Q

        R =hex(R).split('x')[-1]
        G =hex(G).split('x')[-1]
        B =hex(B).split('x')[-1]
        if len(R) == 1:
            R='0'+R
        if len(G) == 1:
            G='0'+G
        if len(B) == 1:
            B='0'+B

        return (R+G+B)

    def trans_0_2_A(self,num):
### Translate protein family from 00N_00N_00N type to A_B_C type ###
        string = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        out = ''
        for n in num.split('_'):
            b = int(n)-1
            decimal = 0
            while b > len(string)-1:
                decimal +=1
                b = b - len(string)
            if decimal != 0:
                if decimal>len(string):
                    out += 'A'+string[(decimal-len(string))]+string[b]+'_'
                else:
                    out += string[decimal]+string[b]+'_'
            else:
                out += string[b]+'_'
        return out.rstrip('_')

    def parseFamily(self,inp):
        out = {}
        prots = inp.split(']]')
        for n in prots[:-1]:
            name,acc,fam,desc,spec = n.split('[')
            fam = self.trans_0_2_A(fam)
            out[acc] = {'name':name,'family':fam,'description':desc,'species':spec,'class':'','ligand':'','type':''}
        return out

#### Get colors based on family trees ####
    def get_colours(self,inp):
        cl = {'A':'Class A','B':'Class B','C':'Class C','D':'Class F','E':'cAMP Receptors','K':'Vomeronasal','L':'Taste'}
        self.styles += '<cryst fill=\'#'+ self.HSV_2_RGB((136,130,225))+'\' stroke=\'#'+self.HSV_2_RGB((136,130,225))+'\'/><none fill=\'#EEE\' stroke=\'#DDD\'/>'
        self.classes = []
        self.ligands = []
        self.families = []
    #### Define number of color groups ####
        for acc in inp.keys():
            self.styles += '<%s fill=\'#EEE\' stroke=\'#DDD\' id=\'%s\' class=\'bgfield\'/>' %(inp[acc]['family'],inp[acc]['family'])
            fam = inp[acc]['family'].split('_')
            if len(fam) == 2 or len(fam)==1:
                inp[acc]['class']=fam[0]
                if not fam[0] in self.classes:
                    self.classes.append(fam[0])
                inp[acc]['ligand']='none'
                inp[acc]['type']='none'
            if len(fam)==3:
                inp[acc]['class']=fam[0]
                if not fam[0] in self.classes:
                    self.classes.append(fam[0])
                inp[acc]['ligand']='_'.join(fam[:2])
                if not '_'.join(fam[:2]) in self.ligands:
                    self.ligands.append('_'.join(fam[:2]))
                inp[acc]['type']='none'
            if len(fam)>=4:
                inp[acc]['class']=fam[0]
                if not fam[0] in self.classes:
                    self.classes.append(fam[0])
                inp[acc]['ligand']='_'.join(fam[:2])
                if not '_'.join(fam[:2]) in self.ligands:
                    self.ligands.append('_'.join(fam[:2]))
                inp[acc]['type']='_'.join(fam[:3])
                if not '_'.join(fam[:3]) in self.families:
                    self.families.append('_'.join(fam[:3]))

        self.data = inp
                
    #### Get colours ####
        nclasses = len(self.classes)
        nligands = len(self.ligands)
        ncolours = len(self.families)
        self.ligcolors = {}
        self.famcolors = {}
        hue_step = 256/nligands
        hue_count = 0
        start_hue = 360-(hue_step/2)
        current_hue = 0
        hue_range = [start_hue,start_hue+hue_step]
        for a in sorted(self.ligands):
            colour = self.HSV_2_RGB((int(current_hue),127,128))
            self.styles+=('<%s fill=\'#%s\' stroke=\'#%s\' class=\'chart\' id=\'%s\'/>' %(a,colour,'000000',a))
            self.defaultColours[a]=['#'+colour,'#000000']
            self.ligcolors[a]=[colour,hue_range,0]
            hue_count +=1
            current_hue = hue_step*hue_count
            hue_range[0] += hue_step*hue_count
            hue_range[1] += hue_step*hue_count
        hue_count = 0
        for a in sorted(self.families):
            n=ncolours
            sat_step = 200/n
            counter =0
            lig = '_'.join(a.split('_')[:2])###FIX
            try:
                hue_range = self.ligcolors[lig][1]
                hue_step = (hue_range[1]-hue_range[0])/(ncolours+1)
                self.ligcolors[lig][2]+=hue_step/2
                if hue_step == 0:
                    hue_step = 1
                colour = self.HSV_2_RGB((int(hue_count*hue_step+self.ligcolors[lig][2]),127,180))
                self.styles+=('<%s fill=\'#%s\' stroke=\'#%s\' class=\'chart\' id=\'%s\'/>' %(a,colour,'000000',a))
                self.defaultColours[a]=['#'+colour,'#000000']
                self.famcolors[a]=colour
                counter +=1
                self.ligcolors[lig][2]+=hue_step
            except KeyError:
                continue
            hue_count +=1
               
        v_step = 180/nclasses
        v_count = 0
        for a in self.classes:
            colour = self.HSV_2_RGB((0,0,10+(v_count*int(v_step))))
            self.styles+=('<%s fill=\'#%s\' stroke=\'#%s\' label=\'%s\' labelStyle=\'sectorHighlightText\' class=\'chart\' id=\'%s\'/>' %(a,colour,colour,cl[a].upper(),a))
            self.defaultColours[a]=['#'+colour,'#'+colour]
            self.styles+="<sectorHighlightText font-family='Verdana' font-size='10' font-weight='bold' fill='#FFFFFF' rotate='90'/>"
            v_count +=1
    
    def drawColorPanel(self):

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
        fillcolors = [['#EEEEEE', '#000000']]
        for key,value in presetColors.items():
            if value not in fillcolors:
                fillcolors.append(value)

        colors = ""
        for color in fillcolors:
            colors += "<div class='pick-color selected' id='pick-"+color[0]+"-"+color[1]+"' style='background-color: "+color[0]+";'>&nbsp;</div>"
            
        output = ("<br>Pick color:" +
            colors )



        return boxstyle+ output

    def treeDo(self,infile,branches,family,famdict=None):
        self.famdict=famdict
        d = '/'.join(infile.split('/')[:-1])
        z = open(infile,'r').read()
        w = open(d+'/rong','w').write(z.replace('-',''))
        raw = open(d+'/raw.xml','w')
        Phylo.convert(d+'/rong','newick',raw,'phyloxml')
        raw.close()
        xml = open(d+'/raw.xml','r').readlines()
        out = open(d+'/out.xml','w')
        self.get_colours(family)
        print(branches)
        charts = ''
        flag = False
        flag2 = ''
        stylesflag=False
        charts += "<family type='binary' thickness='10' bufferInner='0'/>"
        charts += "<ligand type='binary' thickness='10' bufferInner='1'/>"
        charts += "<class type='binary' thickness='10' bufferInner='1'/>"
        charts = '<charts>'+charts+'</charts>'
        for line in xml:
            #print('dupajasia')
            if stylesflag == True:
                out.write("<render>"+charts+"<styles>"+self.styles+"</styles></render>")
                stylesflag = False
            ################# Remove header trash #######################
            if 'phyloxml' in line:
                line = line.split('phyloxml')[0]+'phyloxml>'
            line = line.replace('\"','\'').replace('phy:','')
            ################# Remove forced rooting #####################
            if flag == True:
                if '>1.0<' in line:
                    line=line.replace('>1.0<','>0.0<')
                    flag = False
            if "rooted='false'" in line:
                flag = True
                stylesflag=True
            ################# Force even branch lengths #################
            if branches == True:
                if '<branch_length>' in line:
                    if '<branch_length>0.0</branch_length>' not in line:
                        number = line.split('>')[1].split('<')[0]
                        line = line.replace(number,'0.1')
            ################# Reformat names ############################
            if '<name>' in line:
                name = line.split('<')[1].split('>')[1]
                #recname,fam,desc,spec = family[name]
                chart = '<chart>'
                chart += '<family>%s</family>' %(self.data[name]['type'])
                chart += '<ligand>'+self.data[name]['ligand']+'</ligand>'
#                if len(self.classes)>1:
                chart += '<class>'+self.data[name]['class']+'</class>'
                chart += '</chart>'
                flag2 = [name,chart]
                line = line.replace(name,self.data[name]['name']).replace('<name', "<name bgStyle='%s'" %self.data[name]['family'])
#                if recname in crystalized_ligands:
#                    line = line.replace('<name>','<name bgStyle=\'cryst\'>')
#                else:
#                    line = line.replace('<name>','<name bgStyle=\'none\'>')
            ############## Add annotations and descriptions #############
            if '<branch_length>' in line:
                if flag2 != '':
                    line = line.strip('\n')+' <annotation><desc>'+self.data[flag2[0]]['description']+' ('+self.data[flag2[0]]['species']+')'+'</desc><uri>http://tools.gpcr.org/visualise/protein/'+self.data[flag2[0]]['link']+'</uri> </annotation>'+flag2[1]#flag3[3]
                    flag2=''
            out.write(line.strip('\n'))
        ###SVG legend###
        column = 200
        verse = 20
        self.legend =''
        legend2 = ''
        width1 = 300
        width2 = 300
        if 10+verse*len(self.ligands) > verse*len(self.families)+10:
            height = 10+verse*len(self.ligands)
        else:
            height = verse*len(self.families)+10
        self.legend ='<?xml version="1.0" standalone="no"?><svg class="legend" id="legendSVG" width="'+str(width1+width2+100)+'px" height="'+str(height+70)+'px" version="1.1" xmlns="http://www.w3.org/2000/svg"><g id="leg_group" class="legend"><text x="0" y="30" font-family="Verdana" font-size="20" font-weight="bold">Ligand type (Outer ring)</text><rect x="0" y="70" width="'+str(width1)+'" height="'+str(verse*len(self.ligands)+10)+'" stroke="black" fill="transparent" stroke-width="5"/>'
        count = 0
        for n in sorted(self.ligands):
            try:
                name = self.famdict[n]
            except:
                name = n
            self.legend += '<rect x="10" y="'+str(80+count*20)+'" height="10" width="30" class="chart" id="'+n+'" style="stroke:#000000; fill: #'+self.ligcolors[n][0]+'"/><text x="45" y="'+str((count*20)+90)+'" style="font-family:Verdana; font-size:12;">'+name+'</text>'
            count +=1
        self.legend+='<rect x="'+str(width1+40)+'" y="70" width="'+str(width2+45)+'" height="'+str(verse*len(self.families)+10)+'" stroke="black" fill="transparent" stroke-width="5"/><text x="'+str(width1+40)+'" y="30" font-family="Verdana" font-size="20" font-weight="bold">Receptor type (Inner ring)</text>'
        count2 = 0
        for n in sorted(self.families):
            try:
                name2 = self.famdict[n]
            except:
                name2 = n
            self.legend += '<rect x="'+str(width1+50)+'" y="'+str(80+count2*20)+'" height="10" width="30" class="chart" id="'+n+'" style="stroke:#000000; fill: #'+self.famcolors[n]+'"/><text x="'+str(width1+85)+'" y="'+str((count2*20)+90)+'" style="font-family:Verdana; font-size:12;">'+name2+'</text>'
            count2 +=1
        self.legend += '</g></svg>'
        self.box = self.drawColorPanel()

if __name__ == '__main__':
    tree = PrepareTree()
    tree.treeDo(sys.argv[1],sys.argv[2])

