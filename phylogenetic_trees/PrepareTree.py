#!/usr/bin/env python
import urllib
from Bio import Phylo
import sys,os
from collections import defaultdict


class PrepareTree:
    def __init__(self,build):
        self.rings={'crystal': {'include':'False', 'order':6, 'colours':{'crystal_true':'#6dcde1','crystal_false':'#EEE'}, 'color_type':'single', 'items':[], 'parent':None, 'child': None, 'name':'Crystallized structures'},
            'class': {'include':'False', 'order':0, 'colours':{}, 'items':[], 'color_type':'grayscale', 'parent':[], 'child': ['family,ligand'], 'name':'Class'},
            'family': {'include':'False', 'order':1, 'colours':{}, 'items':[], 'color_type':'spectrum', 'parent':[], 'child': ['ligand'], 'name':'Ligand type'},
            'ligand': {'include':'False', 'order':2, 'colours':{}, 'items':[], 'color_type':'spectrum', 'parent':['family','class'], 'child': [], 'name':'Receptor family'},
            'mutant': {'include':'False', 'order':3, 'colours':{'mutant_true':'#6dcde1','mutant_false':'#EEE'}, 'color_type':'single', 'items':[], 'parent':[], 'child': ['mutant_plus','mutant_minus'], 'name':'Mutated proteins'},
            'mutant_plus': {'include':'False', 'order':4, 'colours':{'mutant_plus_true':'#6dcde1','mutant_plus_false':'#EEE'}, 'color_type':'single', 'items':[], 'parent':'mutant', 'child': [], 'name':'Positive affinity mutants'},
            'mutant_minus': {'include':'False', 'order':5, 'colours':{'mutant_minus_true':'#6dcde1','mutant_minus_false':'#EEE'}, 'color_type':'single', 'items':[], 'parent':'mutant', 'child': [], 'name':'Negative affinity mutants'}
            }
        self.prots = {}
        self.build=build

    def get_grayscale_colours(self, itemlist):
        colours_dict = {}
        if itemlist:
            v_step = 180/len(itemlist)
        else:
            v_step = 0
        v_count = 0
        for a in sorted(itemlist):
            colour = self.HSV_2_RGB((0,0,10+(v_count*int(v_step))))
            v_count +=1
            colours_dict[a]='#'+colour
        return colours_dict


    def get_spectrum_colours(self, itemlist,range=(0,255)):
        colours_dict = {}
        if itemlist:
            v_step = (range[1]-range[0])/len(itemlist)
        else:
            v_step = 0
        v_count = 0
        for a in sorted(itemlist):
            colour = self.HSV_2_RGB((int(range[0]+v_step*v_count),127,200))
            v_count +=1
            colours_dict[a]='#'+colour
        return colours_dict

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
                                    '-': ['#000000', '#000000']
                                    }
        fillcolours = [['#EEEEEE', '#000000']]
        for key,value in presetColors.items():
            if value not in fillcolours:
                fillcolours.append(value)

        colours = ""
        for color in fillcolours:
            colours += "<div class='pick-color selected' id='pick-"+color[0]+"-"+color[1]+"' style='background-color: "+color[0]+";'>&nbsp;</div>"

        output = ("<br><b>Pick color:</b>" +
            colours )

        return boxstyle+ output


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

    def trim_colour(self,m):
        if m > 256:
            return m-256
        else:
            return m

    def get_family_meta(self,inp):
        for entry_name in inp.keys():
            fam = inp[entry_name]['family'].split('_')
            prot_class=fam[0]
            family='_'.join(fam[:2])
            ligand='_'.join(fam[:3])
            prot='_'.join(fam)
            if not prot_class in self.rings['class']['items']:
                self.rings['class']['items'].append(prot_class)
            if not family in self.rings['family']['items']:
                self.rings['family']['items'].append(family)
            if not ligand in self.rings['ligand']['items']:
                self.rings['ligand']['items'].append(ligand)
            self.prots[entry_name]={'class':prot_class,'family':family,'ligand':ligand, 'acc':prot, 'name':inp[entry_name]['name'], 'link':inp[entry_name]['link'], 'species':inp[entry_name]['species'], 'desc':inp[entry_name]['description']}
            self.prots[prot]=entry_name

    def get_tree_data(self, Additional_info):
           for datatype in Additional_info:
               self.rings[datatype]['include']=Additional_info[datatype]['include']
               self.rings[datatype]['items']=Additional_info[datatype]['proteins']
               self.rings[datatype]['colours']=Additional_info[datatype]['colours']


    def get_charts(self):
        charts = '<charts>'
        for ring in sorted(self.rings.items(), key=lambda x: x[1]['order'], reverse=True):
            if ring[1]['include'] == 'True':
                charts += "<%s type='binary' thickness='10' bufferInner='0'/>" %(ring[0])
        charts = charts+'</charts>'
        return charts

    def get_colours(self):
        self.rings['class']['colours']= self.get_grayscale_colours(self.rings['class']['items'])
        self.rings['family']['colours']= self.get_spectrum_colours(self.rings['family']['items'],(0,255))
        self.rings['ligand']['colours']= self.get_spectrum_colours(self.rings['ligand']['items'],(0,255))

    def build_legend(self):
        verse = 20
        self.legend =''
        width = 400
        length = 0
        total_rings = 0
        for ring in self.rings:
            if self.rings[ring]['include']=='True':
                total_rings+=1
                if len(self.rings[ring]['items'])>length:
                   length = len(self.rings[ring]['items'])
        height = verse*length+10
        ring_no = 0
        verse_count = 0
        self.legend ='<?xml version="1.0" standalone="no"?><svg class="legend" id="legendSVG" width="'+str((width+40)*total_rings+100)+'px" height="'+str(height+70)+'px" version="1.1" xmlns="http://www.w3.org/2000/svg"><g id="leg_group" class="legend">'
        for ring in sorted(self.rings.items(), key=lambda x: x[1]['order']):
            if ring[1]['include']=='True':
                self.legend+='<g class="chart'+str(ring[1]['order'])+'"><text x="'+str(ring_no*(width+40))+'" y="30" font-family="Verdana" font-size="20" font-weight="bold">'+self.rings[ring[0]]['name']+'</text><rect x="'+str(ring_no*(width+40))+'" y="70" width="'+str(width)+'" height="'+str(verse*len(self.rings[ring[0]]['items'])+10)+'" stroke="black" fill="transparent" stroke-width="5"/>'
                verse_count = 0
                for item in sorted(ring[1]['items']):
                    if ring[1]['color_type']=='single':
                        self.legend += '<rect x="'+str(ring_no*(width+40)+10)+'" y="'+str(80+verse_count*20)+'" height="10" width="30" class="chart" id="'+item+'" stroke="#000000" fill="'+ring[1]['colours'][ring[0]+'_true']+'"/><text x="'+str(ring_no*(width+40)+45)+'" y="'+str(90+verse_count*20)+'" style="font-family:Verdana; font-size:12;">'+self.famdict['_'.join(item.split('_')[:4])].replace('<sub>','').replace('</sub>','').replace('<i>','').replace('</i>','')+'</text>'
                    else:
                        self.legend += '<rect x="'+str(ring_no*(width+40)+10)+'" y="'+str(80+verse_count*20)+'" height="10" width="30" class="chart" id="'+item+'" stroke="#000000" fill="'+ring[1]['colours'][item]+'"/><text x="'+str(ring_no*(width+40)+45)+'" y="'+str(90+verse_count*20)+'" style="font-family:Verdana; font-size:12;">'+self.famdict['_'.join(item.split('_')[:4])].replace('<sub>','').replace('</sub>','').replace('<i>','').replace('</i>','')+'</text>'
                    verse_count +=1
                ring_no+=1
                self.legend+='</g>'
        self.legend+='</g></svg>'

    def get_styles(self):
        self.styles = ''
        for prot in self.prots:
            if len(prot.split('_')) == 1:
                if self.build != False:
                    if self.prots[prot]['acc'] in self.rings['crystal']['items']:
                        self.styles += '<'+self.prots[prot]['acc']+' fill=\''+self.rings['crystal']['colours']['crystal_true']+'\' stroke=\'#DDD\' id=\''+self.prots[prot]['acc']+'\' class=\'bgfield\'/>'
                    else:
                        self.styles += '<%s fill=\'#EEE\' stroke=\'#DDD\' id=\'%s\' class=\'bgfield\'/>' %(self.prots[prot]['acc'],self.prots[prot]['acc'])
                else:
                    self.styles += '<%s fill=\'#EEE\' stroke=\'#DDD\' id=\'%s\' class=\'bgfield\'/>' %(self.prots[prot]['acc'],self.prots[prot]['acc'])
        for ring in self.rings:
            if self.rings[ring]['include']== 'True':
                if ring == 'class':
                    self.styles+="<sectorHighlightText font-family='Verdana' font-size='10' font-weight='bold' fill='#FFFFFF' rotate='90'/>"
                colour_dict = self.rings[ring]['colours']
                for item in colour_dict:
                    if ring == 'class':
                        self.styles+='<%s fill=\'%s\' stroke=\'%s\' label=\'%s\' labelStyle=\'sectorHighlightText\' class=\'chart0\' id=\'%s\'/>' %(item,colour_dict[item],colour_dict[item],self.famdict[item],item)
                    else:
                        self.styles +='<%s fill=\'%s\' stroke=\'%s\' id=\'%s\' class=\'%s\'/>' %(item,colour_dict[item],colour_dict[item],item, 'chart'+str(self.rings[ring]['order']))




    def treeDo(self,d, infile,branches,family,Additional_info, famdict=None):
        self.famdict=famdict
        d = '/'.join(d.split('/'))
        z = infile
        w = open(d+'/rong','w').write(z.replace('-',''))
        raw = open(d+'/raw.xml','w')
        Phylo.convert(d+'/rong','newick',raw,'phyloxml')
        raw.close()
        xml = open(d+'/raw.xml','r').readlines()
        out = open(d+'/out.xml','w')
        self.get_tree_data(Additional_info)
        if self.build !=False:
            self.rings['class']['include']=False
            #self.rings['ligand']['include']=False
            #self.rings['family']['include']=False
        self.get_family_meta(family)
        charts=self.get_charts()
        self.get_colours()
        self.build_legend()
        self.get_styles()
        flag = False
        flag2 = ''
        stylesflag=False
        for line in xml:
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
                chart = '<chart>'
                for ring in self.rings:
                    if self.rings[ring]['include']=='True':
                        if self.rings[ring]['color_type']=='single':
                            if self.prots[name]['acc'] in self.rings[ring]['items']:
                                chart += '<%s>%s_true</%s>' %(ring,ring,ring)
                            else:
                                chart += '<%s>%s_false</%s>' %(ring,ring,ring)
                        else:
                            chart += '<'+ring+'>'+self.prots[name][ring]+'</'+ring+'>'
                chart += '</chart>'
                flag2 = [name,chart]
                line = line.replace(name,self.prots[name]['name']).replace('<name', "<name bgStyle='%s'" %self.prots[name]['acc'])
            ############## Add annotations and descriptions #############
            if '<branch_length>' in line:
                line=line.replace('>1E05<','>0.00001<').replace('-','')
                if '>0.0<' in line and flag == True:
                    line=line.replace('>0.0<','>0.00001<')
                if flag2 != '':
                    line = line.strip('\n')+' <annotation><desc>'+self.prots[flag2[0]]['desc']+' ('+self.prots[flag2[0]]['species']+')'+'</desc><uri>/protein/'+self.prots[flag2[0]]['link']+'</uri> </annotation>'+flag2[1]
                    flag2=''
            out.write(line)

        self.box = self.drawColorPanel()

if __name__ == '__main__':
    tree = PrepareTree(False)
    tree.treeDo(sys.argv[1],sys.argv[2])
