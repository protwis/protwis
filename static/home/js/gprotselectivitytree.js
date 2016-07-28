
var r = 1200 / 2;
var innerRadius = r - 170 // change inner radius of tree with this argument
var names=0; // indexing for all nodes


var cluster = d3.layout.cluster()
    .size([360, innerRadius])
    .sort(function(a, b) { return (a.value - b.value) || d3.ascending(a.length, b.length); })
    .value(function(d) { return d.length; })
    .children(function(d) { return d.branchset; })
    .separation(function(a, b) { return 1; });

function project(d) {
  var r = d.y, a = (d.x - 90) / 180 * Math.PI;
  return [r * Math.cos(a), r * Math.sin(a)];
}

function cross(a, b) { return a[0] * b[1] - a[1] * b[0]; }
function dot(a, b) { return a[0] * b[0] + a[1] * b[1]; }

function step(d) {
  var s = project(d.source),
      m = project({x: d.target.x, y: d.source.y}),
      t = project(d.target),
      r = d.source.y,
      sweep = d.target.x > d.source.x ? 1 : 0;
  return (
    "M" + s[0] + "," + s[1] +
    "A" + r + "," + r + " 0 0," + sweep + " " + m[0] + "," + m[1] +
    "L" + t[0] + "," + t[1]);
}

var rect_dim = (r-90)* 2 ;
var svg_dim = (r-90)* 2 ;


var wrap = d3.select("#selectivitydata").append("svg")
    .attr("width", svg_dim)
    .attr("height", svg_dim)
    .style("-webkit-backface-visibility", "hidden");

// Catch mouse events in Safari.
wrap.append("rect")
    .attr("width", rect_dim)
    .attr("height", rect_dim)
    .attr("fill", "none");


var vis = wrap.append("g")
    .attr("transform", "translate(" + rect_dim/2 + "," + rect_dim/2 + ")");

var start = null,
    rotate = 0,
    div = document.getElementById("selectivitydata");

//function to catch XY coordinates for mouse
function mouse(e) {

  return [
    e.pageX - div.offsetLeft - r,
    e.pageY - div.offsetTop - r
  ];
}

wrap.on("mousedown", function() {
  wrap.style("cursor", "move");
  start = mouse(d3.event);
  d3.event.preventDefault();
});
d3.select(window)
  .on("mouseup", function() {
    if (start) {
      wrap.style("cursor", "auto");
      var m = mouse(d3.event);
      var delta = Math.atan2(cross(start, m), dot(start, m)) * 180 / Math.PI;
      rotate += delta;
      if (rotate > 360) rotate %= 360;
      else if (rotate < 0) rotate = (360 + rotate) % 360;
      start = null;
      wrap.style("-webkit-transform", null);
      vis
          .attr("transform", "translate(" + rect_dim/2 + "," + rect_dim/2 + ")rotate(" + rotate + ")")
        .selectAll("text")
          .attr("text-anchor", function(d) { return (d.x + rotate) % 360 < 180 ? "start" : "end"; })
          .attr("transform", function(d) {
            return "rotate(" + (d.x - 90) + ")translate(" + (r - 170 + 2) + ")rotate(" + ((d.x + rotate) % 360 < 180 ? 0 : 180) + ")";
          });
    }
  })
  .on("mousemove", function() {
    if (start) {
      var m = mouse(d3.event);
      var delta = Math.atan2(cross(start, m), dot(start, m)) * 180 / Math.PI;
      wrap.style("-webkit-transform", "rotateZ(" + delta + "deg)");
    }
  });

//Branch length function
function phylo(n, offset) {
  if (n.length != null) offset += n.length * 115;
  n.y = offset;
  if (n.children)
    n.children.forEach(function(n) {
      phylo(n, offset);
    });
}

var x = "(((5HT1D:0.14656,5HT1B:0.18681):0.08093,(5HT1F:0.18620,5HT1E:0.20912):0.12352):0.09868,(5HT5A:0.61709,(5HT7R:0.48940,((DRD1:0.13081,DRD5:0.11617):0.37543,((((DRD2:0.15236,DRD3:0.18945):0.24327,DRD4:0.48536):0.05986,((ADA2A:0.13115,ADA2C:0.18228):0.01936,ADA2B:0.15621):0.42124):0.05763,(((ADRB2:0.25547,ADRB1:0.19778):0.08565,ADRB3:0.34080):0.26142,(((ADA1A:0.17733,ADA1B:0.14958):0.02087,ADA1D:0.20860):0.31132,((HRH2:0.52426,(5HT4R:0.52452,((TAAR1:0.26729,(TAAR2:0.22284,TAAR3:0.31017):0.23539):0.22600,(((TAAR6:0.08400,TAAR8:0.19222):0.07454,TAAR9:0.21355):0.22951,TAAR5:0.47985):0.02793):0.23241):0.04297):0.05334,((((5HT2A:0.17106,5HT2C:0.16504):0.04954,5HT2B:0.25906):0.41256,(5HT6R:0.69434,GP119:1.09071):0.07887):0.04207,(((((ACM2:0.10343,ACM4:0.08889):0.11475,((ACM3:0.11552,ACM1:0.18469):0.01030,ACM5:0.09545):0.11701):0.40955,HRH1:0.66818):0.06411,(HRH4:0.47944,HRH3:0.32974):0.35975):0.14632,((((AA1R:0.28657,AA3R:0.32509):0.15582,(AA2AR:0.22538,AA2BR:0.20286):0.13983):0.45766,((((MSHR:0.38244,((MC3R:0.16962,MC5R:0.21034):0.04580,MC4R:0.17196):0.09610):0.10745,ACTHR:0.44663):0.33195,((GPR6:0.26054,GPR12:0.27634):0.00641,GPR3:0.19274):0.56708):0.06336,(((S1PR1:0.27466,(S1PR2:0.39143,(S1PR3:0.36509,(S1PR5:0.36473,S1PR4:0.55222):0.02790):0.00867):0.16091):0.26558,(LPAR1:0.27356,(LPAR2:0.36350,LPAR3:0.29943):0.02502):0.28322):0.13154,(CNR1:0.35087,CNR2:0.35507):0.62427):0.12276):0.19783):0.04243,(((((GP161:0.88256,GP101:0.81561):0.07534,GP135:0.80381):0.04794,(GPR45:0.26599,GPR63:0.34295):0.52601):0.08401,((MTR1A:0.20256,MTR1B:0.19320):0.15610,MTR1L:0.45969):0.54665):0.04891,(((((OPN4:0.56879,OPN5:0.92903):0.06501,((((OPSG:0.02697,OPSR:0.03342):0.53998,OPSB:0.52647):0.01057,OPSD:0.44439):0.19490,OPN3:0.80151):0.19502):0.18268,((RXFP1:0.21942,RXFP2:0.26012):0.66213,(((TSHR:0.16173,LSHR:0.21867):0.04430,FSHR:0.15058):0.48573,((LGR5:0.24949,LGR6:0.45298):0.09865,LGR4:0.35743):0.99322):0.08326):0.46464):0.07709,(((FFAR4:0.92854,GPR22:1.18891):0.07305,GP176:1.12688):0.05328,((((MRGX4:0.10513,((MRGX1:0.04116,MRGX3:0.12617):0.05361,MRGX2:0.27642):0.04263):0.29790,((MRGRE:0.61765,MRGRG:0.73434):0.25048,MAS1L:0.90053):0.03282):0.06336,(MAS:0.77797,(MRGRD:0.50798,MRGRF:0.77014):0.08596):0.11897):0.69213,((((PE2R2:0.41800,PD2R:0.67895):0.10363,PI2R:0.67801):0.24656,PE2R4:0.72942):0.13337,(((PE2R1:0.58381,PF2R:0.56725):0.05786,TA2R:0.56254):0.20544,PE2R3:0.59249):0.28215):0.56311):0.09818):0.03669):0.05237,((((((GPR61:0.53396,GPR62:0.79568):0.45475,(GPR26:0.35671,GPR78:0.55697):0.67563):0.13541,UR2R:0.90308):0.05448,(((((((RL3R2:0.41574,RL3R1:0.37875):0.33106,(AGTR2:0.55722,AGTR1:0.55967):0.13173):0.05873,(APJ:0.67628,(GPR15:0.58010,GPR25:0.61852):0.18176):0.00778):0.02780,(((((CXCR2:0.05388,CXCR1:0.08153):0.42425,(CXCR5:0.50867,CXCR3:0.45350):0.06213):0.01352,CXCR4:0.56322):0.03272,((((CCR8:0.42058,((((CCR2:0.08743,CCR5:0.10207):0.16440,(CCR3:0.18374,CCR1:0.18158):0.07935):0.02277,CCRL2:0.52166):0.07748,CCR4:0.38947):0.05710):0.05044,CX3C1:0.38371):0.06819,(XCR1:0.54471,ACKR2:0.54304):0.04238):0.08031,(CXCR6:0.50090,((ACKR4:0.46379,(CCR7:0.42440,CCR9:0.38556):0.05089):0.04825,CCR6:0.47437):0.00461):0.12806):0.01989):0.14820,((ACKR3:0.52247,GP182:0.59468):0.23813,GPER1:0.88006):0.01691):0.01708):0.02247,(BKRB2:0.51447,BKRB1:0.56345):0.29484):0.05314,(((((HCAR2:0.05555,HCAR3:0.10376):0.18659,HCAR1:0.35694):0.24507,OXER1:0.47038):0.11069,GPR31:0.61743):0.35707,(((((((GPR35:0.72834,(GPR55:0.66501,FFAR1:0.82706):0.21885):0.07208,LPAR5:0.68476):0.04084,GPR20:0.77307):0.07754,((PAR4:0.57467,(((PAR2:0.52345,PAR1:0.57365):0.02890,PAR3:0.54309):0.04601,P2RY8:0.59538):0.03028):0.12099,(FFAR2:0.42973,(FFAR3:0.00538,GPR42:0.00422):0.37178):0.49994):0.05897):0.01895,(((PSYR:0.57559,(GPR4:0.40370,OGR1:0.39889):0.18363):0.09103,GP132:0.70545):0.24656,(((((LPAR6:0.25135,LPAR4:0.32636):0.28952,(PTAFR:0.67270,(GP174:0.36331,P2Y10:0.36053):0.42764):0.08679):0.09420,(CLTR1:0.52254,CLTR2:0.54393):0.15646):0.01173,GPR17:0.65080):0.04675,(((P2RY1:0.51184,P2Y11:0.83945):0.06289,((P2RY2:0.15600,P2RY4:0.32583):0.25068,P2RY6:0.50242):0.14435):0.09686,(OXGR1:0.62196,SUCR1:0.62748):0.09278):0.13068):0.02514):0.04548):0.02914,GP183:0.77503):0.01692,(((GPR34:0.63511,((GPR87:0.35245,(P2Y14:0.36677,(P2Y12:0.26004,P2Y13:0.40774):0.03744):0.08916):0.20002,GP171:0.69543):0.15216):0.06879,GPR82:1.06732):0.09381,(GPR18:0.89407,GP141:1.71148):0.10901):0.06110):0.09217):0.07603):0.06065,((((LT4R1:0.34230,LT4R2:0.56786):0.20423,GP152:1.16379):0.13340,(((C3AR:0.47649,(C5AR1:0.41334,C5AR2:0.63635):0.17172):0.02610,PD2R2:0.67680):0.03026,((CML1:0.43811,GPR1:0.57985):0.07689,((GPR32:0.71839,((FPR2:0.14255,FPR1:0.18677):0.03982,FPR3:0.22440):0.34075):0.07286,GPR33:0.82218):0.02563):0.02168):0.09214):0.14307,((CCR10:0.51827,ACKR1:1.92511):0.18804,(((GP162:0.40624,GP153:0.80663):1.72843,GPR75:1.26731):0.05945,(GPBAR:1.48202,(GP149:2.44392,GP160:1.72540):0.31634):0.18513):0.09889):0.06230):0.04630):0.04141):0.02581,(((((SSR1:0.17885,SSR4:0.18260):0.11696,((SSR3:0.20393,SSR5:0.23343):0.05093,SSR2:0.28178):0.06928):0.12990,(((OPRK:0.12704,(OPRD:0.15715,OPRM:0.15980):0.04678):0.06932,OPRX:0.26870):0.24080,(NPBW2:0.18687,NPBW1:0.20734):0.44035):0.03465):0.10010,(MCHR1:0.52414,MCHR2:0.65206):0.24588):0.12716,(GALR1:0.44912,((GALR2:0.20275,GALR3:0.25153):0.32957,KISSR:0.71877):0.09052):0.14818):0.06828):0.04164,((((((QRFPR:0.68527,((OX2R:0.10301,OX1R:0.13914):0.48857,(NPFF2:0.21754,NPFF1:0.21929):0.38045):0.05377):0.03012,((PKR2:0.02894,PKR1:0.05433):0.94818,(((NPY2R:0.55492,GPR83:0.61282):0.05075,PRLHR:0.52930):0.02798,((NPY4R:0.39508,(NPY1R:0.28705,NPY6R:0.46511):0.11830):0.21474,NPY5R:0.70556):0.04245):0.10298):0.03365):0.02436,((NK1R:0.19827,NK3R:0.17196):0.05848,NK2R:0.28508):0.64876):0.05066,((((EDNRA:0.27381,EDNRB:0.12535):0.49073,(GPR37:0.16443,ETBR2:0.32750):0.72694):0.17684,((GRPR:0.23985,NMBR:0.15452):0.06152,BRS3:0.23636):0.38078):0.25329,((MTLR:0.33253,GHSR:0.28880):0.35765,(((NTR1:0.30166,NTR2:0.48134):0.19543,GPR39:0.79339):0.20210,(NMUR2:0.26599,NMUR1:0.31573):0.37745):0.08651):0.21233):0.04012):0.03969,(((GASR:0.22694,CCKAR:0.25077):0.47959,(((TRFR:0.72308,(GP151:0.98114,GP148:2.65017):0.40020):0.18809,GP146:1.53033):0.20362,(GPR84:1.09609,GPR88:1.55637):0.25109):0.01643):0.05945,((GPR19:0.91434,(GPR21:0.11729,GPR52:0.15469):1.06171):0.05945,(((GP173:0.24912,GPR85:0.24442):0.22342,GPR27:0.44745):0.95429,(GP142:0.39179,GP139:0.29840):1.13574):0.05795):0.09169):0.02103):0.00363,(((NPSR1:0.63239,GP150:1.02736):0.02766,(((V1AR:0.23225,V1BR:0.25766):0.04238,OXYR:0.29802):0.17148,V2R:0.51306):0.25307):0.13452,(GNRHR:0.94238,GNRR2:0.34480):0.53223):0.27855):0.01449):0.03852):0.02177):0.12428):0.07151):0.05842):0.05872):0.01847):0.00164):0.01488):0.01919):0.06475):0.03004):0.08337,5HT1A:0.36804)OROOT;"
 x = newick.parse(x);
var nodes = cluster.nodes(x);

  nodes.forEach(function(n){

    if (n.name==""){
      n.name = names.toString();
      names++;
    }

  });

//Uncomment the line below to show branch length
// phylo(nodes[0], 0);

var link = vis.selectAll("path.link")
    .data(cluster.links(nodes))
  .enter().append("path")
    .attr("class", "link")
    .attr("d", step);


//D3 selection to distinguish between inner and leaf nodes
var node = vis.selectAll("g.node")
    .data(nodes.filter(function(n) { return n.x !== undefined; }))
  .enter().append("g")
    .attr("class", function(n) {
        if (n.children) {
          return "inner node";
        } else {
          //return "leaf node";
            return ('X'+n.name);
        }
      })
    .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })

//node.append("circle")
    //.attr("r", 2.5);

var innernodes= vis.selectAll('g.inner.node')
      .append("circle")
      .attr("r", 2.5);


DrawSelectivity(selectivitydata);


//Color code for Gprot selectivity
function DrawSelectivity(selectivityinfo) {


for (var x in selectivityinfo){
    //console.log(selectivityinfo[x]);

    var spacer = 8


    if(selectivityinfo[x].indexOf("Gs family") >= 0){
    var leafwithname = vis.selectAll('g.X'+x)
        .append("circle")
        .attr("r", 3.25)
        .style("fill", "blue")
        .attr("transform", "translate(" + (23 + spacer) + ",0)");

      }

    if(selectivityinfo[x].indexOf("Gi/Go family") >= 0){
    var leafwithname = vis.selectAll('g.X'+x)
        .append("circle")
        .attr("r", 3.25)
        .style("fill", "red")
        .attr("transform", "translate(" + (23  + 2*spacer) + ",0)");

      }

    if(selectivityinfo[x].indexOf("Gq/G11 family") >= 0){
    var leafwithname = vis.selectAll('g.X'+x)
        .append("circle")
        .attr("r", 3.25)
        .style("fill", "black")
        .attr("transform", "translate(" + (23 + 3*spacer) + ",0)");

      }

    if( selectivityinfo[x].indexOf("G12/G13 family") >= 0){
    var leafwithname = vis.selectAll('g.X'+x)
        .append("circle")
        .attr("r", 3.25)
        .style("fill", "green")
        .attr("transform", "translate(" + (23 + 4*spacer) + ",0)");

      }

    }
}


  var label = vis.selectAll("text")
      .data(nodes.filter(function(d) { return d.x !== undefined && !d.children; }))
    .enter().append("text")
      .attr("dy", ".31em")
      .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
      .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + (r - 170 + 2) + ")rotate(" + (d.x < 180 ? 0 : 180) + ")"; })
      .text(function(d) { return d.name.replace(/_/g, ' '); });


//click count to change clade color selection and unhighlight subtree nodes
  var click_count =0;

//On click event to get subtree of clicked node
    node.on("click", function(d){

      //doubleclade(d);
      singleclade(d);

  }

    );

//Highlight clade when mouse is over node
node.on("mouseover", function(h){

  traverse(h);

});

//Unhighlight clade when mouse is out of node
node.on("mouseout", function(u){

  var unhighlight = vis.selectAll("path.link")
    .style("stroke", "#ccc")
    .style("stroke-width", "1.5");

  // var unhighlight_labels = vis.selectAll("text").style("font-weight", "normal");

});


//iterate over links for (source,target) info to highlight subtree path
function highlightlink(src,tgt){
                      var link = d3.selectAll("path.link")[0].filter(function(d){
                         var flag = (d3.select(d).data()[0].source.name == src && d3.select(d).data()[0].target.name == tgt);
                        return flag;
                      });

                      d3.selectAll(link)
                        .style("stroke", "#000000")
                        .style("stroke-width", "2r");

                    }
//iterate over labels for leaf nodes info to highlight labels
function highlightlabel(tgt){

var highlight_labels = vis.selectAll("text")[0].filter(function(t){
                         var th = (d3.select(t).data()[0].name == tgt);
                        
                        return th;
                      });

                       d3.selectAll(highlight_labels)
                         .style("font-weight", "bold");
}

//Recursive function to highlight all links of a subtree
function traverse(node){


if(node.children){

node.children.forEach(function(d)
{

      
      //highlightlabel(d.name);
      highlightlink(node.name, d.name);

        traverse(d);
});

  }


}

//Recursive function to traverse a subtree and deselect all nodes
function deselect(node){

node.isSelected = false;

if(node.children){

node.children.forEach(function(d)
{
    //Select all children
  d.isSelected = false;

        deselect(d);
});}
    

}


//Recursive function to traverse a node and its subtree
function selectSubtree(node, subtree){

node.isSelected = true;

if(node.children){

node.children.forEach(function(d)
{

  d.isSelected = true;
  highlightlabel(d.name)
  
  if (d.name.length>3){
  subtree.push(d.name);
    }

  selectSubtree(d, subtree);

});

  }
    
return subtree;
}

//Single clade selection
function singleclade(node){

if (click_count>0){

  var deselection = vis.selectAll('g.inner.node').selectAll("circle")
    .attr("r", 2.5)
    .style("fill", "white");

  var unhighlight_labels = vis.selectAll("text").style("font-weight", "normal");


}


      var subtree = [];
      var subtreeArray = selectSubtree(node, subtree);

      for (i = 0; i < subtreeArray.length; i++)
        subtreeArray[i] = subtreeArray[i].toLowerCase()+"_human"


  var subnode = vis.selectAll('g.inner.node')
    .attr("id", function(d){
      if (d.isSelected==true)
      {
        return 'sub1';}
      

    });


    subnode1=vis.selectAll('#sub1').selectAll("circle")
    .attr("r", 2.5)
    .style("fill", "#000000");


      // console.log(subtreeArray);

      makeUL(subtreeArray);

      deselect(node);

  click_count++;

}

//Double clade selection
function doubleclade(node){

  if (click_count%2==0 && click_count>1){

  var deselection = vis.selectAll('g.inner.node').selectAll("circle")
    .attr("r", 2.5)
    .style("fill", "white");

  deselection = vis.selectAll('g.inner.node')
  .attr("id", "x");

}


      var subtree = [];
      var subtreeArray = selectSubtree(node, subtree);

      for (i = 0; i < subtreeArray.length; i++)
        subtreeArray[i] = subtreeArray[i]+"_human"


  var subnode = vis.selectAll('g.inner.node')
    .attr("id", function(d){
      if (d.isSelected==true && click_count%2==0)
      {
        return 'sub1';}
      else if (d.isSelected==true && click_count%2==1)
      {
        return 'sub2';}

    });




if (click_count % 2 == 0){

    subnode1=vis.selectAll('#sub1').selectAll("circle")
    .attr("r", 2.5)
    .style("fill", "red");

    subnode1=vis.selectAll('g#sub1.inner.node')
    .attr("id", "deselect");

      // console.log(subtreeArray);

      makeUL(subtreeArray);


  deselect(node);
  click_count++;

    }


  else if (click_count % 2 == 1) {

  var subtree = [];

    subnode2=vis.selectAll('#sub2').selectAll("circle")
    .attr("r", 2.5)
    .style("fill", "blue");

    subnode2=vis.selectAll('g#sub2.inner.node')
    .attr("id", "deselect");

      // console.log(subtreeArray);
      makeUL(subtreeArray);
      

  deselect(node);
  click_count++;

  }

}

//Fill text area with subtree data
function makeUL(descendents){

var textarea = document.getElementById('input-targets');
textarea.value = descendents.join("\n");


}


