bezier_cache = [];
function wherebezier(p04,p14,p24,step,stop,p34,allow_cache) {
    //https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    pos = 0;
    length = 0;
    var p = p04;
    xy = [0,0];

    // if (stop<0) {
    //     stop = lengthbezier(p04,p14,p24,step,p34)+stop;
    // }
    i = 0;
    while (pos <= 1) {
        if (length>stop) { //stop if it reached the length along the line
            break;
        }
        if (i in bezier_cache && allow_cache) {
            console.log('using cache!');
            xy = bezier_cache[i][0];
            length = bezier_cache[i][1];
            // xy = bezier_high(p04,p14,p24,p34,pos);
            // length += Math.sqrt( Math.pow(xy[0]-p[0],2) + Math.pow(xy[1]-p[1],2) );
        } else {

            if (p34 === undefined) {
                xy = bezier(p04,p14,p24,pos);
            } else {
                xy = bezier_high(p04,p14,p24,p34,pos);
            }
            // console.log(xy,p);
            length += Math.round(100*Math.sqrt( Math.pow(xy[0]-p[0],2) + Math.pow(xy[1]-p[1],2) ))/100;
        }
        p = xy;
        pos += step;
        bezier_cache[i] =  [xy,length]
        i += 1;
        // console.log('pos '+pos+' length '+length+' stop at '+stop);
    }
    return [xy,length];
}

function lengthbezier(p03,p13,p23,step,p33) {
    //https://en.wikipedia.org/wiki/B%C3%A9zier_curve

    pos = 0;
    length = 0;
    p = p03;
    while (pos <= 1) {

        if (p33 === undefined) {
            xy = bezier(p03,p13,p23,pos);
        } else {
            xy = bezier_high(p03,p13,p23,p33,pos);
        }

        length += Math.sqrt( Math.pow(xy[0]-p[0],2) + Math.pow(xy[1]-p[1],2) );
        p = xy;
        pos += step;
    }
    return Math.round(length);
}

function bezier(p02,p12,p22,t){
    //https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    v1x = p12[0]-p02[0];
    v1y = p12[1]-p02[1];

    i12 = [p02[0]+(p12[0]-p02[0])*t,p02[1]+(p12[1]-p02[1])*t];
    i22 = [p12[0]+(p22[0]-p12[0])*t,p12[1]+(p22[1]-p12[1])*t];

    return [i12[0]+(i22[0]-i12[0])*t,i12[1]+(i22[1]-i12[1])*t];
}

function bezier_high(p01,p11,p21,p31,t){
    //https://en.wikipedia.org/wiki/B%C3%A9zier_curve
    i1 = bezier(p01,p11,p21,t);
    i2 = bezier(p11,p21,p31,t);

    return [i1[0]+(i2[0]-i1[0])*t,i1[1]+(i2[1]-i1[1])*t];
}

function add_box(x,y,name,color){
    box_text = "";
    name = name.substring(0, 6);
    n = name.length;
    width = n*7+2;
    box_text += "<rect class='"+name+" construct_custom' x="+(x-width/2)+" y="+(y-10)+" rx=5 ry=5 width='"+width+"' height='20' stroke='"+color+"' fill='white' stroke-width='2' style2='fill:red;stroke:black;stroke-width:5;opacity:0.5'/>";
    box_text += "<text  class='"+name+" construct_custom' x="+(x)+" y="+(y+4)+" text-anchor='middle' font-size=12 font-family='helvetica'>"+name+"</text>";
    return box_text;
}