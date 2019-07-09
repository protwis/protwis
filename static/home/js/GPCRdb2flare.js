function parseGPCRdb2flare(data) {
  if (typeof data == "string") {
    data = JSON.parse(data);
  }

  var dataFlare = {
    edges: [],
    tracks: [
      {
        trackLabel: "Secondary Structure",
        trackProperties: []
      }
    ],
    trees: [
      {
        treeLabel: "Secondary Structure",
        treePaths: []
      }
    ],
    default: {
      edgeColor: "rgba(100,100,100,100)",
      edgeWidth: 1
    }
  };

  var segmentColors = {
    "1": "#A8DDB5",
    "2": "#91C6AD",
    "3": "#7AB0A6",
    "4": "#63999E",
    "5": "#4C8397",
    "6": "#356C8F",
    "7": "#1E5688",
    "8": "#084081",
    "0": "#D0D0D0"
  };

  var segmentGrayColors = {
    "1": "#ccc",
    "2": "#bbb",
    "3": "#aaa",
    "4": "#999",
    "5": "#888",
    "6": "#777",
    "7": "#666",
    "8": "#555",
    "0": "#FCC0C0"
  };

  var segmentRainbowColors = {
    "1": "#1500D6",
    "2": "#006BDB",
    "3": "#00E1D1",
    "4": "#00E74F",
    "5": "#39ED00",
    "6": "#C9F300",
    "7": "#F99100",
    "8": "#FF0000",
    "0": "#EEE"
  };

  var segmentRainbowColors2 = {
    "1": "#736DA7",
    "2": "#5EB7B7",
    "3": "#CE9AC6",
    "4": "#DD7D7E",
    "5": "#E6AF7C",
    "6": "#DEDB75",
    "7": "#80B96F",
    "8": "#000000",
    "0": "#EEE"
  };

  function assignColor(segment) {
    var color = "";
    switch (segment) {
      case "TM1":
        color = segmentColors["1"];
        break;
      case "TM2":
        color = segmentColors["2"];
        break;
      case "TM3":
        color = segmentColors["3"];
        break;
      case "TM4":
        color = segmentColors["4"];
        break;
      case "TM5":
        color = segmentColors["5"];
        break;
      case "TM6":
        color = segmentColors["6"];
        break;
      case "TM7":
        color = segmentColors["7"];
        break;
      case "H8":
        color = segmentColors["8"];
        break;
      default:
        color = segmentColors["0"];
    }
    return color;
  }

  function assignRainbowColor(segment) {
    var color = "";
    switch (segment) {
      case "TM1":
        color = segmentRainbowColors2["1"];
        break;
      case "TM2":
        color = segmentRainbowColors2["2"];
        break;
      case "TM3":
        color = segmentRainbowColors2["3"];
        break;
      case "TM4":
        color = segmentRainbowColors2["4"];
        break;
      case "TM5":
        color = segmentRainbowColors2["5"];
        break;
      case "TM6":
        color = segmentRainbowColors2["6"];
        break;
      case "TM7":
        color = segmentRainbowColors2["7"];
        break;
      case "H8":
        color = segmentRainbowColors2["8"];
        break;
      default:
        color = segmentRainbowColors2["0"];
    }
    return color;
  }

  function separatePair(stringPair) {
    var regex = /([0-9x]+),([0-9x]+)/g;
    var m;

    matches = regex.exec(stringPair);

    return [matches[1], matches[2]];
  }

  // Fill tracks and trees
  Object.keys(data.segment_map).forEach(function(residue) {
    dataFlare.tracks[0].trackProperties.push({
      nodeName: residue,
      color: assignColor(data.segment_map[residue]),
      rainbow: assignRainbowColor(data.segment_map[residue]),
      size: 1,
      segment: data.segment_map[residue]
    });

    dataFlare.trees[0].treePaths.push(
      data.segment_map[residue] + "." + residue
    );
  });

  // Fill edges
  Object.keys(data.interactions).forEach(function(pair) {
    pairResidues = separatePair(pair);

    var values = data.interactions[pair];
    // merge all interactions in data
    // var allInteractions = [];
    // Object.keys(data.interactions[pair]).forEach( key => {
    //   var tmp = [];
    //   for (var i = 0; i < data.interactions[pair][key].length; i++) {
    //     tmp.push(getFriendlyInteractionName(data.interactions[pair][key][i]).replace(' ', '-'));
    //   }
    //   allInteractions.push(Array.from(new Set(tmp)));
    // });
    // allInteractions = [].concat.apply([],allInteractions); // flatten

    // // Essential for multiple x-rays in group - count num X-rays with interaction
    // var interactions = { any : Object.keys(data.interactions[pair]).length };
    // new Set(allInteractions).forEach( i => { interactions[i] = 0; });
    // allInteractions.forEach( i => { interactions[i] = interactions[i] + 1; });

    if ("pdbs2" in data){
        c1 = values['pdbs1'].length;
        c2 = values['pdbs2'].length;
        c3 = c1+c2;
        f1 = c1/data['pdbs1'].length;
        f2 = c2/data['pdbs2'].length;
        f3 = f1 - f2;

        var frequency = [f1, f2, f3];
        var count = [c1, c2, c3];
        dataFlare.edges.push({
          name1: pairResidues[0],
          name2: pairResidues[1],
          frames: [0],
          color: "#A0A0A0", // Default gray coloring of edges
          interactions: values['types'], // For frequency and type coloring
          // split between 1 and 2 groups
          frequency: frequency,
          count: count,
          segment: assignColor(data.segment_map[pairResidues[0]]), // Segment coloring
          rainbow: assignRainbowColor(data.segment_map[pairResidues[0]]), // Rainbow coloring
        });
    } else {
        c1 = values['pdbs'].length;
        f1 = c1/data['pdbs'].length;
        dataFlare.edges.push({
          name1: pairResidues[0],
          name2: pairResidues[1],
          frames: [0],
          color: "#A0A0A0", // Default gray coloring of edges
          interactions: values['types'], // For frequency and type coloring
          frequency: f1,
          count: c1,
          segment: assignColor(data.segment_map[pairResidues[0]]), // Segment coloring
          rainbow: assignRainbowColor(data.segment_map[pairResidues[0]]), // Rainbow coloring
        });
    }
  });

  return dataFlare;
}
