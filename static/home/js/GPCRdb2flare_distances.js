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
      size: 1,
      segment: data.segment_map[residue]
    });

    dataFlare.trees[0].treePaths.push(
      data.segment_map[residue] + "." + residue
    );
  });

  // Fill edges
  Object.keys(data.interactions).reverse().forEach(function(pair) {
    pairResidues = separatePair(pair);

    // merge all interactions in data
    if ("max_dispersion" in data) {
      // Single group has a max_dispersion value
      var dispersion = data.interactions[pair][1];
      var frequency = dispersion / data.max_dispersion;
      var rgb = { r: 255, g: Math.round(255-frequency*255), b: Math.round(255-frequency*255) };
    } else if ("max_diff" in data) {
      // two groups has a max_dispersion value
      var difference = data.interactions[pair][0];
      var frequency = difference / data.max_diff;
      var rgb = getGradientColor(frequency, true);
    }
    // console.log(pair,data.interactions[pair]);
    if ("frequency" in data){
        dataFlare.edges.push({
          name1: pairResidues[0],
          name2: pairResidues[1],
          frames: [0],
          color: "#A0A0A0", // Default gray coloring of edges
          interactions: interactions, // For frequency and type coloring
          // split between 1 and 2 groups
          frequency: data.frequency[pair],
          count: data.count[pair],
          segment: assignColor(data.segment_map[pairResidues[0]]), // Segment coloring
        });
    } else {
        dataFlare.edges.push({
          name1: pairResidues[0],
          name2: pairResidues[1],
          frames: [0],
          color: "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")", // Default gray coloring of edges
          // interactions: interactions, // For frequency and type coloring
          // // split between 1 and 2 groups
          // frequency: interactions["any"]/data.pdbs.length,
          // count: interactions["any"],
          segment: assignColor(data.segment_map[pairResidues[0]]), // Segment coloring
        });
    }
    return
  });

  return dataFlare;
}
