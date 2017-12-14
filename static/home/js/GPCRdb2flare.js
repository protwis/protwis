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
    "1": "#1500D6",
    "2": "#006BDB",
    "3": "#00E1D1",
    "4": "#00E74F",
    "5": "#39ED00",
    "6": "#C9F300",
    "7": "#F99100",
    "8": "#FF0000",
    "0": "#FFFFFF"
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
      size: 1
    });

    dataFlare.trees[0].treePaths.push(
      data.segment_map[residue] + "." + residue
    );
  });

  // Fill edges
  if (data.generic == false) {
    // data.generic == false in case of single crystal

    Object.keys(data.interactions).forEach(function(pair) {
      pairResidues = separatePair(pair);

      dataFlare.edges.push({
        name1: pairResidues[0],
        name2: pairResidues[1],
        frames: [0]
      });
    });
  } else {
    crystals = Object.keys(data.aa_map);

    Object.keys(data.interactions).forEach(function(pair) {
      pairResidues = separatePair(pair);

      var frames = [];

      Object.keys(data.interactions[pair]).forEach(function(crystal) {
        var frame = crystals.indexOf(crystal);
        frames.push(frame);
      });

      dataFlare.edges.push({
        name1: pairResidues[0],
        name2: pairResidues[1],
        frames: frames
      });
    });
  }

  return dataFlare;
}
