function createSchematicPlot(data, containerSelector, options, data1, data2) {
  const config = {
    w: 900, // Width of the circle
    h: 900, // Height of the circle
    margin: {
      top: 20,
      right: 20,
      bottom: 20,
      left: 20,
    },
    type: 'singleCrystal', // ['singleCrystal', 'singleCrystalGroup', 'twoCrystalGroups']
    isContiguousPlot: true,
  };

  // Put all of the options into a variable called config
  if (typeof options !== 'undefined') {
    Object.keys(options).forEach((i) => {
      if (typeof options[i] !== 'undefined') {
        config[i] = options[i];
      }
    });
  }

  let isGeneric = false;

  if (config.type !== 'singleCrystal') {
    isGeneric = true;
  }

  // Process the data
  if (typeof data === 'string') {
    data = JSON.parse(data);
  }

  const interactions = data.interactions;
  const segment_map_full = data.segment_map_full;
  const segment_map_full_gn = data.segment_map_full_gn;
  const sequence_numbers = data.sequence_numbers;
  const aa_map = data.aa_map[Object.keys(data.aa_map)[0]];
  const gen_map = data.generic_map;
  const gen_map_full = data.generic_map_full;
  const num_seq_numbers = Object.keys(data.sequence_numbers).length;

  // Compute segment offsets
  let i;

  const segments = [];

  let seg;
  let prevSeg = isGeneric
    ? segment_map_full_gn[sequence_numbers[0]]
    : segment_map_full[sequence_numbers[0]];
  let seqStart = 0;

  for (i = 0; i < num_seq_numbers; i++) {
    seg = isGeneric
      ? segment_map_full_gn[sequence_numbers[i]]
      : segment_map_full[sequence_numbers[i]];

    if (seg === prevSeg) {
      continue;
    }

    segments.push({
      seg: prevSeg,
      start: seqStart,
      end: i - 1,
    });

    seqStart = i;
    prevSeg = seg;
  }
  // Push last segment
  segments.push({
    seg: prevSeg,
    start: seqStart,
    end: i - 1,
  });

  const segmentList = [
    // 'N-term',
    'TM1',
    // 'ICL1',
    'TM2',
    // 'ECL1',
    'TM3',
    // 'ICL2',
    'TM4',
    // 'ECL2',
    'TM5',
    // 'ICL3',
    'TM6',
    // 'ECL3',
    'TM7',
    'H8',
    // 'C-term',
  ];


  //
  // $(containerSelector).html('')

  // Remove whatever chart with the same id/class was present before
  d3
    .select(containerSelector)
    .select('svg')
    .remove();

  // Initiate the SVG
  const svg = d3
    .select(containerSelector)
    .append('svg')
    .attr('width', config.w + config.margin.left + config.margin.right)
    .attr('height', config.h + config.margin.top + config.margin.bottom)
    .attr('class', 'schematic2d')
    .attr('style','margin: auto; display: block');

  // Draw Reisues
  let oldCol = 0;
  let oldI = 0;

  const colSpace = 120;

  const rectWidth = 30;
  const rectHeight = 14;

  const paths = svg.append('g');
  const residues = svg.append('g');

  const g = residues
    .selectAll('g')
    .data(isGeneric ? Object.keys(segment_map_full_gn) : Object.keys(segment_map_full))
    .enter()
    .append('g')
    .attr('class', d => `node aa-${d}`)
    .attr('data-aa', d => d)
    .attr('data-segment', d => (isGeneric ? segment_map_full_gn[d] : segment_map_full[d]))
    .attr('transform', (d, i) => {
      const col = segmentList.indexOf(isGeneric ? segment_map_full_gn[d] : segment_map_full[d]);

      const height = rectHeight + 1;

      if (oldCol !== col) {
        oldI = i - 1;
      }

      let x;
      let y;

      if (config.isContiguousPlot) {
        x = colSpace * col;

        if (col % 2 === 0) {
          y = (i - oldI) * height;
        } else {
          y = 500 - (i - oldI) * height;
        }
      } else {
        // non-contiguous plot
        switch (col) {
          case 0:
            x = 0;
            y = config.margin.top + 250 + (i - oldI) * height;
            break;
          case 1:
            x = 300;
            y = config.h - config.margin.bottom - 200 - (i - oldI) * height;
            break;
          case 2:
            x = 600;
            y = config.margin.top + 50 + (i - oldI) * height;
            break;
          case 3:
            x = 900;
            y = config.h - config.margin.bottom - 100 - (i - oldI) * height;
            break;
          case 4:
            x = 750;
            y = config.margin.top  + (i - oldI) * height;
            break;
          case 5:
            x = 450;
            y = config.h - config.margin.bottom - 400 - (i - oldI) * height;
            break;
          case 6:
            x = 150;
            y = config.margin.top + 0 + (i - oldI) * height;
            break;
          case 7:
            x = 0 + (i - oldI) * (rectWidth + 1);
            y = config.h - config.margin.bottom;
            break;
          default:
            x = -100;
            y = config.h - config.margin.bottom - (i - oldI) * height;
            break;
        }
      }

      oldCol = col;

      return `translate(${x},${y})`;
    });


  g
    .append('rect')
    .attr('width', rectWidth)
    .attr('height', rectHeight)
    .style('fill', 'white')
    .style('stroke', 'black');

  g
    .append('text')
    .attr('x', rectWidth / 2)
    .attr('y', rectHeight - 3)
    .attr('font-size', 10)
    .attr('text-anchor', 'middle')
    .text(d => (isGeneric ? d : gen_map_full[d]));
  // .text(d => d);

  switch (config.type) {
    case 'singleCrystal':
      // svg.style('background-color', '#f0f0f0');
      renderSchematicSingleCrystal(getInteractionsSingleCrystal());
      createLegendSingleCrystal();
      break;
    case 'singleCrystalGroup':
      // getInteractionsCrystalGroup();
      // svg.style('background-color', '#f0f0f0');
      renderSchematicSingleCrystalGroup();
      createLegendSingleCrystalGroup();
      break;
    case 'twoCrystalGroups':
      // getInteractionsCrystalGroup();
      // svg.style('background-color', '#f0f0f0');
      renderSchematicTwoCrystalGroups();
      createLegendTwoCrystalGroups();
      break;
    default:
      break;
  }



  function getInteractionsSingleCrystal() {
    const interactionsList = [];
    Object.keys(interactions).forEach((interaction) => {
      const pair = separatePair(interaction);
      if (pair[0] in segment_map_full && pair[1] in segment_map_full) {
        if (
          config.isContiguousPlot
            ? isContiguous(segment_map_full[pair[0]], segment_map_full[pair[1]])
            : isNonContiguous(segment_map_full[pair[0]], segment_map_full[pair[1]])
        ) {
          getInteractionTypesFromPdbObject(Object.values(interactions[interaction])).forEach((interactionType) => {
            const d = {
              pair: interaction,
              interactionType,
            };
            interactionsList.push(d);
          });
        }
      }
    });
    return interactionsList;
  }

  function getInteractionsCrystalGroup() {
    const interactionsList = [];
    Object.keys(interactions).forEach((interaction) => {
      const pair = separatePair(interaction);
      if (pair[0] in segment_map_full_gn && pair[1] in segment_map_full_gn) {
        if (isContiguous(segment_map_full_gn[pair[0]], segment_map_full_gn[pair[1]])) {
          getInteractionTypesFromPdbObject(Object.values(interactions[interaction])).forEach((interactionType) => {
            const d = {
              pair: interaction,
              interactionType,
            };
            interactionsList.push(d);
          });
        }
      }
    });
    return interactionsList;
  }

  function renderSchematicSingleCrystal(interactionsList) {
    paths
      .selectAll('path')
      .data(interactionsList)
      .enter()
      .append('path')
      .attr('d', (d) => {
        const coord = getCoordPair(d.pair);
        return `M ${coord.sourceX} ${coord.sourceY} L ${coord.targetX} ${coord.targetY}`;
      })
      .attr('data-source-segment', d => segment_map_full[separatePair(d.pair)[0]])
      .attr('data-target-segment', d => segment_map_full[separatePair(d.pair)[1]])
      .attr('data-interaction-type', d => d.interactionType)
      .style('stroke', (d) => {
        const rgb = getInteractionColor(d.interactionType);
        const hex = rgb2hex(rgb.r, rgb.g, rgb.b);
        return hex;
      })
      .style('stroke-width', '3')
      .style('opacity', '0.6')
      .attr(
        'class',
        d =>
          `${getFriendlyInteractionName(d.interactionType).replace(/ /g, '-')} edge edge-${d.pair}`,
      );
  }

  function renderSchematicSingleCrystalGroup() {
    paths
      .selectAll('path')
      .data(Object.keys(interactions))
      .enter()
      .append('path')
      .filter((d) => {
        const pair = separatePair(d);
        if (pair[0] in segment_map_full_gn && pair[1] in segment_map_full_gn) {
          if (
              config.isContiguousPlot
              ? isContiguous(segment_map_full_gn[pair[0]], segment_map_full_gn[pair[1]])
              : isNonContiguous(segment_map_full_gn[pair[0]], segment_map_full_gn[pair[1]])
             ) {
            return d;
          }
        }
      })
      .attr('d', (d) => {
        const coord = getCoordPair(d);
        return `M ${coord.sourceX} ${coord.sourceY} L ${coord.targetX} ${coord.targetY}`;
      })
      .attr('data-source-segment', d => segment_map_full_gn[separatePair(d)[0]])
      .attr('data-target-segment', d => segment_map_full_gn[separatePair(d)[1]])
      .attr('data-pdbs', d => Object.keys(interactions[d]))
      .attr('data-num-interactions', (d) => {
        const nInteractions = Object.keys(interactions[d]).length;
        return nInteractions;
      })
      .style('stroke', (d) => {
        const nInteractions = Object.keys(interactions[d]).length;
        const frequency = nInteractions / data.pdbs.length;
        return d3.interpolateReds(frequency / 1.5);
      })
      .style('stroke-width', '2')
      .style('opacity', '0.6')
      .attr('class', d => `edge edge-${d}`)
      .on('mouseover', function (d) {
        d3.select(this).classed('highlighted', true);

        const coord = getCoordPair(d);

        const xPosition = (coord.sourceX + coord.targetX) / 2;
        const yPosition = (coord.sourceY + coord.targetY) / 2 + 20;

        // Update the tooltip position and value
        svg
          .append('text')
          .attr('id', 'tooltip')
          .attr('x', xPosition)
          .attr('y', yPosition)
          .attr('text-anchor', 'middle')
          .attr('font-family', 'sans-serif')
          .attr('font-size', '11px')
          .attr('font-weight', 'bold')
          .attr('fill', 'blue')
          .text($(this).data('pdbs'));

        // TODO use jQuery UI
        // $(this).tooltip({
        //   container: containerSelector,
        //   placement: 'top',
        //   delay: 0,
        //   html: true,
        //   title: 'demo',
        // });
      })
      .on('mouseout', function (d) {
        d3.select(this).classed('highlighted', false);
        d3.select('#tooltip').remove();
      });
  }

  function renderSchematicTwoCrystalGroups() {
    paths
      .selectAll('path')
      .data(Object.keys(interactions))
      .enter()
      .append('path')
      .filter((d) => {
        const pair = separatePair(d);
        if (pair[0] in segment_map_full_gn && pair[1] in segment_map_full_gn) {
          if (
              config.isContiguousPlot
              ? isContiguous(segment_map_full_gn[pair[0]], segment_map_full_gn[pair[1]])
              : isNonContiguous(segment_map_full_gn[pair[0]], segment_map_full_gn[pair[1]])
             ) {
            return d;
          }
        }
      })
      .attr('class', d => `edge edge-${d}`)
      .attr('d', (d) => {
        const coord = getCoordPair(d);
        return `M ${coord.sourceX} ${coord.sourceY} L ${coord.targetX} ${coord.targetY}`;
      })
      .attr('data-source-segment', d => segment_map_full_gn[separatePair(d)[0]])
      .attr('data-target-segment', d => segment_map_full_gn[separatePair(d)[1]])
      .attr('data-pdbs', d => Object.keys(interactions[d]))
      .attrs((d) => {
        let n1 = 0;
        let n2 = 0;

        if (d in data1.interactions) {
          n1 = Object.keys(data1.interactions[d]).length;
        }

        if (d in data2.interactions) {
          n2 = Object.keys(data2.interactions[d]).length;
        }

        if (d in data1.interactions || d in data2.interactions) {
          const f1 = n1 / data1.pdbs.length;
          const f2 = n2 / data2.pdbs.length;
          const fDiff = f1 - f2;

          return {
            'data-frequency-diff': fDiff,
            'data-group-1-num-ints': n1,
            'data-group-2-num-ints': n2,
            'data-group-1-num-pdbs': data1.pdbs.length,
            'data-group-2-num-pdbs': data2.pdbs.length,
            'data-group-1-freq': f1.toFixed(2),
            'data-group-2-freq': f2.toFixed(2),
          };
        }
      })

      .style('stroke', (d) => {
        let n1 = 0;
        let n2 = 0;

        if (d in data1.interactions) {
          n1 = Object.keys(data1.interactions[d]).length;
        }

        if (d in data2.interactions) {
          n2 = Object.keys(data2.interactions[d]).length;
        }

        if (d in data1.interactions || d in data2.interactions) {
          const f1 = n1 / data1.pdbs.length;
          const f2 = n2 / data2.pdbs.length;
          const fDiff = f1 - f2;

          let rgb;
          if (fDiff <= 0) {
            // If fDiff is close to -1, we want a red color
            rgb = { r: 255, g: Math.round(255 - 255 * -fDiff), b: Math.round(255 - 255 * -fDiff) };
          } else {
            // If fDiff is close to 1 we want a blue color
            rgb = { r: Math.round(255 - 255 * fDiff), g: Math.round(255 - 255 * fDiff), b: 255 };
          }

          return `rgb(${[rgb.r, rgb.g, rgb.b].join(',')})`;
        }
      })
      .style('stroke-width', '2')
      .on('mouseover', function (d) {
        d3.select(this).classed('highlighted', true);

        const coord = getCoordPair(d);

        const xPosition = (coord.sourceX + coord.targetX) / 2;
        const yPosition = (coord.sourceY + coord.targetY) / 2 + 20;

        // Update the tooltip position and value
        svg
          .append('text')
          .attr('id', 'tooltip')
          .attr('x', xPosition)
          .attr('y', yPosition)
          .attr('text-anchor', 'middle')
          .attr('font-family', 'sans-serif')
          .attr('font-size', '11px')
          .attr('font-weight', 'bold')
          .attr('fill', 'blue')
          .text($(this).data('pdbs'));

        // TODO use jQuery UI
        // $(this).tooltip({
        //   container: containerSelector,
        //   placement: 'top',
        //   delay: 0,
        //   html: true,
        //   title: 'demo',
        // });
      })
      .on('mouseout', function (d) {
        d3.select(this).classed('highlighted', false);
        d3.select('#tooltip').remove();
      });
  }

  function repositionSegment(segment) {
    if (segment !== 'TM1') {
      let run = 0;
      while (run < 10) {
        let gradientSum = 0;
        let contactCount = 0;

        const done_paths = [];
        $(containerSelector + ` path[data-target-segment='${segment}']`).each((i, path) => {
          const d = path.getAttribute('d');

          const regex = /M (.+) (.+) L (.+) (.+)/;
          const matches = regex.exec(d);

          const gradient = (matches[4] - matches[2]) / (matches[3] - matches[1]);
          if (done_paths.includes(d)) {
              return
          }
          done_paths.push(d);
          gradientSum += gradient;
          contactCount += 1;
        });

        const gradientMean = gradientSum / contactCount;
        const shiftY = gradientMean * colSpace*0.8;
        // if shift is small, then stop.
        if (shiftY<5 && shiftY>-5) run = 10;

        // Reposition each node
        $(containerSelector + ` g.node[data-segment='${segment}']`).each((i, g) => {
          const transformValue = g.getAttribute('transform');

          const regex = /\((.+),(.+)\)/;
          const matches = regex.exec(transformValue);

          const x = matches[1];
          const y = matches[2] - shiftY;

          if (Number.isNaN(x) || Number.isNaN(y)) {
            // FIXME: handle error
          } else {
            g.setAttribute('transform', `translate(${x},${y})`);
          }
        });

        // Reposition the edges TERMINATING at the repositioned nodes
        $(containerSelector + ` path[data-target-segment='${segment}']`).each((i, path) => {
          const d = path.getAttribute('d');

          const regex = /M (.+) (.+) L (.+) (.+)/;
          const matches = regex.exec(d);

          const newTargetY = matches[4] - shiftY;
          if (Number.isNaN(matches[1]) || Number.isNaN(matches[2]) || Number.isNaN(matches[3]) || Number.isNaN(newTargetY)) {
            // FIXME: handle error
          } else {
            path.setAttribute('d', `M ${matches[1]} ${matches[2]} L ${matches[3]} ${newTargetY}`);
          }
        });

        // Reposition the edges ORIGINATING at the repositioned nodes
        $(containerSelector + ` path[data-source-segment='${segment}']`).each((i, path) => {
          const d = path.getAttribute('d');

          const regex = /M (.+) (.+) L (.+) (.+)/;
          const matches = regex.exec(d);

          const newSourceY = matches[2] - shiftY;

          if (Number.isNaN(matches[1]) || Number.isNaN(matches[3]) || Number.isNaN(matches[4]) || Number.isNaN(newSourceY)) {
            // FIXME: handle error
          } else {
            path.setAttribute('d', `M ${matches[1]} ${newSourceY} L ${matches[3]} ${matches[4]}`);
          }
        });
        run += 1;
      }
    }
  }


  if (config.isContiguousPlot) {
    // Reposition the helices by minimizing the sum of gradients of contacts
    segmentList.forEach(s => repositionSegment(s));
    // console.log('=================');
  }

  if (config.isContiguousPlot) {
    // console.log('^^^^^^^^^^^^^^^^^');
    // console.log('Gradient Before');
    // console.log('=================');
    // console.log('Gradient After');
    // console.log('vvvvvvvvvvvvvvvvv');
    segmentList.forEach((segment) => {
      if (segment !== 'TM1') {
        let gradientSum = 0;
        let contactCount = 0;
        const done_paths = [];
        $(`path[data-target-segment='${segment}']`).each((i, path) => {
          const d = path.getAttribute('d');

          const regex = /M (.+) (.+) L (.+) (.+)/;
          const matches = regex.exec(d);
          // Do not use the next path if an identical path already been used
          // Applies to interation pairs that have many types.
          if (done_paths.includes(d)) {
            return
          }
          done_paths.push(d);
          const gradient = (matches[4] - matches[2]) / (matches[3] - matches[1]);

          gradientSum += gradient;
          contactCount += 1;
        });

        $(`g.node[data-segment='${segment}']`).each((i, g) => {
          const transformValue = g.getAttribute('transform');

          const regex = /\((.+),(.+)\)/;
          const matches = regex.exec(transformValue);

          const y = parseFloat(matches[2]);

        });

        const gradientMean = gradientSum / contactCount;
        // console.log(segment,contactCount,gradientSum,-gradientMean);
      }
    });
  }

    max_y = 0;
    min_y = 0;
    segmentList.forEach((segment) => {
      $(containerSelector + ` g.node[data-segment='${segment}']`).each((i, g) => {
          const transformValue = g.getAttribute('transform');

          const regex = /\((.+),(.+)\)/;
          const matches = regex.exec(transformValue);

          const y = parseFloat(matches[2]);
          if ( y > max_y ) max_y = y;
          if ( y < min_y ) min_y = y;

        });
    });

    width = config.w + config.margin.left + config.margin.right
    height = max_y-min_y+ config.margin.top + config.margin.bottom

    svg.attr('height', height);
    svg.attr('viewBox', "0 "+(min_y- config.margin.top)+" "+width+" "+height);

  function getCoordPair(pair) {
    const AAs = separatePair(pair);
    const coordSourceAA = getCoordAA(AAs[0]);
    const coordTargetAA = getCoordAA(AAs[1]);

    if (config.isContiguousPlot) {
      return {
        sourceY: parseFloat(coordSourceAA.y) + rectHeight / 2,
        sourceX: parseFloat(coordSourceAA.x) + rectWidth,
        targetX: parseFloat(coordTargetAA.x),
        targetY: parseFloat(coordTargetAA.y) + rectHeight / 2,
      };
    }

    // Deduce side of rect to start line
    sourceXpush = 0;
    targetXpush = 0;
    if (coordSourceAA.x<coordTargetAA.x){
      sourceXpush = rectWidth;
    } else {
      targetXpush = rectWidth;
    }

    return {
      sourceY: parseFloat(coordSourceAA.y) + rectHeight / 2,
      sourceX: parseFloat(coordSourceAA.x) + sourceXpush,
      targetX: parseFloat(coordTargetAA.x) + targetXpush,
      targetY: parseFloat(coordTargetAA.y) + rectHeight / 2,
    };
  }

  function isContiguous(segment1, segment2) {
    const regex = /[TMH]+([0-9])/;

    if (!regex.exec(segment1) || !regex.exec(segment2)) {
      return false;
    }
    const segNo1 = regex.exec(segment1)[1];
    const segNo2 = regex.exec(segment2)[1];

    if (segNo2 - segNo1 === 1) {
      return true;
    }
    return false;
  }

  function isNonContiguous(segment1, segment2) {
    const regex = /[TMH]+([0-9])/;

    if (!regex.exec(segment1) || !regex.exec(segment2)) {
      return false;
    }
    const segNo1 = regex.exec(segment1)[1];
    const segNo2 = regex.exec(segment2)[1];

    const nonContiguousPairs = {
      1: [3, 4, 5, 6, 7, 8],
      2: [4, 5, 6, 7, 8],
      3: [1, 5, 6, 7, 8],
      4: [1, 2, 6, 7, 8],
      5: [1, 2, 3, 7, 8],
      6: [1, 2, 3, 4, 8],
      7: [1, 2, 3, 4, 5],
      8: [1, 2, 3, 4, 5, 6],
    };

    if (nonContiguousPairs[segNo1].indexOf(+segNo2) > -1) {
      return true;
    }

    return false;
  }

  function getCoordAA(aa) {
    const translate = d3.select(containerSelector + ' .aa-'+aa).attr('transform');

    const regex = /(-?[0-9]+),(-?[0-9]+)/;

    const matches = regex.exec(translate);
    return {
      x: matches[1],
      y: matches[2],
    };
  }

  function separatePair(stringPair) {
    const regex = /([0-9x]+),([0-9x]+)/;

    const matches = regex.exec(stringPair);

    return [matches[1], matches[2]];
  }

  function createLegendSingleCrystal() {
    let legendHtml = ""
    /*let interactionTypes = new Set();

    $(`${containerSelector} .edge`).each(function () {
      const friendlyName = getFriendlyInteractionName($(this).data('interaction-type'));
      interactionTypes.add(friendlyName);
    });

    // empty previous one
    $(`${containerSelector} .schematic-legend`).html('');

    // Add interactions color legend
    let legendHtml = '<ul>';

    interactionTypes = Array.from(interactionTypes).sort((i1, i2) => getInteractionStrength(i2) - getInteractionStrength(i1));

    interactionTypes.forEach((i) => {
      const rgb = getInteractionColor(i);
      legendHtml =
        `${legendHtml}<li>` +
        `<div class="color-box" style="background-color: ${rgb2hex(rgb.r, rgb.g, rgb.b)}">` +
        `<input type="checkbox" data-interaction-type="${i.replace(/ /g, '-')}"></input>` +
        `</div><p>${i}</p>` +
        '</li>';
    });
    legendHtml += '</ul>';
    */

    // // Add SVG download button
    // legendHtml +=
    //   `<button onclick="downloadSVG('${containerSelector}schematic', 'interactions.svg')" type="button" class="btn btn-primary pull-right svg-download-button" aria-label="Left Align">` +
    //   '<span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download SVG' +
    //   '</button>';

    // // Add CSV download button
    // legendHtml +=
    //   `<br /><button onclick="downloadSingleCrystalCSV('${containerSelector}schematic', 'interactions.csv')" type="button" class="btn btn-success pull-right csv-download-button" aria-label="Left Align"><span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download CSV` +
    //   '</button>';

    // $(`${containerSelector} .schematic-legend`).html(legendHtml);

    /*$(`${containerSelector} .schematic-legend input[type=checkbox]`).each(function () {
      $(this).prop('checked', true);
      $(this).change(function () {
        const interactionType = $(this).data('interaction-type');
        const paths = $(`${containerSelector} path.${interactionType}`);

        if ($(this).is(':checked')) {
          paths.show();
        } else {
          paths.hide();
        }
      });
    });*/
  }

  function createLegendSingleCrystalGroup() {
    // Populate schematic legend
    // Changed from separate min/max sliders to one range slider - to REMOVE if OK
    /*let legendHtml =
      `${'<h4 class="center">Interaction count</h4>' +
        '<p>From: <span class="min-value">0</span></p>' +
        '<input class="min-interactions-range" type="range" min="0" max="'}${
        data.pdbs.length
      }" value="0" step="1" />` +
      '<div class="temperature-scale">' +
      '<span class="white-to-red"></span>' +
      '</div>' +
      `<p>To: <span class="max-value">${data.pdbs.length}</span></p>` +
      `<input class="max-interactions-range" type="range" min="0" max="${
        data.pdbs.length
      }" value="${data.pdbs.length}" step="1" />` +
      '<div class="temperature-scale">' +
      '<span class="white-to-red"></span>' +
      '</div>';*/

    let legendHtml = "";
    /*let legendHtml = '<h4 class="center">Frequency (#PDBs)</h4>'
          + `<p>Range: <span id="clscg-pdbs-range">1 - ${data.pdbs.length}</span></p>`
          + '<div class="slider-range" data-text-id="clscg-pdbs-range" id="clscg-pdbs-range-slider"></div>'
          + '<div class="temperature-scale">'
          + '<span class="white-to-red"></span>'
          + '</div>';*/

    // // Add SVG download button
    // legendHtml +=
    //   `<button onclick="downloadSVG('${containerSelector} .schematic', 'interactions.svg')" type="button" class="btn btn-primary pull-right svg-download-button" aria-label="Left Align">` +
    //   '<span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download SVG' +
    //   '</button>';

    // // Add CSV download button
    // legendHtml +=
    //   `<br /><button onclick="downloadSingleCrystalGroupCSV('${containerSelector} .schematic', 'interactions.csv')" type="button" class="btn btn-success pull-right csv-download-button" aria-label="Left Align"><span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download CSV` +
    //   '</button>';

    // $(`${containerSelector} .schematic-legend`).html(legendHtml);

    /*
    // Changed from separate min/max sliders to one range slider - to REMOVE if OK
    function getRangeChangeFunction() {
      return function () {
        const tMin = $(`${containerSelector} .schematic-legend .min-interactions-range`).val();
        const tMax = $(`${containerSelector} .schematic-legend .max-interactions-range`).val();

        $(`${containerSelector} .schematic-legend .min-value`).html(tMin);
        $(`${containerSelector} .schematic-legend .max-value`).html(tMax);

        // Hide all below min treshold
        $(`${containerSelector} .edge`).each(function () {
          const n = $(this).data('num-interactions');
          if (n < tMin || tMax < n) {
            $(this).hide();
          } else {
            $(this).show();
          }
        });
      };
    }

    $(`${containerSelector} .schematic-legend .min-interactions-range`).change(getRangeChangeFunction());

    $(`${containerSelector} .schematic-legend .max-interactions-range`).change(getRangeChangeFunction());*/

    /*$( function() {
      $(`${containerSelector} #clscg-pdbs-range-slider`).slider({
        range: true,
        min: 1,
        max: `${data.pdbs.length}`,
        step: 1,
        values: [0, `${data.pdbs.length}`],
        slide: function( event, ui ) {
          $( `${containerSelector} #`+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
          getCLSCGRangeChangeFunction(ui.values[ 0 ], ui.values[ 1 ]);
        }
      });
    } );

    function getCLSCGRangeChangeFunction(tMin, tMax) {
          // Hide all below min or above max treshold
          $(`${containerSelector} .edge`).each(function () {
            const n = $(this).data('num-interactions');
            if (n < tMin || tMax < n) {
              $(this).hide();
            } else {
              $(this).show();
            }
          });
    }*/
  }

  function createLegendTwoCrystalGroups() {
    // Populate heatmap legend
    // Changed from separate min/max sliders to one range slider - to REMOVE if OK
    /*let legendHtml =
      '<h4 class="center">Frequency</h4>' +
      '<p>From: <span class="min-value">-1</span></p>' +
      '<input class="min-interactions-range" type="range" min="-1" max="1" value="-1" step="0.01" />' +
      '<div class="temperature-scale">' +
      '<span class="red-to-white"></span>' +
      '<span class="white-to-blue"></span>' +
      '</div>' +
      '<p>To: <span class="max-value">1</span></p>' +
      '<input class="max-interactions-range" type="range" min="-1" max="1" value="1" step="0.01" />' +
      '<div class="temperature-scale">' +
      '<span class="red-to-white"></span>' +
      '<span class="white-to-blue"></span>' +
      '</div>';*/

    /*let legendHtml = '<h4 class="center">Frequency</h4>'
          + `<p>Group 1 range: <span id="cltcg-group1-range">0 - 1</span></p>`
          + '<div class="slider-range" data-text-id="cltcg-group1-range" id="cltcg-group1-range-slider"></div>'
          + `<p>Group 2 range: <span id="cltcg-group2-range">0 - 1</span></p>`
          + '<div class="slider-range" data-text-id="cltcg-group2-range" id="cltcg-group2-range-slider"></div>'
          + `<p>Freq difference range:<span id="cltcg-diff-range"> -1 - 1</span></p>`
          + '<div class="slider-range" data-text-id="cltcg-diff-range" id="cltcg-diff-range-slider"></div>'
          + '<div class="temperature-scale">'
          + '<span class="red-to-gray"></span>'
          + '<span class="gray-to-blue"></span>'
          + '</div>';*/

    let legendHtml = ""

    // // Add SVG download button
    // legendHtml +=
    //   `<button onclick="downloadSVG('${containerSelector} .schematic', 'interactions.svg')" type="button" class="btn btn-primary pull-right svg-download-button" aria-label="Left Align">` +
    //   '<span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download SVG' +
    //   '</button>';

    // // Add CSV download button
    // legendHtml +=
    //   `<br /><button onclick="downloadTwoCrystalGroupsCSV('${containerSelector} .schematic', 'interactions.csv')" type="button" class="btn btn-success pull-right csv-download-button" aria-label="Left Align"><span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download CSV` +
    //   '</button>';

    // $(`${containerSelector} .schematic-legend`).html(legendHtml);

    // Changed from separate min/max sliders to one range slider - to REMOVE if OK
    /*
    function getRangeChangeFunction() {
      return function () {
        const tMin = $(`${containerSelector} .schematic-legend .min-interactions-range`).val();
        const tMax = $(`${containerSelector} .schematic-legend .max-interactions-range`).val();

        $(`${containerSelector} .schematic-legend .min-value`).html(tMin);
        $(`${containerSelector} .schematic-legend .max-value`).html(tMax);

        // Hide all below min treshold
        $(`${containerSelector} .edge`).each(function () {
          const f = $(this).data('frequency-diff');
          if (f < tMin || f > tMax) {
            $(this).hide();
          } else {
            $(this).show();
          }
        });
      };
    }

    $(`${containerSelector} .schematic-legend .min-interactions-range`).change(getRangeChangeFunction());

    $(`${containerSelector} .schematic-legend .max-interactions-range`).change(getRangeChangeFunction());*/

    /*$( function() {
      $( containerSelector+" #cltcg-group1-range-slider" ).data({ "referenceContainer" : containerSelector });
      $( containerSelector+" #cltcg-group2-range-slider" ).data({ "referenceContainer" : containerSelector });
      $( containerSelector+" #cltcg-diff-range-slider" ).data({ "referenceContainer" : containerSelector });
      $( containerSelector+" #cltcg-group1-range-slider" ).slider({
        range: true,
        min: 0,
        max: 1,
        step: 0.01,
        values: [0,1],
        slide: function( event, ui ) {
          $( $(this).data("referenceContainer")+" #"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
          getCLTCGRangeChangeFunction($(this)[0].id, ui, $(this).data("referenceContainer"));
        }
      });
      $( containerSelector+" #cltcg-group2-range-slider" ).slider({
        range: true,
        min: 0,
        max: 1,
        step: 0.01,
        values: [0,1],
        slide: function( event, ui ) {
          $( $(this).data("referenceContainer")+" #"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
          getCLTCGRangeChangeFunction($(this)[0].id, ui, $(this).data("referenceContainer"));
        }
      });
      $( containerSelector+" #cltcg-diff-range-slider" ).slider({
        range: true,
        min: -1,
        max: 1,
        step: 0.01,
        values: [-1,1],
        slide: function( event, ui ) {
          $( $(this).data("referenceContainer")+" #"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
          getCLTCGRangeChangeFunction($(this)[0].id, ui, $(this).data("referenceContainer"));
        }
      });
    });

    function getCLTCGRangeChangeFunction(origin, ui, container) {
        // grab data from all sliders
        var g1 = $( container+" #cltcg-group1-range-slider").slider("values");
        var g2 = $( container+" #cltcg-group2-range-slider" ).slider("values");
        var d = $( container+" #cltcg-diff-range-slider" ).slider("values");

        // update with current slide values
        if (origin == "cltcg-group1-range-slider")
            g1 = ui.values;
        else if (origin == "cltcg-group2-range-slider")
            g2 = ui.values;
        else if (origin == "cltcg-diff-range-slider")
            d = ui.values;

        // Hide all below min or above max treshold
        $(`${containerSelector} .edge`).each(function () {
          const group1 = $(this).data('group-1Freq');
          const group2 = $(this).data('group-2Freq');
          const diff = $(this).data('frequencyDiff');

          $(this).show();
          // group 1
          if (group1 < g1[0] || group1 > g1[1]) {
              $(this).hide();
          }
          // group 2
          else if (group2 < g2[0] || group2 > g2[1]) {
              $(this).hide();
          }
          // difference
          else if (diff < d[0] || diff > d[1]) {
            $(this).hide();
          }
        });
    }*/
  }
}
