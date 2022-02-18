function GeneFusionPlot(config) {

    // INITIALIZATION

    //var self = this;

    // Plot dimensions & target
    var targetElementID = config.targetElement;
    var cytobandFile = config.cytobandFile;
    var downloadID = config.downloadID;
    var targetElement = document.getElementById(targetElementID) || targetElementID || document.body // Where to append the plot (svg)

    var width = this.width = config.width || targetElement.offsetWidth || 1000;
    var height = this.height = config.height || targetElement.offsetHeight || 700;

    this.svgClasses = "genefusion" + targetElementID;
    var buffer = 150;

    this.cytobandFile = cytobandFile;
    this.cytobandData = config.cytobandData;
    this.buffer = buffer;    
    this.bg_height = 10;
    this.coord_bar_offset = 53;
    this.tick_offset = 10;
    this.bg_offset = 0;
    this.region_offset = this.bg_offset - 3;
    this.coord_margin = 100;
    this.exon_height = 16;
    this.domain_height = 30;
    if (isNaN(config.opacity))
        config.opacity = 0.2;
    else
        config.opacity = parseFloat(config.opacity);
    var has_protein = config.genes.has_protein;
    //this.width = width;


    // INIT SVG
    var svg;
    var topnode;
    if (config.responsive == 'resize') {
        topnode = d3.select(targetElement).append("svg")
            .attr("width", '100%')
            .attr("height", '100%')
            .attr('viewBox', '0 0 ' + Math.min(width) + ' ' + Math.min(height))
            .attr('class', 'brush');
        svg = topnode
            .append("g")
            .attr("class", this.svgClasses)
            .attr("transform", "translate(0,0)");
    } else {

        var svg = d3.select(targetElement).append("svg")
            .attr("width", width)
            .attr("height", height)
            .attr("class", this.svgClasses + " brush");
        topnode = svg;
    }


    // DEFINE SCALES

    minCoord = 999999999;
    maxCoord = 0;
    config.genes.gene1.exons.forEach(function(exon) {        
        if (exon.start_pos < minCoord)
            minCoord = exon.start_pos;
        if (exon.end_pos > maxCoord)
            maxCoord = exon.end_pos;
    });

    coord_margin = 100;
    this.minCoord = minCoord - coord_margin;
    this.maxCoord = maxCoord + coord_margin;
    //alert(this.minCoord + ',' + this.maxCoord);


    var y1 = d3.scale.linear()
        .domain([1, 100])
        .range([0, height])
        .nice();
    this.y = y1;

    var y2 = d3.scale.linear()
        .domain([1, 20])
        .range([0, buffer])
        .nice();

    function adjustCoordRange(coordRange, w) {
        var mid_point = (coordRange[0] + coordRange[1]) / 2;
        coordRange[0] = mid_point - w / 2;
        coordRange[1] = mid_point + w / 2;
        return coordRange;
    }

    function mapPoints(pointList) {
        return pointList.map(function(d) {
            return [d.x, d.y].join(",");
        }).join(" ");
    }

    function drawPolygon(pointList, color, className, opacity) {
        var points = mapPoints(pointList);
        var p = d3.select(className).selectAll()
            .data(["dummy"]).enter()
            .append("polygon")
            .attr("class", className)
            .attr("points", points)
            .attr("fill", color)
            .attr("opacity", opacity);
        return p;
    }

    function getChrShapeCoord(x1_top, x2_top, x1_down, x2_down, y_top, y_down) {
        return [{
            "x": x1_top,
            "y": y_top
        }, {
            "x": x2_top,
            "y": y_top
        }, {
            "x": x1_down,
            "y": y_down
        }, {
            "x": x2_down,
            "y": y_down
        }];
    }

    function getShapeCoord(x1_top, x2_top, x1_down, x2_down, y_top, y_down, bar_height) {
        return [
        {
            "x": x1_top,
            "y": y_top + bar_height
        }, {
            "x": x1_top,
            "y": y_top
        }, {
            "x": x2_top,
            "y": y_top
        }, {
            "x": x2_top,
            "y": y_top + bar_height
        }, {
            "x": x2_down,
            "y": y_down
        }, {
            "x": x2_down,
            "y": y_down + bar_height
        }, {
            "x": x1_down,
            "y": y_down + bar_height
        }, {
            "x": x1_down,
            "y": y_down
        }];
    }

    function writeDownloadLink(){
        //console.log("clicked...");
        var html = d3.select("svg")
                        .attr("title", "svg_title")
                        .attr("version", 1.1)
                        .attr("xmlns", "http://www.w3.org/2000/svg")
                        .node().parentNode.innerHTML;
                    
        d3.select(this)
                        .attr("href-lang", "image/svg+xml")
                        .attr("href", "data:image/svg+xml;base64,\n" + btoa(html))
                        .attr("target", "_blank");
    }

    function getExonPos(strand, exons, start=0) {
        var exon_pos = [];
        var length = start;
        //console.log(start);
        if (strand == '+') {
            for (var i = 0; i < exons.length - 1; i++) {
                if (exons[i].type == 'cds') {
                    length += exons[i].hint.Length;
                    //console.log(i + ':' + exons[i].start_pos + ' ' + exons[i].end_pos);
                    exon_pos.push(length);
                }
            }
        } else {
            for (var i = exons.length - 1; i > 1; i--) {
                if (exons[i].type == 'cds') {
                    length += exons[i].hint.Length;
                    //console.log(i + ':' + exons[i].start_pos + ' ' + exons[i].end_pos);
                    exon_pos.push(length);
                }
            }
        }
        return exon_pos;
    }
    //start the plot...

    gene1 = config.genes.gene1;
    gene2 = config.genes.gene2;
    
    gene1.exons = this.parseIntRegions(gene1.exons);    
    gene2.exons = this.parseIntRegions(gene2.exons);        

    for (var i = 0; i < gene1.exons.length; i++)
        gene1.exons[i].color = gene1.color;
    
    for (var i = 0; i < gene2.exons.length; i++)
        gene2.exons[i].color = gene2.color;

    fusedGene = this.getFusedGene(gene1, gene2);
    
    coordRange1 = this.getCoordRange(gene1.exons, gene1.junction);

    coord_width1 = coordRange1[1] - coordRange1[0] + 1;
    maxExp1 = this.getMaxExp(gene1.exons);
    //alert(maxExp1);

    coordRange2 = this.getCoordRange(gene2.exons, gene2.junction);
    maxExp2 = this.getMaxExp(gene2.exons);
    maxExp = Math.max(maxExp1, maxExp2);
    coord_width2 = coordRange2[1] - coordRange2[0] + 1;
    coordFusedRange = this.getCoordRange(fusedGene.exons, fusedGene.junction);
    coordFused_width = coordFusedRange[1] - coordFusedRange[0] + 1;

    //console.log(JSON.stringify(coordRange1));
    //console.log(JSON.stringify(coordRange2));
    //console.log(JSON.stringify(coordFusedRange));
    var coord_width = Math.max(coord_width1, coord_width2, coordFused_width);
    coordRange1 = adjustCoordRange(coordRange1, coord_width1);
    coordRange2 = adjustCoordRange(coordRange2, coord_width2);
    //coordFusedRange = [0, coord_width];
    coordFusedRange = [0, coordFused_width];
    

    coord_width1 = coordRange1[1] - coordRange1[0] + 1;
    coord_width2 = coordRange2[1] - coordRange2[0] + 1;
    coordFused_width = coordFusedRange[1] - coordFusedRange[0] + 1;
    title_dna_y = 5;
    gene1_y = 20;
    gene2_y = 38;
    fuse_y = 52;
    title_protein_y = 62;
    protein1_y = 70;
    protein2_y = 80;
    protein_fuse_y = 90;
    chr_offset = 13;

    if (!has_protein) {
        title_dna_y = 10;
        gene1_y = 35;
        gene2_y = 65;
        fuse_y = 85;
        chr_offset = 18;
    }
    x_pos = buffer;
    var y = this.y;
    var duration_time = 1000;
    var domain_height = this.domain_height;
    var exon_height = this.exon_height;

    var currentObj = this; 

    $(".d3-tip").remove();
    $(".d3-tip\\:after").remove();
    $(".d3-tip\\.n\\:after").remove();

    var tooltip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([-10, 0])
        .html(function(r) {
                //return "<strong>Frequency:</strong> <span style='color:red'>" + d.frequency + "</span>";
                var html = '<table><tr><td><a href="javascript:closeTooltip();">[close]</a></td><td></td></tr>';
                //var html = '<table>';
                var desc = r.hint;
                for (var prop in desc) {
                    if (desc.hasOwnProperty(prop)) {
                        html += '<tr><td><B>' + prop + ':</B></td><td><font color=red>' + desc[prop] + '</font></td></tr>';
                    }
                }
                html += '</table>';
                return html;
        });

    this.tooltip = tooltip;
    
    svg.call(tooltip);
    svg.on('click', tooltip.hide);
    
    this.drawTitle(svg, "DNA", 10, y(title_dna_y), 100, 30);
    if (has_protein)
        this.drawTitle(svg, "Protein", 10, y(title_protein_y), 100, 30);

	var chrRange1 = currentObj.getChromosomeRange(gene1.chr, this.cytobandData);
    var chrRange2 = currentObj.getChromosomeRange(gene2.chr, this.cytobandData);
    var maxChrCoord = Math.max(chrRange1[0], chrRange2[0]);        
    var range1 = currentObj.drawChromosome(svg, gene1.chr, this.cytobandData, maxChrCoord, x_pos, gene1_y - chr_offset, chrRange1, coordRange1, gene1.junction);
    var range2 = currentObj.drawChromosome(svg, gene2.chr, this.cytobandData, maxChrCoord, x_pos, gene2_y - chr_offset, chrRange2, coordRange2, gene2.junction);

    if (has_protein) {
            var maxProteinCoord = Math.max(gene1.pep_length, gene2.pep_length, config.genes.fused_pep.length);
            var bg1 = [{"start_pos" : 1, "end_pos" : (gene1.pep_length == 0)? 1: gene1.pep_length, "color" : gene1.color, "name" : ""}];
            var bg2 = [{"start_pos" : 1, "end_pos" : (gene2.pep_length == 0)? 1: gene2.pep_length, "color" : gene2.color, "name" : ""}];
            var bg_fuse = [];
            gene1.exon_pos_in_protein = getExonPos(gene1.strand, gene1.exons);            
            gene2.exon_pos_in_protein = getExonPos(gene2.strand, gene2.exons);
            //console.log(JSON.stringify(gene1.exon_pos_in_protein));
            //console.log(JSON.stringify(gene2.exon_pos_in_protein));
            if (gene1.pep_length > 0)
                bg_fuse.push({"start_pos" : 1, "end_pos" : gene1.pep_junction, "color" : gene1.color, "name" : ""});
            if (gene2.pep_length > 0) {
                //console.log("gene1.pep_junction" + gene1.pep_junction);
                bg_fuse.push({"start_pos" : ((gene1.pep_junction < 0)? 1 : gene1.pep_junction) , "end_pos" : config.genes.fused_pep.length, "color" : gene2.color, "name" : ""});
            }
            if (config.genes.domains[gene1.trans] == null)
                config.genes.domains[gene1.trans] = [];
            if (config.genes.domains[gene2.trans] == null)
                config.genes.domains[gene2.trans] = [];            
            if (config.genes.domains["fused"] == null)
                config.genes.domains["fused"] = [];
            config.genes.domains[gene1.trans] = currentObj.parseIntRegions(config.genes.domains[gene1.trans]);
            config.genes.domains[gene2.trans] = currentObj.parseIntRegions(config.genes.domains[gene2.trans]);
            config.genes.domains["fused"] = currentObj.parseIntRegions(config.genes.domains["fused"]);

            for (var i = 0; i < config.genes.domains[gene1.trans].length; i++)
                config.genes.domains[gene1.trans][i].color = "blue";
            for (var i = 0; i < config.genes.domains[gene2.trans].length; i++)
                config.genes.domains[gene2.trans][i].color = "blue";
            for (var i = 0; i < config.genes.domains["fused"].length; i++)
                config.genes.domains["fused"][i].color = "blue";
            protein1 = [0, 0, 0];
            //console.log("fused:" + JSON.stringify(config.genes.domains["fused"]));

            var fused_left_text = gene1.name;
            var fused_right_text = gene2.name;
            //if right gene intact
            //console.log(JSON.stringify(gene1.exon_pos_in_protein));
            fusedGene.exon_pos_in_protein = getExonPos(fusedGene.strand, fusedGene.exons, (gene1.exon_pos_in_protein.length == 0)?gene1.pep_length*3:0);
            //console.log(JSON.stringify(fusedGene.exon_pos_in_protein));
            if (gene1.pep_length == 0) {
                fused_left_text = gene2.name;
                fused_right_text = null;
                fusedGene.exon_pos_in_protein = gene2.exon_pos_in_protein;
                bg_fuse = bg2;
            }
            if (gene2.pep_length == 0) {
                fusedGene.exon_pos_in_protein = gene1.exon_pos_in_protein;                
                fused_right_text = null;
                bg_fuse = bg1;
            }
            fuse_protein = currentObj.drawRegions(svg, "."+ currentObj.svgClasses, config.genes.domains["fused"], "+", null, config.genes.fused_pep.length + 'aa', fused_left_text, fused_right_text, gene1.pep_junction, x_pos, protein_fuse_y, 0, 0.2, domain_height, [1,maxProteinCoord], 0, bg_fuse, false, true, 0, fusedGene.exon_pos_in_protein, 6);
            if (gene1.pep_length > 0) {
                protein1 = currentObj.drawRegions(svg, "."+ currentObj.svgClasses, config.genes.domains[gene1.trans], "+", null, gene1.pep_length + 'aa', gene1.name , null, gene1.pep_junction, x_pos, protein1_y, 0, 0.2, domain_height, [1,maxProteinCoord], 0, bg1, true, true, 0, gene1.exon_pos_in_protein, 6);
                var shapeProteinFrom1 = getShapeCoord(protein1[0], protein1[1], protein1[0], protein1[1], y(protein1_y), y(protein1_y), domain_height);
                var shapeProteinTo1 = getShapeCoord(protein1[0], protein1[1], fuse_protein[0], fuse_protein[1], y(protein1_y), y(protein_fuse_y), domain_height);
                var polyProtein1 = drawPolygon(shapeProteinFrom1, "blue", "."+ currentObj.svgClasses, config.opacity);
                var pointsProteinTo1 = mapPoints(shapeProteinTo1);
                var shapeProtein1Class = targetElementID + "_shapeProtein1";
                polyProtein1.transition().duration(duration_time).attr("points", pointsProteinTo1).attr('pointer-events', 'none').each("end", function() {
                //polyProtein2.transition().duration(duration_time).attr("points", pointsProteinTo2);
                    d3.selectAll('.regionGroup').moveToFront();
                });
            }
            if (gene2.pep_length > 0) {
                protein2 = currentObj.drawRegions(svg, "."+ currentObj.svgClasses, config.genes.domains[gene2.trans], "+", null, gene2.pep_length + 'aa', gene2.name, null, gene2.pep_junction, x_pos, protein2_y, 0, 0.2, domain_height, [1,maxProteinCoord], 0, bg2, false, true, 0, gene2.exon_pos_in_protein, 6);
                var shapeProteinFrom2 = getShapeCoord(protein2[1], protein2[2], protein2[1], protein2[2], y(protein2_y), y(protein2_y), domain_height);
                var shapeProteinTo2 = getShapeCoord(protein2[1], protein2[2], fuse_protein[1], fuse_protein[2], y(protein2_y), y(protein_fuse_y), domain_height);
                var polyProtein2 = drawPolygon(shapeProteinFrom2, "red", "."+ currentObj.svgClasses, config.opacity);
                var pointsProteinTo2 = mapPoints(shapeProteinTo2);
                var shapeProtein2Class = targetElementID + "_shapeProtein2";
                polyProtein2.transition().duration(duration_time).attr("points", pointsProteinTo2).each("end", function() {
                //polyProtein2.transition().duration(duration_time).attr("points", pointsProteinTo2);
                    d3.selectAll('.regionGroup').moveToFront();
                });
            }            
    }

    var shapeChrFrom1 = getChrShapeCoord(range1.selected_start, range1.selected_end, range1.selected_start, range1.selected_end, y(gene1_y - chr_offset) + exon_height, y(gene1_y - chr_offset) + exon_height);
    var shapeChrTo1 = getChrShapeCoord(range1.selected_start, range1.selected_end, range1.min, range1.max, y(gene1_y - chr_offset) + exon_height, y(gene1_y) - currentObj.coord_bar_offset);
    var shapeChrFrom2 = getChrShapeCoord(range2.selected_start, range2.selected_end, range2.selected_start, range2.selected_end, y(gene2_y - chr_offset) + exon_height, y(gene2_y - chr_offset) + exon_height);
    var shapeChrTo2 = getChrShapeCoord(range2.selected_start, range2.selected_end, range2.min, range2.max, y(gene2_y - chr_offset) + exon_height, y(gene2_y) - currentObj.coord_bar_offset);
        
    var shape1Class = targetElementID + "_chr_shape1";
    var shape2Class = targetElementID + "_chr_shape2";
    var poly1 = drawPolygon(shapeChrFrom1, "black", "." + currentObj.svgClasses, config.opacity);
    var poly2 = drawPolygon(shapeChrFrom2, "black", "." + currentObj.svgClasses, config.opacity);
    var pointsTo1 = mapPoints(shapeChrTo1);
    var pointsTo2 = mapPoints(shapeChrTo2);
    poly1.transition().duration(duration_time).attr("points", pointsTo1).attr('pointer-events', 'none').each("end", function() {                
            gene1_points = currentObj.drawRegions(svg, "."+ currentObj.svgClasses, gene1.exons, gene1.strand, gene1.chr, null, gene1.name, null, gene1.junction, x_pos, gene1_y, 0, 1, exon_height, coordRange1, coord_margin, null, true, false, maxExp1);
            poly2.transition().duration(duration_time).attr("points", pointsTo2).attr('pointer-events', 'none').each("end", function() {
                gene2_points = currentObj.drawRegions(svg, "."+ currentObj.svgClasses, gene2.exons, gene2.strand, gene2.chr, null, gene2.name, null, gene2.junction, x_pos, gene2_y, 0, 1, exon_height, coordRange2, coord_margin, null, true, false, maxExp2);
                fuse_points = currentObj.drawRegions(svg, "."+ currentObj.svgClasses, fusedGene.exons, fusedGene.strand, fusedGene.chr, null, fusedGene.name, null, fusedGene.junction, x_pos, fuse_y, 0, 1, exon_height, coordFusedRange, coord_margin, null, true, false, maxExp);
                var gene1_x = (gene1.strand == "+") ? {
                    "from": gene1_points[0],
                    "to": gene1_points[1]
                } : {
                    "from": gene1_points[2],
                    "to": gene1_points[1]
                };
                var gene2_x = (gene2.strand == "+") ? {
                    "from": gene2_points[1],
                    "to": gene2_points[2]
                } : {
                    "from": gene2_points[1],
                    "to": gene2_points[0]
                };

                var shapeFrom1 = getShapeCoord(gene1_x.from, gene1_x.to, gene1_x.from, gene1_x.to, y(gene1_y), y(gene1_y), exon_height);
                var shapeTo1 = getShapeCoord(gene1_x.from, gene1_x.to, fuse_points[0], fuse_points[1], y(gene1_y), y(fuse_y), exon_height);
                var shapeFrom2 = getShapeCoord(gene2_x.from, gene2_x.to, gene2_x.from, gene2_x.to, y(gene2_y), y(gene2_y), exon_height);
                var shapeTo2 = getShapeCoord(gene2_x.from, gene2_x.to, fuse_points[1], fuse_points[2], y(gene2_y), y(fuse_y), exon_height);

                
                var shape1Class = targetElementID + "_shape1";
                var shape2Class = targetElementID + "_shape2";
                var poly1 = drawPolygon(shapeFrom1, "blue", "."+ currentObj.svgClasses, config.opacity);
                var poly2 = drawPolygon(shapeFrom2, "red", "."+ currentObj.svgClasses, config.opacity);


                var pointsTo1 = mapPoints(shapeTo1);
                var pointsTo2 = mapPoints(shapeTo2);



                poly1.transition().duration(duration_time).attr("points", pointsTo1).attr('pointer-events', 'none').each("end", function() {
                    poly2.transition().duration(duration_time).attr('pointer-events', 'none').attr("points", pointsTo2);
                    d3.selectAll('.regionGroup').moveToFront();
                    d3.selectAll('.junction_arrow').moveToFront();
                });

                //console.log("downloadID:" + downloadID);
                //console.log("downloadID2:" + this.downloadID);
                d3.select("#" + downloadID)
                    .on("mouseover", function() {          
                                var html = d3.select("svg")
                                                .attr("title", "svg_title")
                                                .attr("version", 1.1)
                                                .attr("xmlns", "http://www.w3.org/2000/svg")
                                                .node().parentNode.innerHTML;
                                            
                                d3.select(this)
                                                .attr("href-lang", "image/svg+xml")                        
                                                .on('click', function() { 
                                                    var blob = new Blob([html], {type: "image/svg+xml"});
                                                    saveAs(blob, gene1.name + "-" + gene2.name + ".svg");
                                                  });
                        }
                    );
            });
    });            
}

GeneFusionPlot.prototype.drawTitle = function(svg, title, x_pos, y_pos, width, height) {
    //draw title text
    svg.append("rect")
        .attr("x", x_pos)
        .attr("y", y_pos-20)
        .attr("ry", 4)
        .attr("rx", 4)
        .attr("width", width)
        .attr("height", height)
        .style("fill", "orange")
        .attr("opacity", 0.3)
        .style("stroke", d3.rgb("orange").darker());
    svg.append("text")
        .attr("text-anchor", "middle")
        .attr("x", x_pos + width / 2)
        .attr("y", y_pos)
        .text(title)
        .attr("fill", "blue")
        .attr('font-size', 20);
}

GeneFusionPlot.prototype.getChromosomeRange = function(chromosome, cytobandData) {
    var maxCoord = 0;
    var midPoint = 0;
    cytobandData.forEach(function(d) {
        if (d.chromosome == chromosome) {
            var end_pos = parseInt(d.bp_stop);
            if (end_pos > maxCoord)
                maxCoord = end_pos;
            if (d.band.substring(0,1) == "p" && d.stain == "acen")
                midPoint = d.bp_stop;
        }
    });
    return [maxCoord, midPoint];
}

GeneFusionPlot.prototype.drawChromosome = function(svg, chromosome, cytobandData, maxCoord, x_pos, y_pos, chrRange, selectedRange, junction) {
    var colors = {
        "gneg": "#FFF",
        "gpos25": "#C8C8C8",
        "gpos33": "#BBB",
        "gpos50": "#999",
        "gpos66": "#888",
        "gpos75": "#777",
        "gpos100": "#444",
        "acen": "pink"
    };
    var coord_margin = this.coord_margin;
    var minCoord = 0;
    var y = this.y;
    var region_height = this.exon_height;
    var text_x_offset = 120;
    var text_y_offset = 10;

    var chrData = [];
    var first = true;
    cytobandData.forEach(function(d) {
        if (d.chromosome == chromosome) {
            var end_pos = parseInt(d.bp_stop);
            var left_rounded = false;
            var right_rounded = false;
            d.arm = d.band.substring(0,1);            

            if (first) {
                left_rounded = true;
                first = false;
            }
            if (end_pos == maxCoord)
                right_rounded = true;
            if (d.stain == "acen") {
                if (d.arm == "p")
                    right_rounded = true;
                if (d.arm == "q")
                    left_rounded = true;
            }
            var bandData = {
                "start": d.bp_start,
                "end": d.bp_stop,
                "color": colors[d.stain],
                "stain" : d.stain,
                "arm" : d.arm,
                "left_rounded": left_rounded,
                "right_rounded": right_rounded
            };
            chrData.push(bandData);
        }
    });

//console.log(chrData.length);
    chrData[chrData.length - 1].right_rounded = true;

    var x = d3.scale.linear()
        .domain([minCoord, maxCoord])
        .range([x_pos, this.width - x_pos / 2]);


    /*
    x: x-coordinate
    y: y-coordinate
    w: width
    h: height
    r: corner radius
    tl: top_left rounded?
    tr: top_right rounded?
    bl: bottom_left rounded?
    br: bottom_right rounded?
    */

    function rounded_rect(x, y, w, h, r, tl, bl, tr, br) {
        var retval;              
        retval = "M" + (x + r) + "," + y;
        retval += "h" + (w - 2 * r);
        if (tr) {            
            retval += "a" + r + "," + r + " 0 0 1 " + r + "," + r;
        } else {
            retval += "h" + r;
            retval += "v" + r;
        }
        retval += "v" + (h - 2 * r);
        if (br) {
            retval += "a" + r + "," + r + " 0 0 1 " + -r + "," + r;
        } else {
            retval += "v" + r;
            retval += "h" + -r;
        }
        retval += "h" + (2 * r - w);
        if (bl) {
            retval += "a" + r + "," + r + " 0 0 1 " + -r + "," + -r;
        } else {
            retval += "h" + -r;
            retval += "v" + -r;
        }
        retval += "v" + (2 * r - h);
        if (tl) {
            retval += "a" + r + "," + r + " 0 0 1 " + r + "," + -r;
        } else {
            retval += "v" + -r;
            retval += "h" + r;
        }
        retval += "z";
        return retval;
    }

    //draw band
    var band = d3.select("." + this.svgClasses).selectAll()
        .data(chrData)
        .enter()
        .append("g")
        .attr("class", "regionGroup");

    band.append("path")
        .attr("d", function(d) {
            if (d.stain == "acen") {
                if (d.arm == "p")
                    return "M" + x(d.start) + "," + y(y_pos) + "L" + x(d.end) + "," + (y(y_pos) + region_height / 2) + "L" + x(d.start) + "," + (y(y_pos) + region_height);
                if (d.arm == "q")
                    return "M" + x(d.end) + "," + y(y_pos) + "L" + x(d.start) + "," + (y(y_pos) + region_height / 2) + "L" + x(d.end) + "," + (y(y_pos) + region_height);
            } else
              return rounded_rect(x(d.start), y(y_pos), x(d.end) - x(d.start), region_height, 5, d.left_rounded, d.left_rounded, d.right_rounded, d.right_rounded)
        })
        .style("fill", function(d) {
            return d.color
        }).style("stroke", function(d) {
            return "black"
        });

    var selected = d3.select("." + this.svgClasses).selectAll()
        .data("dummy")
        .enter()        
        .append("rect")
        .attr("class", "regionGroup")
        .attr("x", x(selectedRange[0]))
        .attr("y", y(y_pos))
        .attr("width", x(selectedRange[1]) - x(selectedRange[0]))
        .attr("height", region_height)
        .style("stroke", "red");

    //draw chromosome text
    svg.append("text")
        .attr("class", "y-label")
        .attr("text-anchor", "right")
        .attr("x", x_pos - text_x_offset)
        .attr("y", y(y_pos) + text_y_offset)
        .text(chromosome)
        .attr("fill", "red")
        .attr('font-size', 16);

    //draw chromosome junction arrow
    var junctionPoint = d3.select("." + this.svgClasses).selectAll()
        .data(["dummy"]).enter()
        .append("line")
        .attr("x1", x(junction))
        .attr("y1", y(y_pos) - 12)
        .attr("x2", x(junction))
        .attr("y2", y(y_pos) - 2)
        .attr("stroke", "red")
        .attr("stroke-width", "3")
        .style("marker-end", "url(#arrow-down)");
    //band.attr("transform", "translate(-300,-150) rotate(90)");

    /*
        band.append("rect")
                .attr("x", function (r) {
                    return x(r.start);
                })
                .attr("y", y(y_pos) )
                .attr("ry", "0")
                .attr("rx", "0")
                .attr("width", function (r) {
                    return x(r.end) - x(r.start)
                })
                .attr("height", region_height)
                .style("fill", function (r) {
                    return r.color
                })
                .style("stroke", function (data) {
                    return d3.rgb(data.color).darker()
                });              
    */
    return {
        "min": x(minCoord),
        "max": x(maxCoord),
        "selected_start": x(selectedRange[0]),
        "selected_end": x(selectedRange[1])
    };
}

GeneFusionPlot.prototype.getCoordRange = function(exons, junction = -1) {

    var minCoord = 999999999;
    var maxCoord = 0;
    if (Array.isArray(exons)) {
        exons.forEach(function(exon) {        
            minCoord = Math.min(minCoord, exon.start_pos);
            maxCoord = Math.max(maxCoord, exon.end_pos);
            if (junction != -1) {
                minCoord = Math.min(minCoord, junction);
                maxCoord = Math.max(maxCoord, junction);
            }
        });
    }
    return [minCoord, maxCoord];
}

GeneFusionPlot.prototype.parseIntRegions = function(regions) {
    //console.log(JSON.stringify(regions));
    if (Array.isArray(regions)) {
        regions.forEach(function(region) {
            region.start_pos = parseInt(region.start_pos);
            region.end_pos = parseInt(region.end_pos);
            if (region.hasOwnProperty("value"))
                region.value = parseFloat(region.value);
        });
    }
    return regions;
}

GeneFusionPlot.prototype.getMaxExp = function(exons) {
    var maxExp = 0;
    exons.forEach(function(exon) {
        if (exon.hasOwnProperty("value"))
            maxExp = Math.max(exon.value, maxExp);
    });
    return maxExp;
}


GeneFusionPlot.prototype.drawRegions = function(svg, className, region_data, strand, chr, upper_text, left_text, right_text, junction, x_pos, y_pos, round_conor, region_opacity, region_height, coordRange, coord_margin, bg, draw_coord, show_text, maxExp=0, exon_pos_in_protein=null, domain_round_conor=0) {


    var buffer = this.buffer;
    //buffer = 200;
    var margin = 40;
    //var region_height = this.region_height;
    var bg_height = this.bg_height;
    var bg_offset = this.bg_offset;
    var region_offset = this.region_offset;
    var colors = this.colorMap;
    //var coord_margin = this.coord_margin;
    var coord_bar_offset = this.coord_bar_offset;
    var tick_offset = this.tick_offset;
    var text_offset = this.text_offset;
    var y = this.y;
    var tooltip = this.tooltip;

    var below = true;

    minCoord = coordRange[0] - coord_margin;
    maxCoord = coordRange[1] + coord_margin;

    //console.log(coordRange[0] + ' ' + minCoord + ' ' + coordRange[1] + ' ' + maxCoord + ' ' + coord_margin + ' ' + x_pos);
    //console.log(JSON.stringify(bg));
    var x_padding = 20;

    var x = d3.scale.linear()
        .domain([minCoord, maxCoord])
        .range([x_pos, this.width - x_pos / 2]);

    //define arrows
    svg.append("svg:defs")
        .append("svg:marker")
        .attr("id", "arrow-right")
        .attr("refX", 2)
        .attr("refY", 5)
        .attr("markerWidth", 10)
        .attr("markerHeight", 10)
        .attr("orient", "auto")
        .append("svg:path")
        .attr("d", "M2,2 L2,8 L8,5 L2,2")
        .attr("fill", "red");

    svg.append("svg:defs")
        .append("svg:marker")
        .attr("id", "arrow-left")
        .attr("refX", 8)
        .attr("refY", 5)
        .attr("markerWidth", 10)
        .attr("markerHeight", 10)
        .attr("orient", "auto")
        .append("svg:path")
        .attr("d", "M8,2 L8,8 L2,5 L8,2")
        .attr("fill", "red");

    svg.append("svg:defs")
        .append("svg:marker")
        .attr("id", "arrow-down")
        .attr("refX", 5)
        .attr("refY", 4)
        .attr("markerWidth", 10)
        .attr("markerHeight", 10)
        .attr("orient", "auto")
        .append("svg:path")
        .attr("d", "M2,2 L4,4 L2,6 L6,4 L2,2")
        .attr("fill", "red");

    getColor = this.colorScale;

    //var bg_offset = 0;
    //var region_offset = bg_offset-3;
    var text_offset = bg_offset + 5;
    if (below != true) {
        text_offset = bg_offset + 5;
    }

    bg_offset = region_height / 2;
    text_x_offset = 120;
    text_y_offset = 50;
    var num_ticks = (this.width - x_pos / 2) / 200;
    var coords = this.getCoordRange(region_data);  

    //draw background  
    if (bg != null) {
        var regionsBG = d3.select(className).selectAll()
            .data(bg).enter()
            .insert("g", ":first-child")
            .attr("class", "regionsBG")
            .append("rect")
            .attr("x", function(r) {
                return x(r.start_pos);
            })
            .attr("y", y(y_pos))
            .attr("ry", round_conor)
            .attr("rx", round_conor)
            .attr("width", function(r) {                
                var w = x(r.end_pos) - x(r.start_pos);
                if (w < 0)
                    w = 0;
                return w;
            })
            .attr("height", region_height)
            .style("fill", function(r) {
                return r.color;
            })            
            .style("stroke", function(r) {
                return d3.rgb(r.color).darker();
            });        
    } else {  //draw bg line
        var regionsBG = d3.select(className).selectAll()
            .data(["dummy"]).enter()
            .insert("g", ":first-child")
            .attr("class", "regionsBG")
            .append("line")
            .attr("x1", x(minCoord))
            .attr("y1", y(y_pos) + bg_offset)
            .attr("x2", x(maxCoord))
            .attr("y2", y(y_pos) + bg_offset)
            .attr("stroke", "red")
            .attr("stroke-width", "3");
        //append arrow
        if (strand == "+")
            regionsBG.style("marker-end", "url(#arrow-right)")
        else
            regionsBG.style("marker-start", "url(#arrow-left)");
    }

    if (draw_coord) {
        coordBar = d3.svg.axis().scale(x).orient("top");

        svg.append("svg:g")
            .attr("class", "x-axis chr_axis")
            .attr("transform", "translate(0," + (y(y_pos) - coord_bar_offset) + ")")
            .call(coordBar);
        
    }
          
    //draw gene name text

    if (upper_text != null) {
        svg.append("text")
            .attr("class", "y-label")
            .attr("text-anchor", "right")
            .attr("x", x(bg[bg.length-1].end_pos ) - 40)
            .attr("y", y(y_pos) - coord_bar_offset + text_y_offset)
            .text(upper_text)
            .attr('font-size', 12);
    }
    if (left_text != null) {
        svg.append("text")
            .attr("class", "y-label")
            .attr("text-anchor", "right")
            .attr("x", x_pos - text_x_offset)
            .attr("y", y(y_pos) + region_height / 3 * 2)
            //.attr("y", y(y_pos) - coord_bar_offset + text_y_offset)
            .text(left_text)
            .attr('font-size', 15);
    }
    if (right_text != null) {
        svg.append("text")
            .attr("class", "y-label")
            .attr("text-anchor", "right")
            .attr("x", x(bg[bg.length-1].end_pos) + 10)
            .attr("y", y(y_pos) + region_height / 3 * 2)
            .text(right_text)
            .attr('font-size', 15);
    }

    d3.select(".extent")
        .attr("y", y(0) + region_offset - 10);


    //draw exons/protein domains    


    var regions = regionsBG = d3.select(className).selectAll()
        .data(region_data)
        .enter()
        .append("g")
        .attr("class", "regionGroup");

    regions.append("rect")
        .attr("x", function(r) {
            return x(r.start_pos);
        })
        .attr("y", function(r) {            
            if (r.type)
                if (r.type.toLowerCase() == "utr3" || r.type.toLowerCase() == "utr5")
                    return y(y_pos) + region_height/4;
            return y(y_pos);
        })
        .attr("ry", domain_round_conor)
        .attr("rx", domain_round_conor)
        .attr("width", function(r) {
            var w = x(r.end_pos) - x(r.start_pos);
            //if (w < 0)
             //   console.log("neg: " + console.log(JSON.stringify(r)));
            return x(r.end_pos) - x(r.start_pos)
        })
        .attr("height", function(r) {
            if (r.type)
                if (r.type.toLowerCase() == "utr3" || r.type.toLowerCase() == "utr5")
                    return region_height/2;
            return region_height;
        })
        .style("fill", function(r) {
            return r.color
        })
        .attr("opacity", region_opacity)
        .attr('cursor', 'pointer')
        .on('mouseover', function(d){tooltip.show(d); if(tooltip.handle != undefined) clearTimeout( tooltip.handle );}) 
        //.on('mouseout', tooltip.hide)
        .on('mouseout', function(d){tooltip.handle = setTimeout( tooltip.hide, delay_sec );})
        .style("stroke", function(data) {
            return d3.rgb(data.color).darker()
        });

    if (exon_pos_in_protein != null) {
        var junctionPoint = d3.select(className).selectAll()
            .data(exon_pos_in_protein).enter()
            .append("line")
            .attr("class", "exon_tick")
            .attr("x1", function (d) { return x(d/3)})
            .attr("y1", y(y_pos))
            .attr("x2", function (d) { return x(d/3)})
            .attr("y2", y(y_pos)+region_height)
            .attr("stroke", "grey")
            .attr("stroke-width", "1");
    }
    var exp_gap = 3;
    var delay_sec = (show_text)? 2000 : 0;
    if (maxExp > 0) {
        var exp_height = 33;
        var exp_scale = d3.scale.linear()
            .domain([0, maxExp])
            .range([y(y_pos) - exp_gap, y(y_pos) - exp_height - exp_gap]);

        regions.append("rect")
        .attr("x", function(r) {
            return x(r.start_pos);
        })
        .attr("y", function(r) {
            return exp_scale(r.value);
        })
        .attr("width", function(r) {
            return x(r.end_pos) - x(r.start_pos)
        })
        .attr("height", function(r) {
            return exp_scale(0) - exp_scale(r.value);
        })
        .style("fill", function(r) {
            return "lightgrey";
        })
        .attr("opacity", region_opacity)
        .attr('cursor', 'pointer')
        .on('mouseover', tooltip.show) 
        .on('mouseout', tooltip.hide)
        //.on('mouseout', function(d){setTimeout( tooltip.hide, delay_sec );})
        .style("stroke", function(data) {
            return d3.rgb(data.color).darker()
        });

        xAxis = d3.svg.axis().scale(x).orient("bottom");

        svg.append("svg:g")
            .attr("class", "x-axis axis")
            .attr("transform", "translate(0," + (y(y_pos) - exp_gap) + ")")
            .call(xAxis).selectAll("text").remove();

        yAxis = d3.svg.axis().scale(exp_scale).orient("left").ticks(2);

        svg.append("svg:g")
            .attr("class", "y-axis axis")
            .attr("transform", "translate(" + (x(minCoord) - 3 )  + ",0)")
            .call(yAxis).append("text");

        svg.append("text")
        .attr("class", "y-label")
        .attr("text-anchor", "right")
        .attr("x", x(minCoord) - 40)
        .attr("y", exp_scale(maxExp) - exp_gap)
        .text("TPM")
        .attr('font-size', 11);
        
        //regions.append("line").attr("x1", x(minCoord) - 5).attr("x2", x(minCoord) - 5).attr("y1", y(y_pos) - exp_scale(maxExp) - exp_gap).attr("y2", y(y_pos) - exp_gap).attr("stroke", "red").attr("stroke-width", "3");;

    }
    //regions.attr('cursor', 'pointer').attr('pointer-events', 'all');        

    //show_text = false;

    if (show_text) {
        regions.append("text")
            .attr("text-anchor", "middle")
            .attr("fill", "black")
            .attr("opacity", 0.8)
            .attr("x", function (r) {
                return x(r.start_pos) + (x(r.end_pos) - x(r.start_pos)) / 2;
            })
            .attr("y", function(r) {return y(y_pos) + (region_height) / 2;} )
            .attr("dy", "0.35em")
            .style("font-size", "14px")
            //.style("text-decoration", "bold")
            .text(function (r) {
                if (getTextWidth(r.name, 14) < (x(r.end_pos) - x(r.start_pos)))
                    return r.name;
            });
        regions.attr('cursor', 'pointer');
        //.on('click', tooltip.show);
        //.on('mouseover', tooltip.show)
        //.on('mouseout', function(d){setTimeout(3000);} );
    }    
    //draw junction arrow
    if (junction < 1)
        junction = 1;
    var junctionPoint = d3.select(className).selectAll()
        .data(["dummy"]).enter()
        .append("line")
        .attr("class", "junction_arrow")
        .attr("x1", x(junction))
        .attr("y1", y(y_pos) - 12)
        .attr("x2", x(junction))
        .attr("y2", y(y_pos) - 2)
        .attr("stroke", "red")
        .attr("stroke-width", "3")
        .style("marker-end", "url(#arrow-down)");    
    
    /*
    var divLine = d3.select(className).selectAll()
        .data(["dummy"]).enter()
        .append("line")
        .attr("x1", x_pos - text_x_offset)
        .attr("y1", y(y_pos) + 30)
        .attr("x2", x(maxCoord))
        .attr("y2", y(y_pos) + 30)
        .style("stroke", "lightgrey");
    */
    //divLine.moveToFront();                
            
    if (show_text) 
        var coords = this.getCoordRange(bg);    
    return [x(coords[0]), x(junction), x(coords[1])];

};

GeneFusionPlot.prototype.calculateCoord = function(exons, ref_pos, junction, positive, left) {
        newExons = [];
        //var gene_start_pos = (positive)? exons[0].start_pos : exons[exons.length - 1].end_pos;
        //if (left) {

        if (left) {
            if (positive)
                gene_start_pos = Math.min(exons[0].start_pos, junction);
            else
                gene_start_pos = Math.max(exons[exons.length - 1].end_pos, junction);
        }
        else
            gene_start_pos = junction;
        
        /*} else
            if (positive)
                gene_start_pos = Math.min(exons[0].start_pos, junction);
            else
                gene_start_pos = Math.max(exons[exons.length - 1].end_pos, junction);
            if (junction < gene_start_pos)
                 gene_start_pos = junction;
        if (this.XOR(positive, left)) {
            gene_start_pos = ;
            if (gene_start_pos > junction)
                gene_start_pos = junction;
        }*/
        for (var i = 0; i < exons.length; i++) {            
            var exon = exons[i];
            var start = exon.start_pos; //start: start pos of exon in the fused part
            var end = exon.end_pos;     //end: end pos of exon in the fused part
            var adj_ref_pos = ref_pos;            
            //postive, left or negative, right
            if (!this.XOR(positive, left)) {
                //if exon is in the right part, skip
                if (exon.start_pos > junction) break;
                if (exon.end_pos >= junction) end = junction;
                adj_ref_pos = ref_pos;
            //negative, left or positive, right
            } else {
                //if exon is in the left part, skip
                if (exon.end_pos < junction) continue;
                if (exon.start_pos <= junction) start = junction;
            }

            if (!left) {
                //if (start > gene_start_pos)
                    adj_ref_pos = ref_pos - Math.abs(start - gene_start_pos);
                //else
                //    adj_ref_pos = ref_pos + Math.abs(junction - gene_start_pos);
            }
            //reverse coord if not positive
            if (!positive)
                end = [start, start = end][0];
            //var new_start_pos = Math.abs(start - gene_start_pos) + adj_ref_pos;
            //var new_end_pos = Math.abs(end - gene_start_pos) + adj_ref_pos;            
            var new_start_pos = Math.abs(start - gene_start_pos) + ref_pos;
            var new_end_pos = Math.abs(end - gene_start_pos) + ref_pos;

            //if (!left) {
            //    console.log("start:" + start + " end:" + end + " gene_start_pos:" + gene_start_pos + " ref_pos:"  + ref_pos);
            //    console.log("new_start_pos:" + new_start_pos + " new_end_pos:" + new_end_pos);
            //}
            if (new_start_pos > new_end_pos) {
                new_end_pos = [new_start_pos, new_start_pos = new_end_pos][0];
            }
            if (new_start_pos < 0 || new_end_pos < 0)
                continue;
            var newExon = {"start_pos" : new_start_pos, "end_pos" : new_end_pos, "type" : exon.type, "value" : exon.value, "hint" : {"Type" : exon.type, "Coordiante" : new_start_pos + ' - ' + new_end_pos, "Length" : (new_end_pos - new_start_pos), "Expression" : exon.value}, "color": exon.color};
            newExons.push(newExon);
        }
        return newExons;
}

GeneFusionPlot.prototype.getFusedGene = function (gene1, gene2) {
        var fused = [];
        var adj_coord1 = this.calculateCoord(gene1.exons, 0, gene1.junction, (gene1.strand == "+"), true);
        if (gene1.strand == "-")
            adj_coord1.reverse();
        //console.log("adj_coord1");
        //console.log(JSON.stringify(adj_coord1));
        if (adj_coord1.length == 0)
            fusion_point = 0;//(gene1.strand == "-") ? gene1.exons[0].start_pos : gene1.exons[gene1.exons.length-1].end_pos;
        else    
            fusion_point = adj_coord1[adj_coord1.length - 1].end_pos;
        var adj_coord2 = this.calculateCoord(gene2.exons, fusion_point, gene2.junction, (gene2.strand == "+"), false);
        //console.log("adj_coord2");
        //console.log(JSON.stringify(adj_coord2));
        if (gene2.strand == "-")
            adj_coord2.reverse();
        var fusedCoord = adj_coord1.concat(adj_coord2);
        //console.log("fused");
        //console.log(JSON.stringify(fusedCoord));
        return {
            "name": "FusedGene",
            "chr": gene1.chr + "+" + gene2.chr,
            "strand": "+",
            "junction": fusion_point,
            "exons": fusedCoord
        };

}

GeneFusionPlot.prototype.XOR = function(a, b) {
        return (a || b) && !(a && b);
}

GeneFusionPlot.prototype.drawAxes = function(svg) {

    var y = this.y;
    var x = this.x;

    xAxis = d3.svg.axis().scale(x).orient("bottom");

    svg.append("svg:g")
        .attr("class", "x-axis axis")
        .attr("transform", "translate(0," + (this.height - this.buffer) + ")")
        .call(xAxis);

    yAxis = d3.svg.axis().scale(y).orient("left");


    /*svg.append("svg:g")
      .attr("class", "y-axis axis")
      .attr("transform", "translate(" + (this.buffer * 1.5 )  + ",0)")
      .call(yAxis);
*/
    // appearance for x and y legend
    d3.selectAll(".axis path")
        .attr('fill', 'none');
    d3.selectAll(".domain")
        .attr('stroke', 'black')
        .attr('stroke-width', 1);

    svg.append("text")
        .attr("class", "y-label")
        .attr("text-anchor", "middle")
        .attr("transform", "translate(" + (this.buffer / 3) + "," + (this.height / 2) + "), rotate(-90)")
        .text(this.legends.y)
        .attr('font-size', 12);

    /*svg.append("text")
          .attr("class", "x-label")
          .attr("text-anchor", "middle")
          .attr("transform", "translate(" + (this.width / 2) + "," + (this.height - this.buffer / 3) + ")")
          .text(this.legends.x)
        .attr('font-weight', 'bold')
        .attr('font-size', 10);*/

};

function closeTooltip() {
    nodel = d3.select('.d3-tip');
    nodel.style({ opacity: 0, 'pointer-events': 'none' });
}