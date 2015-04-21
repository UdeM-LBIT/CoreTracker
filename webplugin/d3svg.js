

function update_dot(data, svg, line, x, y, color, radius) {

    radius = radius || 4;

    data.forEach(function(d) {
        d.filtered = +d.filtered;
        d.global = +d.global;
    });

    svg.selectAll(".dot").remove()
    
    svg.selectAll(".dot")
    .data(data)
    .enter().append("circle")
    .attr("class", "dot")
    .attr("r", radius)
    .attr("cx", function(d) {
        return x(d.global);
    })
    .attr("cy", function(d) {
        return y(d.filtered);
    })
    .style("fill", function(d) {
        species = d.species.split('_');
        if (species.length >1) {            
            var grad = svg.append("defs")
                .append("linearGradient").attr("id", "grad"+d.species)
                .attr("x1", "0%").attr("x2", "0%").attr("y1", "100%").attr("y2", "0%");

            grad.append("stop").attr("offset", "50%").style("stop-color", color[species[0]]);
            grad.append("stop").attr("offset", "50%").style("stop-color", color[species[1]]);
            return "url(#grad"+d.species+")";
    }

        return color[d.species];
    })

    .append("title")
    .text(function(d) {
        return d.species;
    });

    svg.append('svg:path')
    .attr("d", line)
    .attr("class", "link")
    .attr("fill", "none");
}

function get_color(datalist, saturation) {

    saturation = saturation || 0.5;
    lightnessMin = 0.15;
    lightnessMax = 0.85;
    N = datalist.length;
    var base = 4;
    var lightnessDecay = 30;
    var tmp = "";
    var hue, lightness;

    color = {};
    for (var i = 0; i < N; i++) {

        tmp = i.toString(base).split("").reverse().join("");
        //console.log(tmp);
        hue = 360 * parseInt(tmp, base) / Math.pow(base, tmp.length);
        lightness = lightnessMin + (lightnessMax - lightnessMin) * (1 - Math.exp(-i / lightnessDecay));

        color[datalist[i]] = d3.hsl(hue, saturation +(0.2/i) ,lightness);

    }
    return color;
}


function aa_set_svg(svg, xAxis, yAxis, width, height, title) {

    svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(0," + height + ")")
    .call(xAxis)
    .append("text")
    .attr("class", "label")
    .attr("x", width)
    .attr("y", -6)
    .style("text-anchor", "end")
    .text("Global");
    svg.append("g")
    .attr("class", "y axis")
    .call(yAxis)
    .append("text")
    .attr("class", "label")
    .attr("transform", "rotate(-90)")
    .attr("y", 6)
    .attr("dy", ".71em")
    .style("text-anchor", "end")
    .text("Filtered");

    svg.append("text")
    .attr("x", (width / 2))
    .attr("y", height+35)
    .attr("text-anchor", "middle")
    .style("font-size", "13px")
    .style("text-decoration", "underline")
    .text(title);

}


function firstline(aafrequency, aause) {


    var margin = {
        top: 30,
        right: 30,
        bottom: 50,
        left: 40
    },
    width = 500 - margin.left - margin.right,
    height = 400 - margin.top - margin.bottom;

    var aasvg = d3.select("#aafreq").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


    var p_aasvg = d3.select("#aause").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    queue()
    .defer(d3.json, aafrequency)
    .defer(d3.json, aause)
    .await(readyall);

    function readyall(error, freqjson, usejson){

        if (error) return console.warn(error);

        var x1 = d3.scale.linear()
        .range([0, width]);
        var y1 = d3.scale.linear()
        .range([height, 0]);
        var xAxis1 = d3.svg.axis()
        .scale(x1)
        .orient("bottom");
        var yAxis1 = d3.svg.axis()
        .scale(y1)
        .orient("left");

        // set second axis
        var x2 = d3.scale.linear()
        .range([0, width]);
        var y2 = d3.scale.linear()
        .range([height, 0]);
        var xAxis2 = d3.svg.axis()
        .scale(x2)
        .orient("bottom");
        var yAxis2 = d3.svg.axis()
        .scale(y2)
        .orient("left");  


        var linefunc1 = d3.svg.line()
        .x(function(d) {
            return x1(d.x);
        })
        .y(function(d) {
            return y1(d.y);
        })
        .interpolate('linear');

        var linefunc2 = d3.svg.line()
        .x(function(d) {
            return x2(d.x);
        })
        .y(function(d) {
            return y2(d.y);
        })
        .interpolate('linear');       // set json data


        injson = freqjson['AA'];
        expected = freqjson['EXP'];
        max_val = +freqjson['MAX'] + 0.5;
        keys = d3.keys(usejson); // those keys are aa
        default_key = keys[0];
        species_list = []
        for (key in injson[default_key]){
            species_list.push(injson[default_key][key].species)
        }

        //setting colors
        //var bezInterpolator = chroma.interpolate.bezier(['lightyellow', 'pink',  'orange', 'red', 'blue']);

        // var chromacolor = chroma.scale(bezInterpolator).correctLightness(true).domain([0, species_list.length]);
        // dynamic bit...
        //var color = d3.scale.category10();
        //var color = chromacolor //d3.scale.ordinal().domain(species_list).range(chromacolor);
        color = get_color(species_list)
        
        var lineData1 = [{
            x: 0,
            y: 0
        }, {
            x: max_val,
            y: max_val
        }];


        var lineData2 = [{
            x: 0,
            y: 0
        }, {
            x: 1.1,
            y: 1.1
        }];

        var expectedData = [{
            x: 1,
            y: 0
        }, {
            x: 1,
            y: max_val
        }];

        x1.domain([0, max_val]).nice();
        y1.domain([0, max_val]).nice();

        x2.domain([0, 1.1]).nice();
        y2.domain([0, 1.1]).nice();

        // configure svg
        aa_set_svg(aasvg, xAxis1, yAxis1, width, height, "Amino acid usage (Observed-Freq / Expected-Freq from alignment)");
        aa_set_svg(p_aasvg, xAxis2, yAxis2, width, height, "Amino acid usage for each pair");

        aasvg.append('svg:path')
        .attr("d", linefunc1(expectedData))
        .attr("class", "link")
        .style("stroke-dasharray", ("2, 2"))
        .style("stroke", "#333")
        .attr("fill", "none")
        .append("title")
        .text("Expected frequence");

        update_dot(injson[default_key], aasvg, linefunc1(lineData1), x1, y1, color);
        update_dot(usejson[default_key], p_aasvg, linefunc2(lineData2), x2, y2, color);


        $("#aa").autocomplete({
            source: keys,
            autoFocus: true,
            select: function(a, b) {
                $(this).val(b.item.value);
                update_dot(injson[b.item.value], aasvg, linefunc1(lineData1), x1, y1, color);
                update_dot(usejson[b.item.value], p_aasvg, linefunc2(lineData2), x2, y2, color);
            }
        }).val(default_key).data('autocomplete');

    }
}

function secondline(similarity){



    var margin = {
        top: 30,
        right: 30,
        bottom: 50,
        left: 40
    },
    width = 600 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

      var idsvg = d3.select("#identity").append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      // sequence identity
      d3.json(similarity, function(error, json) {
        var x = d3.scale.linear()
          .range([0, width]);
        var y = d3.scale.linear()
          .range([height, 0]);
        var xAxis = d3.svg.axis()
          .scale(x)
          .orient("bottom");
        var yAxis = d3.svg.axis()
          .scale(y)
          .orient("left");
        var linefunc = d3.svg.line()
          .x(function(d) {
            return x(d.x);
          })
          .y(function(d) {
            return y(d.y);
          })
          .interpolate('linear');
        if (error) return console.warn(error);

        keys = d3.keys(json)
        data = json[keys[0]]
        var lineData = [{
          x: 0,
          y: 0
        }, {
          x: 100,
          y: 100
        }]
        x.domain([0, 110]).nice();
        y.domain([0, 110]).nice();

        species_list = []
        for (key in data){
            species_list.push(data[key].species)
        }
        color = get_color(species_list)

        aa_set_svg(idsvg, xAxis, yAxis, width, height, "Pairwise identity in global and filtered alignment (%)");
    
          
        update_dot(json[keys[0]], idsvg, linefunc(lineData), x, y, color);

        $("#gene").autocomplete({
          source: keys,
          autoFocus: true,
          select: function(a, b) {
            $(this).val(b.item.value);
            update_dot(json[b.item.value], idsvg, linefunc(lineData), x, y, color);

          }
        }).val(keys[0]).data('autocomplete');
      });
}