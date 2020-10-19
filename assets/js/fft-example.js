// !preview r2d3 data=dt, dependencies = "d3-jetpack"
//

var { svg, margin, height, width } = d3.conventions({
  sel: svg,
  totalWidth: width,
  totalHeight: height,
  margin: { top: 30, right: 10, bottom: 40, left: 40 }
});

var shape = d3.scaleOrdinal(data.map(d => d.type), d3.symbols.map(s => d3.symbol().type(s)()));

var x = d3.scaleLinear()
      .domain([0, 22])
      .range([margin.left, width - margin.right]);
var y = d3.scaleLinear()
      .domain([0, 1])
      .range([height - margin.bottom, margin.top]);
var y2 = d3.scaleLinear()
      .domain(d3.extent(data, d => d.fft12_mod)).nice()
      .range([height - margin.bottom, margin.top]);

var phase = theta => d3.interpolateRainbow((theta + Math.PI) / (2 * Math.PI));

var xAxis = g => g
    .attr("transform", `translate(0,${height - margin.bottom})`)
    .call(d3.axisBottom(x).ticks(width / 80))
    .call(g => g.select(".domain").remove())
    .call(g => g.append("text")
        .attr("x", width)
        .attr("y", margin.bottom - 4)
        .attr("fill", "currentColor")
        .attr("text-anchor", "end")
        .text("Time"));
        
yAxis = g => g
    .attr("transform", `translate(${margin.left},0)`)
    .call(d3.axisLeft(y))
    .call(g => g.select(".domain").remove())
    .call(g => g.select("text.axis-label").remove())
    .call(g => g.append("text")
        .attr("x", -margin.left)
        .attr("y", 10)
        .attr("class", "axis-label")
        .attr("fill", "currentColor")
        .attr("text-anchor", "start")
        .text("Probability density function"));
        
yAxis2 = g => g
    .attr("transform", `translate(${margin.left},0)`)
    .call(d3.axisLeft(y2))
    .call(g => g.select(".domain").remove())
    .call(g => g.select("text.axis-label").remove())
    .call(g => g.append("text")
        .attr("x", -margin.left)
        .attr("y", 10)
        .attr("class", "axis-label")
        .attr("fill", "currentColor")
        .attr("text-anchor", "start")
        .text("|FFT(f(t))|"));

grid = g => g
    .attr("stroke", "currentColor")
    .attr("stroke-opacity", 0.1)
    .call(g => g.append("g")
      .selectAll("line")
      .data(x.ticks())
      .join("line")
        .attr("x1", d => 0.5 + x(d))
        .attr("x2", d => 0.5 + x(d))
        .attr("y1", margin.top)
        .attr("y2", height - margin.bottom))
    .call(g => g.append("g")
      .selectAll("line")
      .data(y.ticks())
      .join("line")
        .attr("y1", d => 0.5 + y(d))
        .attr("y2", d => 0.5 + y(d))
        .attr("x1", margin.left)
        .attr("x2", width - margin.right));

var points = svg.append("g")
    .attr("id", "plot-data")
    .attr("stroke-width", 1.5)
    .attr("font-family", "sans-serif")
    .attr("font-size", 10)
  .selectAll("path");
  
  
var nextButton = svg.append("g")
    .attr("transform", `translate(${width - margin.right - 20},0)`)
  .append("path")
    .attr("fill", "#808080")
    .attr("cursor", "pointer")
    .attr("d", "M0,7.5 L10,7.5 L10,0 L20,15 L10,30 L10,22.5 L0,22.5 Z");

var tick = function(step) {
  switch(parseInt(step)) {
    case 1:
      points = points
        .data(data.filter(d => d.type == "Model"), d => d.t)
        .join(
          enter => enter.append("path")
              .attr("d", d => shape(d.type))
              .attr("fill-opacity", 0)
              .attr("transform", d => `translate(${x(d.t)},${y(d.f1)})`)
              .call(path => path.transition()
                .delay((d, i) => d.t * 120)
                .attr("fill-opacity", 1)
              ),
          update => update.attr("fill", "#808080"),
          exit => exit.remove()
        );
      break;
    case 2:
      points = points.data(data, d => d.t)
        .join(
          enter => enter.append("path")
              .attr("transform", d => `translate(${x(d.t)},${y(d.f1)})`)
              .attr("fill-opacity", 0)
              .attr("d", d => shape(d.type))
              .call(path => path
                .transition()
                  .delay((d, i) => (d.t - 10) * 120)
                  .attr("fill-opacity", 1))
        );
      break;
    case 3:
      // data.forEach(d => console.log(JSON.stringify(d)));
      points = points.data(data, d => d.t)
        .join(
          enter => enter.append("path")
              .attr("d", d => shape(d.type))
              .attr("transform", d => `translate(${x(d.t)},${y2(d.fft_mod)})`)
              .attr("fill-opacity", 1),
          update => update
              .call(path => path.transition()
                .ease(d3.easeSin)
                .duration(1500)
                .attr("fill", d => phase(d.fft_arg))
                .attr("transform", d => `translate(${x(d.t)},${y2(d.fft_mod)})`)),
          exit => exit.remove()
        );
        svg.select("#y-axis").call(yAxis2);
  
      break;
    case 4:
      svg.select("#plot-data")
        .selectAll("rect")
        .data(data, d => d.i + 64)
        .join(
          enter => enter.append("rect")
            .attr("x", d => x(d.t) - 3)
            .attr("y", d => y2(d.fft2_mod) - 3)
            .attr("width", 6)
            .attr("height", 6)
            .attr("fill", d => phase(d.fft2_arg))
            .attr("fill-opacity", 0)
            .call(rect => rect.transition()
              .delay((d, i) => i * 10)
              .attr("fill-opacity", 1)),
          update => update,
          exit => exit.remove()
        );
      break;
    case 5:
      svg.select("#plot-data")
        .selectAll("path")
        .data(data, d => d.t)
        // Update
        .transition()
          .ease(d3.easeSin)
          .duration(1500)
          .attr("fill", d => phase(d.fft12_arg))
          .attr("transform", d => `translate(${x(d.t)},${y2(d.fft12_mod)})`);
      svg.select("#plot-data")
        .selectAll("rect")
        .remove();
      break;
    case 6:
      svg.select("#plot-data")
        .selectAll("path")
        .data(data.filter(d => d.type == "Model"), d => d.t)
        .join(
          enter => enter.append("path"),
          update => update.call(
            path => path.transition()
            .ease(d3.easeSin)
            .duration(1500)
            .attr("fill", "#000000")
            .attr("transform", d => `translate(${x(d.t)},${y(d.f12)})`)
          ),
          exit => exit.remove()
        );
      svg.select("#y-axis").call(yAxis);
  }
};
  
nextButton.attr("data-step", 1);
nextButton.on("click", () => {
  var step = nextButton.attr("data-step");
  tick(step);
  nextButton.attr("data-step", parseInt(step) + 1);
});

svg.append("g")
    .call(xAxis);
svg.append("g")
    .attr("id", "y-axis")
    .call(yAxis);
svg.append("g")
    .call(grid);
