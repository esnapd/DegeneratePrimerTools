var prms, tre, clst, rt,lnk, lnkext;

HTMLWidgets.widget({

  name: 'phylogenetictree',

  type: 'output',

  initialize: function(el, width, height) {
    return {};
  },

  drawGraphic: function(el, params, width, height) {
  // remove existing children
  while (el.firstChild)
  el.removeChild(el.firstChild);

  console.log(params);
  prms = params;

  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  //Basic Margin Information
  var margin = { top: params.top, right: params.right,
                 bottom: params.bottom, left: params.left };
  width  = width - margin.left - margin.right;
  height = height - margin.top - margin.bottom;

  tre = params.tree;

  // add checkbox
  d3.select("#" + el.id).append("label")
    .attr("id", "show-length")
    .text("Show Branch length")
      .append("input")
      .attr("type", "checkbox");


  var outerRadius = params.outerradius,
      innerRadius = params.innerradius;

  var color = d3.scale.category10()
      .domain(["Bacteria", "Eukaryota", "Archaea"]);

  var cluster = d3.layout.cluster()
      .size([360, innerRadius])
      .children(function(d) { return d.branchset; })
      .value(function(d) { return 1; })
      .sort(function(a, b) { return (a.value - b.value) || d3.ascending(a.length, b.length); })
      .separation(function(a, b) { return 1; });

  clt = cluster;

  var svg = d3.select("#" + el.id).append("svg")
      .attr("width", outerRadius * 2)
      .attr("height", outerRadius * 2);

  var legend = svg.append("g")
      .attr("class", "legend")
    .selectAll("g")
      .data(color.domain())
    .enter().append("g")
      .attr("transform", function(d, i) { return "translate(" + (outerRadius * 2 - 10) + "," + (i * 20 + 10) + ")"; });

  legend.append("rect")
      .attr("x", -18)
      .attr("width", 18)
      .attr("height", 18)
      .style("fill", color);

  legend.append("text")
      .attr("x", -24)
      .attr("y", 9)
      .attr("dy", ".35em")
      .style("text-anchor", "end")
      .text(function(d) { return d; });

  var chart = svg.append("g")
      .attr("transform", "translate(" + outerRadius + "," + outerRadius + ")");

  var root = parseNewick(params.tree),
      nodes = cluster.nodes(root),
      links = cluster.links(nodes),
      input = d3.select("#show-length input").on("change", changed),
      timeout = setTimeout(function() { input.property("checked", true).each(changed); }, 2000);

  rt = root;

  setRadius(root, root.length = 0, innerRadius / maxLength(root));
  setColor(root);

  var linkExtension = chart.append("g")
      .attr("class", "link-extensions")
    .selectAll("path")
      .data(links.filter(function(d) { return !d.target.children; }))
    .enter().append("path")
      .each(function(d) { d.target.linkExtensionNode = this; })
      .attr("d", function(d) { return step(d.target.x, d.target.y, d.target.x, innerRadius); });

  var link = chart.append("g")
      .attr("class", "links")
    .selectAll("path")
      .data(links)
    .enter().append("path")
      .each(function(d) { d.target.linkNode = this; })
      .attr("d", function(d) { return step(d.source.x, d.source.y, d.target.x, d.target.y) })
      .style("stroke", function(d) { return d.target.color; });

  lnk = link;
  lnkext = linkExtension;

  chart.append("g")
      .attr("class", "labels")
    .selectAll("text")
      .data(nodes.filter(function(d) { return !d.children; }))
    .enter().append("text")
      .attr("dy", ".31em")
      .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + (innerRadius + 4) + ",0)" + (d.x < 180 ? "" : "rotate(180)"); })
      .style("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
      .text(function(d) { return d.name.replace(/_/g, " ");  })
      .on("mouseover", mouseovered(true))
      .on("mouseout", mouseovered(false));

  function changed() {
    clearTimeout(timeout);
    var checked = this.checked;
    d3.transition().duration(750).each(function() {
      linkExtension.transition().attr("d", function(d) { return step(d.target.x, checked ? d.target.radius : d.target.y, d.target.x, innerRadius); });
      link.transition().attr("d", function(d) { return step(d.source.x, checked ? d.source.radius : d.source.y, d.target.x, checked ? d.target.radius : d.target.y) });
    });
  }

  function mouseovered(active) {
    return function(d) {
      d3.select(this).classed("label--active", active);
      d3.select(d.linkExtensionNode).classed("link-extension--active", active).each(moveToFront);
      do d3.select(d.linkNode).classed("link--active", active).each(moveToFront); while (d = d.parent);
    };
  }

  function moveToFront() {
    this.parentNode.appendChild(this);
  }

  // Compute the maximum cumulative length of any node in the tree.
  function maxLength(d) {
    return d.length + (d.children ? d3.max(d.children, maxLength) : 0);
  }

  // Set the radius of each node by recursively summing and scaling the distance from the root.
  function setRadius(d, y0, k) {
    d.radius = (y0 += d.length) * k;
    if (d.children) d.children.forEach(function(d) { setRadius(d, y0, k); });
  }

  // Set the color of each node by recursively inheriting.
  function setColor(d) {
    d.color = color.domain().indexOf(d.name) >= 0 ? color(d.name) : d.parent ? d.parent.color : null;
    if (d.children) d.children.forEach(setColor);
  }

  // Like d3.svg.diagonal.radial, but with square corners.
  function step(startAngle, startRadius, endAngle, endRadius) {
    var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
        s0 = Math.sin(startAngle),
        c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
        s1 = Math.sin(endAngle);
    return "M" + startRadius * c0 + "," + startRadius * s0
        + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
        + "L" + endRadius * c1 + "," + endRadius * s1;
  }
  },

  renderValue: function(el, params, instance) {
    instance.params = params;
    this.drawGraphic(el, params, el.offsetWidth, el.offsetHeight);
  },

  resize: function(el, width, height, instance) {
  }

});
