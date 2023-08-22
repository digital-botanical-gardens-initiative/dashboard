// Wait until amCharts 5 library is fully loaded
am5.ready(function() {
      
    // Create root elements for charts
    // Create a root element for the first chart (chartdiv will be its container)
    // https://www.amcharts.com/docs/v5/getting-started/#Root_element
    var root = am5.Root.new("chartdiv");
    // Create a root element for the second chart (chartdiv2 will be its container)
    var root2 = am5.Root.new("chartdiv2");
      
    // Set themes
    // https://www.amcharts.com/docs/v5/concepts/themes/
    root.setThemes([
      am5themes_Animated.new(root)
    ]);
    root2.setThemes([
      am5themes_Animated.new(root2)
    ]);

    // Create radar charts
    // https://www.amcharts.com/docs/v5/charts/radar-chart/
    var chart = root.container.children.push(am5radar.RadarChart.new(root, {
      panX: false,
      panY: false,
      wheelX: "panX",
      wheelY: "zoomX",
      innerRadius: am5.percent(20),
      startAngle: -90,
      endAngle: 180
    }));

    var chart2 = root2.container.children.push(am5radar.RadarChart.new(root2, {
      panX: false,
      panY: false,
      wheelX: "panX",
      wheelY: "zoomX",
      innerRadius: am5.percent(20),
      startAngle: -90,
      endAngle: 180
    }));


    // Define the data for the charts
    // Data for the first chart (relating to all taxons)
    var data = [{
      category: "Species",
      value: CollectedSpecies/TotalSpecies * 100,
      full: 100,
      columnSettings: {
        fill: chart.get("colors").getIndex(0)
      }
    }, {
      category: "Genus",
      value: CollectedGenus/TotalGenus * 100,
      full: 100,
      columnSettings: {
        fill: chart.get("colors").getIndex(1)
      }
    }, {
      category: "Family",
      value: CollectedFamily/TotalFamily *100,
      full: 100,
      columnSettings: {
        fill: chart.get("colors").getIndex(2)
      }
    }, {
      category: "Order",
      value: CollectedOrder/TotalOrder * 100,
      full: 100,
      columnSettings: {
        fill: chart.get("colors").getIndex(3)
      }
    },
    {
      category: "Phylum",
      value: CollectedKingdom/TotalKingdom * 100,
      full: 100,
      columnSettings: {
        fill: chart.get("colors").getIndex(5)
      }
    }];

    // Data for the second chart (specific to plants)
    var data2 = [{
      category: "Species",
      value: CollectedSpecies/TotalPlantSpecies * 100,
      full: 100,
      columnSettings: {
        fill: chart2.get("colors").getIndex(0)
      }
    }, {
      category: "Genus",
      value: CollectedGenus/TotalPlantGenus * 100,
      full: 100,
      columnSettings: {
        fill: chart2.get("colors").getIndex(1)
      }
    }, {
      category: "Family",
      value: CollectedFamily/TotalPlantFamily *100,
      full: 100,
      columnSettings: {
        fill: chart2.get("colors").getIndex(2)
      }
    }, {
      category: "Order",
      value: CollectedOrder/TotalPlantOrder * 100,
      full: 100,
      columnSettings: {
        fill: chart2.get("colors").getIndex(3)
      }
    }
    ];

    // Add chart titles
    chart.children.unshift(am5.Label.new(root, {
      text: "Progress of Sampling by All Taxon Rank",
      fontSize: 25,
      fontWeight: "500",
      textAlign: "center",
      x: am5.percent(50),
      centerX: am5.percent(50),
      y: 0,
      paddingTop: 0,
      paddingBottom: 0
    }));
    chart2.children.unshift(am5.Label.new(root2, {
      text: "Progress of Sampling by Plant Taxon Rank",
      fontSize: 25,
      fontWeight: "500",
      textAlign: "center",
      x: am5.percent(50),
      centerX: am5.percent(50),
      y: 0,
      paddingTop: 0,
      paddingBottom: 0
    }));

    // Configure cursor for zooming    
    // https://www.amcharts.com/docs/v5/charts/radar-chart/#Cursor
    var cursor = chart.set("cursor", am5radar.RadarCursor.new(root, {
      behavior: "zoomX"
    }));
    cursor.lineY.set("visible", false);
    var cursor2 = chart2.set("cursor", am5radar.RadarCursor.new(root2, {
      behavior: "zoomX"
    }));
    cursor2.lineY.set("visible", false);

    // Create and set up axes
    // https://www.amcharts.com/docs/v5/charts/radar-chart/#Adding_axes
    var xRenderer = am5radar.AxisRendererCircular.new(root, {
      //minGridDistance: 50
    });
    var xRenderer2 = am5radar.AxisRendererCircular.new(root2, {
      //minGridDistance: 50
    });

    // Modify appearance of x-axis labels and grid    
    xRenderer.labels.template.setAll({
      radius: 10
    });
    xRenderer2.labels.template.setAll({
      radius: 10
    });
    xRenderer.grid.template.setAll({
      forceHidden: true  
    });
    xRenderer2.grid.template.setAll({
      forceHidden: true  
    });

    // Create x-axes for the charts
    var xAxis = chart.xAxes.push(am5xy.ValueAxis.new(root, {
      renderer: xRenderer,
      min: 0,
      max: 100,
      strictMinMax: true,
      numberFormat: "#'%'",
      tooltip: am5.Tooltip.new(root, {})
    }));
    var xAxis2 = chart2.xAxes.push(am5xy.ValueAxis.new(root2, {
      renderer: xRenderer2,
      min: 0,
      max: 100,
      strictMinMax: true,
      numberFormat: "#'%'",
      tooltip: am5.Tooltip.new(root2, {})
    }));

    // Set up y-axis renderers
    var yRenderer = am5radar.AxisRendererRadial.new(root, {
      minGridDistance: 20
    });
    var yRenderer2 = am5radar.AxisRendererRadial.new(root2, {
      minGridDistance: 20
    });

    // Modify appearance of y-axis labels and grid
    yRenderer.labels.template.setAll({
      centerX: am5.p100,
      fontWeight: "500",
      fontSize: 18,
      templateField: "columnSettings"
    });
    yRenderer2.labels.template.setAll({
      centerX: am5.p100,
      fontWeight: "500",
      fontSize: 18,
      templateField: "columnSettings"
    });
    yRenderer.grid.template.setAll({
      forceHidden: true
    });
    yRenderer2.grid.template.setAll({
      forceHidden: true
    });

    // Create y-axes for the charts
    var yAxis = chart.yAxes.push(am5xy.CategoryAxis.new(root, {
      categoryField: "category",
      renderer: yRenderer
    }));
    var yAxis2 = chart2.yAxes.push(am5xy.CategoryAxis.new(root2, {
      categoryField: "category",
      renderer: yRenderer2
    }));

    // Assign data to y-axes
    yAxis.data.setAll(data);
    yAxis2.data.setAll(data2);


    // Create series for both charts
    // Series for background (100% bars)
    // https://www.amcharts.com/docs/v5/charts/radar-chart/#Adding_series
    var series1 = chart.series.push(am5radar.RadarColumnSeries.new(root, {
      xAxis: xAxis,
      yAxis: yAxis,
      clustered: false,
      valueXField: "full",
      categoryYField: "category",
      fill: root.interfaceColors.get("alternativeBackground")
    }));
    var series1_2 = chart2.series.push(am5radar.RadarColumnSeries.new(root2, {
      xAxis: xAxis2,
      yAxis: yAxis2,
      clustered: false,
      valueXField: "full",
      categoryYField: "category",
      fill: root2.interfaceColors.get("alternativeBackground")
    }));

    // Configure appearance of the background series columns
    series1.columns.template.setAll({
      width: am5.p100,
      fillOpacity: 0.08,
      strokeOpacity: 0,
      cornerRadius: 20
    });
    series1_2.columns.template.setAll({
      width: am5.p100,
      fillOpacity: 0.08,
      strokeOpacity: 0,
      cornerRadius: 20
    });

    // Assign data to background series
    series1.data.setAll(data);
    series1_2.data.setAll(data2);

    // Series for data values
    var series2 = chart.series.push(am5radar.RadarColumnSeries.new(root, {
      xAxis: xAxis,
      yAxis: yAxis,
      clustered: false,
      valueXField: "value",
      categoryYField: "category"
    }));
    var series2_2 = chart2.series.push(am5radar.RadarColumnSeries.new(root2, {
      xAxis: xAxis2,
      yAxis: yAxis2,
      clustered: false,
      valueXField: "value",
      categoryYField: "category"
    }));

    // Configure appearance of the data series columns
    series2.columns.template.setAll({
      width: am5.p100,
      strokeOpacity: 0,
      tooltipText: "{category}: {valueX}%",
      cornerRadius: 20,
      templateField: "columnSettings"
    });
    series2_2.columns.template.setAll({
      width: am5.p100,
      strokeOpacity: 0,
      tooltipText: "{category}: {valueX}%",
      cornerRadius: 20,
      templateField: "columnSettings"
    });


    // Assign data to the data series
    series2.data.setAll(data);
    series2_2.data.setAll(data2);


    // Apply initial animations to the charts
    // https://www.amcharts.com/docs/v5/concepts/animations/#Initial_animation
    series1.appear(1000);
    series2.appear(1000);
    chart.appear(1000, 100);

    series1_2.appear(1000);
    series2_2.appear(1000);
    chart2.appear(1000, 100);

}); // end am5.ready()