var builder = Project.buildRun("model.js").setStartTime(0).setStepCount(10);
builder.configureModel("model.js").setSketch({
  element: "Flow 1",
  start: 5,
  values: [5, 2, 4, 4]
});
var run = builder.run();