var EPSILON = Math.pow(2, -52);

function __Mash() {
  this.n = 0xefc8249d;
}
__Mash.prototype.mash = function (data) {
  data = data.toString();
  var n = this.n;
  for (var i = 0; i < data.length; i++) {
    n += data.charCodeAt(i);
    var h = 0.02519603282416938 * n;
    n = h >>> 0;
    h -= n;
    h *= n;
    n = h >>> 0;
    h -= n;
    n += h * 0x100000000; // 2^32
  }
  this.n = n;
  return (n >>> 0) * 2.3283064365386963e-10; // 2^-32
};

function __Alea(seed) {
  var mash = new __Mash();

  // Apply the seeding algorithm from Baagoe.
  var s0 = mash.mash(" ");
  var s1 = mash.mash(" ");
  var s2 = mash.mash(" ");
  s0 -= mash.mash(seed);
  if (s0 < 0) {
    s0 += 1;
  }
  s1 -= mash.mash(seed);
  if (s1 < 0) {
    s1 += 1;
  }
  s2 -= mash.mash(seed);
  if (s2 < 0) {
    s2 += 1;
  }

  this.c = 1;
  this.s0 = s0;
  this.s1 = s1;
  this.s2 = s2;
}
__Alea.prototype.sample = function () {
  var t = 2091639 * this.s0 + this.c * 2.3283064365386963e-10; // 2^-32
  this.s0 = this.s1;
  this.s1 = this.s2;
  return (this.s2 = t - (this.c = t | 0));
};

function __Model(tickCount, elementCount, state) {
  this.bufSize = tickCount;
  this.outputs = [];
  for (var i = 0; i < elementCount; ++i) {
    this.outputs.push(new Float64Array(this.bufSize));
  }
  this.state = state;
}
__Model.prototype.getSketchValues = function (id) {
  var data = new Float64Array(this.bufSize);

  for (var i = 0; i < this.bufSize; ++i) {
    data[i] = this.computeSketch(this.state, null, i, id, NaN);
  }

  return data;
};
__Model.prototype.setInputByOffset = function (offset, value) {
  this.state.values[offset] = value;
  this.state.calc[offset] = 0;
};
__Model.prototype.getLookupValue = function (id, x) {
  return __lookup(this.state.lk[id], x);
};
__Model.prototype.getLookupValues = function (id, min, max, steps) {
  var data = new Float64Array(steps + 1);
  var lookup = this.state.lk[id];
  let dx = (max - min) / steps;
  for (var i = 0; i <= steps; ++i) {
    data[i] = __lookup(lookup, min + dx * i);
  }
  return data;
};
__Model.prototype.setDenseLookupById = function (id, table, min, max, steps) {
  this.state.lk[id] = __createDenseLookupFromSparse(table, min, max, steps);
};
__Model.prototype.setLookupById = function (id, table) {
  this.state.lk[id] = __createSparseLookup(table);
};
__Model.prototype.step = function () {
  for (var i = 0; i < this.state.res; ++i) {
    this.state.t += 1;
    this.calc();
  }
};
__Model.prototype.runToEnd = function () {
  if (this.state.t === 0) {
    this.calc();
  }

  while (this.state.t < this.bufSize) {
    this.step();
  }
};

function __delay(tick, resolution, closure, offset, defaultValue, legacy) {
  if (legacy) {
    var calcTick = tick - Math.round(offset * resolution);
    if (calcTick < 0) {
      if (defaultValue === null) {
        calcTick = 0;
      } else {
        return defaultValue;
      }
    }
    return closure(calcTick);
  } else {
    if (offset <= 0) {
      return defaultValue;
    } else {
      var calcTick = tick - Math.round(offset * resolution);
      return calcTick < 0 ? defaultValue : closure(calcTick);
    }
  }
}

function __denseLookup(table, min, max, step, x) {
  if (x <= min) {
    return table[0];
  } else if (x >= max) {
    return table[table.length - 1];
  } else {
    var delta = (x - min) / step;
    var left = Math.floor(delta);
    var right = Math.ceil(delta);
    if (left === right) {
      return table[left];
    } else {
      var frac = delta - left;
      return table[left] + frac * (table[right] - table[left]);
    }
  }
}

function __buildLookups(lookups, modelId) {
  var lookupMap = {};

  if (!lookups) {
    return lookupMap;
  }

  var modelLookups = lookups[modelId];

  if (modelLookups != null) {
    for (var key in modelLookups) {
      lookupMap[key] = __createSparseLookup(modelLookups[key]);
    }
  }

  return lookupMap;
}

function __buildSketches(sketches, modelId, containingModelId, submodelId) {
  var sketchMap = {};

  var modelSketches = sketches[modelId];
  if (modelSketches != null && modelSketches.elements != null) {
    for (var key in modelSketches.elements) {
      sketchMap[key] = modelSketches.elements[key];
    }
  }

  if (containingModelId != null && submodelId != null) {
    var containingModelSketches = sketches[containingModelId];
    if (
      containingModelSketches != null &&
      containingModelSketches.submodels != null
    ) {
      var overlaySketches = containingModelSketches.submodels[submodelId];
      if (overlaySketches != null) {
        for (var key in overlaySketches) {
          sketchMap[key] = overlaySketches[key];
        }
      }
    }
  }

  return sketchMap;
}

function __sparseLookup(xs, ys, x) {
  if (x <= xs[0]) {
    return ys[0];
  } else if (x >= xs[xs.length - 1]) {
    return ys[ys.length - 1];
  } else {
    var leftX = xs[0];
    for (var i = 1; i < xs.length; ++i) {
      var rightX = xs[i];
      if (x == rightX) {
        return ys[i];
      } else if (x < rightX) {
        var delta = ys[i] - ys[i - 1];
        var frac = (x - leftX) / (rightX - leftX);
        return ys[i - 1] + frac * delta;
      } else {
        leftX = rightX;
      }
    }
  }
}

function __createSparseLookup(table) {
  var vs = [];
  for (var key in table) {
    vs.push({ x: parseFloat(key), y: +table[key] });
  }
  vs.sort(function (a, b) {
    return a.x - b.x;
  });
  if (vs.length === 0) {
    vs.push({ x: 0, y: 0 });
  }
  var xs = new Float64Array(vs.length);
  var ys = new Float64Array(vs.length);
  vs.forEach(function (p, i) {
    xs[i] = p.x;
    ys[i] = p.y;
  });
  return { sparse: true, xs: xs, ys: ys };
}

function __createDenseLookupFromSparse(table, min, max, steps) {
  if (max <= min) {
    throw new Error("Lookup max must be strictly less than min.");
  }
  if (steps < 0) {
    throw new Error("Lookup steps must be greater than zero.");
  }
  var sparse = __createSparseLookup(table);
  let values = new Float64Array(steps + 1);
  let dx = (max - min) / steps;
  for (var i = 0; i <= steps; ++i) {
    values[i] = __sparseLookup(sparse.xs, sparse.ys, min + i * dx);
  }
  return { sparse: false, values: values, step: dx, min: min, max: max };
}

function __lookup(table, x) {
  if (table === undefined) {
    return 0;
  } else if (table.sparse) {
    return __sparseLookup(table.xs, table.ys, x);
  } else {
    return __denseLookup(table.values, table.min, table.max, table.step, x);
  }
}

function __sketchLookup(sketch, t, blankValue) {
  if (sketch === undefined) {
    return blankValue;
  }

  var index = t - sketch.start;
  if (index < 0 || index > sketch.values.length - 1) {
    return blankValue;
  }

  var leftIdx = Math.floor(index);
  var left = sketch.values[leftIdx];
  if (isNaN(left)) {
    return blankValue;
  }
  if (leftIdx === index) {
    return left;
  }
  var rightIdx = Math.ceil(index);
  if (rightIdx >= sketch.values.length) {
    return blankValue;
  }
  var right = sketch.values[rightIdx];
  if (isNaN(right)) {
    return blankValue;
  }
  return left + (index - leftIdx) * (right - left);
}

function __getTableColumn(tables, tableId, columnId) {
  if (!tables) return;
  var table = tables[tableId];
  if (!table) return;
  return table[columnId];
}

var TOLERANCE = 3 * EPSILON;

function __eq(a, b) {
  return a === b || Math.abs(a - b) < TOLERANCE;
}

function __leq(a, b) {
  return a < b || __eq(a, b);
}

function __geq(a, b) {
  return a > b || __eq(a, b);
}

function RunBuilder(project, modelName) {
  this.project = project;
  this.model = project.byName[modelName];
  this.modelConfigs = {};
  this.startTime = 0;
  this.stepCount = 100;
  this.resolution = 1;
}

RunBuilder.prototype.configureModel = function (modelName) {
  if (!this.modelConfigs[modelName]) {
    this.modelConfigs[modelName] = new ModelConfiguration(this, modelName);
  }

  return this.modelConfigs[modelName];
};

RunBuilder.prototype.setStartTime = function (startTime) {
  this.startTime = startTime;
  return this;
};

RunBuilder.prototype.setStepCount = function (stepCount) {
  this.stepCount = stepCount;
  return this;
};

RunBuilder.prototype.setResolution = function (resolution) {
  this.resolution = resolution;
  return this;
};

RunBuilder.prototype.start = function () {
  var sketches = {};
  var lookups = {};
  var tables = {};

  for (var name in this.project.byName) {
    var modelconstructor = this.project.byName[name];

    var config = this.configureModel(name);
    var modelId = null;

    for (var id in this.project.byId) {
      if (this.project.byId[id] === modelconstructor) {
        modelId = +id;
        break;
      }
    }

    sketches[modelId] = config.sketches;
    lookups[modelId] = config.lookups;
    tables[modelId] = config.tables;
  }

  return new this.model(
    this.startTime,
    this.resolution,
    this.stepCount,
    sketches,
    lookups,
    tables,
  );
};

RunBuilder.prototype.run = function () {
  var run = this.start();
  run.runToEnd();
  return run;
};

function ModelConfiguration(builder, modelName) {
  this.builder = builder;
  this.project = builder.project;
  var model = (this.model = this.project.byName[modelName]);
  this.sketches = { elements: {}, submodels: {} };
  this.lookups = {};
  var tables = (this.tables = {});

  for (var elementId in model.elementMetadata) {
    var element = model.elementMetadata[elementId];

    if (element.type === "DataTable") {
      tables[element.id] = {};
    }
  }

  this.rootPanel = model.getPanels().find(function (panel) {
    return panel.root === true;
  }).name;
}

function findElement(model, panel, element) {
  if (typeof element === "number") {
    return model.elementMetadata[element];
  } else {
    var id = model.getId(panel, element);
    if (id === undefined) {
      throw new Error('Could not find element "' + element + '".');
    }
    return model.elementMetadata[id];
  }
}

ModelConfiguration.prototype.setSketch = function (config) {
  var element = findElement(
    this.model,
    config.panel || this.rootPanel,
    config.element,
  );
  this.sketches.elements[element.id] = {
    start: config.start === undefined ? this.builder.startTime : config.start,
    values: config.values,
  };

  return this;
};

ModelConfiguration.prototype.setSubmodelSketch = function (config) {
  var element = findElement(
    this.model,
    config.panel || this.rootPanel,
    config.element,
  );

  if (this.sketches.submodels[element.id] == null) {
    this.sketches.submodels[element.id] = {};
  }

  var sketchMap = this.sketches.submodels[element.id];

  var targetModelBuilder = this.builder.configureModel(element.model);

  var targetElement = findElement(
    targetModelBuilder.model,
    config.targetPanel || targetModelBuilder.rootPanel,
    config.targetElement,
  );

  sketchMap[targetElement.id] = {
    start: config.start === undefined ? this.builder.startTime : config.start,
    values: config.values,
  };

  return this;
};

ModelConfiguration.prototype.setLookup = function (config) {
  var element = findElement(
    this.model,
    config.panel || this.rootPanel,
    config.element,
  );
  this.lookups[element.id] = config.table;

  return this;
};

ModelConfiguration.prototype.setDataset = function (config) {
  var table = {};

  var element = findElement(
    this.model,
    config.panel || this.rootPanel,
    config.element,
  );

  for (var columnName in config.columns) {
    var columnData = config.columns[columnName];

    var column = element.columns.find(function (col) {
      return col.name === columnName;
    });

    table[column.id] = {
      start:
        columnData.start === undefined
          ? this.builder.startTime
          : columnData.start,
      values: columnData.values,
    };
  }

  this.tables[element.id] = table;

  return this;
};
function _abs(x) {return Math.abs(x)};function _ceil(x) {return Math.ceil(x)};function _exp(x) {return Math.exp(x)};function _floor(x) {return Math.floor(x)};function _sqrt(x) {return Math.sqrt(x)};function _div(x, y, z) {return y == 0 ? z : x / y};function _frac(x) {return x - _int(x)};function _int(x) {return x < 0 ? Math.ceil(x) : Math.floor(x)};function _ln(x) {return x <= 0 ? 0 : Math.log(x)};function _log(x, base) {return x <= 0 ? 0 : (Math.log(x) / Math.log(base))};function _max(xs) {return xs.length === 0 ? 0 : Math.max.apply(Math, xs)};function _min(xs) {return xs.length === 0 ? 0 : Math.min.apply(Math, xs)};function _bound(x, min, max) {return x < min ? min : x > max ? max : x};function _power(x, e) {return Math.pow(x, e);};function _root(x, n) {if (x <= 0) {
  return 0;
} else {
  return Math.pow(Math.E, Math.log(x) / n);
}
};function _sign(x) {return x < 0 ? -1 : x > 0 ? 1 : 0};function _round2(x, places) {var y = places === 0 ? x : Math.pow(10, places) * x;
var i = Math.floor(y);
var f = y - i;

if (f === 0.5) {
  y = i % 2 === 0 ? i : i + 1;
} else {
  y = Math.round(y);
}

if (places === 0) {
  return y;
} else {
  return y / Math.pow(10, places);
}
};function _round(x, places) {if (places === 0) {
  return Math.round(x);
} else {
  var m = Math.pow(10, places);
  return Math.round(m * x) / m;
}
};function _random(state, max) {return state.sample() * max};function _step(tick, value, start) {return tick < start ? 0 : value};function _pulse(t, resolution, value, start, interval) {if (t < start) {
  return 0;
} else if (interval === 0) {
  return t === start ? value : 0;
} else {
  return (t - start) % Math.round(interval * resolution) === 0 ? value : 0;
}
};function _ramp(t, resolution, rate, start) {return t <= start ? 0 : (t - start) * (rate / resolution)};function _rand_bernoulli(state, p) {return state.sample() < p ? 1 : 0};function _rand_binomial(state, n, p) {if (n >= 40) {
  var np = n * p;
  if (np < 5 && p < 0.05) {
    return _rand_poisson(state, np);
  } else {
    return Math.min(
      n,
      Math.round(_rand_gaussian(state, np, Math.sqrt(np * (1 - p)))),
    );
  }
}

var x = 0;
for (var i = 0; i < n; ++i) {
  x += _rand_bernoulli(state, p);
}
return x;
};function _rand_discrete(state, min, max) {return min + Math.floor(state.sample() * ((max + 1) - min))};function _rand_uniform(state, min, max) {return min + (state.sample() * (max - min))};function _rand_exp(state, lambda) {return _rand_exp_percentile(lambda, state.sample())};function _rand_exp_percentile(lambda, perc) {return -Math.log(1 - perc) / lambda};function _rand_gaussian(state, mu, sigma) {do {
  var u1 = state.sample();
  var u2 = state.sample();
} while (u1 <= Number.MIN_VALUE);
var z0 = Math.sqrt(-2 * Math.log(u1)) * Math.cos(Math.PI * 2 * u2);
return z0 * sigma + mu;
};function _rand_gaussian_percentile(mu, sigma, q) {var x = 2 * q - 1;
var a = 0.140012;
var logtop = Math.log(1 - x * x);
var logterm = logtop / 2;
var piterm = 2 / (Math.PI * a) + logterm;
var inverf =
  _sign(x) * Math.sqrt(Math.sqrt(piterm * piterm - logtop / a) - piterm);

return mu + sigma * Math.SQRT2 * inverf;
};function _rand_pareto(state, xm, alpha) {return xm / Math.pow(1 - state.sample(), 1 / alpha)};function _rand_poisson(state, lambda) {return _rand_poisson_percentile(lambda, state.sample())};function _rand_poisson_percentile(lambda, q) {if (lambda === 0) return 0;

// Knuth-style naive simulation
if (lambda < 20) {
  // Cap Q to avoid infinite loops.
  if (q > 0.999999) {
    q = 0.999999;
  }
  var f = Math.exp(-lambda);
  var F = f;
  var k = 1;
  while (F < q) {
    f = (lambda / k) * f;
    F = F + f;
    k = k + 1;
  }
  return k - 1;
}

var q = _rand_gaussian_percentile(0, 1, q);
var t1 = q * Math.pow(lambda, 0.5);
var t2 = (q * q - 1) / 6;
var t3 = (Math.pow(lambda, -0.5) * (q - q * q * q)) / 72;

return Math.max(0, Math.floor(t1 + t2 + t3 + lambda + 0.5));
};function _rand_poisson_binomial(state, ps) {var x = 0,
  n = ps.length;
for (var i = 0; i < n; ++i) {
  x += _rand_bernoulli(state, ps[i]);
}
return x;
};function _rand_student(state, dof) {var w;
do {
  var u1 = 2 * state.sample() - 1;
  var u2 = 2 * state.sample() - 1;
} while ((w = u1 * u1 + u2 * u2) > 1);

return u1 * Math.sqrt((dof * (Math.exp((-2 / dof) * Math.log(w)) - 1)) / w);
};
function Model160(startTime, resolution, stepCount, sketches, lookups, tables, containingModelId, submodelId) {
  __Model.call(this, 1 + stepCount * resolution, 21, new ModelState160(startTime, resolution, stepCount, sketches, lookups, tables, containingModelId, submodelId));
}
Model160.prototype = Object.create(__Model.prototype);
Model160.prototype.getOffsetById = function getOffsetById(elementId) {
  switch (elementId) {
    case 161:
      return 0;
    case 162:
      return 1;
    case 170:
      return 2;
    case 171:
      return 3;
    case 172:
      return 4;
    case 174:
      return 5;
    case 175:
      return 6;
    case 176:
      return 7;
    case 178:
      return 8;
    case 184:
      return 9;
    case 185:
      return 10;
    case 186:
      return 11;
    case 187:
      return 12;
    case 189:
      return 13;
    case 190:
      return 14;
    case 341:
      return 15;
    case 351:
      return 16;
    case 359:
      return 17;
    case 360:
      return 18;
    case 376:
      return 19;
    case 378:
      return 20;
  }
};
Model160.prototype.getSubmodelById = function getSubmodelById(elementId) {
  switch (elementId) {
  }
};
Model160.prototype.computeSketch = function computeSketch($$, $o, $t, id, $blank) {
  var self = this;
  var $res = $$.res;
  var $dt = 1 / $res;
  var $v = $$.values;
  var $c = $$.calc;
  switch (id) {
    default:
      return __sketchLookup($$.sk[id], $t / $res, $blank);
  }
};
Model160.prototype.calc1 = function calc1() {
  var self = this;
  var $o = this.outputs;
  var $$ = this.state;
  var $res = $$.res;
  var $t = $$.t;
  var $dt = 1 / $res;
  var $v = $$.values;
  var $c = $$.calc;
  if ($t === 0 && $c[4] === 1) $v[4] = 0;
  if ($c[6] === 1) $v[6] = 5;
  if ($c[8] === 1) $v[8] = 14;
  if ($c[10] === 1) $v[10] = 100;
  if ($c[11] === 1) $v[11] = 5;
  if ($c[12] === 1) $v[12] = 2;
  if ($c[14] === 1) $v[14] = __sketchLookup(__getTableColumn($$.tables, 360, 361), $t / $res, 0);
  if ($c[16] === 1) $v[16] = __sketchLookup(__getTableColumn($$.tables, 360, 362), $t / $res, 0);
  if ($c[17] === 1) $v[17] = __sketchLookup(__getTableColumn($$.tables, 360, 363), $t / $res, 0);
  if ($c[19] === 1) $v[19] = __sketchLookup(__getTableColumn($$.tables, 360, 373), $t / $res, 0);
  if ($c[20] === 1) $v[20] = self.computeSketch($$, $o, $t, 378, 0);
  if ($t === 0 && $c[2] === 1) $v[2] = $v[14] - $v[16];
  if ($t === 0 && $c[3] === 1) $v[3] = $v[16];
  if ($c[9] === 1) $v[9] = _div($v[4], $v[14], 0) * 100;
  if ($c[1] === 1) $v[1] = _div($v[3], $v[8], 0);
  $$._171_acc -= $v[1] * $dt;
  $$._172_acc += $v[1] * $dt;
  if ($c[13] === 1) $v[13] = _div($v[3], $v[14], 0) * 100;
  if ($c[7] === 1) $v[7] = $v[13] < $v[10] ? $v[11] : $v[12];
  if ($c[15] === 1) $v[15] = $v[7] * $v[3];
  if ($c[5] === 1) $v[5] = _max([0, _div($v[15] * $v[2], $v[14], 0)]);
  if ($c[0] === 1) $v[0] = $v[17] + _min([$v[2], $v[5] * $v[6] / 100]) * $v[20];
  $$._170_acc -= $v[0] * $dt;
  $$._171_acc += $v[0] * $dt;
  for (var i = 0; i < 21; ++i) $o[i][$t] = $v[i];
  $v[4] += $$._172_acc;
  $$._172_acc = 0;
  $v[2] += $$._170_acc;
  $$._170_acc = 0;
  $v[3] += $$._171_acc;
  $$._171_acc = 0;
};
Model160.prototype.calc = function calc() {
  this.calc1();
};
function ModelState160(startTime, resolution, stepCount, sketches, lookups, tables, containingModelId, submodelId) {
  this._172_acc = 0;
  this._170_acc = 0;
  this._171_acc = 0;
  this.st = startTime;
  this.scenario = null;
  this.zt = -startTime * resolution;
  this.res = resolution;
  this.sk = __buildSketches(sketches, 160, containingModelId, submodelId);
  this.tables = tables[160];
  this.lk = __buildLookups(lookups, 160);
  this.t = 0;
  this.values = new Float64Array(21);
  this.calc = new Int8Array(21);
  for (var i = 0; i < this.calc.length; ++i) this.calc[i] = 1;
}
Model160.prototype.getValues = function getValues(panelName, elementName) {
  var offset = this.getOffsetById(Model160.getId(panelName, elementName));
  return offset === undefined ? undefined : this.outputs[offset];
};
Model160.prototype.setInput = function setInput(panelName, elementName, value) {
  var offset = this.getOffsetById(Model160.getId(panelName, elementName));
  if (offset !== undefined) {
    this.setInputByOffset(offset, value);
  }
};
Model160.getId = function getId(panelName, elementName) {
  switch (panelName) {
    case "Simplest \'SIR\' model":
      switch (elementName) {
        case "Newly infected per day":
          return 161;
        case "Recover per day":
          return 162;
        case "Susceptible (uninfected)":
          return 170;
        case "Infected people":
          return 171;
        case "Resistant people":
          return 172;
        case "contacts/day between Infected and Susceptible people":
          return 174;
        case "% of contact events pass infection":
          return 175;
        case "contacts/day per person":
          return 176;
        case "Recovery time days":
          return 178;
        case "% resistant people":
          return 184;
        case "SET % infected population at which to cut contact rate":
          return 185;
        case "Normal contacts/day per person":
          return 186;
        case "SET constrained contacts/day per person":
          return 187;
        case "Infected % of population":
          return 189;
        case "Initial population":
          return 190;
        case "total contacts/day by Infected people":
          return 341;
        case "Initial number infected":
          return 351;
        case "in-bound infections/day":
          return 359;
        case "initial and inbound data":
          return 360;
        case "Infections per day data":
          return 376;
        case "Reduction in contacts per day %":
          return 378;
      }
      return;
  }
};
Model160.getPanels = function getPanels() {
  return [{
    name: "Simplest \'SIR\' model",
    root: true
  }];
};
Model160.elementMetadata = {
  161: {
    id: 161,
    panel: "Simplest \'SIR\' model",
    name: "Newly infected per day",
    type: "Flow"
  },
  162: {
    id: 162,
    panel: "Simplest \'SIR\' model",
    name: "Recover per day",
    type: "Flow"
  },
  170: {
    id: 170,
    panel: "Simplest \'SIR\' model",
    name: "Susceptible (uninfected)",
    type: "Stock"
  },
  171: {
    id: 171,
    panel: "Simplest \'SIR\' model",
    name: "Infected people",
    type: "Stock"
  },
  172: {
    id: 172,
    panel: "Simplest \'SIR\' model",
    name: "Resistant people",
    type: "Stock"
  },
  174: {
    id: 174,
    panel: "Simplest \'SIR\' model",
    name: "contacts/day between Infected and Susceptible people",
    type: "Variable"
  },
  175: {
    id: 175,
    panel: "Simplest \'SIR\' model",
    name: "% of contact events pass infection",
    type: "Variable"
  },
  176: {
    id: 176,
    panel: "Simplest \'SIR\' model",
    name: "contacts/day per person",
    type: "Variable"
  },
  178: {
    id: 178,
    panel: "Simplest \'SIR\' model",
    name: "Recovery time days",
    type: "Variable"
  },
  184: {
    id: 184,
    panel: "Simplest \'SIR\' model",
    name: "% resistant people",
    type: "Variable"
  },
  185: {
    id: 185,
    panel: "Simplest \'SIR\' model",
    name: "SET % infected population at which to cut contact rate",
    type: "Variable"
  },
  186: {
    id: 186,
    panel: "Simplest \'SIR\' model",
    name: "Normal contacts/day per person",
    type: "Variable"
  },
  187: {
    id: 187,
    panel: "Simplest \'SIR\' model",
    name: "SET constrained contacts/day per person",
    type: "Variable"
  },
  189: {
    id: 189,
    panel: "Simplest \'SIR\' model",
    name: "Infected % of population",
    type: "Variable"
  },
  190: {
    id: 190,
    panel: "Simplest \'SIR\' model",
    name: "Initial population",
    type: "Variable"
  },
  341: {
    id: 341,
    panel: "Simplest \'SIR\' model",
    name: "total contacts/day by Infected people",
    type: "Variable"
  },
  351: {
    id: 351,
    panel: "Simplest \'SIR\' model",
    name: "Initial number infected",
    type: "Variable"
  },
  359: {
    id: 359,
    panel: "Simplest \'SIR\' model",
    name: "in-bound infections/day",
    type: "Variable"
  },
  360: {
    id: 360,
    panel: "Simplest \'SIR\' model",
    name: "initial and inbound data",
    type: "DataTable",
    columns: [{
      id: 361,
      name: "init population"
    }, {
      id: 362,
      name: "init infected"
    }, {
      id: 363,
      name: "inbound infections"
    }, {
      id: 373,
      name: "infections per day data"
    }, {
      id: 374,
      name: "deaths per day data"
    }]
  },
  376: {
    id: 376,
    panel: "Simplest \'SIR\' model",
    name: "Infections per day data",
    type: "Variable"
  },
  378: {
    id: 378,
    panel: "Simplest \'SIR\' model",
    name: "Reduction in contacts per day %",
    type: "Variable"
  }
};
function Model371(startTime, resolution, stepCount, sketches, lookups, tables, containingModelId, submodelId) {
  __Model.call(this, 1 + stepCount * resolution, 0, new ModelState371(startTime, resolution, stepCount, sketches, lookups, tables, containingModelId, submodelId));
}
Model371.prototype = Object.create(__Model.prototype);
Model371.prototype.getOffsetById = function getOffsetById(elementId) {
  switch (elementId) {
  }
};
Model371.prototype.getSubmodelById = function getSubmodelById(elementId) {
  switch (elementId) {
  }
};
Model371.prototype.computeSketch = function computeSketch($$, $o, $t, id, $blank) {
  var self = this;
  var $res = $$.res;
  var $dt = 1 / $res;
  var $v = $$.values;
  var $c = $$.calc;
  switch (id) {
    default:
      return __sketchLookup($$.sk[id], $t / $res, $blank);
  }
};
Model371.prototype.calc1 = function calc1() {
  var self = this;
  var $o = this.outputs;
  var $$ = this.state;
  var $res = $$.res;
  var $t = $$.t;
  var $dt = 1 / $res;
  var $v = $$.values;
  var $c = $$.calc;
  for (var i = 0; i < 0; ++i) $o[i][$t] = $v[i];
};
Model371.prototype.calc = function calc() {
  this.calc1();
};
function ModelState371(startTime, resolution, stepCount, sketches, lookups, tables, containingModelId, submodelId) {
  this.st = startTime;
  this.scenario = null;
  this.zt = -startTime * resolution;
  this.res = resolution;
  this.sk = __buildSketches(sketches, 371, containingModelId, submodelId);
  this.tables = tables[371];
  this.lk = __buildLookups(lookups, 371);
  this.t = 0;
  this.values = new Float64Array(0);
  this.calc = new Int8Array(0);
  for (var i = 0; i < this.calc.length; ++i) this.calc[i] = 1;
}
Model371.prototype.getValues = function getValues(panelName, elementName) {
  var offset = this.getOffsetById(Model371.getId(panelName, elementName));
  return offset === undefined ? undefined : this.outputs[offset];
};
Model371.prototype.setInput = function setInput(panelName, elementName, value) {
  var offset = this.getOffsetById(Model371.getId(panelName, elementName));
  if (offset !== undefined) {
    this.setInputByOffset(offset, value);
  }
};
Model371.getId = function getId(panelName, elementName) {
  switch (panelName) {
    case "Model 1":
      switch (elementName) {
      }
      return;
  }
};
Model371.getPanels = function getPanels() {
  return [{
    name: "Model 1",
    root: true
  }];
};
Model371.elementMetadata = {};
var Project = {
  byName: {
    "Simplest \'SIR\' model": Model160,
    "Model 1": Model371
  },
  byId: {
    160: Model160,
    371: Model371
  },
  buildRun: function (modelName) {
    return new RunBuilder(Project, modelName);
  }
};

const export_to_json = (data) => {
    var fs = require("fs");
    fs.writeFile ("output.json", JSON.stringify(data), (err) => {
        if (err) throw err;
        console.log("complete");
    })
}
const range = (start, stop, step) => {
    var a = [start], b = start;
    while (b < stop) {
        a.push(b += step || 1);
    }
    return a;
}

const loss = (arr1, arr2) => {
    var res = 0;
    for(let i = 0; i<arr1.length; i++) {
        res += (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]);
    }
    return res
}

const optimize_simple = (range1, range2) => {
    var param1_array = range(range1[0], range1[2], 0.001);
    var param2_array = range(range2[0], range2[1], 0.001);
    var cur_loss = Infinity;
    var ans = Array(2);

    for(let i = 0; i<param1_array.length; i++) {
        for(let j=0; j<param2_array.length; j++) {
            var new_loss = loss(simple_test(param1_array[i], param2_array[j]), infections_per_day_data);

            if (new_loss < cur_loss) {
                console.log("lose:", new_loss);
                cur_loss = new_loss;
                
                ans[0] = param1_array[i];
                ans[1] = param2_array[j];
            }
        }
    }

    return ans;
}


const optimize = (range1, range2) => {
    var param1_array = [5] //range(range1[0], range1[2], 1);
    var param2_array = range(range2[0], range2[1], 0.001);
    var cur_loss = Infinity;
    var ans = [];
    var ans1 = Array(1);
    var ans2 = Array(120);
    for(let i = 0; i<param1_array.length; i++) {
        for (let z = 0; z<infections_per_day_data.length; z++) {
            cur_loss = Infinity;
            for(let j = 0; j<param2_array.length; j++) {
                var new_loss = loss([test(param1_array[i], param2_array[j], z)], [infections_per_day_data[z]]);
                
                if (new_loss < cur_loss) {
                    console.log("loss:", new_loss);
                    cur_loss = new_loss;
                    
                    ans2[z] = param2_array[j];
                    ans1[0] = param1_array[i];
                }

            }
        }
    }
    return {
        "param2": ans2,
    };
}

const simple_test = (val1, val2) => {
    var run = builder.start();
    run.setInput("Simplest 'SIR' model", "Normal contacts/day per person", val1); //change parameter 1
    run.setInput("Simplest 'SIR' model", "Reduction in contacts per day %", val2); //change parameter 2
    run.runToEnd();
    var ans = run.getValues("Simplest 'SIR' model", "Newly infected per day");
    return ans;   
}

const test = (val1, val2, index) => {
    var run = builder.start();
    run.setInput("Simplest 'SIR' model", "Normal contacts/day per person", val1); //change parameter 1
    run.setInput("Simplest 'SIR' model", "Reduction in contacts per day %", val2); //change parameter 2
    run.runToEnd();
    var ans = run.getValues("Simplest 'SIR' model", "Newly infected per day");
    return ans[index];
}

/* PARSE input.JSON containing model inputs*/
const data = JSON.parse(require('fs').readFileSync("./input.JSON", 'utf8')); // reads the JSON file

var init_population = data["init population"];
var init_infected = data["init infected"];
var inbound_infections = data["inbound infections"];
var infections_per_day_data = data["infections per day data"]

/* RUN THE MODEL */
var builder = Project.buildRun("Simplest 'SIR' model").setStartTime(0).setStepCount(init_population.length-1);

var modelConfig = builder.configureModel("Simplest 'SIR' model");
modelConfig.setDataset({
  element: "initial and inbound data",
  columns: {
    "init population": {
      start: 0,
      values: init_population,
    },
    "init infected": {
      start: 0,
      values: init_infected,
    },
    "inbound infections": {
      start: 0,
      values: inbound_infections,
    },
    "infections per day data": {
        start: 0,
        values: infections_per_day_data,
    }
  },
});

/* CHANGE THE CODE BELOW THIS LINE */
/************************************************************************************/

/* EXAMPLE 1: Use flat tuning to find optimal hyperparameters between (4,6) and (0,1) */
optimization_result = optimize_simple([4, 6], [0, 1]); //calls optimization function
export_to_json(optimization_result); //store results of function in output.json

