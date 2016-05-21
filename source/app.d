import std.stdio;
import mir.ndslice;
import std.algorithm;
import std.conv;
import std.array;
import std.range;
import std.string;
import std.math;
import gsl_qrng;
import gsl_rng;
import gsl_randist;

void main(string[] args)
{
  // execute with: tcbuilder <input_file> <output_file>

  if (args.length < 3){
    writeln("usage: tcbuilder <input_file> <output_file>");
    return;
  }
  auto inputFile = File(args[1]);
  auto outputFile = File(args[2], "wt");
  
  // parse input file
  // network-level params
  auto netParams = inputFile
    .readln()
    .chomp()
    .splitter('\t')
    .filter!(a => a != "")
    .take(5)
    .map!(to!int)
    .array();
  int numClasses = netParams[0];
  int numUnits = netParams[1];
  int numConns = netParams[2];
  int numLayers = netParams[3];
  int maxDelay = netParams[4];
  
  // layer map object (currently only a configurable distance between layers, in um)
  // TODO: mappable inter-layer plane geometry
  auto interLayerDelays = inputFile
    .readln()
    .chomp()
    .splitter('\t')
    .filter!(a => a != "")
    .map!(a => a.to!double/1000)
    .take(numClasses - 1)
    .array();

  // class-level params
  string[] className;  // name
  real[] abundance; // relative abundance
  int[] layer;   // soma layer (Z) - following neuro conventions, deeper layers have higher indices
  char[] xDistro;  // X distribution type
  real[] xMin;  // X min coordinate
  real[] xMax;  // X max coordinate
  char[] yDistro;  // Y distribution type
  real[] yMin;  // Y min coordinate
  real[] yMax;  // Y max coordinate
  real[] dyn_C;  // C (membrane conductance)
  real[] dyn_k;  // k
  real[] dyn_vr;  // vr (rheobase?)
  real[] dyn_vt;  // vt (threshold?)
  real[] dyn_peak;  // spike peak voltage
  real[] dyn_a;  // a 
  real[] dyn_b;  // b
  real[] dyn_bhyp;  // b under hyperpolarized conditions
  real[] dyn_c;  // c
  real[] dyn_d;  // d
  real[] dyn_umax;  // umax
  real[] dyn_caInact;  // calcium inactivation current
  real[] stdp_Aplus;  // A plus (STDP positive intercept)
  real[] stdp_Aminus;  // A minus (STDP negative intercept)
  real[] stdp_tauPlus;  // tau plus (STDP positive decay constant)
  real[] stdp_tauMinus;  // tau minus (STDP negative decay constant)
  bool[] plastic;  // plastic post-synapses?
  real[] maxWeight;  // max synaptic weight
  bool[] record;  // record this class of neurons?
  string[] target;  // postsynaptic target
  double[] dend_totLen;  // total dendritic length (mm)
  double[] axon_totLen;  // total axonal length (mm)
  // per-layer dendritic X,Y extents (mm)
  auto dend_classLayerCov = (0.0).repeat(numClasses*numLayers*4)
    .array.sliced(numClasses, numLayers, 2, 2).pack!2;
  // per-layer axonal X,Y extents (mm)
  auto axon_classLayerCov = (0.0).repeat(numClasses*numLayers*4)
    .array.sliced(numClasses, numLayers, 2, 2).pack!2;
  // per-layer dendritic Z abundance/prior
  auto dend_zPrior = new double[numClasses*numLayers].sliced(numClasses, numLayers);
  // per-layer axonal Z abundance/prior
  auto axon_zPrior = new double[numClasses*numLayers].sliced(numClasses, numLayers);

  string[] classParams;
  int i,j,l;
  for (i = 0; i < numClasses; i++){
    classParams = inputFile
      .readln
      .chomp
      .splitter('\t')
      .array();
    j = 0;
    // name
    className ~= classParams[j++];
    // relative abundance
    abundance ~= classParams[j++].to!real;
    // soma layer (Z)
    layer ~= classParams[j++].to!int;
    // X distribution type
    xDistro ~= classParams[j][0];
    // X min coordinate
    auto splitToken = classParams[j++]
      .splitter(',')
      .array();
    xMin ~= splitToken[0][2..$].to!real;
    // X max coordinate
    xMax ~= splitToken[1][0..$-1].to!real;
    // Y distribution type
    yDistro ~= classParams[j][0];
    // Y min coordinate
    splitToken = classParams[j++]
      .splitter(',')
      .array();
    yMin ~= splitToken[0][2..$].to!real;
    // Y max coordinate
    yMax ~= splitToken[1][0..$-1].to!real;
    // C (membrane conductance)
    dyn_C ~= classParams[j++].to!real;
    // k
    dyn_k ~= classParams[j++].to!real;
    // vr rheobase?
    dyn_vr ~= classParams[j++].to!real;
    // vt threshold?
    dyn_vt ~= classParams[j++].to!real;
    // spike peak voltage
    dyn_peak ~= classParams[j++].to!real;
    // a
    dyn_a ~= classParams[j++].to!real;
    // b
    dyn_b ~= classParams[j++].to!real;
    // bhyp
    dyn_bhyp ~= classParams[j++].to!real;
    // c
    dyn_c ~= classParams[j++].to!real;
    // d
    dyn_d ~= classParams[j++].to!real;
    // umax
    dyn_umax ~= classParams[j++].to!real;
    // calcium inactivation current
    dyn_caInact ~= classParams[j++].to!real;
    // A plus (STDP potentiating intercept)
    stdp_Aplus ~= classParams[j++].to!real;
    // A minus (STDP depressive intercept)
    stdp_Aminus ~= classParams[j++].to!real;
    // tau plus (STDP potentiating decay constant)
    stdp_tauPlus ~= classParams[j++].to!real;
    // tau minus (STDP depressive decay constant)
    stdp_tauMinus ~= classParams[j++].to!real;
    // plastic?
    plastic ~= classParams[j++].to!int != 0;
    // max synaptic weight
    maxWeight ~= classParams[j++].to!real;
    // record?
    record ~= classParams[j++].to!int != 0;
    // postsynaptic target
    target ~= classParams[j++];
    // total dendritic length (mm)
    dend_totLen ~= classParams[j++].to!real/1000;
    // total axonal length (mm)
    axon_totLen ~= classParams[j++].to!real/1000;
    // per-layer dendritic X,Y extents
    for (l = 0; l < numLayers; l++){
      dend_classLayerCov[i,l][0,0] = classParams[j].to!real/1000;
      dend_classLayerCov[i,l][1,1] = classParams[j++].to!real/1000;
    }
    // per-layer axonal X,Y extents
    for (l = 0; l < numLayers; l++){
      axon_classLayerCov[i,l][0,0] = classParams[j].to!real/1000;
      axon_classLayerCov[i,l][1,1] = classParams[j++].to!real/1000;
    }
    // per-layer dendritic Z abundance/prior
    for (l = 0; l < numLayers; l++){
      dend_zPrior[i,l] = classParams[j++].to!real;
    }
    // per-layer axonal Z abundance/prior
    for (l = 0; l < numLayers; l++){
      axon_zPrior[i,l] = classParams[j++].to!real;
    }
  }

  //roundabout way to avoid rounding errors
  // for some reason recurrence! wasn't working, so zip!ping sequence!s instead
  auto prefixSum = sequence!
    ((a,n) => cast(int)(sum(abundance[0..n])*numUnits))(0)
    .take(numClasses + 1);
  int[] unitsPerClass = zip(prefixSum.dropOne, prefixSum.take(numClasses))
    .map!("a[0] - a[1]")
    .array();
  auto unitClass = unitsPerClass
    .enumerate
    .map!("a.index.repeat(a.value)")
    .joiner;
  version(imperative){
    auto unitClassArr = unitClass.array;
  }

  //normalize process lengths as a probability
  auto dend_classPrior = dend_totLen.dup;
  dend_classPrior[] /= dend_totLen.sum;
  auto axon_classPrior = axon_totLen.dup;
  axon_classPrior[] /= axon_totLen.sum;
  
  // place neurons
  Slice!(2,double*)[] locations;
  Slice!(2,double*) delegate() choose;

  const(gsl_rng_type)* rngType;
  gsl_rng* rng;
  gsl_rng_env_setup();
  rngType = gsl_rng_default;
  rng = gsl_rng_alloc(rngType);
  gsl_qrng* qrng;
  
  foreach(classNum; numClasses.iota){
    j = 0;
    switch(xDistro[classNum]){ // ignoring y distro (i.e. assuming both equal) for expediency for now.
    case 'U': // uniform random
      choose = delegate Slice!(2,double*)() {
	auto pair = new double[2].sliced(2,1);
	pair[0,0] = gsl_rng_uniform(rng) * (xMax[classNum] - xMin[classNum]) + xMin[classNum];
	pair[1,0] = gsl_rng_uniform(rng) * (yMax[classNum] - yMin[classNum]) + yMin[classNum];
	return pair;
      };
      break;
    case 'L': // in a line
      choose = delegate Slice!(2,double*)() {
      	auto pair = new double[2].sliced(2,1);
      	pair[0,0] = xMin[classNum]
	+ ((xMax[classNum] - xMin[classNum])/(unitsPerClass[classNum] - 1)) * j;
	pair[1,0] = yMin[classNum]
	+ ((yMax[classNum] - yMin[classNum])/(unitsPerClass[classNum] - 1)) * j++;
	return pair;
      };
      break;
    case 'H': // Halton sequence
      if (qrng is null) {
	qrng = gsl_qrng_alloc(gsl_qrng_halton, 2);
      } else if (qrng.type != gsl_qrng_halton){
	gsl_qrng_free(qrng);
	qrng = gsl_qrng_alloc(gsl_qrng_halton, 2);	
      }
      choose = delegate Slice!(2,double*)() {
	double longer = max(xMax[classNum] - xMin[classNum], yMax[classNum] - yMin[classNum]);
	double[] x = new double[](2);
	auto pair = x.sliced(2,1);
	// draw while not in bounds
	// strategy is to scale both dimensions equally and throw out
	// out-of-bounds draws, so that we don't get crowding along a shorter
	// dimension.
	// double initialized to NaN, so enter loop with check for NaN in first condition
	while (x[0].isNaN || x[0] < xMin[classNum] || x[0] > xMax[classNum]
	       || x[1] < yMin[classNum] || x[1] > yMax[classNum]){ 
	  gsl_qrng_get(qrng, &x[0]);
	  pair[] *= longer;
	}
	return pair;
      };
      break;
    case 'N': // Niederreiter sequence
      if (qrng is null) {
	qrng = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
      } else if (qrng.type != gsl_qrng_niederreiter_2){
	gsl_qrng_free(qrng);
	qrng = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);	
      }
      choose = delegate Slice!(2,double*)() {
	double longer = max(xMax[classNum] - xMin[classNum], yMax[classNum] - yMin[classNum]);
	double[] x = new double[](2);
	auto pair = x.sliced(2,1);
	// same as delegate above
	while (x[0].isNaN || x[0] < xMin[classNum] || x[0] > xMax[classNum]
	       || x[1] < yMin[classNum] || x[1] > yMax[classNum]){ 
	  gsl_qrng_get(qrng, &x[0]);
	  pair[] *= longer;
	}
	return pair;
      };
      break;
    default:
      break;
    }
    locations ~= generate!choose.take(unitsPerClass[classNum]).array.sort!("a[0,0] < b[0,0]").array;
  }

  // calculate connection probabilities for each pair of neurons
  // 1/(sqrt(det(2pi*(S1+S2)))) * exp ^ (-1/2 * (m1 - m2)^T*(S1 + S2)^-1*(m1-m2))
  version (functional){
    // warning: functional code not updated to handle zero variance convention
    auto scales = unitClass.enumerate
      .map!(pre => unitClass.enumerate
	    .map!(post => numLayers.iota
		  .map!(layer =>
			(1 / (sqrt(det(elMul(msum(axon_classLayerCov[pre.value, layer],
						  dend_classLayerCov[post.value, layer]),
					     (2*PI))
				       ))))
			* exp(-0.5*(msum(locations[pre.index], locations[post.index].elMul(-1.0))
				    .transposed
				    .mmul(inv(msum(axon_classLayerCov[pre.value, layer],
						   dend_classLayerCov[post.value, layer]))
					  .mmul(msum(locations[pre.index],
						     locations[post.index].elMul(-1.0)))
					  ))[0,0])
			)));
    //writeln(scales);
  }
  version (imperative){ // seemingly faster
    auto scales = new double[](numUnits * numUnits * numLayers)
      .sliced(numUnits, numUnits, numLayers);
    for (auto elems = scales.byElement; !elems.empty; elems.popFront){
      auto S1 = axon_classLayerCov[unitClassArr[elems.index[0]], elems.index[2]].slice;
      auto S2 = dend_classLayerCov[unitClassArr[elems.index[1]], elems.index[2]].slice;
      // input has zero listed as variance if no probability mass for that class in that layer.
      if (S1[0,0] != 0 && S2[0,0] != 0){
	S2[] += S1;
	auto S1p2tau = S2.slice;
	S1p2tau[] *= 2*PI;
	double coeff = 1/sqrt(det(S1p2tau));
	auto mdiff = locations[elems.index[0]].slice;
	mdiff[] -= locations[elems.index[1]];
	double exponent = mdiff.transposed.mmul(S2.inv).mmul(mdiff)[0,0] * (-0.5);
	elems.front = coeff * exp(exponent);
      } else {
	//zero
	elems.front = 0.0;
      }
    }
    //writeln(scales);
  }
  // need to store these for later sampling of synapse locations
  // S3 = S1*(S1+S2)^-1*S2
  auto productCovs = new double[](numUnits * numUnits * numLayers * 2 * 2)
    .sliced(numUnits, numUnits, numLayers, 2, 2)
    .pack!2;
  // m3 = S2*(S1+S2)^-1*m1 + S1*(S1+S2)^-1*m2
  // see: http://math.stackexchange.com/questions/157172 \
  //   /product-of-two-multivariate-gaussians-distributions
  auto productMeans = new double[](numUnits * numUnits * numLayers * 2)
    .sliced(numUnits, numUnits, numLayers, 2, 1)
    .pack!2;
  //could do this with byElement, but can't easily combine combinable computations
  for (i = 0; i < numUnits; i++){
    for (j = 0; j < numUnits; j++){
      for (l = 0; l < numLayers; l++){
	auto s1 = axon_classLayerCov[unitClassArr[i], l].slice;
	auto s2 = dend_classLayerCov[unitClassArr[j], l].slice;
	auto s12inv = s1.slice;
	s12inv[] += s2;
	s12inv = s12inv.inv;
	productCovs[i,j,l][] = s1.mmul(s12inv).mmul(s2);
	productMeans[i,j,l][] = s2.mmul(s12inv).mmul(locations[i]);
	productMeans[i,j,l][] += s1.mmul(s12inv).mmul(locations[j]);
      }
    }
  }

  // these probability calculations don't take into account layer boundaries,
  // but maybe I don't care anymore...?
  
  // [presynaptic neuron, connection number, postsynaptic neuron (idx 0) and layer (idx 1)]
  auto connections = new int[](numUnits * numConns * 2).sliced(numUnits, numConns, 2);
  // choose connections
  for (i = 0; i < numUnits; i++){
    // TODO: normalize against max sum (dump extra probability into null self-connections)
    auto unitLayerConnProbs = scales[i,0..$,0..$].joiner.array;
    unitLayerConnProbs[] *= unitClass
      .map!(post => numLayers.iota
	    .map!(layer => axon_classPrior[unitClassArr[i]] //relative length of presynaptic axon
		  * axon_zPrior[unitClassArr[i], layer] //relative amount in given layer
		  * dend_classPrior[post] // relative length of postsynaptic dendrite
		  * dend_zPrior[post, layer] // relative amount in given layer
		  ))
      .joiner
      .array[];
    gsl_ran_discrete_t* sampler =
      gsl_ran_discrete_preproc(numUnits * numLayers, &unitLayerConnProbs[0]); 
    connections[i][] +=
      generate!( () => gsl_ran_discrete(rng, sampler))
      .take(numConns)
      .map!(a => [a/numLayers, a % numLayers])
      .array;
  }

  // choose synapse locations for each connection
  auto synapseLocs = new double[](numUnits * numConns * 2)
    .sliced(numUnits, numConns, 2, 1)
    .pack!2;
  for (i = 0; i < numUnits; i++){
    for (j = 0; j < numConns; j++){
      gsl_ran_bivariate_gaussian(rng,
				 productCovs[i, //presynaptic unit
					     connections[i,j,0], // postsynaptic unit
					     connections[i,j,1] // layer selected for synapse location
					     ][0,0], // x
				 productCovs[i, //ditto
					     connections[i,j,0], //ditto
					     connections[i,j,1] //ditto
					     ][1,1], // y
				 0, // rho; just assuming independence for now - can calculate
				 &synapseLocs[i,j][0,0], // x location
				 &synapseLocs[i,j][1,0]  // y location
				 );				 
      // gsl returns a point from a standard distribution about the origin
      // need to shift that to the mean of the distribution of synapse locations for this connection
      synapseLocs[i,j][] += productMeans[i, //pre
					 connections[i,j,0], //post
					 connections[i,j,1] //layer of connection
					 ][];
    }
  }
  
  // calculate distances (presynaptic to synapse + synapse to postsynaptic)
  auto distances = new double[](numUnits * numConns).sliced(numUnits, numConns);
  for (i = 0; i < numUnits; i++){
    for (j = 0; j < numConns; j++){
      // pre (x,y) to synapse (x,y)
      double horiz = sqrt((synapseLocs[i, j][0,0] - locations[i][0,0])^^2 +
			     (synapseLocs[i, j][1,0] - locations[i][1,0])^^2);
      int top = min(layer[unitClassArr[i]], connections[i,j,1]);
      int bottom = max(layer[unitClassArr[i]], connections[i,j,1]);
      double vert = top == bottom ? 0 : interLayerDelays[top..bottom].sum;
      double dist1 = sqrt(horiz^^2 + vert^^2);
      // synapse (x,y) to post (x,y)
      horiz = sqrt((locations[connections[i,j,0]][0,0] - synapseLocs[i, j][0,0])^^2 +
		   (locations[connections[i,j,0]][1,0] - synapseLocs[i, j][1,0])^^2);
      top = min(layer[unitClassArr[connections[i,j,0]]], connections[i,j,1]);
      bottom = max(layer[unitClassArr[connections[i,j,0]]], connections[i,j,1]);
      vert = top == bottom ? 0: interLayerDelays[top..bottom].sum;
      distances[i,j] = sqrt(horiz^^2 + vert^^2) + dist1;
    }
  }

  // write network-level params
  outputFile.writeln(numClasses, ",", numUnits, ",", numConns, ",", maxDelay);
  // write class-level params
  foreach (cl; numClasses.iota){
    outputFile.writeln(unitsPerClass[cl], ",",
		       dyn_C[cl], ",",
		       dyn_k[cl], ",",
		       dyn_vr[cl], ",",
		       dyn_vt[cl], ",",
		       dyn_peak[cl], ",",
		       dyn_a[cl], ",",
		       dyn_b[cl], ",",
		       dyn_bhyp[cl], ",",
		       dyn_c[cl], ",",
		       dyn_d[cl], ",",
		       dyn_umax[cl], ",",
		       dyn_caInact[cl], ",",
		       stdp_Aplus[cl], ",",
		       stdp_Aminus[cl], ",",
		       stdp_tauPlus[cl], ",",
		       stdp_tauMinus[cl], ",",
		       plastic[cl] ? 1 : 0, ",",
		       maxWeight[cl], ",",
		       record[cl] ? 1 : 0
		       );
  }
  // write neuron-level params
  foreach (ucl; unitClass){
    outputFile.writeln(ucl);
  }
  // write synapse-level params
  for (i = 0; i < numUnits; i++){
    for (j = 0; j < numConns; j++){
      outputFile.write(connections[i,j,0], ",",
		       "<init wt>", ",", //TODO: calculate initial weight
		       0, ",", //initial change in weight (zero for new network)
		       "<delay>", ";"); //TODO: calculate delay from distance
    }
    outputFile.writeln();
  }
}

// determinant; only works for 2x2 matrix, which is all I need...
private double det(Slice!(2,double*) m){
  return m[0,0]*m[1,1]-m[0,1]*m[1,0];
}

// matrix inverse; only works for 2x2 diagonal matrices (also all I need)
private auto inv(Slice!(2,double*) m) {
  auto n = (cast(double[])[1, 0, 0, 1]).sliced(2,2);
  n.diagonal[] /= m.diagonal;
  return n;
}

//...and a naive matrix multiply. THAT'S all I need.
private auto mmul(ulong M, ulong N)(Slice!(M,double*) a, Slice!(N,double*) b){
  static assert(M == 2 && N == 2);
  assert (a.shape[1] == b.shape[0]);
  auto c = new double[](a.shape[0] * b.shape[1])
    .sliced(a.shape[0], b.shape[1]);
  for (auto els = c.byElement; !els.empty; els.popFront){
    els.front =
      zip(a[els.index[0], 0..$], b[0..$, els.index[1]])
      .map!("a[0] * a[1]")
      .sum;
  }
  return c;
}

// ...and this lamp.
private auto msum(ulong N)(Slice!(N,double*) a, Slice!(N,double*) b){
  assert(a.shape == b.shape);
  static if( N > 1) {
    return zip(a.joiner, b.joiner).map!("a[0] + a[1]").array.sliced(a.shape);
  } else {
    return zip(a,b).map!("a[0] + a[1]").array.sliced(a.shape);
  }
}

private auto elMul(ulong N)(Slice!(N,double*) a, double k){
  static if( N > 1) {
    return zip(a.joiner, k.repeat(a.elementsCount)).map!("a[0] * a[1]").array.sliced(a.shape);
  } else {
    return zip(a, k.repeat(a.elementsCount)).map!("a[0] * a[1]").array.sliced(a.shape);
  }
}

unittest {
  auto a = ((4.0).iota.array).sliced(2,2);
  auto b = ((6.0).iota.array).sliced(2,3);
  auto c = (cast(double[])[2, 0, 0, 2]).sliced(2,2);
  assert(msum(a, c) == [2.0, 1.0, 2.0, 5.0].sliced(2,2));
  assert(elMul(b, 2.0) == [0.0, 2.0, 4.0,
			   6.0, 8.0, 10.0].sliced(2,3));
  assert(a.mmul(b) ==
	 [[3, 4, 5],
	  [9, 14, 19]]
	 );
  assert(a.det == -2);
  assert(c.inv ==
	 [[0.5, 0],
	  [0, 0.5]]
	 );
}
