import std.stdio;
import std.experimental.ndslice;
import std.algorithm.iteration;
import std.algorithm.searching;
import std.conv;
import std.array;
import std.range;
import std.string;
import gsl_qrng;
import gsl_rng;
import gsl_randist;

struct pair{
  real x;
  real y;
}

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
    .map!(to!int)
    .take(numClasses - 1)
    .array();

  // class-level params
  string[] className;  // name
  real[] abundance; // relative abundance
  int[] layer;   // soma layer (Z)
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
  real[] dend_totLen;  // total dendritic length (mm)
  real[] axon_totLen;  // total axonal length (mm)
  // per-layer dendritic X,Y extents (mm)
  auto dend_xyLayerSd = new real[numClasses*numLayers].sliced(numClasses, numLayers);
  // per-layer axonal X,Y extents (mm)
  auto axon_xyLayerSd = new real[numClasses*numLayers].sliced(numClasses, numLayers);
  // per-layer dendritic Z abundance/prior
  auto dend_zPrior = new real[numClasses*numLayers].sliced(numClasses, numLayers);
  // per-layer axonal Z abundance/prior
  auto axon_zPrior = new real[numClasses*numLayers].sliced(numClasses, numLayers);

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
      dend_xyLayerSd[i,l] = classParams[j++].to!real/1000;
    }
    // per-layer axonal X,Y extents
    for (l = 0; l < numLayers; l++){
      axon_xyLayerSd[i,l] = classParams[j++].to!real/1000;
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
  auto prefixSum = sequence!
    ((a,n) => cast(int)(sum(abundance[0..n])*numUnits))(0)
    .take(numClasses + 1);
  int[] unitsPerClass = zip(prefixSum.dropOne, prefixSum.take(numClasses))
    .map!("a[0] - a[1]")
    .array();
  
  // place neurons
  Slice!(1,double*)[] locations;
  Slice!(1,double*) delegate() choose;

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
      choose = delegate Slice!(1,double*)() {
	auto pair = new double[2].sliced(2);
	// TODO: bounds
	pair[0] = gsl_rng_uniform(rng);
	pair[1] = gsl_rng_uniform(rng);
	return pair;
      };
      break;
    case 'L': // in a line
      choose = delegate Slice!(1,double*)() {
      	auto pair = new double[2].sliced(2);
      	pair[0] = xMin[classNum] + ((xMax[classNum] - xMin[classNum])/(unitsPerClass[classNum] - 1)) * j;
	pair[1] = yMin[classNum] + ((yMax[classNum] - yMin[classNum])/(unitsPerClass[classNum] - 1)) * j++;
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
      choose = delegate Slice!(1,double*)() {
	// TODO: bounds
	double[] x = new double[](2);
	auto pair = x.sliced(2);
	gsl_qrng_get(qrng, &x[0]);
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
      choose = delegate Slice!(1,double*)() {
	// TODO: bounds
	double[] x = new double[](2);
	auto pair = x.sliced(2);
	gsl_qrng_get(qrng, &x[0]);
	return pair;
      };
      break;
    default:
      break;
    }
    locations ~= generate!choose.take(unitsPerClass[classNum]).array();
  }

  // calculate connection probabilities for each pair of neurons
  auto scales = new double[](numUnits*numUnits)
    .sliced(numUnits, numUnits);
  // need to store these for later sampling of synapse locations
  auto productCovs = new double[](numUnits * numUnits * 2 * 2)
    .sliced(numUnits, numUnits, 2, 2)
    .pack!2;
  // ditto
  auto productMeans = new double[](numUnits * numUnits)
    .sliced(numUnits, numUnits);


  // choose connections
  // choose synapse locations for each connection
  // calculate distances (presynaptic to synapse + synapse to postsynaptic)
  // write network-level params
  // write class-level params
  // write neuron-level params
  // write synapse-level params
	
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
private auto mm(Slice!(2,double*) a, Slice!(2,double*) b){
  assert(a.shape[1] == b.shape[0]);
  auto c = new double[](a.shape[0] * b.shape[1])
    .sliced(a.shape[0], b.shape[1]);
  for (auto els = c.byElement; !els.empty; els.popFront){
    // should I have to allocate here?
    auto temp = a[els.index[0], 0..$].array;
    temp[] *= b[0..$, els.index[1]].array[];
    //    writeln(temp);
    els.front = sum(temp[]);
  }
  return c;
}

// ...and this lamp.
unittest {
  auto a = ((4.0).iota.array).sliced(2,2);
  auto b = ((6.0).iota.array).sliced(2,3);
  auto c = (cast(double[])[2, 0, 0, 2]).sliced(2,2);
  assert(a.mm(b) ==
	 [[3, 4, 5],
	  [9, 14, 19]]
	 );
  assert(a.det == -2);
  assert(c.inv ==
	 [[0.5, 0],
	  [0, 0.5]]
	 );
}
