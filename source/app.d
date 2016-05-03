import std.stdio;
import std.experimental.ndslice;
import std.algorithm.iteration;
import std.conv;
import std.array;
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

  auto inputFile = File(args[1]);
  auto outputFile = File(args[2]);
  
  // parse input file
  // network-level params
  auto netParams = inputFile
    .readln()
    .splitter('\t')
    .filter!(a => a != "")
    .array()[0..5]
    .map!(to!int);
  int numClasses = netParams[0];
  int numUnits = netParams[1];
  int numConns = netParams[2];
  int numLayers = netParams[3];
  int maxDelay = netParams[4];
  
  // layer map object (currently only a configurable delay between layers)
  // TODO: mappable inter-layer plane geometry
  auto interLayerDelays = inputFile
    .readln()
    .splitter('\t')
    .filter!(a => a != "")
    .map!(to!int)
    .array()[0..numClasses - 1];

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
      dend_zPrior[i,l] = classParams[j++].to!real/1000;
    }
    // per-layer axonal Z abundance/prior
    for (l = 0; l < numLayers; l++){
      axon_zPrior[i,l] = classParams[j++].to!real/1000;
    }
  }

  // place neurons
  
  // calculate connection probabilities for each pair of neurons
  // choose connections
  // choose synapse locations for each connection
  // calculate distances (presynaptic to synapse + synapse to postsynaptic)
  // write network-level params
  // write class-level params
  // write neuron-level params
  // write synapse-level params
	
}

