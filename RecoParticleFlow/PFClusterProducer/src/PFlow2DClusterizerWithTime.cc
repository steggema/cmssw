#include "PFlow2DClusterizerWithTime.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "Math/GenVector/VectorUtil.h"

#include "vdt/vdtMath.h"

#include "TMath.h"

#include <iterator>

#ifdef PFLOW_DEBUG
#define LOGVERB(x) edm::LogVerbatim(x)
#define LOGWARN(x) edm::LogWarning(x)
#define LOGERR(x) edm::LogError(x)
#define LOGDRESSED(x) edm::LogInfo(x)
#else
#define LOGVERB(x) LogTrace(
#define LOGWARN(x) edm::LogWarning(x)
#define LOGERR(x) edm::LogError(x)
#define LOGDRESSED(x) LogDebug(x)
#endif

PFlow2DClusterizerWithTime::
PFlow2DClusterizerWithTime(const edm::ParameterSet& conf) :
    PFClusterBuilderBase(conf),
    _maxIterations(conf.getParameter<unsigned>("maxIterations")),
    _stoppingTolerance(conf.getParameter<double>("stoppingTolerance")),
    _showerSigma2(std::pow(conf.getParameter<double>("showerSigma"),2.0)),
    _timeSigma_eb(std::pow(conf.getParameter<double>("timeSigmaEB"),2.0)),
    _timeSigma_ee(std::pow(conf.getParameter<double>("timeSigmaEE"),2.0)),
    _excludeOtherSeeds(conf.getParameter<bool>("excludeOtherSeeds")),
    _minFracTot(conf.getParameter<double>("minFracTot")),
    _applyChi2ProbCut(conf.getParameter<bool>("applyChi2ProbCut")),
    _useConstantTimeResolution(conf.getParameter<bool>("useConstantTimeResolution")),

    

    _layerMap({ {"PS2",(int)PFLayer::PS2},
	        {"PS1",(int)PFLayer::PS1},
	        {"ECAL_ENDCAP",(int)PFLayer::ECAL_ENDCAP},
	        {"ECAL_BARREL",(int)PFLayer::ECAL_BARREL},
	        {"NONE",(int)PFLayer::NONE},
	        {"HCAL_BARREL1",(int)PFLayer::HCAL_BARREL1},
	        {"HCAL_BARREL2_RING0",(int)PFLayer::HCAL_BARREL2},
		{"HCAL_BARREL2_RING1",100*(int)PFLayer::HCAL_BARREL2},
	        {"HCAL_ENDCAP",(int)PFLayer::HCAL_ENDCAP},
	        {"HF_EM",(int)PFLayer::HF_EM},
		{"HF_HAD",(int)PFLayer::HF_HAD} }) { 

  const std::vector<edm::ParameterSet>& thresholds =
    conf.getParameterSetVector("recHitEnergyNorms");
  for( const auto& pset : thresholds ) {
    const std::string& det = pset.getParameter<std::string>("detector");
    const double& rhE_norm = pset.getParameter<double>("recHitEnergyNorm");    
    auto entry = _layerMap.find(det);
    if( entry == _layerMap.end() ) {
      throw cms::Exception("InvalidDetectorLayer")
	<< "Detector layer : " << det << " is not in the list of recognized"
	<< " detector layers!";
    }
    _recHitEnergyNorms.emplace(_layerMap.find(det)->second,rhE_norm);
  }
  
  _allCellsPosCalc.reset(NULL);
  if( conf.exists("allCellsPositionCalc") ) {
    const edm::ParameterSet& acConf = 
      conf.getParameterSet("allCellsPositionCalc");
    const std::string& algoac = 
      acConf.getParameter<std::string>("algoName");
    PosCalc* accalc = 
      PFCPositionCalculatorFactory::get()->create(algoac, acConf);
    _allCellsPosCalc.reset(accalc);
  }
  // if necessary a third pos calc for convergence testing
  _convergencePosCalc.reset(NULL);
  if( conf.exists("positionCalcForConvergence") ) {
    const edm::ParameterSet& convConf = 
      conf.getParameterSet("positionCalcForConvergence");
    const std::string& algoconv = 
      convConf.getParameter<std::string>("algoName");
    PosCalc* convcalc = 
      PFCPositionCalculatorFactory::get()->create(algoconv, convConf);
    _convergencePosCalc.reset(convcalc);
  }
}

double timeResolution(double energy, bool isEB)
{
  // JAN - ECAL time cleaning based on sigmas
  // this should go to another class that gets the
  // parameters from some global config
  
  double N_EB = 26.4021428571; // # 27.45 <- EGM PAS
  double N_EE = 40.8921428571; // # 41.5937142857 # 36.08 <- EGM PAS
  double N_EE_10 = 49.4773571429; // # For energy < 1 GeV, larger time error
  double N_L_EB = 0.042; //
  double N_L_EE = 0.14; //
  double C_EB = 0.428192; // # 0.27 <- EGM PAS
  double C_EE = 0.; // 0.15528 # 0.18 <- EGM PAS
  double N_EB_05 = 31.4007142857; //
  double P2_EB_05 = 0.0510871; //


  double res = 100.;

  if (isEB)
  {
    if (energy < 0.5 && energy > 0.)
      res = N_EB_05 * N_L_EB/energy + P2_EB_05/pow(energy, 2);
    else if (energy < 5.)
      res = sqrt(pow((N_EB * N_L_EB / energy), 2) + pow(C_EB, 2));
    else if (energy >= 5.)
      res = sqrt(pow((N_EB * N_L_EB / 5.), 2) + pow(C_EB, 2));
  }
  else
  {
      if (energy < 1. && energy > 0.)
        res = sqrt(N_EE_10 * N_L_EE/energy + pow(C_EE, 2));
      else if (energy < 10.)
        res = sqrt(pow((N_EE * N_L_EE / energy), 2) + pow(C_EE, 2));
      else if (energy >= 5.)
        res = sqrt(pow((N_EE * N_L_EE / 10.), 2) + pow(C_EE, 2));
  }
  if (res > 100.)
    res = 100.;
  return res;
  
}

std::pair<double, double> clusterTimeResolution(const reco::PFCluster& cluster)
{
  double sumTimeSigma2 = 0.;
  double sumSigma2 = 0.;

  if (cluster.recHitFractions().size() == 0)
    return std::make_pair(0., 999999.);

  for (unsigned ic = 0; ic < cluster.recHitFractions().size(); ++ic)
  {
    const reco::PFRecHit& rh = *(cluster.recHitFractions()[ic].recHitRef());

    double fraction = cluster.recHitFractions()[ic].fraction();

    bool isEB = (rh.layer() == PFLayer::ECAL_BARREL);

    double res = timeResolution(rh.energy()*fraction, isEB);
    double time = rh.time();
    sumTimeSigma2 += time/res/res;
    sumSigma2 += 1./res/res;
  }

  double clusterRes = sqrt(1./sumSigma2);
  double clusterTime = sumTimeSigma2/sumSigma2;

  return std::make_pair(clusterTime, clusterRes);
}

double clusterChi2Prob(const reco::PFCluster& cluster, const reco::PFRecHit& rhNew, double threshold)
{ 
  if (cluster.recHitFractions().size() == 0)
    return 1.;

  if (rhNew.layer() != PFLayer::ECAL_BARREL && rhNew.layer() != PFLayer::ECAL_ENDCAP)
    return 1.;

  double sumTimeSigma2 = 0.;
  double sumSigma2 = 0.;

  // New hit
  double resNew = timeResolution(rhNew.energy(), (rhNew.layer() == PFLayer::ECAL_BARREL));
  sumTimeSigma2 += rhNew.time()/resNew/resNew;
  sumSigma2 += 1./resNew/resNew;

  for (unsigned ic = 0; ic < cluster.recHitFractions().size(); ++ic)
  {
    const reco::PFRecHit& rh = *(cluster.recHitFractions()[ic].recHitRef());

    // Only look up to the current rechit
    if (rh.detId() == rhNew.detId())
      break;

    double fraction = cluster.recHitFractions()[ic].fraction();


    bool isEB = (rh.layer() == PFLayer::ECAL_BARREL);

    double res = timeResolution(rh.energy()*fraction, isEB);
    double time = rh.time();
    sumTimeSigma2 += time/res/res;
    sumSigma2 += 1./res/res;
  }

  double clusterTime = sumTimeSigma2/sumSigma2;

  // Start chi2 calculation
  double chi2 = 0.;
  // Don't count the first hit (no ndof), but count the outlier below
  unsigned ndof = 0;

  for (unsigned ic = 0; ic < cluster.recHitFractions().size(); ++ic)
  {
    const reco::PFRecHit& rh = *(cluster.recHitFractions()[ic].recHitRef());

    // Only look up to the current rechit
    if (rh.detId() == rhNew.detId())
      break;

    double fraction = cluster.recHitFractions()[ic].fraction();

    if (fraction == 0.)
      continue;

    // A valid hit. 
    ndof += 1;

    bool isEB = (rh.layer() == PFLayer::ECAL_BARREL);
    double res = timeResolution(rh.energy()*fraction, isEB);
    double dtime = rh.time() - clusterTime;

    // Here it's a bit tricky. We know the previous hits were allowed
    // in the cluster only if the chi2 prob > threshold, so we skip
    // those where it's not the case. 
    if (ndof > 0){ 
      if (TMath::Prob(chi2, ndof) > threshold)
        chi2 += dtime*dtime/res/res;
      else
        ndof -= 1;
    }
  }

  // do stuff for the new hit
  double dtime = rhNew.time() - clusterTime;
  chi2 += dtime*dtime/resNew/resNew;

  return TMath::Prob(chi2, ndof);
}



void PFlow2DClusterizerWithTime::
buildClusters(const reco::PFClusterCollection& input,
	      const std::vector<bool>& seedable,
	      reco::PFClusterCollection& output) {
  reco::PFClusterCollection clustersInTopo;
  for( const auto& topocluster : input ) {
    clustersInTopo.clear();
    seedPFClustersFromTopo(topocluster,seedable,clustersInTopo);
    const unsigned tolScal = 
      std::pow(std::max(1.0,clustersInTopo.size()-1.0),2.0);
    growPFClusters(topocluster,seedable,tolScal,0,tolScal,clustersInTopo);
    // step added by Josh Bendavid, removes low-fraction clusters
    // did not impact position resolution with fraction cut of 1e-7
    // decreases the size of each pf cluster considerably
    prunePFClusters(clustersInTopo);
    // recalculate the positions of the pruned clusters
    if( _convergencePosCalc ) { 
      // if defined, use the special position calculation for convergence tests
      _convergencePosCalc->calculateAndSetPositions(clustersInTopo);
    } else {
      if( clustersInTopo.size() == 1 && _allCellsPosCalc ) {
	_allCellsPosCalc->calculateAndSetPosition(clustersInTopo.back());
      } else {
	_positionCalc->calculateAndSetPositions(clustersInTopo);
      }   
    }
    for( auto& clusterout : clustersInTopo ) {
      output.insert(output.end(),std::move(clusterout));
    }
  }
}

void PFlow2DClusterizerWithTime::
seedPFClustersFromTopo(const reco::PFCluster& topo,
		       const std::vector<bool>& seedable,
		       reco::PFClusterCollection& initialPFClusters) const {
  const auto& recHitFractions = topo.recHitFractions();
  for( const auto& rhf : recHitFractions ) {
    if( !seedable[rhf.recHitRef().key()] ) continue;
    initialPFClusters.push_back(reco::PFCluster());
    reco::PFCluster& current = initialPFClusters.back();
    current.recHitFractions().push_back(rhf);
    current.setSeed(rhf.recHitRef()->detId());   
    if( _convergencePosCalc ) {
      _convergencePosCalc->calculateAndSetPosition(current);
    } else {
      _positionCalc->calculateAndSetPosition(current);
    }
  }
}

void PFlow2DClusterizerWithTime::
growPFClusters(const reco::PFCluster& topo,
	       const std::vector<bool>& seedable,
	       const unsigned toleranceScaling,
	       const unsigned iter,
	       double diff,
	       reco::PFClusterCollection& clusters) const {
  if( iter >= _maxIterations ) {
    LOGDRESSED("PFlow2DClusterizerWithTime:growAndStabilizePFClusters")
      <<"reached " << _maxIterations << " iterations, terminated position "
      << "fit with diff = " << diff;
  }      
  if( iter >= _maxIterations || 
      diff <= _stoppingTolerance*toleranceScaling) return;
  // reset the rechits in this cluster, keeping the previous position    
  std::vector<reco::PFCluster::REPPoint> clus_prev_pos;  
  std::vector<double> prevTimes;
  std::vector<double> prevTimeResolutions;
  for( auto& cluster : clusters) {
    const std::pair<double, double> clusterTimeRes = clusterTimeResolution(cluster);
    prevTimes.push_back(clusterTimeRes.first);
    prevTimeResolutions.push_back(clusterTimeRes.second);
    // cluster.setTime(clusterTimeResolution(cluster).first);
    const reco::PFCluster::REPPoint& repp = cluster.positionREP();
    clus_prev_pos.emplace_back(repp.rho(),repp.eta(),repp.phi());
    if( _convergencePosCalc ) {
      if( clusters.size() == 1 && _allCellsPosCalc ) {
	_allCellsPosCalc->calculateAndSetPosition(cluster);
      } else {
	_positionCalc->calculateAndSetPosition(cluster);
      }
    }
    cluster.recHitFractions().clear();
  }
  // loop over topo cluster and grow current PFCluster hypothesis 
  std::vector<double> dist2, frac;


  double fractot = 0, fraction = 0;
  for( const reco::PFRecHitFraction& rhf : topo.recHitFractions() ) {
    const reco::PFRecHitRef& refhit = rhf.recHitRef();
    int cell_layer = (int)refhit->layer();
    if( cell_layer == PFLayer::HCAL_BARREL2 && 
	std::abs(refhit->positionREP().eta()) > 0.34 ) {
      cell_layer *= 100;
    }  
    const double recHitEnergyNorm = 
      _recHitEnergyNorms.find(cell_layer)->second; 
    const math::XYZPoint& topocellpos_xyz = refhit->position();
    dist2.clear(); frac.clear(); fractot = 0;
    // add rechits to clusters, calculating fraction based on distance
    // for( auto& cluster : clusters ) {      
    for( unsigned ic = 0; ic < clusters.size(); ++ic ){
      const reco::PFCluster& cluster = clusters[ic];
      const math::XYZPoint& clusterpos_xyz = cluster.position();
      fraction = 0.0;
      const math::XYZVector deltav = clusterpos_xyz - topocellpos_xyz;
      double d2 = deltav.Mag2()/_showerSigma2;
      const double t2M =(cluster.time()-refhit->time())*(cluster.time()-refhit->time());
      
      double clusterTime = prevTimes[ic];
      double clusterTimeRes = prevTimeResolutions[ic];

      const double t2 = (clusterTime - refhit->time())*(clusterTime - refhit->time());

      
      if (cell_layer == PFLayer::HCAL_BARREL1 ||
	  cell_layer == PFLayer::HCAL_BARREL2 ||
	  cell_layer == PFLayer::ECAL_BARREL) {
        if (_useConstantTimeResolution)
	       d2=d2+t2M/_timeSigma_eb;
       else
       {
          const double tres2 = sqrt(clusterTimeRes*clusterTimeRes + pow(timeResolution(refhit->energy(), true), 2));
          d2 = d2 + t2/tres2;
        }
        
      }

      else if (cell_layer == PFLayer::HCAL_ENDCAP ||
	       cell_layer == PFLayer::HF_EM ||
	       cell_layer == PFLayer::HF_HAD ||
         cell_layer == PFLayer::ECAL_ENDCAP) { // JAN - why did michalis miss the ECAL_ENDCAP?
        if (_useConstantTimeResolution)
          d2=d2+t2M/_timeSigma_ee;
        else
        {
          const double tres2 = sqrt(clusterTimeRes*clusterTimeRes + pow(timeResolution(refhit->energy(), false), 2));
          d2 = d2 + t2/tres2;
        }

      }
      dist2.emplace_back( d2);

      if( d2 > 100 ) {
	LOGDRESSED("PFlow2DClusterizerWithTime:growAndStabilizePFClusters")
	  << "Warning! :: pfcluster-topocell distance is too large! d= "
	  << d2;
      }
      // fraction assignment logic
      if( refhit->detId() == cluster.seed() && _excludeOtherSeeds ) {
	fraction = 1.0;	
      } else if ( seedable[refhit.key()] && _excludeOtherSeeds ) {
	fraction = 0.0;
      } else {
        if (_applyChi2ProbCut && clusterChi2Prob(cluster, *refhit, 0.01) < 0.01)
        {
          fraction = 0.;
        }
        else
          fraction = cluster.energy()/recHitEnergyNorm * vdt::fast_expf( -0.5*d2 );
      }      
      fractot += fraction;
      frac.emplace_back(fraction);
    }
    for( unsigned i = 0; i < clusters.size(); ++i ) {      
      if( fractot > _minFracTot || 
	  ( refhit->detId() == clusters[i].seed() && fractot > 0.0 ) ) {
	frac[i]/=fractot;
      } else {
	continue;
      }
      // if the fraction has been set to 0, the cell 
      // is now added to the cluster - careful ! (PJ, 19/07/08)
      // BUT KEEP ONLY CLOSE CELLS OTHERWISE MEMORY JUST EXPLOSES
      // (PJ, 15/09/08 <- similar to what existed before the 
      // previous bug fix, but keeps the close seeds inside, 
      // even if their fraction was set to zero.)
      // Also add a protection to keep the seed in the cluster 
      // when the latter gets far from the former. These cases
      // (about 1% of the clusters) need to be studied, as 
      // they create fake photons, in general.
      // (PJ, 16/09/08) 
      if( dist2[i] < 100.0 || frac[i] > 0.9999 ) {	
	clusters[i].recHitFractions().emplace_back(refhit, frac[i]);
      }
    }
  }
  // recalculate positions and calculate convergence parameter
  double diff2 = 0.0;  
  for( unsigned i = 0; i < clusters.size(); ++i ) {
    if( _convergencePosCalc ) {
      _convergencePosCalc->calculateAndSetPosition(clusters[i]);
    } else {
      if( clusters.size() == 1 && _allCellsPosCalc ) {
	_allCellsPosCalc->calculateAndSetPosition(clusters[i]);
      } else {
	_positionCalc->calculateAndSetPosition(clusters[i]);
      }
    }
    const double delta2 = 
      reco::deltaR2(clusters[i].positionREP(),clus_prev_pos[i]);    
    if( delta2 > diff2 ) diff2 = delta2;
  }
  diff = std::sqrt(diff2);
  dist2.clear(); frac.clear(); clus_prev_pos.clear();// avoid badness
  growPFClusters(topo,seedable,toleranceScaling,iter+1,diff,clusters);
}

void PFlow2DClusterizerWithTime::
prunePFClusters(reco::PFClusterCollection& clusters) const {
  for( auto& cluster : clusters ) {
    std::vector<reco::PFRecHitFraction> allFracs = 
      std::move(cluster.recHitFractions());
    std::vector<reco::PFRecHitFraction> prunedFracs;
    prunedFracs.reserve(cluster.recHitFractions().size());
    for( const auto& frac : allFracs ) {
      if( frac.fraction() > _minFractionToKeep ) {
	prunedFracs.push_back(std::move(frac));
	cluster.addHitAndFraction( prunedFracs.back().recHitRef()->detId(),
				   prunedFracs.back().fraction()            );
      }
    }
    prunedFracs.shrink_to_fit();
    cluster.recHitFractions() = std::move(prunedFracs);
  }
}


