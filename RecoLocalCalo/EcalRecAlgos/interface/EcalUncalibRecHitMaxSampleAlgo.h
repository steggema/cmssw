#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMaxSampleAlgo_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMaxSampleAlgo_HH

/** \class EcalUncalibRecHitMaxSampleAlgo
  *  Amplitude reconstucted by the difference MAX_adc - min_adc
  *  jitter is sample number of MAX_adc, pedestal is min_adc
  *
  *  \author G. Franzoni, E. Di Marco
  */

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"

template<class C> class EcalUncalibRecHitMaxSampleAlgo : public EcalUncalibRecHitRecAbsAlgo<C>
{
  
 public:
  
  virtual ~EcalUncalibRecHitMaxSampleAlgo<C>() { };
  virtual EcalUncalibratedRecHit makeRecHit(const C& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain);
  virtual EcalUncalibratedRecHit makeRecHit(const C& dataFrame, const double* pedestals,
                                            const double* gainRatios,
                                            const EcalWeightSet::EcalWeightMatrix** weights, 
                                            const EcalWeightSet::EcalChi2WeightMatrix** chi2Matrix) { return EcalUncalibratedRecHit(); }

 private:

};

/// compute rechits
template<class C> EcalUncalibratedRecHit  
EcalUncalibRecHitMaxSampleAlgo<C>::makeRecHit(const C& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain) {

  double maxamplitude = -std::numeric_limits<double>::max();
  double maxpedestal  = 4095;
  double maxjitter    = -1;
  double maxchi2      = -1;
  //bool isSaturated = 0;
  uint32_t maxflags = 0;
  for(int16_t iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
    
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    
    double amplitude = 0.;
    int gainId = sample.gainId();
    
    double pedestal = 0.;
    double gainratio = 1.;
    
    uint32_t flags = 0;
    
    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      gainratio = 1.;
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      gainratio = aGain->gain12Over6();
    }

    amplitude = double(((double)(sample.adc()) - pedestal) * gainratio);
    
    if (gainId == 0) {
      flags = EcalUncalibratedRecHit::kSaturated;
      amplitude = double((4095. - pedestal) * gainratio);
    }
    
    if (amplitude>maxamplitude) {
      maxamplitude = amplitude;
      maxpedestal = pedestal;
      maxjitter = (iSample-5);
      maxflags = flags;
    }
    
    

  }// loop on samples
      
      
  return EcalUncalibratedRecHit( dataFrame.id(), maxamplitude , maxpedestal, maxjitter, maxchi2, maxflags );
}

#endif
