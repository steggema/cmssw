/** \class EcalMaxSampleUncalibRecHitProducer
 *   produce ECAL uncalibrated rechits from dataframes 
 *
 *  \author G. Franzoni, E. Di Marco
 *
 */
#include "RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMaxSample.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

EcalUncalibRecHitWorkerMaxSample::EcalUncalibRecHitWorkerMaxSample(const edm::ParameterSet& ps, edm::ConsumesCollector& c) :
  EcalUncalibRecHitWorkerBaseClass( ps ,c)
{
}


void
EcalUncalibRecHitWorkerMaxSample::set(const edm::EventSetup& es)
{
  es.get<EcalGainRatiosRcd>().get(gains);
  es.get<EcalPedestalsRcd>().get(peds);
}

bool
EcalUncalibRecHitWorkerMaxSample::run( const edm::Event & evt, 
                const EcalDigiCollection::const_iterator & itdg, 
                EcalUncalibratedRecHitCollection & result )
{
        DetId detid(itdg->id());

        if ( detid.subdetId() == EcalBarrel ) {
                const EcalPedestals::Item * aped = &peds->barrel(EBDetId(detid).hashedIndex());
                const EcalMGPAGainRatio * aGain  = &gains->barrel(EBDetId(detid).hashedIndex());
                result.push_back( ebAlgo_.makeRecHit(*itdg, aped, aGain) );
        } else {
                const EcalPedestals::Item * aped = &peds->endcap(EEDetId(detid).hashedIndex());
                const EcalMGPAGainRatio * aGain  = &gains->endcap(EEDetId(detid).hashedIndex());          
                result.push_back( eeAlgo_.makeRecHit(*itdg, aped, aGain) );
        }
            
        return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalUncalibRecHitWorkerFactory, EcalUncalibRecHitWorkerMaxSample, "EcalUncalibRecHitWorkerMaxSample" );
