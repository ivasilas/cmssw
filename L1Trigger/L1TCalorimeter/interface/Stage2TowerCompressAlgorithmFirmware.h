///
/// Description: Firmware headers
///
/// Implementation:
///    Concrete firmware implementations
///
/// \author: Jim Brooke - University of Bristol
///

//
//

#ifndef Stage2TowerCompressAlgorithmFirmware_H
#define Stage2TowerCompressAlgorithmFirmware_H

#include "L1Trigger/L1TCalorimeter/interface/Stage2TowerCompressAlgorithm.h"
#include "CondFormats/L1TObjects/interface/CaloParams.h"

namespace l1t {

  // Imp1 is for v1 and v2
  class Stage2TowerCompressAlgorithmFirmwareImp1 : public Stage2TowerCompressAlgorithm {
  public:
    Stage2TowerCompressAlgorithmFirmwareImp1(CaloParams* params);
    virtual ~Stage2TowerCompressAlgorithmFirmwareImp1();
    virtual void processEvent(const std::vector<l1t::CaloTower> & inTowers,
			      std::vector<l1t::CaloTower> & outTowers);
  private:
    CaloParams* params_;

  };
  
}

#endif
