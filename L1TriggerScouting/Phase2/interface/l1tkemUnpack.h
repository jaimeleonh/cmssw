#ifndef L1TriggerScouting_Phase2_l1tkemUnpack_h
#define L1TriggerScouting_Phase2_l1tkemUnpack_h
#include <cstdint>
#include <cmath>

// Section 6.3.3.2 https://www.overleaf.com/project/5fbecbed9fa10637655869e2
// 0 valid
// 16-1 pt
// 29-17 phi
// 43-30 eta
// 47-44 quality
// 58-48 isolation
// 59 charge (for electrons, 1 is negative)
// 69-60 z0 (for electrons, LSB is 0.05cm)
// 77-70 idScore (for electrons, GT encoding, LSB is 1/128)

namespace l1tkemUnpack {
  inline void readshared(const uint64_t datalow,
                         const uint32_t datahigh,
                         uint16_t &pt,
                         int16_t &eta,
                         int16_t &phi,
                         uint8_t &quality,
                         uint16_t &isolation) {  //int
    // bit 0 is the valid bit but it's always true so not worth unpacking
    pt = ((datalow >> 16) & 1) ? ((datalow >> 1) | (-0x8000)) : ((datalow >> 1) & (0xFFFF));         // 16 bits
    phi = ((datalow >> 29) & 1) ? ((datalow >> 17) | (-0x1000)) : ((datalow >> 17) & (0x1FFF));      // 13 bits
    eta = ((datalow >> 43) & 1) ? ((datalow >> 30) | (-0x2000)) : ((datalow >> 30) & (0x3FFF));      // 14 bits
    quality = ((datalow >> 47) & 1) ? ((datalow >> 44) | (-0x8)) : ((datalow >> 44) & (0xF));        // 4 bits
    isolation = ((datalow >> 58) & 1) ? ((datalow >> 48) | (-0x400)) : ((datalow >> 48) & (0x7FF));  // 11 bits
  }
  inline void readshared(const uint64_t datalow,
                         const uint32_t datahigh,
                         float &pt,
                         float &eta,
                         float &phi,
                         uint8_t &quality,
                         float &isolation) {  //float

    // bit 0 is the valid bit but it's always true so not worth unpacking

    uint16_t ptint = ((datalow >> 16) & 1) ? ((datalow >> 1) | (-0x8000)) : ((datalow >> 1) & (0xFFFF));
    pt = ptint * 0.03125f;

    int phiint = ((datalow >> 29) & 1) ? ((datalow >> 17) | (-0x1000)) : ((datalow >> 17) & (0x1FFF));
    phi = phiint * float(M_PI / 4096.);

    int etaint = ((datalow >> 43) & 1) ? ((datalow >> 30) | (-0x2000)) : ((datalow >> 30) & (0x3FFF));
    eta = etaint * float(M_PI / 4096.);

    quality = ((datalow >> 47) & 1) ? ((datalow >> 44) | (-0x8)) : ((datalow >> 44) & (0xF));

    int isolationint = ((datalow >> 58) & 1) ? ((datalow >> 48) | (-0x400)) : ((datalow >> 48) & (0x7FF));
    isolation = isolationint * 0.25f;
  }
  inline void readele(
      const uint64_t datalow, const uint32_t datahigh, int8_t &charge, int16_t &z0, uint8_t &idScore) {  //int
    charge = (datalow & (1llu << 59)) ? -1 : +1;

    uint16_t z0raw = ((datahigh & 0x3f) << 4) | (datalow >> 60);    // 6 bits from high, 4 from low
    z0 = (z0raw & 0x200) ? (z0raw | (-0x200)) : (z0raw & (0x3FF));  // 10 bits

    idScore = (datahigh >> 6) & 0xff;  // 8 bits, unsigned
  }
  inline void readele(
      const uint64_t datalow, const uint32_t datahigh, int8_t &charge, float &z0, float &idScore) {  //float
    int16_t z0int;
    uint8_t idScoreInt;
    readele(datalow, datahigh, charge, z0int, idScoreInt);
    z0 = z0int * 0.05f;
    idScore = idScoreInt * float(1 / 128.f);
  }
}  // namespace l1tkemUnpack

#endif
