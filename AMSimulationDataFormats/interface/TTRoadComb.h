#ifndef AMSimulationDataFormats_TTRoadComb_h_
#define AMSimulationDataFormats_TTRoadComb_h_

#include <vector>
#include <iosfwd>
#include <bitset>

#include <iostream>


namespace slhcl1tt {

struct TTRoadComb {
  unsigned roadRef;
  unsigned combRef;
  unsigned patternRef;
  float ptSegment;
  unsigned hitBits;

  std::vector<unsigned> stubRefs;
  std::vector<float> stubs_r;
  std::vector<float> stubs_phi;
  std::vector<float> stubs_z;
  std::vector<bool> stubs_bool;
  std::vector<std::string> stubs_bitString;

  int64_t stubs_phi_int(const size_t index) const { return convert<18, 1>(index, 1); }  // 18 bits starting from position 1
  int64_t stubs_r_int(const size_t index) const { return convert<18, 1>(index, 1 + 18); }
  int64_t stubs_z_int(const size_t index) const { return convert<18, 1>(index, 1 + 18 + 18); }
  int64_t bend_int(const size_t index) const { return convert<4, 1>(index, 1 + 18 + 18 + 18); }
  int64_t strip_int(const size_t index) const { return convert<7, 0>(index, 1 + 18 + 18 + 18 + 4); }

 private:
  template<int N, int S>
  // N=number of bits, S=signed or unsigned
  int64_t convert(const size_t index, const size_t pos = 0, const size_t n = N) const
  {
    std::bitset<N> bits(stubs_bitString.at(index), pos, n);
    int64_t ret = static_cast<int64_t>(bits.to_ulong());
    if (S > 0 && bits.test(N - 1)) {  // is signed and is negative
      static const uint64_t ffffffff = -1;
      ret |= (ffffffff << N);
    }
    return ret;
  }
};

// _____________________________________________________________________________
// Output streams
std::ostream& operator<<(std::ostream& o, const TTRoadComb& comb);

}  // namespace slhcl1tt

#endif
