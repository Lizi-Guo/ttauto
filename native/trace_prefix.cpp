// SPDX-License-Identifier: GPL-3.0-or-later
// Migrated from output/disc_even_20261026/trace_prefix.cpp. Exact, bounded
// Newton jobs with a complete deferred frontier; no server or classifier code.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <vector>
using I=std::int64_t; using U=std::uint64_t; using W=__int128;
int g,horizon; std::vector<I> bound,c,tr;
std::vector<std::unordered_set<U>> allowed;
bool restricted=false,stopped=false; U budget,nodes=0,full=0,pruned=0,kept=0,deferred=0;
std::ofstream polys,frontier;
I checked(W a) {
  if(a>std::numeric_limits<I>::max() || a<std::numeric_limits<I>::min())
    throw std::overflow_error("exact int64 overflow; use the Python integer source");
  return I(a);
}
W floorq(W a,I b) { return a/b-((a%b)<0); }
void postpone(int k) {
  ++deferred; frontier<<k-1;
  for(int j=1;j<k;++j) frontier<<' '<<c[j];
  frontier<<'\n';
}
void visit(int k,U bits) {
  if(stopped || nodes>=budget) { stopped=true; postpone(k); return; }
  ++nodes;
  if(k<=g) {
    W t=0; for(int j=1;j<k;++j) t-=(W)c[j]*tr[k-j];
    I lo=checked(-floorq(-(t-bound[k]),k)),hi=checked(floorq(t+bound[k],k));
    for(I a=lo;a<=hi;++a) {
      U next=(bits<<1)|U((a%2+2)%2);
      if(restricted && !allowed[k].count(next)) continue;
      c[k]=a; tr[k]=checked(t-(W)k*a); visit(k+1,next);
      if(a==std::numeric_limits<I>::max()) break;
    }
    return;
  }
  ++full; for(int j=g+1;j<=2*g;++j) c[j]=c[2*g-j];
  for(int j=g+1;j<=horizon;++j) {
    W t=0; for(int i=1;i<j && i<=2*g;++i) t-=(W)c[i]*tr[j-i];
    if(j<=2*g) t-=(W)j*c[j];
    if(t < -bound[j] || t > bound[j]) { ++pruned; return; }
    tr[j]=checked(t);
  }
  ++kept; for(int j=0;j<=2*g;++j) polys<<(j?" ":"")<<c[j];
  polys<<'\n'; if(!polys) throw std::runtime_error("polynomial output failed");
}
int main(int argc,char**argv) { try {
  if(argc!=5) throw std::runtime_error("INPUT POLYNOMIALS REPORT FRONTIER");
  std::ifstream input(argv[1]); int count,depth;
  if(!(input>>g>>horizon>>count) || g<1 || g>24 || horizon<g || horizon>192 || count<0 || count>1000000)
    throw std::runtime_error("invalid native dimensions; use the Python integer source");
  bound.resize(horizon+1); tr.resize(horizon+1); c.resize(2*g+1); allowed.resize(g+1);
  c[0]=1; tr[0]=2*g;
  // Stored traces <=1e15 and coefficients <=int64; a sum of <=192 products
  // fits signed 128-bit arithmetic. Overflowing coefficients fail explicitly.
  for(int k=1;k<=horizon;++k)
    if(!(input>>bound[k]) || bound[k]<0 || bound[k]>1000000000000000LL)
      throw std::runtime_error("trace bound exceeds native exact range; use Python integers");
  restricted=count>0;
  for(int i=0;i<count;++i) {
    U b; if(!(input>>b) || b>=(U(1)<<(2*g+1)) || !(b>>(2*g))) throw std::runtime_error("invalid parity mask");
    for(int k=0;k<=g;++k) allowed[k].insert(b>>(2*g-k));
  }
  if(!(input>>depth) || depth<0 || depth>g) throw std::runtime_error("invalid prefix length");
  U bits=1; bool rejected=false;
  for(int k=1;k<=depth;++k) {
    if(!(input>>c[k])) throw std::runtime_error("missing prefix");
    if(rejected) continue;
    W t=-(W)k*c[k]; for(int j=1;j<k;++j) t-=(W)c[j]*tr[k-j];
    tr[k]=checked(t); bits=(bits<<1)|U((c[k]%2+2)%2);
    if(t < -bound[k] || t > bound[k] || (restricted && !allowed[k].count(bits))) rejected=true;
  }
  if(!(input>>budget) || budget<1) throw std::runtime_error("invalid node budget");
  std::string rest; if(input>>rest) throw std::runtime_error("trailing input");
  polys.open(argv[2]); frontier.open(argv[4]);
  if(!polys || !frontier) throw std::runtime_error("output open failed");
  if(!rejected) visit(depth+1,bits);
  polys.close(); frontier.close();
  if(!polys || !frontier) throw std::runtime_error("output close failed");
  std::ofstream report(argv[3]);
  report<<"{\"completed\":"<<(deferred?"false":"true")<<",\"prefix_rejected\":"<<(rejected?"true":"false")
        <<",\"nodes\":"<<nodes<<",\"parity_admissible_completions\":"<<full<<",\"trace_pruned\":"<<pruned
        <<",\"retained\":"<<kept<<",\"deferred_prefixes\":"<<deferred<<"}\n";
  if(!report) throw std::runtime_error("report write failed");
} catch(const std::exception&e) { std::cerr<<e.what()<<'\n'; return 1; } }
