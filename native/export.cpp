// SPDX-License-Identifier: GPL-3.0-or-later
// JSON bridge to the in-tree fold_with_map algorithm. Signed maps are preserved.
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "traintracks/traintrack.hpp"
#include "traintracks/build.hpp"
#include "traintracks/fold_map.hpp"
#include "ttauto/ttfoldgraph.hpp"

template<class V> void array(std::ostream& out, const V& values)
{
  out << '[';
  bool comma = false;
  for (const auto value : values) { if (comma) out << ','; out << value; comma = true; }
  out << ']';
}

int main(int argc, char** argv)
{
  try {
    if (argc != 4) throw std::invalid_argument("usage: ttauto_pipeline_export N STRATUM OUTPUT.json");
    const int n = std::stoi(argv[1]), index = std::stoi(argv[2]);
    if (n < 3 || index < 0) throw std::invalid_argument("invalid n or stratum index");
    auto seeds = traintracks::build_traintrack_list(n);
    if (index >= static_cast<int>(seeds.size())) throw std::invalid_argument("stratum index out of range");
    ttauto::ttfoldgraph<traintracks::traintrack> graph(seeds[index]);
    std::ofstream out(argv[3]);
    if (!out) throw std::runtime_error("cannot open output");
    out << "{\"schema\":\"ttauto-native-export-v1\",\"punctures\":" << n << ",\"stratum_index\":" << index << ",\"states\":[";
    for (int v=0; v<graph.vertices(); ++v) {
      if (v) out << ',';
      const auto& source = graph.traintrack(v);
      const auto numbering = source.numbering();
      out << "{\"main_edges\":" << numbering.nedges() << ",\"fans\":[";
      for (int q=0; q<numbering.nprongs(); ++q) { if (q) out << ','; array(out,numbering.prong_letters[q]); }
      out << "],\"side_to\":[";
      for (int q=0; q<numbering.nprongs(); ++q) { if (q) out << ','; out << numbering.side_to(q); }
      out << "],\"punctured\":[";
      for (int q=0; q<numbering.nprongs(); ++q) { if (q) out << ','; out << (numbering.prong[q].punctured ? "true" : "false"); }
      out << "],\"arcs\":[";
      int branch=0;
      for (int f=0; f<source.foldings(); ++f) {
        auto target=source;
        traintracks::fold_map_data fm;
        if (!target.fold_with_map(f,fm)) continue;
        if (branch >= graph.foldings(v) || !(target == graph.traintrack(graph.target_vertex(v,branch))))
          throw std::logic_error("native fold target disagrees with automaton");
        if (branch) out << ',';
        out << "{\"target\":" << graph.target_vertex(v,branch) << ",\"vertex_images\":";
        array(out,fm.prong_image);
        out << ",\"edge_images\":[";
        const auto words=fm.to_freeauto();
        for (int e=1; e<=fm.nmain+fm.nsides; ++e) { if (e>1) out << ','; array(out,words[e]); }
        out << "]}";
        ++branch;
      }
      if (branch != graph.foldings(v)) throw std::logic_error("native arc count disagrees");
      out << "]}";
    }
    out << "]}\n";
    if (!out) throw std::runtime_error("output write failed");
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
