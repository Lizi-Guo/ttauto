// SPDX-License-Identifier: GPL-3.0-or-later
// The three-puncture builder uses a third singularity-enumeration slot.
// Run this test under AddressSanitizer to detect the former out-of-bounds write.
#include "traintracks/build.hpp"
#include "traintracks/traintrack.hpp"
#include "ttauto/ttfoldgraph.hpp"
#include "check.hpp"

int main()
{
  for (int repetition=0; repetition<32; ++repetition) {
    const auto seeds = traintracks::build_traintrack_list(3);
    CHECK(seeds.size() == 1);
    ttauto::ttfoldgraph<traintracks::traintrack> graph(seeds[0]);
    CHECK(graph.vertices() == 1);
    CHECK(graph.foldings() == 2);
    CHECK(graph.edges() == 2);
  }
}
