%{
#include "mlsa/rtl_bp.h"

namespace ord {
// Define in OpenRoad.i
mlsa::BlobPlacer*
getBlobPlacer();
}

using ord::getBlobPlacer;
%}

%include "../../Exception.i"

%inline %{

namespace mlsa {

bool rtl_blob_placer_cmd(const float tolerance,
                         const int max_num_level,
                         const float coarsening_ratio,
                         const int num_bundled_ios,
                         const int large_net_threshould,
                         const int signature_net_threshold,
                         const int intended_num_clusters,
                         const int num_tiers,
                         const int tier_resolution,
                         const float wire_length_weight,
                         const float vertical_weight,
                         const float area_upper_bound,
                         const float swap_prob,
                         const float init_prob,
                         const int max_steps,
                         const int num_perturb_per_step,
                         const float timing_penalty,
                         const unsigned seed) {
    
    auto blob_placer = getBlobPlacer();
    return blob_placer->place(tolerance,
                              max_num_level,
                              coarsening_ratio,
                              num_bundled_ios,
                              large_net_threshould,
                              signature_net_threshold,
                              intended_num_clusters,
                              num_tiers,
                              tier_resolution,
                              wire_length_weight,
                              vertical_weight,
                              area_upper_bound,
                              swap_prob,
                              init_prob,
                              max_steps,
                              num_perturb_per_step,
                              timing_penalty,
                              seed);
}


} // namespace mlsa

%}