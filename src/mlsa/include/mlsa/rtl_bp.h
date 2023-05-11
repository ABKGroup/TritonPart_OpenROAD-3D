#pragma once

#include <memory>

namespace odb {
    class dbDatabase;
} // namespace odb

namespace sta {
    class dbNetwork;
    class dbSta;
} // namespace sta

namespace utl {
    class Logger;
}

namespace par {
    class PartitionMgr;
}

namespace mlsa {
    class HierBlobPlacement;

    class BlobPlacer {
        public:
            BlobPlacer();
            ~BlobPlacer();

            void init(sta::dbNetwork* network,
                      odb::dbDatabase* db,
                      sta::dbSta* sta,
                      utl::Logger* logger,
                      par::PartitionMgr* tritonpart);
            
            bool place(const float tolerance,
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
                       const unsigned seed);
        
        private:
            std::unique_ptr<HierBlobPlacement> hier_blob_;
    };
    
} // namespace mlsa