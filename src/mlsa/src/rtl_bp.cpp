#include "mlsa/rtl_bp.h"
#include "hier_mlsa.h"
#include "object.h"
#include "utl/Logger.h"

namespace mlsa {
    using odb::dbDatabase;
    using std::string;
    using std::unordered_map;
    using utl::Logger;

    BlobPlacer::BlobPlacer() = default;
    BlobPlacer::~BlobPlacer() = default;

    void BlobPlacer::init(sta::dbNetwork* network,
                      odb::dbDatabase* db,
                      sta::dbSta* sta,
                      utl::Logger* logger,
                      par::PartitionMgr* tritonpart) {
        hier_blob_
            = std::make_unique<HierBlobPlacement>(network, db, sta, logger, tritonpart);
    }

    bool BlobPlacer::place(const float tolerance,
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
        hier_blob_->SetClusterSizeTolerance(tolerance);
        hier_blob_->SetMaxNumLevel(max_num_level);
        hier_blob_->SetClusterSizeRatioPerLevel(coarsening_ratio);
        hier_blob_->SetNumBundledIOsPerBoundary(num_bundled_ios);
        hier_blob_->SetLargeNetThreshold(large_net_threshould);
        hier_blob_->SetSignatureNetThreshold(signature_net_threshold);
        hier_blob_->SetNumClustersBase(intended_num_clusters);
        hier_blob_->SetNumTiers(num_tiers);
        hier_blob_->SetTierResolution(tier_resolution);
        hier_blob_->SetWireLengthWeight(wire_length_weight);
        hier_blob_->SetVerticalWeight(vertical_weight);
        hier_blob_->SetAreaUpperBound(area_upper_bound);
        hier_blob_->SetSwapProb(swap_prob);
        hier_blob_->SetInitProb(init_prob);
        hier_blob_->SetMaxSteps(max_steps);
        hier_blob_->SetNumPerturbPerStep(num_perturb_per_step);
        hier_blob_->SetTimingPenalty(timing_penalty);
        hier_blob_->SetSeed(seed);

        // Run placer
        hier_blob_->BlobPlacement();

        return true;
    }
} // namespace mlsa