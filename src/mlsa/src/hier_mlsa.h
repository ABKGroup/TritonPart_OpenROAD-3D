#pragma once

#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include "object.h"

#include "db_sta/dbReadVerilog.hh"
#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "sta/Bfs.hh"
#include "sta/Graph.hh"
#include "sta/Liberty.hh"
#include "sta/Sta.hh"

namespace odb {
    class dbBlock;
    class dbBTerm;
    class dbDatabase;
    class dbITerm;
    class dbInst;
    class dbModule;
    class dbNet;
} // namespace odb

namespace sta {
    class dbNetwork;
    class dbSta;
    //class Instances;
    //class NetworkReader;
    //class Library;
    //class Port;
    //class Net;
} // namespace sta

namespace utl {
    class Logger;
}

namespace par {
    class PartitionMgr;
}

namespace mlsa {
    // Define all needed structures
    struct BundledNet;
    class Cluster;
    class Metrics;

    // ***************************************************************************
    // * Hierarchical Blob Placement - 3DIC
    // ***************************************************************************
    /**
     * @author Zhiang Wang, Jiantao Liu
     * @date  04/03/2023
     * @brief Support Multi-level Hierarchical clustering
     *        Support design with IO Pads
     *        Important factors:
     *        1) RTL information is identified based on logical hierarchy, connection 
     *           signature and dataflow.
     *        2) Connection signature is defined as the connection topology respect to 
     *           other clusters.
     *        3) Dataflow is defined based on sequential graph.
     *        4) SA is called based on multi-threaded go-with-the-winner manner
     *        5) Results refinement is based on the slidding window manner
     *        6) odb::dbIntProperty::create() is used to add cluster_id to all components
     */
    class HierBlobPlacement {
        public:
            HierBlobPlacement(sta::dbNetwork* network,
                              odb::dbDatabase* db,
                              sta::dbSta* sta,
                              utl::Logger* logger,
                              par::PartitionMgr* tritonpart
            );
            ~HierBlobPlacement();

            // Top Level Interface Function
            // Interface function for calling HierBlobPlacement --> Multi-level Clustering + GWW SA
            // Steps for the function
            // Step 1:
            //        Traverse the logical hierarchy, get all the statistics of each logical module
            // Step 2:
            //        Create Bundled pins and treat each bundled pin as a cluster without area
            // Step 3:
            //        Create physical hierarchy tree in a DFS manner, large cluster will be break by
            //        doing a bottom-up clustering.
            // Step 4:
            //        SA will be called for different levels in the physical hierarchy tree,
            //        with a multi-threaded go-with-the-winner method
            void BlobPlacement();

            // ***************************************************************************
            // * Hierarchical Clustering
            // ***************************************************************************
            void SetNumBundledIOsPerBoundary(int num_bundled_ios);
            void SetNumClustersBase(int num_clusters_base);
            void SetClusterSize(int max_num_inst,
                                int min_num_inst);
            void SetClusterSizeTolerance(float tolerance);
            void SetMaxNumLevel(int max_num_level);
            void SetClusterSizeRatioPerLevel(float coarsening_ratio);
            void SetLargeNetThreshold(int large_net_threshold);
            void SetSignatureNetThreshold(int signature_net_threshold);
            void SetTopKTimingPaths(int num_paths);
            void PrintThresholdParameters();

            // Set Parameters for Simulated Annealing
            void SetNumTiers(const int num_tiers);
            void SetTierResolution(const int tier_resolution);
            void SetWireLengthWeight(const float wire_length_weight);
            void SetVerticalWeight(const float vertical_weight);
            void SetAreaUpperBound(const float area_upper_bound);
            void SetSwapProb(const float swap_prob);
            void SetInitProb(const float init_prob);
            void SetMaxSteps(const int max_steps);
            void SetNumPerturbPerStep(const int num_perturb_per_step);
            void SetTimingPenalty(const float timing_penalty);
            void SetSeed(const unsigned seed);

        private:
            void SetDefaultThresholds();
            void CreateDataFlow();
            // void UpdateDataFlow();
            void DataFlowDFSIOPin(int parent,
                                  int idx,
                                  std::vector<std::set<odb::dbInst*>>& insts,
                                  std::map<int, odb::dbBTerm*>& io_pin_vertex,
                                  std::map<int, odb::dbInst*>& std_cell_vertex,
                                  std::vector<bool>& stop_flag_vec,
                                  std::vector<bool>& visited,
                                  std::vector<std::vector<int>>& vertices,
                                  std::vector<std::vector<int>>& hyperedges,
                                  bool backward_flag);

            // Traverse the logical hierarchy
            Metrics* ComputeMetrics(odb::dbModule* module);
            // Calculate metrics for the specified cluster
            // For now, we only supports StdCellCluster
            void SetClusterMetrics(Cluster* cluster);
            // Map IOs to Pads (if any)
            void MapIOPads();
            // Create bundled IOs as clusters (for designs with IO pins or pads)
            void CreateBundledIOs();

            // Create the timing paths for gathering the timing information
            // Get all critical timing paths
            void FindTimingPaths(Hypergraph& hypergraph,
                                 std::vector<TimingPath>& timing_paths,
                                 const int parent_cluster_id);
            
            // Create global timing information 
            void AddGlobalTimingInfo();
            
            // Set the corresponding timing factor
            void SetTimingFactor(float intend_timing_factor);

            // Update the cluster_id property of insts in the cluster
            // Need to call this before each level of SA runs
            // The physical hierarchy tree will be created in a pre-order DFS manner
            void MultiLevelCluster(Cluster* parent);
            void SetInstProperty(Cluster* cluster);
            void SetInstProperty(odb::dbModule* module,
                                 int cluster_id);
            void SetInstLocProperty(Cluster* cluster, float x_loc, float y_loc, float z_loc);
            void BreakCluster(Cluster* parent);
            void MergeClusters(std::vector<Cluster*>& candidate_clusters);
            void CalculateConnection();
            void CalculateCliqueConnection();
            void PrintConnection();
            void PrintClusters();
            void UpdateSubTree(Cluster* parent);

            // Simulated Annealing related functions
            void ActivateClusters(int level);
            void CleanSlidingWindowVariable();
            void ClearBundledNetMap();
            void CondenseBundledNet();
            void ExtractClusters(int level, Bbox bbox);
            void ExtractBoth(Bbox bbox, int level);
            
            // Large flat clusters(no chilren) will be broken by a bottom-up clustering
            void BreakLargeFlatCluster_TP(Cluster* cluster);
            void BreakLargeFlatCluster_BU(Cluster* cluster);
            void BestChoiceClustering(Hypergraph& hypergraph, 
                                      std::vector<Group>& groups,
                                      int desired_num_groups,
                                      float area_limit);
            void InitializePriorityQueue(Hypergraph& hypergraph,
                                         AccessiblePriorityQueue<PQElement>& storing_queue,
                                         std::vector<Group>& groups);
            void MergeTwoGroups(AccessiblePriorityQueue<PQElement>& input_queue,
                                std::set<int>& existing_groups,
                                std::vector<Group>& groups,
                                Hypergraph& hypergraph,
                                float area_limit);
            void UpdatePriorityQueue(AccessiblePriorityQueue<PQElement>& input_queue,
                                     int group_id,
                                     std::vector<Group>& groups,
                                     float area_limit,
                                     Hypergraph& hypergraph,
                                     std::set<int>& existing_groups
                                     );
            void RefineGrouping(std::vector<Group>& groups, 
                                std::set<int>& unique_groups, 
                                float limit, 
                                float number_fixed_vertices);
            void PrintPriorityQueue(AccessiblePriorityQueue<PQElement> storing_queue);
            void PrintGroupInfo(Group sel_group);
            void PrintTimingPathInfo(TimingPath sel_timing_path);

            // Print the Physical hierarchical tree in a DFS manner
            void PrintPhysicalHierarchyTree(Cluster* parent, int level);

            // This function writes out the final mapping between cells and the corresponding location
            void GenerateFinalOutputFile(std::string output_file_name, int intended_level);
            void WriteClusterCellInfo(std::ofstream& out_file, Cluster* sel_cluster);
            void WriteClusterCellInfo(std::ofstream& out_file, float X, float Y, float Z, odb::dbModule* module);
            
            // Global status checker
            void CalGlobalCost();
            
            // Force-directed Placer to generate guides
            // This maybe usefull if large amount of clusters in the partition is
            // obtained from a bottom-up clustering. We can use this to generate the 
            // initial placement result.
            // void FDInitialPlacement();

            // Variable Definition
            sta::dbNetwork* network_ = nullptr;
            odb::dbDatabase* db_ = nullptr;
            odb::dbBlock* block_ = nullptr;
            sta::dbSta* sta_ = nullptr;
            utl::Logger* logger_ = nullptr;
            par::PartitionMgr* tritonpart_ = nullptr;

            // Timing related variables
            int num_critical_timing_paths_ = 1000;    // Top k timing paths

            // Technology-related variables
            float dbu_ = 0.0;

            // Target utilization of the design
            float target_util_ = 0.60;

            // Design related variables
            float floorplan_lx_ = 0.0;
            float floorplan_ly_ = 0.0;
            float floorplan_ux_ = 0.0;
            float floorplan_uy_ = 0.0;

            // Dataflow parameters
            // Maximum number of FF distances in between the start and end
            int max_num_ff_dist_ = 5; 
            float dataflow_factor_ = 2.0;
            float dataflow_weight_ = 1;
            std::vector<std::pair<odb::dbBTerm*, std::vector<std::set<odb::dbInst*>>>>
                io_ffs_conn_map_;
            // TODO: Shall I use map or vector for the normal dbInst?
            std::vector<std::pair<odb::dbInst*, std::vector<std::set<odb::dbInst*>>>>
                non_io_ffs_conn_map_;

            // Store the statistics of the design
            Metrics* metrics_ = nullptr;
            // Store the metric for each hierarchical logical module
            std::map<const odb::dbModule*, Metrics*> logical_module_map_;

            // the virtual weight between std cell part and corresponding IO pins
            // to force them stay together
            float virtual_weight_ = 10.0;

            // User defined parameters
            // Number of bundled IOS on each boundary, also works for port propagation in SA
            int num_bundled_IOs_ = 4;
            int num_clusters_base_ = 0;
            int max_num_inst_base_ = 0;
            int min_num_inst_base_ = 0;

            // Improve the robustness of the clustering engine
            float tolerance_ = 0.1;
            int max_num_inst_ = 0;
            int min_num_inst_ = 0;

            // Multi-level support
            int max_num_level_ = 2;
            int level_ = 0;
            float coarsening_ratio_ = 10.0;

            // Area upper bound to ease the SA process
            float level_1_cluster_area_upper_bound_ = 0;

            // Timing parameters
            float timing_factor_ = 1.0;       // Timing factor used in the net weight calculation
            float maximum_clock_period_ = 0.0;

            // Variables related to connection signature
            // Minimum threshould for identifying two connected clusters
            int signature_net_threshold_ = 20;
            // Global nets are ignored
            int large_net_threshold_ = 100;

            // Maximum retry times during the clustering process
            int max_retry_numer_ = 15;

            // The conversion from database unit to micrometer
            // during the calculation will lose some accuracy.
            const float conversion_tolerance_ = 0.01;

            // At each level we will rebuild the net -> bundled net index map
            // This shall be cleared each time we build the clique based bundled net connection
            std::map<odb::dbNet*, std::vector<int>> net_bundled_net_map_;

            // Physical hierarchical tree
            // root cluster does not correspond to top design
            // it's a special node in the physical hierachy tree
            // Our physical hierarchy tree is as following:
            //             root_cluster_ (top design)
            //            *    *    *    *    *   *
            //           L    T     M1   M2    B    T  (L, T, B, T are clusters for
            //           bundled IOs)
            //                    *   *
            //                   M3    M4     (M1, M2, M3 and M4 are clusters for logical
            //                   modules
            int cluster_id_ = 0;
            Cluster* root_cluster_ = nullptr;       // root_cluster_id = 0
            std::map<int, Cluster*> cluster_map_;   // cluster_id, cluster

            // At each level the bundled IOs are children of root_cluster
            // Bundled IO (Pads)
            // In the bundled IO clusters, we don't store the ios in thier corresponding
            // clusters. However, we store the instances in thier corresponding clusters
            // Map IO pins to Pads (if any)
            std::map<odb::dbBTerm*, odb::dbInst*> io_pad_map_; // Bterm, instance

            // Sliding window related variables
            std::vector<BundledNet> sparse_nets_; // 1 to 1 mapping from dbNet to bundled Net
            std::vector<BundledNet> dense_nets_;             // Compressed Bundled Nets
            std::map<int, int> sparse_dense_map_;            // Map from sparse to dense bundled net connections
            std::vector<Vertex> movable_vertices_;           // Extracted vertices that are in the specified level and Bbox
            std::vector<Edge> selected_dense_nets_;          // Nets that will be traversed in the SA runs
            std::vector<Vertex> fixed_IOs_;                  // IOs that should be fixed in the SA runs
            
            // Sliding window structure
            int kernel_size_ = 2; // The 2D window equals to kernel_size * kernel_size

            // SA related parameters
            int num_tiers_ = 0;
            int tier_resolutions_ = 0;
            float wire_length_weight_ = 1.0;
            float vertical_connection_weight_ = 0.5;
            float area_upper_bound_ = 0.8;     // Ratio of allowed capacities in different partitions
            float swap_prob_ = 0.5;            // Prob for swapping two vertices in different partitions
            float init_prob_ = 0.8;            // Prob for calculating init_temperature
            int max_steps_ = 1000000;           // Num of steps per SA run
            int num_perturb_per_step_ = 10000; // Num perturbs per step
            float timing_penalty_ = 0;         // Timing penalty for cross-tier connections
            unsigned seed_;

            // Number of optimization passes
            int num_sa_passes_ = 1;

    };

} // namespace mlsa