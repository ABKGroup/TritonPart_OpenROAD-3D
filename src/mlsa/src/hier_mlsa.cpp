#include "hier_mlsa.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <queue>
#include <thread>
#include <utility>
#include <random>
#include <string>

// #include "sta/ArcDelayCalc.hh"
// #include "sta/Bfs.hh"
// #include "sta/Corner.hh"
// #include "sta/DcalcAnalysisPt.hh"
// #include "sta/ExceptionPath.hh"
// #include "sta/FuncExpr.hh"
// #include "sta/Graph.hh"
// #include "sta/GraphDelayCalc.hh"
#include "sta/Liberty.hh"
// #include "sta/Network.hh"
// #include "sta/PathAnalysisPt.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
// #include "sta/PathRef.hh"
// #include "sta/PatternMatch.hh"
// #include "sta/PortDirection.hh"
#include "sta/Sdc.hh"
#include "sta/Search.hh"
// #include "sta/SearchPred.hh"
// #include "sta/Sequential.hh"
// #include "sta/Sta.hh"
// #include "sta/Units.hh"
#include "db_sta/dbNetwork.hh"
#include "object.h"
#include "MultilevelSACore.h"
#include "MultilevelSACore_New.h"
#include "odb/db.h"
#include "par/PartitionMgr.h"
#include "db_sta/dbNetwork.hh"
#include "par/PartitionMgr.h"
#include "utl/Logger.h"

namespace mlsa {
    // ***************************************************************************
    // * Hierarchical Blob Placement - 3DIC
    // ***************************************************************************
    HierBlobPlacement::~HierBlobPlacement() = default;

    // Constructors
    // TODO: Add timing extraction part
    HierBlobPlacement::HierBlobPlacement(sta::dbNetwork* network,
                                         odb::dbDatabase* db,
                                         sta::dbSta* sta,
                                         utl::Logger* logger,
                                         par::PartitionMgr* tritonpart
    ) {
        network_ = network;
        db_ = db;
        sta_ = sta;
        logger_ = logger;
        tritonpart_ = tritonpart;
    }

    // ***************************************************************************
    // * Hierarchical Clustering - Set parameters
    // ***************************************************************************
    void HierBlobPlacement::SetNumBundledIOsPerBoundary(int num_bundled_ios){
        num_bundled_IOs_ = num_bundled_ios;
    }

    void HierBlobPlacement::SetNumClustersBase(int num_clusters_base) {
        num_clusters_base_ = num_clusters_base;
    }

    void HierBlobPlacement::SetClusterSize(int max_num_inst,
                                           int min_num_inst
    ) {
        max_num_inst_ = max_num_inst;
        min_num_inst_ = min_num_inst;
    }

    void HierBlobPlacement::SetClusterSizeTolerance(float tolerance) {
        tolerance_ = tolerance;
    }

    void HierBlobPlacement::SetTopKTimingPaths(int num_paths) {
        num_critical_timing_paths_ = num_paths;
    }

    void HierBlobPlacement::SetMaxNumLevel(int max_num_level) {
        max_num_level_ = max_num_level;
    }

    void HierBlobPlacement::SetClusterSizeRatioPerLevel(float coarsening_ratio) {
        coarsening_ratio_ = coarsening_ratio;
    }

    void HierBlobPlacement::SetLargeNetThreshold(int large_net_threshold) {
        large_net_threshold_ = large_net_threshold;
    }

    void HierBlobPlacement::SetSignatureNetThreshold(int signature_net_threshold) {
        signature_net_threshold_ = signature_net_threshold;
    }

    void HierBlobPlacement::SetTimingFactor(float intend_timing_factor) {
        timing_factor_ = intend_timing_factor;
    }

    void HierBlobPlacement::PrintThresholdParameters() {
        logger_->report(
            "**************************Parameters*****************************\n"
        );
        logger_->report("Num_clusters_base: {}", num_clusters_base_);
        logger_->report("Max_num_inst_base: {}", max_num_inst_base_);
        logger_->report("Min_num_inst_base: {}", min_num_inst_base_);
        logger_->report("Max_num_level: {}", max_num_level_);
        logger_->report("Large_net_threshold: {}", large_net_threshold_);
        logger_->report("Signature_net_threshold: {}", signature_net_threshold_);
        logger_->report("Max_num_inst: {}", max_num_inst_);
        logger_->report("Min_num_inst: {}", min_num_inst_);
        logger_->report("Level 1 cluster area upper bound : {}", level_1_cluster_area_upper_bound_);
        logger_->report(
            "*****************************************************************\n"
        );
    }

    // Determine the default threshold parameters, if not specified by the user
    // Set the values for max and min number insts for
    // clusters in the first level
    // TODO: Do we want a number constraint or an area constraint?
    void HierBlobPlacement::SetDefaultThresholds() {
        // If user didn't specify the settings
        // The root node and base metrics shall be set up before
        // calling this function
        if (max_num_inst_base_ <= 0 || min_num_inst_base_ <= 0
            || num_clusters_base_ <= 0) {
            // 1) Default num cluster base value set to 500, if not set
            if (num_clusters_base_ <= 0) {
                num_clusters_base_ = 500;
            }
            
            // 2) Get the constraints in the max level first
            min_num_inst_base_ = 
                std::floor(metrics_->GetNumStdCell()
                           / (num_clusters_base_ * std::pow(coarsening_ratio_, max_num_level_ - 1)));
            
            // TODO: Need to specify better lower bound requirements
            if (min_num_inst_base_ <= 10) {
                min_num_inst_base_ = 10; // Lower bound of the minimum number of instances in one cluster
            }

            // Max num constraints is set to (1 + tolerance) times the minimum num constraints
            max_num_inst_base_ = min_num_inst_base_ * (1 + tolerance_);
        }

        // 
        // Set sizes for the first level based on coarsening_factor and the number 
        // of physical hierarchy levels
        //
        unsigned coarsening_factor = std::pow(coarsening_ratio_, max_num_level_ - 1);
        max_num_inst_base_ = max_num_inst_base_ * coarsening_factor;
        min_num_inst_base_ = min_num_inst_base_ * coarsening_factor;

        // Testing: Set area upper bound for merging
        level_1_cluster_area_upper_bound_ = 
            (((floorplan_ux_ - floorplan_lx_) * (floorplan_uy_ - floorplan_ly_)) / num_clusters_base_);
        
        // Testing
        // Print initial threshold settings --> parameters in the max level except the number of clusters
        PrintThresholdParameters();

        // Report all parameters
        logger_->report(
            "Max num level: {}, max_inst_level_1: {}, min_inst_level_1: {}",
            max_num_level_,
            max_num_inst_base_,
            min_num_inst_base_
        );
    }

    // ***************************************************************************
    // * TOP LEVEL INTERFACE
    // ***************************************************************************
    void HierBlobPlacement::BlobPlacement() {
        // Get the dabase information
        block_ = db_->getChip()->getBlock();
        dbu_ = db_->getTech()->getDbUnitsPerMicron();

        // Get Floorplan info
        // Die Boundary
        odb::Rect die_box = block_->getDieArea();
        floorplan_lx_ = DbuToMicron(die_box.xMin(), dbu_);
        floorplan_ly_ = DbuToMicron(die_box.yMin(), dbu_);
        floorplan_ux_ = DbuToMicron(die_box.xMax(), dbu_);
        floorplan_uy_ = DbuToMicron(die_box.yMax(), dbu_); 

        // Core Boundary
        odb::Rect core_box = block_->getCoreArea();
        float core_lx = DbuToMicron(core_box.xMin(), dbu_);
        float core_ly = DbuToMicron(core_box.yMin(), dbu_);
        float core_ux = DbuToMicron(core_box.xMax(), dbu_);
        float core_uy = DbuToMicron(core_box.yMax(), dbu_);

        logger_->report(
            "Floorplan Outline: ({}, {}) ({}, {}),  Core Outline: ({}, {}) ({}, {})",
            floorplan_lx_,
            floorplan_ly_,
            floorplan_ux_,
            floorplan_uy_,
            core_lx,
            core_ly,
            core_ux,
            core_uy);

        // Calculate the metrics for dbModules
        metrics_ = ComputeMetrics(block_->getTopModule());
        float core_area = (core_ux - core_lx) * (core_uy - core_ly);
        float util
            = (metrics_->GetStdCellArea() / core_area);
        
        // Report Design information
        logger_->report(
            "Finished Traversing logical hierarchy\n"
            "\tNumber of std cell instances : {}\n"
            "\tArea of std cell instances : {:.2f}\n"
            "\tDesign utilization : {:.2f}\n",
            metrics_->GetNumStdCell(),
            metrics_->GetStdCellArea(),
            util
        );

        // Set default thresholds before running 
        SetDefaultThresholds();
        
        // Initialize the physical hierarchy tree
        // Create the root cluster
        cluster_id_ = 0;
        root_cluster_ = new Cluster(cluster_id_, std::string("root"), logger_);
        // Add the top module design metrics to the root cluster
        root_cluster_->AddDbModule(block_->getTopModule());
        root_cluster_->SetMetrics(*metrics_);
        cluster_map_[cluster_id_++] = root_cluster_;
        // Add cluster_id as the attribute to each instance
        for (auto inst: block_->getInsts()) {
            odb::dbIntProperty::create(inst, "cluster_id", cluster_id_);
        }

        // Model each bundled IO as a cluster under the root node
        // cluster_id: 1 to num_bundled_IOs x 4 are reserved for bundled IOs
        // Sequence order: Left -> Top -> Right -> Bottom
        // For multi-level SA port clusters will be recreated and propagate to 
        // the new canvas.
        // IOs will be mapped to pads, if any
        MapIOPads();
        CreateBundledIOs();

        // Testing
        // logger_->report("[DEBUG] Print Physical Hierarchy after generating the IO clusters");
        PrintPhysicalHierarchyTree(root_cluster_, 0);

        // Create Dataflow Information!
        // Only needs to be created once in the entire execution of HierBlobPlacement
        // logger_->report("\nCreate Dataflow Graph");
        // CreateDataFlow();

        // Create the needed physical hierarchy tree in a pre-order DFS manner
        logger_->report("\nPerform Hierarchical Clustering");

        // For now there is no macros allowed in the design pure logic
        if (metrics_->GetNumStdCell() == 0) {
            // ERROR
            logger_->report("[ERROR] This design has no standard cells! WRONG!");
        } else {
            // Conduct the multilevel clustersing 
            MultiLevelCluster(root_cluster_);
        }

        // Print the physical hierarchy tree
        logger_->report("\n Print Physical Hierarchy Tree");
        PrintPhysicalHierarchyTree(root_cluster_, 0);

        // Steps to successfully call SA
        int current_level = 0;                 // We start from the root
        
        // Define the sliding window structure
        while (current_level < max_num_level_) {
            // Increase the level counter
            current_level++;

            // Clean all variables at the beginning
            CleanSlidingWindowVariable();

            // Step 1 : Activate clusters in the selected level
            ActivateClusters(current_level);

            // Step 2 : Calculate Clique Connections (1 to 1 mapping)
            CalculateCliqueConnection();

            // Step 3 : Add timing information to the sparse nets
            AddGlobalTimingInfo();
            
            // Step 4 : Condense nets into bundled nets as input to SA
            CondenseBundledNet();

            // First level is different
            if (current_level == 1) {
                // Clean the storing structure
                movable_vertices_.clear();
                selected_dense_nets_.clear();

                // Get the bounding box
                Bbox bbox;
                bbox.lx = floorplan_lx_;
                bbox.ly = floorplan_ly_;
                bbox.ux = floorplan_ux_;
                bbox.uy = floorplan_uy_;

                // Extract all vertices
                // ExtractClusters(current_level, bbox);

                // New extractor
                ExtractBoth(bbox, current_level);

                logger_->report("[INFO]: {} vertices extracted from the physical hierarchy tree", movable_vertices_.size());

                // Testing
                // for (auto sel_vertex : movable_vertices_) {
                //     sel_vertex.PrintVertexInfo();
                // }

                // Conduct SA runs
                // Create the needed SA worker
                // SimulatedAnnealingCore* sa = new SimulatedAnnealingCore(
                //     movable_vertices_,
                //     area_upper_bound_,
                //     num_tiers_,
                //     tier_resolutions_,
                //     core_ux - core_lx,
                //     core_uy - core_ly,
                //     std::make_pair(core_lx, core_ly),
                //     wire_length_weight_,
                //     vertical_connection_weight_,
                //     timing_penalty_,
                //     swap_prob_,
                //     init_prob_,
                //     max_steps_,
                //     num_perturb_per_step_,
                //     seed_
                // );

                // New SA worker
                SimulatedAnnealingCore_New* sa = new SimulatedAnnealingCore_New(
                    movable_vertices_,
                    selected_dense_nets_,
                    area_upper_bound_,
                    num_tiers_,
                    tier_resolutions_,
                    core_ux - core_lx,
                    core_uy - core_ly,
                    std::make_pair(core_lx, core_ly),
                    wire_length_weight_,
                    vertical_connection_weight_,
                    timing_penalty_,
                    swap_prob_,
                    init_prob_,
                    max_steps_,
                    num_perturb_per_step_,
                    seed_
                );

                // Initialize the SA
                sa->InitializeSA(false);
                sa->RunSA();

                // Get the results from SA
                movable_vertices_.clear();
                sa->GetVertices(movable_vertices_);

                // TODO: Add a function to calculate the overall Cost
                // This is used to validate whether the sliding window is optimizing the global cost
                // or make the results worse

                // Update all location for related clusters
                for (auto sel_vertex : movable_vertices_) {
                    int sel_cluster_id = sel_vertex.GetVertexClusterId();

                    cluster_map_[sel_cluster_id]->SetX(sel_vertex.GetXCord());
                    cluster_map_[sel_cluster_id]->SetY(sel_vertex.GetYCord());
                    cluster_map_[sel_cluster_id]->SetZ(sel_vertex.GetZCord());
                }

                // Global Status
                CalGlobalCost();

            } else {
                for (int i = 0; i < num_sa_passes_; i++) {

                    // Testing
                    logger_->report("[DEBUG] Core Bbox: {}, {}, {}, {}",
                                    core_lx,
                                    core_ly,
                                    core_ux,
                                    core_uy);

                    // Define the sliding window structure
                    for (int i = 0; i < (tier_resolutions_ - kernel_size_ + 1); i++) {
                        for (int j = 0; j < (tier_resolutions_ - kernel_size_ + 1); j++) {
                            // Clear the storing structure
                            movable_vertices_.clear();
                            selected_dense_nets_.clear();

                            // Get the bounding box
                            Bbox bbox;
                            bbox.ux = ((core_ux - core_lx) / tier_resolutions_) * (j + kernel_size_) + core_lx;
                            bbox.uy = (core_uy - core_ly) - ((core_uy - core_ly) / tier_resolutions_) * i + core_ly;
                            bbox.lx = ((core_ux - core_lx) / tier_resolutions_) * j + core_lx;
                            bbox.ly = (core_uy - core_ly) - ((core_uy - core_ly) / tier_resolutions_) * (i + kernel_size_) + core_ly;

                            // Testing
                            logger_->report("[DEBUG] Bbox: ({}, {}, {}, {})",
                                            bbox.lx,
                                            bbox.ly,
                                            bbox.ux,
                                            bbox.uy);

                            // Extract all vertices
                            // ExtractClusters(current_level, bbox);

                            // New Extractor
                            ExtractBoth(bbox, current_level);

                            // Testing
                            logger_->report("\t[DEBUG] {} vertices extracted.", movable_vertices_.size());

                            // SimulatedAnnealingCore* sa = new SimulatedAnnealingCore(
                            //     movable_vertices_,
                            //     area_upper_bound_,
                            //     num_tiers_,
                            //     tier_resolutions_,
                            //     bbox.ux - bbox.lx,
                            //     bbox.uy - bbox.ly,
                            //     std::make_pair(bbox.lx, bbox.ly),
                            //     wire_length_weight_,
                            //     vertical_connection_weight_,
                            //     timing_penalty_,
                            //     swap_prob_,
                            //     init_prob_,
                            //     max_steps_,
                            //     num_perturb_per_step_,
                            //     seed_
                            // );

                            // New SA core
                            SimulatedAnnealingCore_New* sa = new SimulatedAnnealingCore_New(
                                movable_vertices_,
                                selected_dense_nets_,
                                area_upper_bound_,
                                num_tiers_,
                                tier_resolutions_,
                                bbox.ux - bbox.lx,
                                bbox.uy - bbox.ly,
                                std::make_pair(bbox.lx, bbox.ly),
                                wire_length_weight_,
                                vertical_connection_weight_,
                                timing_penalty_,
                                swap_prob_,
                                init_prob_,
                                max_steps_,
                                num_perturb_per_step_,
                                seed_
                            );

                            // Initialize the SA
                            sa->InitializeSA(false);
                            sa->RunSA();

                            // Get the results from SA
                            movable_vertices_.clear();
                            sa->GetVertices(movable_vertices_);

                            // Update all location for related clusters
                            for (auto sel_vertex : movable_vertices_) {
                                int sel_cluster_id = sel_vertex.GetVertexClusterId();

                                cluster_map_[sel_cluster_id]->SetX(sel_vertex.GetXCord());
                                cluster_map_[sel_cluster_id]->SetY(sel_vertex.GetYCord());
                                cluster_map_[sel_cluster_id]->SetZ(sel_vertex.GetZCord());
                            }

                            // Calculate the global cost
                            CalGlobalCost();
                        }
                    } 
                }
            }
        }

        // Clean the memory to avoid memory leakage
        // Release all the pointers

        // Print final physical hierarchy tree
        PrintPhysicalHierarchyTree(root_cluster_, 0);

        // Write out the all the results
        GenerateFinalOutputFile("final_output", max_num_level_);

    }

    // ***************************************************************************
    // * DATAFLOW
    // ***************************************************************************
    // Create Dataflow Information
    // Model each std cell instances, IO pin as vertices
    // For now we only support StdCellClusters
    // TODO: Store connection information to map timing information
    void HierBlobPlacement::CreateDataFlow() {
        // Check whether the dataflow creation is needed or not
        if (max_num_ff_dist_ <= 0) {
            return;
        }

        // Create vertex id property for std cell and IO pin
        // Macro is not supported for now
        std::map<int, odb::dbBTerm*> io_pin_vertex;
        std::map<int, odb::dbInst*> std_cell_vertex;

        std::vector<bool> stop_flag_vec;

        // Assign vertex_id property of each Block Terminal
        // All boundary terminals are marked as sequential stopping points
        for (odb::dbBTerm* term : block_->getBTerms()) {
            odb::dbIntProperty::create(term, "vertex_id", stop_flag_vec.size());
            io_pin_vertex[stop_flag_vec.size()] = term;
            stop_flag_vec.push_back(true);
        }

        // Assign vertex_id for different instances
        for (auto inst : block_->getInsts()) {
            const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
            // Skip if physical only
            if (liberty_cell == nullptr) {
                continue;
            }
            odb::dbMaster* master = inst->getMaster();
            // Check if the instance is a Pad, Cover or a block
            // At this step we skip these instances
            if (master->isPad() || master->isCover() || master->isBlock()) {
                continue;
            }

            // Mark all sequential instances
            odb::dbIntProperty::create(inst, "vertex_id", stop_flag_vec.size());
            std_cell_vertex[stop_flag_vec.size()] = inst;
            if (liberty_cell->hasSequentials()) {
                stop_flag_vec.push_back(true);
            } else {
                stop_flag_vec.push_back(false);
            }
        }

        // Number of vertices will be the number of boundary pins + 
        // number of logical std cells
        // Testing
        // std::cout << "[DEBUG] Number of vertices: " << stop_flag_vec.size() << std::endl;

        // Create Directed Hypergraph
        std::vector<std::vector<int>> vertices(stop_flag_vec.size());
        std::vector<std::vector<int>> backward_vertices(stop_flag_vec.size());
        std::vector<std::vector<int>> hyperedges;  // dircted hypergraph

        // Traverse the netlist
        for (odb::dbNet* net: block_->getNets()) {
            // Ignore all power nets, not modeled for the initial timing driven process
            if (net->getSigType().isSupply()) {
                continue;
            }

            // Driver vertex id
            int driver_id = -1;
            std::set<int> loads_id; // load vertices ids
            bool pad_flag = false;
            // Check the connected instances
            for (odb::dbITerm* iterm : net->getITerms()) {
                odb::dbInst* inst = iterm->getInst();
                
                // Create the liberty cell
                const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
                if (liberty_cell == nullptr) {
                    continue;
                }
                // Get the master cell of the instance 
                odb::dbMaster* master = inst->getMaster();

                // Check if the instance is a Pad, Cover or empty block
                // Nets connecting all above instance types are ignored
                if (master->isPad() || master->isCover()) {
                    pad_flag = true;
                    break;
                }

                int vertex_id = -1;
                // For now we assume there is no macros in the design,
                // all blocks inside the top block will be ignored for now,
                // they shouldn't exist for current design settings
                vertex_id = odb::dbIntProperty::find(inst, "vertex_id")->getValue();

                // If this is output --> driver found!
                if (iterm->getIoType() == odb::dbIoType::OUTPUT) {
                    driver_id = vertex_id;
                } else {
                    loads_id.insert(vertex_id);
                }
            }

            if (pad_flag) {
                continue; // Nets with pads shall be ignored
            }

            // Check the connected IO pins of the net
            for (odb::dbBTerm* bterm : net->getBTerms()) {
                const int vertex_id
                        = odb::dbIntProperty::find(bterm, "vertex_id")->getValue();

                if (bterm->getIoType() == odb::dbIoType::INPUT) {
                    driver_id = vertex_id;
                } else {
                    loads_id.insert(vertex_id);
                }
            }

            // Skip high fanout nets or nets that do not have valid driver or loads
            if (driver_id < 0 || loads_id.size() < 1 || loads_id.size() > large_net_threshold_) {
                continue;
            }

            // Create the hyperedge
            std::vector<int> hyperedge{driver_id};
            for (auto& load: loads_id) {
                if (load != driver_id) {
                    hyperedge.push_back(load);
                }
            }
            // Indirect mapping, the driver_id will be mapped to the corresponding
            // hyperedge(represented by the size of hyperedges)
            vertices[driver_id].push_back(hyperedges.size());
            for (int i = 1; i < hyperedge.size(); i++) {
                // Mapping from load to source
                backward_vertices[hyperedge[i]].push_back(hyperedges.size());
            }
            // Add the hyperedge to the hyperedges structure
            hyperedges.push_back(hyperedge);
        } // End net traversal

        // Testing
        // std::cout << "[DEBUG] Directed Hypergraph Created" << std::endl;

        // Traverse hypergraph to build dataflow
        // Start with IO pins
        for (auto [src, src_pin] : io_pin_vertex) {
            int idx = 0;
            std::vector<bool> visited(vertices.size(), false);
            std::vector<std::set<odb::dbInst*>> insts(max_num_ff_dist_);
            DataFlowDFSIOPin(src,
                             idx,
                             insts,
                             io_pin_vertex,
                             std_cell_vertex,
                             stop_flag_vec,
                             visited,
                             vertices,
                             hyperedges,
                             false
                             );
            DataFlowDFSIOPin(src,
                             idx,
                             insts,
                             io_pin_vertex,
                             std_cell_vertex,
                             stop_flag_vec,
                             visited,
                             vertices,
                             hyperedges,
                             true);
            io_ffs_conn_map_.push_back(
                std::pair<odb::dbBTerm*, std::vector<std::set<odb::dbInst*>>>(src_pin, insts)
            );
        } 
    }

    // Forward or Backward DFS search to find sequential paths from/to IO pins based on hop
    // counts to IO pins
    void HierBlobPlacement::DataFlowDFSIOPin(int parent,
                                             int idx,
                                             std::vector<std::set<odb::dbInst*>>& insts,
                                             std::map<int, odb::dbBTerm*>& io_pin_vertex,
                                             std::map<int, odb::dbInst*>& std_cell_vertex,
                                             std::vector<bool>& stop_flag_vec,
                                             std::vector<bool>& visited,
                                             std::vector<std::vector<int>>& vertices,
                                             std::vector<std::vector<int>>& hyperedges,
                                             bool backward_flag) {
        visited[parent] = true;
        // If this is a sequential stop element
        if (stop_flag_vec[parent]) {
            if (parent < io_pin_vertex.size()) {
                ; // IO Pin to IO Pin connections are not considered for now
            } else if (parent < io_pin_vertex.size() + std_cell_vertex.size()) {
                insts[idx].insert(std_cell_vertex[parent]);
            }
            // Currently we don't support Macro pins
            idx++;
        }

        if (idx >= max_num_ff_dist_) {
            // Reached maximum threshold
            return;
        }

        if (!backward_flag) {
            for (auto& hyperedge : vertices[parent]) {
                for (auto& vertex : hyperedges[hyperedge]) {
                    // Pin to Pin connection are ignored
                    if (visited[vertex] || vertex < io_pin_vertex.size()) {
                        continue;
                    }
                    DataFlowDFSIOPin(vertex,
                                     idx,
                                     insts,
                                     io_pin_vertex,
                                     std_cell_vertex,
                                     stop_flag_vec,
                                     visited,
                                     vertices,
                                     hyperedges,
                                     backward_flag
                                     );

                }
            } // Finish hyperedge
        } else {
            for (auto& hyperedge : vertices[parent]) {
                // Get the driver
                const int vertex = hyperedges[hyperedge][0];
                // Skip pin to pin connection
                if (visited[vertex] || vertex < io_pin_vertex.size()) {
                    continue;
                }
                DataFlowDFSIOPin(vertex,
                                 idx,
                                 insts,
                                 io_pin_vertex,
                                 std_cell_vertex,
                                 stop_flag_vec,
                                 visited,
                                 vertices,
                                 hyperedges,
                                 backward_flag);
            } // Finish the hyperedge
        } // Finish current vertex
    }


    // ***************************************************************************
    // * Hierarchical Clustering - Main Functionalities (private)
    // ***************************************************************************
    // Handle IOs, map IOs to Pads for designs with IO pads. If no
    // IO pads exist, do nothing
    void HierBlobPlacement::MapIOPads() {
        // Check if the design has IO pads
        bool is_pad_design = false;
        for (auto inst : block_->getInsts()) {
            if (inst->getMaster()->isPad()) {
                is_pad_design = true;
                break;
            }
        }
        
        // If there is no pad in the design
        if (!is_pad_design) {
            // Testing
            std::cout << "[DEBUG] No Pad found, quite IO -- PAD mapping process" << std::endl;

            return;
        }

        for (odb::dbNet* net : block_->getNets()) {
            // Check whether the net is connecting to the block terminal
            if (net->getBTerms().size() == 0) {
                continue;
            }

            // If the IO pads exists, and the selected net
            // contains the IO pin and pad
            for (odb::dbBTerm* bterm : net->getBTerms()) {
                for (odb::dbITerm* iterm : net->getITerms()) {
                    odb::dbInst* inst = iterm->getInst();
                    io_pad_map_[bterm] = inst;
                }
            }
        }
    }

    // Create bundled IOs in the following order: L, T, R, B
    // num_bundled_IOs_ x 4 clusters will be created, all pins
    // in the corresponding region will be grouped
    void HierBlobPlacement::CreateBundledIOs() {
        // Convert from micron to dbu
        floorplan_lx_ = MicronToDbu(floorplan_lx_, dbu_);
        floorplan_ly_ = MicronToDbu(floorplan_ly_, dbu_);
        floorplan_ux_ = MicronToDbu(floorplan_ux_, dbu_);
        floorplan_uy_ = MicronToDbu(floorplan_uy_, dbu_);

        // Divide the canvas boundary into different regions
        odb::Rect die_box = block_->getCoreArea();
        int core_lx = die_box.xMin();
        int core_ly = die_box.yMin();
        int core_ux = die_box.xMax();
        int core_uy = die_box.yMax();
        const int x_base = (floorplan_ux_ - floorplan_lx_) / num_bundled_IOs_;
        const int y_base = (floorplan_uy_ - floorplan_ly_) / num_bundled_IOs_;

        // The base cluster_id shall be 1 at the beginning
        int cluster_id_base = cluster_id_;

        // Create all bundled IO clusters
        std::vector<std::string> prefix_vec;
        prefix_vec.push_back(std::string("L"));
        prefix_vec.push_back(std::string("T"));
        prefix_vec.push_back(std::string("R"));
        prefix_vec.push_back(std::string("B"));
        std::map<int, bool> cluster_io_map;
        for (int i = 0; i < 4; i++) {
            // Four boundaries (Left, Top, Right, Bottom)
            for (int j = 0; j < num_bundled_IOs_; j++) {
                // Create the IO cluster
                std::string cluster_name = prefix_vec[i] + std::to_string(j);
                Cluster* cluster = new Cluster(cluster_id_, cluster_name, logger_);

                // Add the IO cluster to the root cluster
                root_cluster_->AddChild(cluster);
                cluster->SetParent(root_cluster_);
                cluster_io_map[cluster_id_] = false;
                cluster_map_[cluster_id_++] = cluster;

                // Calculate coordinates, conversion needed (from dbu to microns)
                // x, y will be set to the center location of the corresponding segment
                // We assume all IOs are in the bottom tier
                int x = 0;
                int y = 0;
                int z = 0;
                int width = 0;
                int height = 0;
                if (i == 0) {
                    // Left Boundary
                    x = floorplan_lx_;
                    y = floorplan_ly_ + (y_base / 2) + (y_base * j);
                    height = y_base;
                } else if (i == 1) {
                    // Top Boundary
                    x = floorplan_lx_ + (x_base / 2) + (x_base * j);
                    y = floorplan_uy_;
                    width = x_base;
                } else if (i == 2) {
                    // Right Boundary
                    x = floorplan_ux_;
                    y = floorplan_uy_ - ((y_base / 2) + (y_base * j));
                    height = y_base;
                } else {
                    // Bottom Boundary
                    x = floorplan_ux_ - ((x_base / 2) + (x_base * j));
                    y = floorplan_ly_;
                    width = x_base;
                }

                // Set the IO cluster flag
                cluster->SetIOClusterFlag(
                    DbuToMicron(x, dbu_),
                    DbuToMicron(y, dbu_),
                    DbuToMicron(z, dbu_),
                    DbuToMicron(width, dbu_),
                    DbuToMicron(height, dbu_)
                );
            }
        }

        // Map all block terminals to bundled IOs
        for (auto term: block_->getBTerms()) {
            int lx = std::numeric_limits<int>::max();
            int ly = std::numeric_limits<int>::max();
            int ux = 0;
            int uy = 0;

            // If the design has IO pads, these block terminals
            // will not have block pins.
            // Otherwise, the design will have IO pins
            for (const auto pin : term->getBPins()) {
                for (const auto box : pin->getBoxes()) {
                    lx = std::min(lx, box->xMin());
                    ly = std::min(ly, box->yMin());
                    ux = std::max(ux, box->xMax());
                    uy = std::max(uy, box->yMax());
                }
            }

            // Remove all power pins
            if (term->getSigType().isSupply()) {
                continue;
            }

            // If the term has a connected pad, get the bbox from the pad inst
            if (io_pad_map_.find(term) != io_pad_map_.end()) {
                lx = io_pad_map_[term]->getBBox()->xMin();
                ly = io_pad_map_[term]->getBBox()->yMin();
                ux = io_pad_map_[term]->getBBox()->xMax();
                uy = io_pad_map_[term]->getBBox()->yMax();

                // Constrain the location to be inside the canvas
                if (lx <= core_lx) {
                    lx = floorplan_lx_;
                }
                if (ly <= core_ly) {
                    ly = floorplan_ly_;
                }
                if (ux >= core_ux) {
                    ux = floorplan_ux_;
                }
                if (uy >= core_uy) {
                    uy = floorplan_uy_;
                }
            }

            // Calculate the cluster id based on the location of the IO pins / pads
            int cluster_id = -1;
            if (lx <= floorplan_lx_) {
                // The IO is on the left boundary
                cluster_id = cluster_id_base
                             + std::floor(((ly + uy) / 2.0 - floorplan_ly_) / y_base);
            } else if (uy >= floorplan_uy_) {
                // The IO is on the top boundary
                cluster_id = cluster_id_base + num_bundled_IOs_
                             + std::floor(((lx + ux) / 2.0 - floorplan_lx_) / x_base);
            } else if (ux >= floorplan_ux_) {
                // The IO is on the right boundary
                cluster_id = cluster_id_base + num_bundled_IOs_ * 2
                             + std::floor((floorplan_uy_ - (ly + uy) / 2.0) / y_base);
            } else if (ly <= floorplan_ly_) {
                // The IO is on the bottom boundary
                cluster_id = cluster_id_base + num_bundled_IOs_ * 3
                             + std::floor((floorplan_ux_ - (lx + ux) / 2.0) / x_base);
            }

            // Check the validity of the obtained IO cluster_id
            if (cluster_id == -1) {
                // Testing
                std::cout << "[ERROR] Floorplan may not be initialized, IO Location error for "
                          << term->getName() << std::endl;

                // In this case, put them in the middle of the left boundary
                logger_->report("[DEBUG] Add the non-initialized pin to the middle of the left boundary");
                cluster_id = 2;
                odb::dbIntProperty::create(term, "cluster_id", cluster_id);
            } else {
                // The location of the IO pin is valid
                odb::dbIntProperty::create(term, "cluster_id", cluster_id);
            }

            // Add the cluster to map
            cluster_io_map[cluster_id] = true;
        }

        // Convert from dbu to micron
        floorplan_lx_ = DbuToMicron(floorplan_lx_, dbu_);
        floorplan_ly_ = DbuToMicron(floorplan_ly_, dbu_);
        floorplan_ux_ = DbuToMicron(floorplan_ux_, dbu_);
        floorplan_uy_ = DbuToMicron(floorplan_uy_, dbu_);

        // Delete the IO clusters that do not have any pins assigned to terminals (empty IO cluster)
        for (auto& [cluster_id, flag] : cluster_io_map) {
            if (!flag) {
                // Testing output
                std::cout << "[WARNING] Remove IO Cluster with no pins: "
                          << cluster_map_[cluster_id]->GetName() << ", id: "
                          << cluster_id << std::endl;
                cluster_map_[cluster_id]->GetParent()->RemoveChild(
                    cluster_map_[cluster_id]
                );
                
                // Memory management
                delete cluster_map_[cluster_id];
                cluster_map_.erase(cluster_id);
            }
        }
    }

    // Create the needed physical hierarchy tree in a pre-order DFS manner
    // Recursive call for creating the physical hierarchy tree.
    // Important: all leaf std cells in the parent cluster will be passed to 
    // child clusters as for the multi-level blob placement
    void HierBlobPlacement::MultiLevelCluster(Cluster* parent) {
        // Variable Defintion
        bool force_split = false;
        if (level_ == 0) {
            // Check whether the root cluster is below the maxsize of the leaf cluster
            // Force the creation of child clusters.
            const int leaf_cluster_size
                = max_num_inst_base_ / std::pow(coarsening_ratio_, max_num_level_ - 1);
            if (parent->GetNumStdCell() < leaf_cluster_size) force_split = true;

            // Testing
            // std::cout << "[DEBUG] Force split : root cluster size: "
            //           << parent->GetNumStdCell() << ", leaf cluster size: "
            //           << leaf_cluster_size << std::endl;
        }

        if (level_ >= max_num_level_) {
            // Testing
            return;
        }

        level_++;

        // Testing
        // std::cout << "[DEBUG] Parent: " << parent->GetName() << ", level: "
        //           << level_ << ", num std cells: " << parent->GetNumStdCell() << std::endl;

        // Compute the size constraints for the current level
        max_num_inst_ = max_num_inst_base_ / std::pow(coarsening_ratio_, level_ - 1);
        min_num_inst_ = min_num_inst_base_ / std::pow(coarsening_ratio_, level_ - 1);

        // Due to the accuracy loss during conversion
        // TODO: Separate the tolerance for clustering and set default
        max_num_inst_ = max_num_inst_ * (1 + tolerance_);
        min_num_inst_ = min_num_inst_ * (1 - tolerance_);

        if (min_num_inst_ <= 0) {
            min_num_inst_ = 100;
            max_num_inst_ = min_num_inst_ * (1 + tolerance_);
        }

        // Testing
        std::cout << "[DEBUG] level: "
                  << level_ << ", min_num_inst_: " << min_num_inst_
                  << ", max_num_inst_: " << max_num_inst_ << std::endl;

        if (force_split || (parent->GetNumStdCell() > max_num_inst_)) {
            // Break the parent cluster into children clusters
            BreakCluster(parent); 
            // Update the subtree to the Physical Hierarchy Tree
            UpdateSubTree(parent);
            
            for (auto& child : parent->GetChildren()) {
                SetInstProperty(child);
            }

            // Recursively call the clustering function
            for (auto& child : parent->GetChildren()) {
                // Testing, Print the basic information of the children cluster
                std::cout << "[DEBUG] Child Cluster: " << child->GetName()<< std::endl;

                MultiLevelCluster(child);
            }
        } else {
            MultiLevelCluster(parent);
        }

        // Back to the last level, reset cluster_id property after building 
        // the physical hierarchy tree
        SetInstProperty(parent);
        level_--;
    }

    // Update the cluster_id property of insts in the cluster
    void HierBlobPlacement::SetInstProperty(Cluster* cluster){
        int cluster_id = cluster->GetId();
        
        // For now we only support StdCellCluster
        ClusterType cluster_type = cluster->GetClusterType();

        // Set property for all std cell instances
        for (auto& inst : cluster->GetLeafStdCells()) {
            odb::dbIntProperty::find(inst, "cluster_id")->setValue(cluster_id);
        }

        // Set int property for submodules
        for (auto& module : cluster->GetDbModules()) {
            SetInstProperty(module, cluster_id);
        }
    }

    // Update the cluster_id property of instances in a logical module
    void HierBlobPlacement::SetInstProperty(odb::dbModule* module,
                                            int cluster_id) {
        // For now we don't include macros in the design
        for (odb::dbInst* inst : module->getInsts()) {
            odb::dbMaster* master = inst->getMaster();

            // Check whether the instance is a PAD, Cover(Filler) or empty block
            if (master->isPad() || master->isCover() || master->isBlock()) {
                continue;
            }
            odb::dbIntProperty::find(inst, "cluster_id")->setValue(cluster_id);
        }

        // Set properties for all submodules
        for (odb::dbModInst* inst : module->getChildren()) {
            SetInstProperty(inst->getMaster(), cluster_id);
        }
    }

    // Update the location of a cluster and all its child clusters
    void HierBlobPlacement::SetInstLocProperty(Cluster* cluster, float x_loc, float y_loc, float z_loc) {
        cluster->SetX(x_loc);
        cluster->SetY(y_loc);
        cluster->SetZ(z_loc);

        // Set the location of all included child clusters
        if (cluster->GetChildren().size() == 0) {
            return;
        } else {
            for (auto child_cluster : cluster->GetChildren()) {
                SetInstLocProperty(child_cluster, x_loc, y_loc, z_loc);
            }
        }
    }

    // Break the parent cluster into multiple children clusters
    // The parent cluster will be expaneded to a subtree based on the logical hierarchy (if exist)
    // in a DFS manner. During the expansion process,
    // Small clusters will be merged if in the same logical hierarchy.
    // We provide two ways to break the large flat cluster(no child in the logical hierarchy)
    // 1) Recursivly call TritonPart to do bi-partition
    // 2) Directly do a bottom up clustering to generate the desired clusters
    void HierBlobPlacement::BreakCluster(Cluster* parent) {
        // Testing
        std::cout << "[DEBUG] Dissolve Cluster: " << parent->GetName() << std::endl;

        // Step 1: Break the parent cluster
        // Possible scenarios
        // (a) parent is an empty cluster
        // (b) parent is a single logical module
        // (c) parent is a cluster with both logical modules and leaf cells
        if ((parent->GetLeafStdCells().size() == 0)
             && (parent->GetDbModules().size() == 0)) {
            // (a) Empty cluster encountered
            // This should not happen

            // Testing
            logger_->report("[ERROR] Break Large Cluster case (a) : Empty cluster");

            return;
        } else if ((parent->GetLeafStdCells().size() == 0)
                    && (parent->GetDbModules().size() == 1)) {
            // Testing
            logger_->report("[DEBUG] Break Large Cluster case (b) : Parent is a single logical module");
            
            // (b) Only one logical in the selected cluster
            // E.X. the root cluster
            odb::dbModule* module = parent->GetDbModules()[0];

            // Check all child logical modules
            // (b.I) if the logical module has no child logical module
            // (a flat cluster, or a leaf logical module)
            // There are two ways to break the cluster
            // 1) Call TP to partition this large flat cluster
            // 2) Conduct a bottom up clustering
            if (module->getChildren().size() == 0) {
                // Testing
                logger_->report("\t\tCase b.I: No child logical module detected --> Flat cluster or leaf logical cluster");

                for (odb::dbInst* inst : module->getInsts()) {
                    const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
                    // Physical only?
                    if (liberty_cell == nullptr) continue;

                    odb::dbMaster* master = inst->getMaster();
                    // Check whether the instance is a Pad, Cover or empty block
                    if (master->isPad() || master->isCover()) {
                        continue;
                    } else if (master->isBlock()) {
                        // Shouldn't happen, we don't support macros for now
                        // Testing
                        std::cout << "[ERROR] There shouldn't be any macros in the design!" << std::endl;
                    } else {
                        parent->AddLeafStdCell(inst);
                    }
                }

                // Check whether this is a flat design (without any hierarchy from the beginning)
                if (parent->GetId() == 0) {
                    // Testing
                    logger_->report("\t\t[Warning] This is a flat design without any hierarchy!");

                    std::string cluster_name
                        = std::string("(") + parent->GetName() + ")_flat";

                    Cluster* cluster = new Cluster(cluster_id_, cluster_name, logger_);
                    for (auto& inst : parent->GetLeafStdCells()) {
                        cluster->AddLeafStdCell(inst);
                    }

                    SetInstProperty(cluster);
                    SetClusterMetrics(cluster);
                    cluster_map_[cluster_id_++] = cluster;
                    // Modify the physical hierarchy tree
                    cluster->SetParent(parent);
                    parent->AddChild(cluster);
                }

                // Testing
                logger_->report("\t\t\t{} cells added to cluster: {}",
                                parent->GetNumStdCell(),
                                parent->GetName());

                // Remove the single module from the parent cluster
                parent->ClearDbModules();
                SetInstProperty(parent);
                return;
            }

            // Testing
            logger_->report("\t\tCase b.II: {} Child logical modules detected",
                            module->getChildren().size());

            // (b.II) if the selected logical module has child logical modules
            // Each of them will be modeled as a cluster
            for (odb::dbModInst* child : module->getChildren()) {
                std::string cluster_name = child->getMaster()->getHierarchicalName();
                Cluster* cluster = new Cluster(cluster_id_, cluster_name, logger_);
                cluster->AddDbModule(child->getMaster());
                SetInstProperty(cluster);
                SetClusterMetrics(cluster);
                cluster_map_[cluster_id_++] = cluster;
                // Testing
                logger_->report("\t\t\t[DEBUG] Add new cluster {} to {}", cluster_name, parent->GetName());

                logger_->report("\t\t\t\tNew Cluster: num_std_cells: {}, std_cell_area: {}", cluster->GetNumStdCell(), cluster->GetStdCellArea());

                // Modify the physical hierarchy tree
                cluster->SetParent(parent);
                parent->AddChild(cluster);
            }

            // Check for existing glue logics(leaf std cells)
            std::string cluster_name
                = std::string("(") + parent->GetName() + ")_glue_logic";
            Cluster* cluster = new Cluster(cluster_id_, cluster_name, logger_);
            // Get all contained leaf_std_cells
            for (odb::dbInst* inst : module->getInsts()) {
                const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
                // Physical only cell?
                if (liberty_cell == nullptr) continue;

                odb::dbMaster* master = inst->getMaster();
                // Check if the instance is a Pad, Cover or empty block
                // Macros are not supported for now
                if (master->isPad() || master->isCover()) {
                    continue;
                } else if (master->isBlock()) {
                    // Not supported for now
                    continue;
                } else {
                    cluster->AddLeafStdCell(inst);
                }
            }

            // If the module has no glue instances
            if (cluster->GetLeafStdCells().size() == 0) {
                delete cluster;
            } else {
                SetInstProperty(cluster);
                SetClusterMetrics(cluster);
                cluster_map_[cluster_id_++] = cluster;
                // Modify the physical hierarchy tree
                cluster->SetParent(parent);
                parent->AddChild(cluster);

                // Testing
                logger_->report("\t\t\t[DEBUG] Glue logic cluster {} added to cluster {} : Num std cells {}", cluster->GetName(), parent->GetName(), cluster->GetNumStdCell());
            }
        } else {
            // Testing
            logger_->report("[DEBUG] Break Large Cluster case (c) : Parent has both logical modules and leaf cells");

            // (c) parent is a cluster with both logical modules and leaf std cells
            // parent cluster has few logic modules or many glue insts(many std insts)
            for (auto& module : parent->GetDbModules()) {
                std::string cluster_name = module->getHierarchicalName();
                Cluster* cluster = new Cluster(cluster_id_, cluster_name, logger_);
                cluster->AddDbModule(module);
                SetInstProperty(cluster);
                SetClusterMetrics(cluster);
                cluster_map_[cluster_id_++] = cluster;
                // Modify the physical hierarchy tree
                cluster->SetParent(parent);
                parent->AddChild(cluster);
            }
            // Check Contained LeafStdCells
            if (parent->GetLeafStdCells().size() > 0) {
                std::string cluster_name
                    = std::string("(") + parent->GetName() + ")_glue_logic";
                
                Cluster* cluster = new Cluster(cluster_id_, cluster_name, logger_);
                for (auto& inst : parent->GetLeafStdCells()) {
                    cluster->AddLeafStdCell(inst);
                }

                SetInstProperty(cluster);
                SetClusterMetrics(cluster);
                cluster_map_[cluster_id_++] = cluster;
                // Modify the physical hierarchy tree
                cluster->SetParent(parent);
                parent->AddChild(cluster);
            }
        }

        // Step 2: Break all large child clusters with logical modules
        // Recursively break down large clusters with logical modules
        // Large flat clusters will be broken in the UpdateSubTree function
        for (auto& child : parent->GetChildren()) {
            if (child->GetDbModules().size() > 0) {
                if (child->GetNumStdCell() > max_num_inst_) {
                    BreakCluster(child);
                }
            }
        }

        // Step 3: Merge small clusters
        std::vector<Cluster*> candidate_clusters;
        for (auto& cluster : parent->GetChildren()) {
            if (!cluster->GetIOClusterFlag() && (cluster->GetNumStdCell() < min_num_inst_)) {
                candidate_clusters.push_back(cluster);
            }
        }

        MergeClusters(candidate_clusters);

        // Update the cluster_id
        // Important to maintain the clustering results
        // Set back to the start point
        SetInstProperty(parent);
    }

    // Merge small clusters with the same parent cluster
    // Recursively merge clusters
    // Here is an example process based on connection signature
    // Iter1 :  A, B, C, D, E, F
    // Iter2 :  A + C,  B + D,  E, F
    // Iter3 :  A + C + F, B + D, E
    // End if there is no same connection signature
    // During the merging process, we support two types of merging
    // Type 1: merging small clusters to their closely connected clusters
    //         For example, if a small cluster A is closely connected to a
    //         well-formed cluster B, (there are also other well-formed clusters
    //         C, D), A is only connected to B and A has no connection with C, D
    // Type 2: merging small clusters with the same connection signature
    //         For example, if we merge small clusters A and B,  A and B will have
    //         exactly the same connections relative to all other clusters (both
    //         small clusters and well-formed clusters). In this case, if A and B
    //         have the same connection signature, A and C have the same connection
    //         signature, then B and C also have the same connection signature.
    // Type 3: mergeing all small clusters that are left after the first two merging
    //         processes.
    // Note in both types, we only merge clusters with the same parent cluster
    void HierBlobPlacement::MergeClusters(std::vector<Cluster*>& candidate_clusters) {
        if (candidate_clusters.size() <= 0) {
            return;
        }

        int merge_iter = 0;

        // Testing
        std::cout << "[DEBUG] Merge Cluster Iter: " << merge_iter++ << std::endl;
        for (auto& cluster : candidate_clusters) {
            std::cout << "[DEBUG] Cluster: " << cluster->GetName()
                      << ", num std cells: " << cluster->GetNumStdCell() << std::endl;
        }

        int num_candidate_clusters = candidate_clusters.size();

        while (true) {
            // Update the connection between clusters
            CalculateConnection();

            // Merge flags
            std::vector<int> cluster_class(num_candidate_clusters, -1);
            // Store cluster id
            std::vector<int> candidate_clusters_id;
            for (auto& cluster : candidate_clusters) {
                candidate_clusters_id.push_back(cluster->GetId());
            }

            // 1) Perform Type 1 merge, if possible
            for (int i = 0; i < num_candidate_clusters; i++) {
                const int cluster_id = candidate_clusters[i]->GetCloseCluster(
                        candidate_clusters_id, signature_net_threshold_);
                
                // Testing
                // std::cout << "Cadnidate cluster: " << candidate_clusters[i]->GetName() << " - "
                //     << (cluster_id != -1 ? cluster_map_[cluster_id]->GetName() : "      ");
                
                // logger_->report("[DEBUG] Candidate cluster id {}", cluster_id);

                if (cluster_id != -1 && !cluster_map_[cluster_id]->GetIOClusterFlag() &&
                    (cluster_map_[cluster_id]->GetStdCellArea() +  candidate_clusters[i]->GetStdCellArea() < level_1_cluster_area_upper_bound_)) {
                    // This close connection cluster exists
                    Cluster* cluster = cluster_map_[cluster_id];
                    bool delete_flag = false;
                    if (cluster->MergeCluster(*candidate_clusters[i], delete_flag)) {
                        if (delete_flag) {
                            cluster_map_.erase(candidate_clusters[i]->GetId());
                            delete candidate_clusters[i];
                        }
                        SetInstProperty(cluster);
                        SetClusterMetrics(cluster);
                        cluster_class[i] = cluster->GetId();
                    }
                }
            }

            // Step to Type 2 merge
            std::vector<Cluster*> new_candidate_clusters;
            for (int i = 0; i < num_candidate_clusters; i++) {
                if (cluster_class[i] == -1) {
                    // This cluster is not merged
                    for (int j = i + 1; j < num_candidate_clusters; j++) {
                        if (cluster_class[j] != -1) {
                            continue;
                        }
                        bool flag = candidate_clusters[i]->IsSameConnectionSignature(
                            *candidate_clusters[j], signature_net_threshold_
                        );

                        // If they have the same signature
                        if (flag &&
                            (candidate_clusters[j]->GetStdCellArea() +  candidate_clusters[i]->GetStdCellArea() < level_1_cluster_area_upper_bound_)) {
                            cluster_class[j] = i;
                            bool delete_flag = false;
                            if (candidate_clusters[i]->MergeCluster(*candidate_clusters[j], delete_flag)) {
                                if (delete_flag) {
                                    cluster_map_.erase(candidate_clusters[j]->GetId());
                                    delete candidate_clusters[j];
                                }
                                SetInstProperty(candidate_clusters[i]);
                                SetClusterMetrics(candidate_clusters[i]);
                            }
                        }
                    }
                }
            }

            // Perform the last type of mergeing, merge all dust clusters
            const int dust_cluster_std_cell = 10;
            for (int i = 0; i < num_candidate_clusters; i++) {
                if (cluster_class[i] == -1) {
                    // Cluster is not merged
                    new_candidate_clusters.push_back(candidate_clusters[i]);
                    if (candidate_clusters[i]->GetNumStdCell() <= dust_cluster_std_cell) {
                        for (int j = i + 1; j < num_candidate_clusters; j++) {
                            if (cluster_class[j] != -1
                                || candidate_clusters[j]->GetNumStdCell() > dust_cluster_std_cell
                                || (candidate_clusters[j]->GetStdCellArea() +  candidate_clusters[i]->GetStdCellArea() > level_1_cluster_area_upper_bound_)) {
                                    continue;
                            }
                            
                            cluster_class[j] = i;
                            bool delete_flag =false;
                            if (candidate_clusters[i]->MergeCluster(*candidate_clusters[j], delete_flag)) {
                                if (delete_flag) {
                                    cluster_map_.erase(candidate_clusters[j]->GetId());
                                    delete candidate_clusters[j];
                                }
                                SetInstProperty(candidate_clusters[i]);
                                SetClusterMetrics(candidate_clusters[i]);
                            }
                        }
                    }
                }
            }

            // Update the candidate clusters
            // Some clusters have become well-formed clusters
            candidate_clusters.clear();
            for (auto& cluster : new_candidate_clusters) {
                if (cluster->GetNumStdCell() < min_num_inst_) {
                    candidate_clusters.push_back(cluster);
                }
            }

            // If nothing happened. Break
            if (num_candidate_clusters == new_candidate_clusters.size()) {
                break;
            }

            num_candidate_clusters = candidate_clusters.size();

            // Testing
            // std::cout << "[DEBUG] Merge Cluster Iter : " << merge_iter++ << std::endl;
            // for (auto& cluster: candidate_clusters) {
            //     std::cout << "[DEBUG] Cluster: " << cluster->GetName() << std::endl;
            // }

            // Break if not needed
            if (candidate_clusters.size() == 0) {
                break;
            }

            // TODO: Shall we directly merge all left clusters?
            // Will this cause any problems? Probably it will
            // Testing
            // std::cout << "[DEBUG] " << candidate_clusters.size() 
            //           << " left unmerged, and will be kept in this level of the physical hierarchy tree." << std::endl;
        }

        // Testing
        std::cout << "[DEBUG] Finish merging clusters" << std::endl;
    }

    // This function has two purposes:
    // 1) Remove all the internal clusters between parent and leaf clusters in its subtree
    // 2) Break Large flat clusters with TP or Bottom-up grouping
    void HierBlobPlacement::UpdateSubTree(Cluster* parent) {
        std::vector<Cluster*> children_clusters;
        std::vector<Cluster*> internal_clusters;
        std::queue<Cluster*> wavefront;  // Define a FIFO

        // Testing
        logger_->report("[DEBUG] Update subtree for cluster : {}", parent->GetName());

        for (auto child : parent->GetChildren()) {
            wavefront.push(child);
        }

        // Testing
        logger_->report("\t{} Child clusters found for parent: {}",
                        parent->GetChildren().size(),
                        parent->GetName());

        while (!wavefront.empty()) {
            // Get the first child
            Cluster* cluster = wavefront.front();
            wavefront.pop();

            if (cluster->GetChildren().size() == 0) {
                // Leaf cluster found
                children_clusters.push_back(cluster);
            } else {
                internal_clusters.push_back(cluster);
                for (auto child : cluster->GetChildren()) {
                    wavefront.push(child);
                }
            }
        }

        // Delete all internal clusters
        for (auto& cluster : internal_clusters) {
            cluster_map_.erase(cluster->GetId());
            delete cluster;
        }

        parent->RemoveChildren();
        parent->AddChildren(children_clusters);
        for (auto& cluster : children_clusters) {
            cluster->SetParent(parent);
            if (cluster->GetNumStdCell() > max_num_inst_) {
                // Here we provide two ways of breaking the large cluster
                // 1) Call TritonPart to do iteratively partitioning
                // 2) Do a timing aware bottom-up clustering
                // Testing
                logger_->report("\tBreaking child cluster: {}", cluster->GetName());

                BreakLargeFlatCluster_TP(cluster);
                // BreakLargeFlatCluster_BU(cluster);
            }
        }
    }

    // Break large flat clusters with TritonPart
    // A large flat cluster does not have a logical module
    void HierBlobPlacement::BreakLargeFlatCluster_TP(Cluster* parent) {
        // Check if the cluster is a large flat cluster
        if (parent->GetDbModules().size() > 0
            || parent->GetLeafStdCells().size() < max_num_inst_) {
                return;
        }

        // Set the instance property for all contained instances in the cluster
        SetInstProperty(parent);
        // Variables for hypergraph
        std::map<int, int> cluster_vertex_id_map;
        std::map<odb::dbInst*, int> inst_vertex_id_map;
        const int parent_cluster_id = parent->GetId();
        std::vector<odb::dbInst*> std_cells = parent->GetLeafStdCells();
        std::vector<std::vector<int>> hyperedges;
        std::vector<float> vertex_weight;

        // Vertices
        // Other clusters behave like fixed vertices
        // We do not consider vertices only between fixed vertices
        // Build the hypergraph for existing clusters
        int vertex_id = 0;
        for (auto& [cluster_id, cluster] : cluster_map_) {
            cluster_vertex_id_map[cluster_id] = vertex_id++;
            vertex_weight.push_back(0.0f);
        }

        // Existing clusters won't be touched during the breakup process
        int num_fixed_vertices = vertex_id;

        for (auto& std_cell : std_cells) {
            inst_vertex_id_map[std_cell] = vertex_id++;
            const sta::LibertyCell* liberty_cell = network_->libertyCell(std_cell);
            vertex_weight.push_back(liberty_cell->area());
        }

        // Create the corresponding hyperedge
        for (odb::dbNet* net : block_->getNets()) {
            // Ignore all power nets
            if (net->getSigType().isSupply()) {
                continue;
            }

            int driver_id = -1;
            std::set<int> loads_id;
            bool pad_flag = false;

            // Check the connected instances
            for (odb::dbITerm* iterm : net->getITerms()) {
                odb::dbInst* inst = iterm->getInst();
                const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
                if (liberty_cell == nullptr) continue;

                odb::dbMaster* master = inst->getMaster();
                // Pad, Cover or Block(macro)?
                if (master->isPad() || master->isCover()) {
                    pad_flag = true;
                    break; // No need to consider this net
                }

                const int cluster_id
                    = odb::dbIntProperty::find(inst, "cluster_id")->getValue();
                
                int vertex_id = (cluster_id != parent_cluster_id) ? cluster_vertex_id_map[cluster_id] : inst_vertex_id_map[inst];
                if (iterm->getIoType() == odb::dbIoType::OUTPUT) {
                    driver_id = vertex_id;
                } else {
                    loads_id.insert(vertex_id);
                }
            }

            // Ignore all nets with IO pads
            if (pad_flag) continue;
            // Check all IO pins
            for (odb::dbBTerm* bterm : net->getBTerms()) {
                const int cluster_id
                    = odb::dbIntProperty::find(bterm, "cluster_id")->getValue();
                
                if (bterm->getIoType() == odb::dbIoType::INPUT) {
                    driver_id = cluster_vertex_id_map[cluster_id];
                } else {
                    loads_id.insert(cluster_vertex_id_map[cluster_id]);
                }
            }
            loads_id.insert(driver_id);
            // Add the net as a hyperedge
            if (driver_id != -1 && loads_id.size() > 1
                && loads_id.size() < large_net_threshold_) {
                std::vector<int> hyperedge;
                hyperedge.insert(hyperedge.end(), loads_id.begin(), loads_id.end());
                hyperedges.push_back(hyperedge);
            }
        }

        // Parameters definition for TritonPart
        const int seed = 0;
        const float balance_constraint = 1.0;
        const int num_parts = 2; // We use two-way partitioning here
        const int num_vertices = static_cast<int>(vertex_weight.size());
        std::vector<float> hyperedge_weights(hyperedges.size(), 1.0f);

        std::vector<int> part = tritonpart_->PartitionKWaySimpleMode(num_parts,
                                                            balance_constraint,
                                                            seed,
                                                            hyperedges,
                                                            vertex_weight,
                                                            hyperedge_weights
        );

        // Create cluster based on partitioning solutions
        // Note that all the std cells are stored in the leaf_std_cells for a flat cluster
        parent->ClearLeafStdCells();
        // we follow binary coding method to differentiate different parts
        // of the cluster
        // cluster_name_0, cluster_name_1
        // cluster_name_0_0, cluster_name_0_1, cluster_name_1_0, cluster_name_1_1
        const std::string cluster_name = parent->GetName();
        // Set the parent cluster for part 0
        // Update the name of the parent cluster
        parent->SetName(cluster_name + std::string("_0"));
        // Create a new cluster for part 1
        Cluster* cluster_part_1
            = new Cluster(cluster_id_, cluster_name + std::string("_1"), logger_);
        // We do not need to touch all those fixed vertices
        for (int i = num_fixed_vertices; i < num_vertices; i++) {
            if (part[i] == 0) {
                parent->AddLeafStdCell(std_cells[i - num_fixed_vertices]);
            } else {
                cluster_part_1->AddLeafStdCell(std_cells[i - num_fixed_vertices]);
            }
        }

        // Update the property of parent cluster
        SetInstProperty(parent);
        SetClusterMetrics(parent);
        // Update the property of cluster_part_1
        SetInstProperty(cluster_part_1);
        SetClusterMetrics(cluster_part_1);
        cluster_map_[cluster_id_++] = cluster_part_1;
        cluster_part_1->SetParent(parent->GetParent());
        parent->GetParent()->AddChild(cluster_part_1);
        // Recursive break the cluster
        // Until the size of the cluster is less than max_num_inst_
        BreakLargeFlatCluster_TP(parent);
        BreakLargeFlatCluster_TP(cluster_part_1);
    }

    // Break large flat clusters with a Bottom up clustering
    // A large flat cluster does not have a logical module
    // We will conduct a bottom-up clustering to break this large flat cluster
    // The Clustering is done in a best-choice manner based on the clustering score
    // Before conducting the clustering, we will construct the corresponding hypergraph
    // All clusters outside the flat large cluster will be modeled as fixed vertices and 
    // won't be touched during the grouping process
    // The clustering scrore between two vertices u, v in the hypergraph is defined as
    // d(u, v) = 1/|e| + timing_weight
    // IMPORTANT!! For all the clustering process we are using a star model instead of a clique model
    void HierBlobPlacement::BreakLargeFlatCluster_BU(Cluster* parent) {
        // Check if the cluster is a large flat cluster
        if (parent->GetDbModules().size() > 0
            || parent->GetLeafStdCells().size() < max_num_inst_) {
                return;
        }

        // Testing
        logger_->report("[Bottom-UP] Breaking Large Flat Cluster {}; Num std cells: {}, Max threshold: {}",
                        parent->GetName(),
                        parent->GetLeafStdCells().size(),
                        max_num_inst_);

        // Step 1: Build the Hypergraph structure
        // We only cares about cells inside this large flat cluster
        // Set the instance property
        SetInstProperty(parent);

        // Define the hypergraph structure
        Hypergraph temp_hypergraph;

        // Testing
        // logger_->report("\t[DEBUG] Starting building hypergraph");

        // Only std cells in the large flat cluster will be considered
        // other clusters will just behave like fixed vertices,
        // their vertex weight is set to 0
        int vertex_id = 0;
        int hyperedge_id = 0;
        int group_id = 0;

        // Grouping parameters
        int desired_num_groups = std::ceil(parent->GetLeafStdCells().size() / max_num_inst_);
        float balanced_area_per_group = parent->GetStdCellArea() / desired_num_groups;

        // Testing
        logger_->report("\tDesired number groups {}, balanced area per group {}",
                        desired_num_groups,
                        balanced_area_per_group * (1 - tolerance_));

        // Sub graph information
        const int parent_cluster_id = parent->GetId();
        std::vector<odb::dbInst*> std_cells = parent->GetLeafStdCells();

        // Store timing information
        std::vector<TimingPath> temp_timing_paths;

        // Store grouping information
        std::vector<Group> temp_grouping_info;

        for (auto& [cluster_id, cluster] : cluster_map_) {
            temp_hypergraph.cluster_vertex_id_map[cluster_id] = vertex_id++;
            temp_hypergraph.vertex_weight.push_back(0.0f);
        }

        // All clusters outside this large flat cluster behave like fixed clusters
        int num_fixed_vertices = vertex_id;

        // Add all std cells to the hypergraph
        // Group will be created for each standard cells in the selected cluster
        for (auto& std_cell : std_cells) {
            // Skip Physical only cells, shouldn't happen, wrong netlist read in
            if (std_cell->getName().find("FILLER") != std::string::npos) {
                // Testing
                // logger_->report("[ERROR] {} detected, Please use netlist after synthesis!", std_cell->getName());

                continue;
            }

            // Construct the hypergraph storing structure
            temp_hypergraph.inst_vertex_id_map[std_cell] = vertex_id++;
            temp_hypergraph.vertex_id_instance_map[vertex_id - 1] = std_cell;
            const sta::LibertyCell* liberty_cell = network_->libertyCell(std_cell);
            temp_hypergraph.vertex_weight.push_back(liberty_cell->area());
            temp_hypergraph.vertex_id_group_id_map[vertex_id - 1] = vertex_id - 1 - num_fixed_vertices;

            // Create the initial group information
            Group temp_group;
            temp_group.included_vertices.push_back(vertex_id - 1);
            temp_group.group_id = vertex_id - 1 - num_fixed_vertices;
            temp_group.group_size = liberty_cell->area();
            temp_grouping_info.push_back(temp_group);
        }

        // Testing
        // logger_->report("\t[INFO] {} Physical only cell removed from the grouping process.", std_cells.size() - temp_grouping_info.size());

        // Traverse all nets to create hyperedges - star model
        for (odb::dbNet* net : block_->getNets()) {
            // Ignore all power nets
            if (net->getSigType().isSupply()) {
                continue;
            }

            // vertex_id of the driver instance
            int driver_id = -1;
            // vertex_id of the sink instances
            std::set<int> loads_id;
            bool pad_flag = false;

            // Check the connected instances
            for (odb::dbITerm* iterm : net->getITerms()) {
                odb::dbInst* inst = iterm->getInst();
                const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
                if (liberty_cell == nullptr) continue;
                odb::dbMaster* master = inst->getMaster();
                // Avoid pads or cover
                if (master->isPad() || master->isCover()) {
                    pad_flag = true;
                    break;
                }

                const int cluster_id = odb::dbIntProperty::find(inst, "cluster_id")->getValue();
                // Get the corresponding vertex id
                int vertex_id = (cluster_id != parent_cluster_id) ? temp_hypergraph.cluster_vertex_id_map[cluster_id] : temp_hypergraph.inst_vertex_id_map[inst];

                if (iterm->getIoType() == odb::dbIoType::OUTPUT) {
                    driver_id = vertex_id;
                } else {
                    loads_id.insert(vertex_id);
                }
            }

            // Ignore nets with IO pads
            if (pad_flag) {
                continue;
            }

            // Check connected IO pins
            for (odb::dbBTerm* bterm : net->getBTerms()) {
                const int cluster_id = odb::dbIntProperty::find(bterm, "cluster_id")->getValue();

                // Check whether it's driver or 
                if (bterm->getIoType() == odb::dbIoType::INPUT) {
                    driver_id = temp_hypergraph.cluster_vertex_id_map[cluster_id];
                } else {
                    loads_id.insert(temp_hypergraph.cluster_vertex_id_map[cluster_id]);
                }
            }
            loads_id.insert(driver_id);

            // Add the net as a hyperedge
            // The initial weight of the hyper edge will be set to 1 (base weight)
            if (driver_id != -1 && loads_id.size() > 1 && loads_id.size() < large_net_threshold_) {
                std::vector<int> hyperedge(loads_id.begin(), loads_id.end());

                // Update all storing structure
                temp_hypergraph.hyperedges.push_back(hyperedge);
                temp_hypergraph.net_hyperedge_id_map[net] = hyperedge_id++;
                temp_hypergraph.hyperedge_weight.push_back(1.0f);
            }
        }

        // Traverse all the hyperedge to create the group storing structure
        // Fixed vertices have no corresponding grouping structure, they are not
        // considered during the bottom-up clustering process
        int hyperedge_counter = 0;
        for (auto hyperedge : temp_hypergraph.hyperedges) {

            // Testing
            // logger_->report("\tHyperedge {}, size {}", hyperedge_counter, hyperedge.size());

            // For all vertices in the hyperedge
            for (auto vertex_id : hyperedge) {
                // Creat the corresponding set variable if this is a movable vertex
                if (vertex_id >= num_fixed_vertices) {
                    for (auto second_vertex_id : hyperedge) {
                        // Add the neighbor only if the second vertex is movable
                        if (second_vertex_id >= num_fixed_vertices && (second_vertex_id != vertex_id)) {
                            temp_grouping_info[vertex_id - num_fixed_vertices].neighbor_vetex_id.insert(second_vertex_id);
                        }
                    }
                    
                    // Insert the hyperedge id to the group
                    temp_grouping_info[vertex_id - num_fixed_vertices].connected_hyperedge_id.insert(hyperedge_counter);
                }

                // Testing
                // logger_->report("\t{} ", vertex_id);
            }

            // Increate the counter value
            hyperedge_counter++;
        }

        // Testing
        // logger_->report("[DEBUG] Print info for all groups generated from cluster {}", parent->GetName());
        // for (auto sel_group : temp_grouping_info) {
        //     PrintGroupInfo(sel_group);
        // }

        // Testing
        // logger_->report("\t[STATUS] Finished building the hypergraph");
        // logger_->report("\t\tNumber of fixed vertices: {}", num_fixed_vertices);
        // logger_->report("\t\tNumber of critical_timing_paths: {}", num_critical_timing_paths_);
        // logger_->report("\t\tNumber of created groups: {}", temp_grouping_info.size());
        // logger_->report("\t\tNumber of created hyperedges: {}", hyperedge_counter);

        // Step 2: Extract all needed timing information
        FindTimingPaths(temp_hypergraph, temp_timing_paths, parent_cluster_id);

        // Testing
        // logger_->report("\t[STATUS] Finished building timing paths");
        // logger_->report("\t[INFO] {} Created timing paths", temp_timing_paths.size());

        // Step 3: Conduct the best choice clustering
        BestChoiceClustering(temp_hypergraph, temp_grouping_info, desired_num_groups, balanced_area_per_group * (1 - tolerance_));

        // Testing
        logger_->report("\t[STATUS] Finished area balanced best-choice clustering for {}", parent_cluster_id);

        // Step 4: Store the clustering results as clusters in the physical hierarchy tree
        // All mapping information is stored in the hypergraph vertex_id_group_mapping
        // The original parent cluster will be set as cluster_name_0, other cluster name will be set as 
        // cluster_name_1, cluster_name_2 ....
        std::set<int> unique_group_ids;
        int cluster_num_counter = 0;

        // Insert the parent cluster as the first node
        parent->ClearLeafStdCells();
        const std::string cluster_name = parent->GetName();
        parent->SetName(cluster_name + std::string("_0"));
        cluster_num_counter++;
        
        // Get all unique group_ids
        for (auto sel_element : temp_hypergraph.vertex_id_group_id_map) {
            if (unique_group_ids.find(sel_element.second) == unique_group_ids.end()) {
                unique_group_ids.insert(sel_element.second);
            }
        }

        // Testing
        logger_->report("\t[DEBUG] Num unique group_ids before refining: {}", unique_group_ids.size());
        // std::cout << "\t[DEBUG] Print all unique group ids " ;
        // for (auto sel_group_id : unique_group_ids) {
        //     std::cout << sel_group_id << " " ;
        // }
        // std::cout << std::endl;

        // Refine the clusturing process, remove those small groups
        RefineGrouping(temp_grouping_info, unique_group_ids, balanced_area_per_group, num_fixed_vertices);

        // Testing
        logger_->report("\t[DEBUG] Num unique group_ids after refining: {}", unique_group_ids.size());
        // std::cout << "\t[DEBUG] Print all unique group ids " ;
        // for (auto sel_group_id : unique_group_ids) {
        //     std::cout << sel_group_id << " " ;
        // }
        // std::cout << std::endl;

        // Create all needed new cluster
        for (auto sel_group_id : unique_group_ids) {
            if (cluster_num_counter == 1) {
                for (auto sel_vertex : temp_grouping_info[sel_group_id].included_vertices) {
                    // parent->AddLeafStdCell(std_cells[sel_vertex - num_fixed_vertices]);
                    parent->AddLeafStdCell(temp_hypergraph.vertex_id_instance_map[sel_vertex]);
                }

                // Testing
                // logger_->report("[DEBUG] Print info for unique group {}", sel_group_id);
                // PrintGroupInfo(temp_grouping_info[sel_group_id]);

                // Update the property
                SetInstProperty(parent);
                SetClusterMetrics(parent);
                cluster_num_counter++;
            } else {
                // New cluster needed
                Cluster* temp_cluster_instance 
                    = new Cluster(cluster_id_, cluster_name + std::string("_") + std::to_string(cluster_num_counter), logger_);

                for (auto sel_vertex: temp_grouping_info[sel_group_id].included_vertices) {
                    // temp_cluster_instance->AddLeafStdCell(std_cells[sel_vertex - num_fixed_vertices]);
                    temp_cluster_instance->AddLeafStdCell(temp_hypergraph.vertex_id_instance_map[sel_vertex]);
                }

                // Testing
                // logger_->report("[DEBUG] Print info for unique group {}", sel_group_id);
                // PrintGroupInfo(temp_grouping_info[sel_group_id]);

                // Update the property
                SetInstProperty(temp_cluster_instance);
                SetClusterMetrics(temp_cluster_instance);

                // Update all the counters
                cluster_num_counter++;
                cluster_map_[cluster_id_++] = temp_cluster_instance;
                temp_cluster_instance->SetParent(parent->GetParent());
                parent->GetParent()->AddChild(temp_cluster_instance);
            }
        }

        // Testing
        logger_->report("[STATUS] Finished Best-Choice Bottom Up Clustering, {} new clusters created", cluster_num_counter - 1);
    }

    // Extract information from all critical timing paths
    // This function is inspired from the implementation in TritonPart.cpp
    // and is similar to gui/src/staGui.cpp
    // Please refer to sta/Search/ReportPath.cc for how to check the timing paths
    // All extracted slack value will be normalized to the clock period
    // **************************** Definition for the findPathend function **********
    // PathEndSeq *findPathEnds(ExceptionFrom *from,
    //              ExceptionThruSeq *thrus,
    //              ExceptionTo *to,
    //              bool unconstrained,
    //              const Corner *corner,
    //              const MinMaxAll *min_max,
    //              int group_count,
    //              int endpoint_count,
    //              bool unique_pins,
    //              float slack_min,
    //              float slack_max,
    //              bool sort_by_slack,
    //              PathGroupNameSet *group_names,
    //              bool setup,
    //              bool hold,
    //              bool recovery,
    //              bool removal,
    //              bool clk_gating_setup,
    //              bool clk_gating_hold);
    // PathEnds represent search endpoints that are either unconstrained or
    // constrained by a timing check, output delay, data check, or path delay.
    // ************************************************************
    // Here we need to use the corresponding group_id.
    void HierBlobPlacement::FindTimingPaths(Hypergraph& hypergraph,
                                            std::vector<TimingPath>& timing_paths,
                                            const int parent_cluster_id) {
        // Check the top_n value
        if (num_critical_timing_paths_ <= 0) {
            // Skip timing calculation
            return;
        }

        // Ensure the timing paths has been built
        sta_->ensureGraph();
        // Make the graph and find delays
        sta_->searchPreamble();
        sta_->ensureLevelized();
        logger_->report("\t[WARNING] We normalized the slack based on the maximum clock period");

        // Step 1: Find the top_n timing paths
        sta::ExceptionFrom* e_from = nullptr;
        sta::ExceptionThruSeq* e_thrus = nullptr;
        sta::ExceptionTo* e_to = nullptr;
        bool include_unconstrained = false;
        // Max is used for setup check , and min used for hold check
        bool get_max = true;

        // Timing paths are grouped into path groups according to the clock associated with
        // the endpoint of the path.
        int group_count = num_critical_timing_paths_;
        // The number of paths to report for each endpoint
        int endpoint_count = 1;

        // Step 2: Get all those paths
        sta::PathEndSeq path_ends = sta_->search()->findPathEnds(
            e_from,  // Return paths from a list of clocks/instances/ports/register clock pins or latch data pins
            e_thrus, // Return paths through a list of instances/ports/nets
            e_to,    // Return paths to a list of clocks/instances/ports or pins
            include_unconstrained, // Return unconstrained paths
            // corner, min_max
            sta_->cmdCorner(),     // Return paths for a process corner
            get_max ? sta::MinMaxAll::max() : sta::MinMaxAll::min(), // Return max/min path checks
            group_count,    // path_count used by GUI
            endpoint_count, // path_count used by GUI
            true,
            -sta::INF,      // Slack min
            sta::INF,       // Slack max
            true,           // Sort by slack
            nullptr,        // group_names
            // Setup, Hold, Recovery, Removal
            get_max,
            !get_max,
            false,
            false,
            // Clk_gating_setup, clk_gating_hold
            false,
            false
        );

        // Step 3: Check all the timing paths and store all needed information
        // Path slack will be added to to all hyperedge that it passes
        // The edge weight will be 1 + net timing weight + path timing weight
        for (auto& path_end : path_ends) {
            // Testing, print timing paths to logger
            // sta_->reportPathEnd(path_end);
            
            auto* path = path_end->path();
            // Initialize the timing path structure
            TimingPath timing_path;
            // Get the corresponding slack information
            const float slack = path_end->slack(sta_);
            
            // Testing
            // logger_->report("[INFO] Slack = {}", slack);

            // Normalize the slack according to the maximum clock period
            // TODO: potentially buggy, please fix the max_timing_calculation
            const float clock_period = path_end->targetClk(sta_)->period();

            // Testing
            // logger_->report("[INFO] Tatget Clk Period = {}", clock_period);

            maximum_clock_period_ = std::max(maximum_clock_period_, clock_period);
            timing_path.slack = slack;

            // Testing
            // logger_->report("[INFO] Timing Path Slack = {}", timing_path.slack);

            sta::PathExpanded expand(path, sta_);
            expand.path(expand.size() - 1);

            // Check all vertices on the path
            for (size_t i = 0; i < expand.size(); i++) {
                // PathRef is the reference to a path vertex
                sta::PathRef* ref = expand.path(i);
                sta::Pin* pin = ref->vertex(sta_)->pin();
                // Get nets connecting the pin
                auto net = network_->net(pin);

                // Check whether this pin is connected to any net
                if (net == nullptr) continue;
                
                // Check whether this is an IO pin
                if (network_->isTopLevelPort(pin) == true) {
                    auto bterm = block_->findBTerm(network_->pathName(pin));
                    int cluster_id = odb::dbIntProperty::find(bterm, "cluster_id")->getValue();
                    int vertex_id = 0;

                    if (cluster_id == -1) {
                        continue; // This should never happen
                    } else {
                        int vertex_id = hypergraph.cluster_vertex_id_map[cluster_id];
                    }

                    // Check whether this vertex is in the Timing path or not
                    if (std::find(
                            timing_path.path.begin(), timing_path.path.end(), vertex_id)
                        == timing_path.path.end()) {
                        timing_path.path.push_back(vertex_id);
                    }
                } else {
                    auto inst = network_->instance(pin);
                    auto db_inst = block_->findInst(network_->pathName(inst));

                    // Get the cluster_id first
                    int cluster_id = odb::dbIntProperty::find(db_inst, "cluster_id")->getValue();
                    int vertex_id = 0;

                    // Check where is the cluster
                    if (cluster_id != parent_cluster_id) {
                        vertex_id = hypergraph.cluster_vertex_id_map[cluster_id];
                    } else {
                        vertex_id = hypergraph.inst_vertex_id_map[db_inst];
                    }

                    // Add the vertex to the path
                    if (std::find(
                            timing_path.path.begin(), timing_path.path.end(), vertex_id)
                        == timing_path.path.end()) {
                        timing_path.path.push_back(vertex_id);
                    }
                }

                // Converting sta::net to dbNet
                auto db_net = block_->findNet(network_->pathName(net));

                // Check whether this edge is important or not
                // TODO: Move this checking to the front to avoid useless computations.
                auto checker = hypergraph.net_hyperedge_id_map.find(db_net);

                if (checker == hypergraph.net_hyperedge_id_map.end()) {
                    continue; // This net is not in the created hypergraph
                }

                int hyperedge_id = hypergraph.net_hyperedge_id_map[db_net];

                // Check whether the hyperedge is in the timing path or not
                if (std::find(timing_path.arcs.begin(), timing_path.arcs.end(), hyperedge_id)
                    == timing_path.arcs.end()) {
                    timing_path.arcs.push_back(hyperedge_id);
                }
            }

            // Add timing path to the overall structure
            timing_paths.push_back(timing_path);
        }

        // Testing
        logger_->report("\t[STATUS] Finished Travsersing timing paths. Now updating path timing weight to the corresponding hypergraph");

        // Step 4: Update the path timing information to all 
        for (auto &sel_timing_path : timing_paths) {
            // Normalize the slack
            sel_timing_path.slack = sel_timing_path.slack / maximum_clock_period_;

            // Testing
            // logger_->report("[INFO] Timing Path Slack = {}", sel_timing_path.slack);

            // Update the path timing weight for each hyperedge included in the timing path
            for (auto edge_id : sel_timing_path.arcs) {
                hypergraph.hyperedge_weight[edge_id] += std::pow((1 - sel_timing_path.slack), timing_factor_);
            }
        }

        // Step 5: Check the slack on each net, Update the hyperedge weight with the net slack weight
        // Add the timing weight to the existing hyperedge weight
        logger_->report("\t[INFO] clock_period: {} second", maximum_clock_period_);

        // Traverse all nets in the block
        for (auto db_net : block_->getNets()) {
            // Check the validity of the net
            auto checker = hypergraph.net_hyperedge_id_map.find(db_net);

            if (checker == hypergraph.net_hyperedge_id_map.end()) {
                continue; // This net is not in the created hypergraph
            }            

            const int hyperedge_id = hypergraph.net_hyperedge_id_map[db_net];

            sta::Net *net = network_->dbToSta(db_net);
            const float slack = sta_->netSlack(net, sta::MinMax::max());
            const float normalized_slack = slack / maximum_clock_period_;

            // Update the slack information
            hypergraph.hyperedge_id_net_slack_weight_map[hyperedge_id] = normalized_slack;
            hypergraph.hyperedge_weight[hyperedge_id] += std::pow((1 - normalized_slack), timing_factor_);
        }

        // Testing
        // logger_->report("[INFO] All timing paths information: ");
        // for (auto sel_timing_path : timing_paths) {
        //     PrintTimingPathInfo(sel_timing_path);
        // }

        // Testing
        // int edge_counter = 0;
        // logger_->report("[INFO] All hyperedge weight information: ");
        // for (auto sel_hyperedge_weight : hypergraph.hyperedge_weight) {
        //     std::cout << "[INFO] Hyperedge(" << edge_counter << ") weight : " << sel_hyperedge_weight << std::endl;
        //     edge_counter++;
        // }

        // Testing
        logger_->report("\t[STATUS] Finished extracting timing information");
    }

    // This function will add net timing weight to the corresponding net instances
    void HierBlobPlacement::AddGlobalTimingInfo() {
        // Check the top_n value
        if (num_critical_timing_paths_ <= 0) {
            return;
        }

        // Testing
        logger_->report("[STATUS] Add global net timing weight");

        // Ensure the timing paths has been built
        sta_->ensureGraph();
        // Make the graph and find delays
        sta_->searchPreamble();
        sta_->ensureLevelized();
        logger_->report("[WARNING] We normalized the slack based on the maximum clock period");

        // Find the top_n timing paths
        sta::ExceptionFrom* e_from = nullptr;
        sta::ExceptionThruSeq* e_thrus = nullptr;
        sta::ExceptionTo* e_to = nullptr;
        bool include_unconstrained = false;
        // Max is used for setup check , and min used for hold check
        bool get_max = true;

        // Timing paths are grouped into path groups according to the clock associated with
        // the endpoint of the path.
        int group_count = num_critical_timing_paths_;
        // The number of paths to report for each endpoint
        int endpoint_count = 1;

        // Step 2: Get all those paths
        sta::PathEndSeq path_ends = sta_->search()->findPathEnds(
            e_from,  // Return paths from a list of clocks/instances/ports/register clock pins or latch data pins
            e_thrus, // Return paths through a list of instances/ports/nets
            e_to,    // Return paths to a list of clocks/instances/ports or pins
            include_unconstrained, // Return unconstrained paths
            // corner, min_max
            sta_->cmdCorner(),     // Return paths for a process corner
            get_max ? sta::MinMaxAll::max() : sta::MinMaxAll::min(), // Return max/min path checks
            group_count,    // path_count used by GUI
            endpoint_count, // path_count used by GUI
            true,
            -sta::INF,      // Slack min
            sta::INF,       // Slack max
            true,           // Sort by slack
            nullptr,        // group_names
            // Setup, Hold, Recovery, Removal
            get_max,
            !get_max,
            false,
            false,
            // Clk_gating_setup, clk_gating_hold
            false,
            false
        );

        // Find the maximum clock period
        for (auto& path_end : path_ends) {
            auto* path = path_end->path();
            // Initialize the timing path structure
            TimingPath timing_path;
            // Get the corresponding slack information
            const float slack = path_end->slack(sta_);

            // Normalize the slack according to the maximum clock period
            // TODO: potentially buggy, please fix the max_timing_calculation
            const float clock_period = path_end->targetClk(sta_)->period();

            maximum_clock_period_ = std::max(maximum_clock_period_, clock_period);
        }

        // Add net timing weight to the created bundled net
        logger_->report("[INFO] clock_period: {} second", maximum_clock_period_);

        // Traverse all nets in the block
        for (auto db_net : block_->getNets()) {
            // Get the related bundled net indexes
            std::vector<int> influenced_nets(net_bundled_net_map_[db_net]);

            sta::Net *net = network_->dbToSta(db_net);
            const float slack = sta_->netSlack(net, sta::MinMax::max());
            const float normalized_slack = slack / maximum_clock_period_;

            float timing_weight = std::pow((1 - normalized_slack), timing_factor_);

            for (auto sel_index : influenced_nets) {
                // sparse_nets_[sel_index].timing_weight = (timing_weight > 1.0) ? timing_weight : 1.0;
                sparse_nets_[sel_index].timing_weight = (timing_weight > 0.0) ? timing_weight : 1.0;
            }
        }

        // Testing
        logger_->report("[STATUS] Finished extracting global timing net");
    }

    // This function is used to facilitate the bottom-up clustering process
    // After building the needed hypergraph, we will conduct the best-choice
    // clustering
    void HierBlobPlacement::BestChoiceClustering(Hypergraph& hypergraph, 
                              std::vector<Group>& groups, 
                              int desired_num_groups,
                              float area_limit) {
        // Define all needed variables
        int initial_group_counter = groups.size();
        std::set<int> existing_groups;

        // Second cluster size upper bound
        int max_insts_per_group = std::ceil(initial_group_counter / desired_num_groups);

        // Testing
        logger_->report("\t[INFO] {} Groups exist", initial_group_counter);

        // Add existing group_ids
        for (int i = 0; i < initial_group_counter; i++) {
            existing_groups.insert(i);
        }

        // Testing
        logger_->report("\t[STATUS] Start Best Choice Clustering");

        // Define the needed priority queue for best-choice clustering
        AccessiblePriorityQueue<PQElement> bu_pri_queue;

        // Step 1: Initialize the needed priority queue
        InitializePriorityQueue(hypergraph, bu_pri_queue, groups);

        // Deadlock counter, check whether it's possible to find a suitable merging group for a selected group
        int dead_lock_counter = 0;

        // Step 2: Conduct the contraction, if the number of groups is not met
        // while ((existing_groups.size() > desired_num_groups) && (!bu_pri_queue.empty()))
        while ((!bu_pri_queue.empty())) {
            // Check the validity of the first element
            if (bu_pri_queue.top().valid) {
                // Valid entry found, check whether it satisfies the area constraints
                std::pair<int, int> group_id_pair = bu_pri_queue.top().group_pair;
                float area_sum = groups[group_id_pair.first].group_size + groups[group_id_pair.second].group_size;
                int member_count = groups[group_id_pair.first].included_vertices.size() + groups[group_id_pair.second].included_vertices.size();
                
                // Testing
                // if (area_sum <= area_limit && member_count <= max_insts_per_group)
                if (area_sum <= area_limit) {
                    // Testing
                    // logger_->report("[INFO] TOP GROUP_1 : {}; TOP GROUP_2 : {}", group_id_pair.first, group_id_pair.second);

                    // Valid merging for the selected two groups
                    MergeTwoGroups(bu_pri_queue, existing_groups, groups, hypergraph, area_limit);
                    dead_lock_counter = 0;
                    continue;
                } else {
                    // Testing 
                    // logger_->report("\t[INFO] Group Pair ({}, {}) exceeds the area or size limit, now updating PQ for group {}",
                    //                 group_id_pair.first,
                    //                 group_id_pair.second,
                    //                 group_id_pair.first);

                    bu_pri_queue.pop();

                    // Increate the deadlock counter
                    dead_lock_counter++;

                    // Check the deadlock counter
                    if (dead_lock_counter > max_retry_numer_) {

                        // Testing
                        // logger_->report("\t[INFO] REMOVE Group Pair ({}, {})",
                        //                 group_id_pair.first,
                        //                 group_id_pair.second);
                        continue;
                    }

                    UpdatePriorityQueue(bu_pri_queue, group_id_pair.first, groups, area_limit, hypergraph, existing_groups);
                }
            } else {
                // First element not valid, Update the priority queue and continue
                std::pair<int, int> group_id_pair = bu_pri_queue.top().group_pair;

                // Testing
                // logger_->report("\t[INFO] Invalid top priority queue element");

                bu_pri_queue.pop();

                // Check whether the group still exists or not
                if (existing_groups.find(group_id_pair.first) == existing_groups.end()) {
                    // Testing
                    // logger_->report("\t\t[INFO] Group {} has already been removed", group_id_pair.first);

                    continue;
                }

                UpdatePriorityQueue(bu_pri_queue, group_id_pair.first, groups, area_limit, hypergraph, existing_groups);
                continue;
            }
        }
    }

    // Merge two selected groups, we will keep the first group and remove the second one
    // Following structure has to be modified as for conducting the merging process
    // For the selected two groups: <group_1, group_2>
    // After the merging process, only group_1 will be kept
    // We store vertex id in each group, original group id = original_vertex_id - num_fixed_vertex
    void HierBlobPlacement::MergeTwoGroups(
        AccessiblePriorityQueue<PQElement>& input_queue,
        std::set<int>& existing_groups,
        std::vector<Group>& groups,
        Hypergraph& hypergraph,
        float area_limit
    ) { 
        // Variable definition
        int num_fixed_vertices = hypergraph.cluster_vertex_id_map.size();

        // Step 0: Get the top element
        PQElement top_element = input_queue.top();
        std::pair<int, int> top_element_pair = top_element.group_pair;
        const int group_1 = top_element_pair.first;
        const int group_2 = top_element_pair.second;
        const int vertex_id_1 = group_1 + num_fixed_vertices;
        const int vertex_id_2 = group_2 + num_fixed_vertices;

        // Testing
        // logger_->report("[INFO] Merge group {} with group {}",
        //                 group_1,
        //                 group_2);
        
        // logger_->report("\t[STATUS] Change all neighboring information for group {}", group_2);
        // std::string test_string = "";
        // test_string += "\tGroup " + std::to_string(group_1) + " included vertices : ";
        // for (auto sel_vertex : groups[group_1].included_vertices) {
        //     test_string += std::to_string(sel_vertex);
        //     test_string += " ";
        // }
        // logger_->report(test_string);
        // test_string = "";
        // test_string += "\tGroup " + std::to_string(group_2) + " included vertices : ";
        // for (auto sel_vertex : groups[group_2].included_vertices) {
        //     test_string += std::to_string(sel_vertex);
        //     test_string += " ";
        // }
        // logger_->report(test_string);

        // Testing
        // test_string = "";
        // test_string += "\tGroup " + std::to_string(group_1) + " neighbor vertices(before) : ";
        // for (auto sel_vertex : groups[group_1].neighbor_vetex_id) {
        //     test_string += std::to_string(sel_vertex);
        //     test_string += " ";
        // }
        // logger_->report(test_string);
        // test_string = "";
        // test_string += "\tGroup " + std::to_string(group_2) + " neighbor vertices(before) : ";
        // for (auto sel_vertex : groups[group_2].neighbor_vetex_id) {
        //     test_string += std::to_string(sel_vertex);
        //     test_string += " ";
        // }
        // logger_->report(test_string);

        // Step 1: Change all negibor vertex information for groups that are related with group 2
        for (auto sel_vertex_id : groups[group_2].neighbor_vetex_id) {
            // Get the real group id
            const int sel_group_id = sel_vertex_id - num_fixed_vertices;

            // Testing
            // logger_->report("\t[DEBUG] Remove vertex {} from group {}", vertex_id_2, sel_group_id);

            // Replace group 2 id with group 1 id
            groups[sel_group_id].neighbor_vetex_id.erase(vertex_id_2);

            // If the selected group equals to the first group, skip adding
            if (sel_group_id == group_1) continue;
            groups[sel_group_id].neighbor_vetex_id.insert(vertex_id_1);

            // Testing
            // logger_->report("\t[DEBUG] Insert vertex {} to group {}", vertex_id_1, sel_group_id);
        }

        // Testing
        // test_string = "";
        // test_string += "\tGroup " + std::to_string(group_1) + " neighbor vertices(after) : ";
        // for (auto sel_vertex : groups[group_1].neighbor_vetex_id) {
        //     test_string += std::to_string(sel_vertex);
        //     test_string += " ";
        // }
        // logger_->report(test_string);
        // test_string = "";
        // test_string += "\tGroup " + std::to_string(group_2) + " neighbor vertices(after) : ";
        // for (auto sel_vertex : groups[group_2].neighbor_vetex_id) {
        //     test_string += std::to_string(sel_vertex);
        //     test_string += " ";
        // }
        // logger_->report(test_string);


        // Step 2: Add included vertices from group 2 to group 1
        // Add the included vertices first
        groups[group_1].included_vertices.insert(groups[group_1].included_vertices.end(), 
                                                 groups[group_2].included_vertices.begin(),
                                                 groups[group_2].included_vertices.end());
        
        // Added the neighbor group_ids from group 2 to group 1
        for (auto sel_vertex_id : groups[group_2].neighbor_vetex_id) {
            // The sel_vertex_id shouldn't be included in the target group.
            if (std::find(groups[group_1].included_vertices.begin(), groups[group_1].included_vertices.end(), sel_vertex_id) != groups[group_1].included_vertices.end()) {
                continue;
            }

            groups[group_1].neighbor_vetex_id.insert(sel_vertex_id);

        }
        
        // Step 3: Add group 2 area size to group 1
        groups[group_1].group_size += groups[group_2].group_size;

        // Step 4: Add included hyperedge from group 2 to group 1
        groups[group_1].connected_hyperedge_id.insert(groups[group_2].connected_hyperedge_id.begin(),
                                                      groups[group_2].connected_hyperedge_id.end());

        // Step 5: Change the mapping of vertices in group 2 to group 1
        for (auto vertex_id : groups[group_2].included_vertices) {
            hypergraph.vertex_id_group_id_map[vertex_id] = group_1;
        }

        // Testing
        // logger_->report("\t[STATUS] Finished Changing, Remove group {} from existing_groups", group_2);

        // Step 6: Remove group 2 from existing group set
        existing_groups.erase(group_2);

        // Testing
        // logger_->report("\t[INFO] After merging: Num existing groups: {}", existing_groups.size());

        // Step 7: Invalidate all entries in the priority queue that contains group 1 or group 2
        input_queue.pop(); // Remove the first element in the priority queue
        // Access the underlying container of the priority queue
        auto& container = input_queue.get_container();

        // Iterate through the elements in the underlying container
        for (auto& elem : container) {
            if (elem.group_pair.first  == group_1 ||
                elem.group_pair.second == group_1 ||
                elem.group_pair.first  == group_2 ||
                elem.group_pair.second == group_2) {
                // Invalidate the entry
                elem.valid = false;
            }
        }

        // Testing
        // logger_->report("\t[INFO] PQ Status");
        // PrintPriorityQueue(input_queue);

        // Step 8: Update the priority with the new pair for group 1, if possible.
        UpdatePriorityQueue(input_queue, group_1, groups, area_limit, hypergraph, existing_groups); 
    }
    

    // Update the value for a specific group_id in the priority queue
    // While satisfying the area constraint balance
    void HierBlobPlacement::UpdatePriorityQueue(AccessiblePriorityQueue<PQElement>& input_queue,
                                                int group_id,
                                                std::vector<Group>& groups,
                                                float area_limit,
                                                Hypergraph& hypergraph,
                                                std::set<int>& existing_groups
                                                ) {
        // Define needed variables
        float max_value = -1;
        std::pair<int, int> max_value_key;

        // Testing
        // logger_->report("\t[INFO] Updating Priority queue for group {}", group_id);

        // Define the variables
        int num_fixed_vertices = hypergraph.cluster_vertex_id_map.size();

        // Define the needed map structure
        std::map<std::pair<int, int>, float> group_pair_clustering_score_map;
        std::map<std::pair<int, int>, float> group_pair_max_net_slack_map;

        // Check all related hyperedges
        for (const auto hyperedge_id : groups[group_id].connected_hyperedge_id) {
            std::set<int> unique_group_id;
            // Obtain all unique group_ids from the selected hyperedge
            for (const auto sel_vertex_id : hypergraph.hyperedges[hyperedge_id]) {
                if (sel_vertex_id >= num_fixed_vertices) {
                    int temp_group_id = hypergraph.vertex_id_group_id_map[sel_vertex_id];
                    if (temp_group_id != group_id) unique_group_id.insert(temp_group_id);
                }
            }

            if (unique_group_id.size() == 0) {

                // Testing
                // logger_->report("\t[WARNING] No unique neighbor found for group {}", group_id);

                continue;
            } else {
                // Check the area constraints of the two groups
                for (auto sel_group : unique_group_id) {
                    float group_area_1 = groups[group_id].group_size;
                    float group_area_2 = groups[sel_group].group_size;

                    // Check counter limit
                    int group_size_1 = groups[group_id].included_vertices.size();
                    int group_size_2 = groups[sel_group].included_vertices.size();

                    // TODO: change size upper bound here
                    if ((group_area_1 + group_area_2 <= area_limit)) {
                        // Selected pair satisfy the area limit constraint
                        std::pair<int, int> final_group_pair(group_id, sel_group);

                        group_pair_clustering_score_map[final_group_pair] 
                                += hypergraph.hyperedge_weight[hyperedge_id] / hypergraph.hyperedges[hyperedge_id].size();
                        group_pair_max_net_slack_map[final_group_pair] 
                                = std::max(group_pair_max_net_slack_map[final_group_pair], 
                                           hypergraph.hyperedge_id_net_slack_weight_map[hyperedge_id] / hypergraph.hyperedges[hyperedge_id].size());
                    }
                }
            }
            
        }

        // Check whether there is a feasible solution or not
        if (group_pair_clustering_score_map.size() == 0) {
            return;
        } else {
            // Valid pair found
            for (auto& item : group_pair_clustering_score_map) {
                // Get the key
                std::pair<int, int> group_pair_key(item.first);

                // Update the clustering score
                item.second += group_pair_max_net_slack_map[group_pair_key];

                // Find the maximum value
                if (item.second > max_value) {
                    max_value = item.second;
                    max_value_key = item.first;
                }
            }
        }
        
        // Insert new element
        if (max_value > 0) {
            
            // Testing
            // logger_->report("\t\t[INFO] New pair ({}, {}), clustering score: {}", max_value_key.first, max_value_key.second,  max_value);

            input_queue.push(PQElement(max_value_key, max_value, true));
        }
    }

    // This function is used to refine the best-choice clustering process
    // After the best choice clustering process, some groups maybe left due to the connected groups
    // already exceeds the area or number limit, we can merge those groups together.
    void HierBlobPlacement::RefineGrouping(std::vector<Group>& groups, 
                                           std::set<int>& unique_groups, 
                                           float limit,
                                           float number_fixed_vertices) {
        // Variable definiton
        std::vector<int> small_groups;
        int deadlock_counter = 0;
        std::random_device rd;
        std::mt19937 rand_gen(seed_);
        int small_group_threshold = limit * (1 - tolerance_);

        // Find all small groups
        for (auto sel_group_id : unique_groups) {
            if (groups[sel_group_id].group_size < small_group_threshold) {
                small_groups.push_back(sel_group_id);
            }
        }

        // Check the necessity
        if (small_groups.size() < max_retry_numer_) return;

        // Testing
        // logger_->report("[INFO] {} small groups found", small_groups.size());
        // std::cout << "\tSmall groups : ";
        // for (auto sel_group_id : small_groups) {
        //     std::cout << sel_group_id << " ";
        // }
        // std::cout << std::endl;

        // Start merging small groups
        while(deadlock_counter < 50) {
            std::uniform_int_distribution<> distribution(0, small_groups.size() - 1);

            // Check the size of the group
            if (small_groups.size() <= 1) break;

            // Randomly select two groups from the detected small groups
            int index_1 = distribution(rand_gen);
            int index_2 = distribution(rand_gen);
            while (index_1 == index_2) {
                index_2 = distribution(rand_gen);
            }

            int group_1 = small_groups[index_1];
            int group_2 = small_groups[index_2];

            // Check whether this two group can be merged
            if (groups[group_1].group_size + groups[group_2].group_size < small_group_threshold) {
                // Testing
                // logger_->report("[INFO] Merge group {} with group {}",
                //                 group_1,
                //                 group_2);

                // Variable definition
                const int vertex_id_1 = group_1 + number_fixed_vertices;
                const int vertex_id_2 = group_2 + number_fixed_vertices;

                // Clear the counter
                deadlock_counter = 0;

                // Merge the selected two groups
                // Step 1: Change all negibor vertex information for groups that are related with group 2
                for (auto sel_vertex_id : groups[group_2].neighbor_vetex_id) {
                    // Get the real group id
                    const int sel_group_id = sel_vertex_id - number_fixed_vertices;

                    // Replace group 2 id with group 1 id
                    groups[sel_group_id].neighbor_vetex_id.erase(vertex_id_2);

                    // If the selected group equals to the first group, skip adding
                    if (sel_group_id == group_1) continue;
                    groups[sel_group_id].neighbor_vetex_id.insert(vertex_id_1);
                }

                // Step 2: Add included vertices from group 2 to group 1
                // Add the included vertices first
                groups[group_1].included_vertices.insert(groups[group_1].included_vertices.end(), 
                                                 groups[group_2].included_vertices.begin(),
                                                 groups[group_2].included_vertices.end());
                
                // Added the neighbor group_ids from group 2 to group 1
                for (auto sel_vertex_id : groups[group_2].neighbor_vetex_id) {
                    // The sel_vertex_id shouldn't be included in the target group.
                    if (std::find(groups[group_1].included_vertices.begin(), groups[group_1].included_vertices.end(), sel_vertex_id) != groups[group_1].included_vertices.end()) {
                        continue;
                    }

                    groups[group_1].neighbor_vetex_id.insert(sel_vertex_id);
                }

                // Step 3: Add group 2 area size to group 1
                groups[group_1].group_size += groups[group_2].group_size;

                // Step 4: Add included hyperedge from group 2 to group 1
                groups[group_1].connected_hyperedge_id.insert(groups[group_2].connected_hyperedge_id.begin(),
                                                      groups[group_2].connected_hyperedge_id.end());

                // Step 5: Remove group 2 from existing group set
                unique_groups.erase(group_2);
                auto it = std::remove(small_groups.begin(), small_groups.end(), group_2);
                small_groups.erase(it, small_groups.end());

                // Testing
                // logger_->report("\t[INFO] After merging: Num existing groups: {}", unique_groups.size());
            }

            // Increase the counter
            deadlock_counter++;
        }

    }

    // This function is used to initialize the priority queue used for bottom-up
    // clustering
    void HierBlobPlacement::InitializePriorityQueue(
        Hypergraph& hypergraph,
        AccessiblePriorityQueue<PQElement>& storing_queue,
        std::vector<Group>& groups
    ) { 
        // Testing
        logger_->report("[STATUS] Initialize the priority queue");

        // Define the variables
        int num_fixed_vertices = hypergraph.cluster_vertex_id_map.size();

        // Testing
        // logger_->report("\tNumber Fixed vertices: {}", num_fixed_vertices);

        // Update the clustering score for all included vertices
        for (const auto inst_vertex_id_pair : hypergraph.inst_vertex_id_map) {
            // Get the vertex_id and group_id
            int vertex_id = inst_vertex_id_pair.second;
            int group_index = vertex_id - num_fixed_vertices;

            // Testing
            // logger_->report("\t[INFO] Check entry for vertex : {}, group : {}", vertex_id, group_index);

            // Define needed map structures
            std::map<std::pair<int, int>, float> group_pair_clustering_score_map;
            std::map<std::pair<int, int>, float> group_pair_max_net_slack_map;

            // Check all hyperedges that include the selected vertex
            for (const auto sel_hyperedge_id : groups[group_index].connected_hyperedge_id) {
                // Check all vertices in the sel_hyperedge
                for (const auto second_vertex_id : hypergraph.hyperedges[sel_hyperedge_id]) {
                    if (second_vertex_id != vertex_id && (second_vertex_id >= num_fixed_vertices)) {
                        // Vertex_pair found
                        int group_index_2 = second_vertex_id - num_fixed_vertices;
                        std::pair<int, int> final_group_pair(group_index, group_index_2);

                        // TODO: Formalize the initialization process for those map variables
                        group_pair_clustering_score_map[final_group_pair] 
                                += hypergraph.hyperedge_weight[sel_hyperedge_id] / hypergraph.hyperedges[sel_hyperedge_id].size();
                        group_pair_max_net_slack_map[final_group_pair] 
                                = std::max(group_pair_max_net_slack_map[final_group_pair], 
                                           hypergraph.hyperedge_id_net_slack_weight_map[sel_hyperedge_id] / hypergraph.hyperedges[sel_hyperedge_id].size());

                        // Testing
                        // logger_->report("\t\t(Group 1, Group 2): ({}, {}), cluster_score: {}, max_net_slack_score: {}",
                        //                 group_index,
                        //                 group_index_2,
                        //                 group_pair_clustering_score_map[final_group_pair],
                        //                 group_pair_max_net_slack_map[final_group_pair]);
                    }
                }
            }

            // Get the complete clustering score and find the maximum one
            float max_value = 0.0;
            std::pair<int, int> max_value_key;

            for (auto& item : group_pair_clustering_score_map) {
                // Get the key
                std::pair<int, int> group_pair_key(item.first);

                // Update the clustering score
                item.second += group_pair_max_net_slack_map[group_pair_key];

                // Find the maximum value
                if (item.second > max_value) {
                    max_value = item.second;
                    max_value_key = item.first;
                }                
            }

            // Testing 
            // logger_->report("\t\t[INFO] Max pair ({}, {}), clustering score: {}", max_value_key.first, max_value_key.second,  max_value);

            // Testing
            // PQElement temp_element(max_value_key, max_value, true);
            // logger_->report("\t\t[INFO] Max value: {}", max_value);
            // logger_->report("\t\t[INFO] Insert priority queue element, ({}, {}), clustering score: {}, Validity: {}", temp_element.group_pair.first, temp_element.group_pair.second, temp_element.cluster_score, temp_element.valid);

            // Insert the key_value pair to priority queue
            storing_queue.push(PQElement(max_value_key, max_value, true));
        }

        // Testing
        // logger_->report("[INFO] PQ Status");
        // PrintPriorityQueue(storing_queue);
        
        // Testing
        logger_->report("[STATUS] Finished initializing the priority queue");
    }

    // Print all elements in the priority queue
    void HierBlobPlacement::PrintPriorityQueue(AccessiblePriorityQueue<PQElement> storing_queue) {
        logger_->report("[INFO] Priority Queue Contents");

        // Iterate througn the priority queue
        while (!storing_queue.empty()) {
            PQElement temp_element = storing_queue.top();
            logger_->report("   Group_pair: <{}, {}>; Clustering Score: {}; Valid: {}",
                            temp_element.group_pair.first,
                            temp_element.group_pair.second,
                            temp_element.cluster_score,
                            temp_element.valid);
            
            storing_queue.pop();
        }

        logger_->report("[INFO] End of the Priority Queue");
    }


    // Calculate Connections between clusters - Star Model used to break multi-pin nets
    void HierBlobPlacement::CalculateConnection() {
        // Initialize the connections of all clusters
        // Clear the connection map
        for (auto& [cluster_id, cluster] : cluster_map_) {
            cluster->InitConnection();
        }

        // Traverse all nets through OpenDB
        for (odb::dbNet* net : block_->getNets()) {
            // Ignore all power net
            if (net->getSigType().isSupply()) {
                continue;
            }
            int driver_id = -1; 
            // TODO: Use set instead of vector in the debug stage
            std::vector<int> loads_id;
            bool pad_flag = false;

            // Check all connected instances
            for (odb::dbITerm* iterm : net->getITerms()) {
                odb::dbInst* inst = iterm->getInst();
                const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
                // Physical only?
                if (liberty_cell == nullptr) continue;

                odb::dbMaster* master = inst->getMaster();
                // Check whether the inst is a Pad, Cover.
                if (master->isPad() || master->isCover()) {
                    pad_flag = true;
                    break;
                }

                // Testing
                // logger_->report("[DEBUG] Find cluster_id for instance {}", inst->getName());

                const int cluster_id
                    = odb::dbIntProperty::find(inst, "cluster_id")->getValue();

                // Testing
                // logger_->report("\t[DEBUG] Cluster ID found for {}", inst->getName());

                if (iterm->getIoType() == odb::dbIoType::OUTPUT) {
                    driver_id = cluster_id;
                } else {
                    loads_id.push_back(cluster_id);
                }
            }
            if (pad_flag) {
                continue;
            }
            bool io_flag = false;
            // Check the connected IO pins
            for (odb::dbBTerm* bterm : net->getBTerms()) {
                // Testing
                // logger_->report("[DEBUG] Find cluster_id for block terminal {}", bterm->getName());

                const int cluster_id
                    = odb::dbIntProperty::find(bterm, "cluster_id")->getValue();

                // Testing
                // logger_->report("\t[DEBUG] Cluster ID found for {}", bterm->getName());

                io_flag = true;
                if (bterm->getIoType() == odb::dbIoType::INPUT) {
                    driver_id = cluster_id;
                } else {
                    loads_id.push_back(cluster_id);
                }
            }

            // Add the net to connections between clusters
            if (driver_id != -1 && loads_id.size() > 0
                    && loads_id.size() < large_net_threshold_) {
                const float weight = (io_flag == true) ? virtual_weight_ : 1.0; 
                for (int i = 0; i < loads_id.size(); i++) {
                    if (loads_id[i] != driver_id) {
                        // All connections are undirected edges
                        cluster_map_[driver_id]->AddConnection(loads_id[i], weight);
                        cluster_map_[loads_id[i]]->AddConnection(driver_id, weight);
                    }
                } // End Sink updates
            }     // End Adding current net
        }         // End Net Traversal
    }

    // Calculate Connections between clusters - Clique Model used to break multi-pin nets
    // This function will store all the information into the BundledNet
    // Two mapping schemes are used
    // 1. Cluster level: map from neighbor cluster_id to the corresponding bundled net index
    // 2. Global: dbNet to bundled net index map
    void HierBlobPlacement::CalculateCliqueConnection() {
        // Clean the map
        ClearBundledNetMap();

        // Clear the bundled net connection
        for (auto& [cluster_id, cluster] : cluster_map_) {
            cluster->InitBundledNetConnection();
        }

        // Testing
        logger_->report("[DEBUG] Start calculating connections between different clusters based on the clique model");

        // Traverse all nets throught OpenDB
        // All multipin net will be broken by a clique model
        for (odb::dbNet* net : block_->getNets()) {
            // Ignore all power net
            if (net->getSigType().isSupply()) {
                continue;
            }

            int driver_id = -1;
            // Define the vector to store all pin locations
            std::vector<int> pins_id;
            bool pad_flag = false;

            // Check all connected instances
            for (odb::dbITerm* iterm : net->getITerms()) {
                odb::dbInst* inst = iterm->getInst();
                const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
                // Physical only?
                if (liberty_cell == nullptr) continue;

                odb::dbMaster* master = inst->getMaster();
                // Check whether the inst is a Pad or Cover
                if (master->isPad() || master->isCover()) {
                    pad_flag = true;
                    break;
                }

                const int cluster_id
                    = odb::dbIntProperty::find(inst, "cluster_id")->getValue();
                // Add the corresponding cluster to the string vector, if not in it
                if (std::find(pins_id.begin(), pins_id.end(), cluster_id) == pins_id.end()) {
                    pins_id.push_back(cluster_id);
                }
            }

            if (pad_flag) {
                    continue;
            }
            // Check the connected IO pins
            for (odb::dbBTerm* bterm : net->getBTerms()) {
                const int cluster_id
                    = odb::dbIntProperty::find(bterm, "cluster_id")->getValue();

                // Check the validity of the selected cluster
                if (std::find(pins_id.begin(), pins_id.end(), cluster_id) == pins_id.end()) {
                    pins_id.push_back(cluster_id);
                }
            }

            // Add Connections to all records
            if (pins_id.size() > 1 && pins_id.size() < large_net_threshold_) {
                // Create the all bundled net
                std::vector<int> net_index_list;
                for (int i = 0; i < pins_id.size(); ++i) {
                    for (int j = i + 1; j < pins_id.size(); ++j) {
                        // Create the needed Bundled net structure and store the indexing structure
                        BundledNet temp_net(pins_id[i], pins_id[j], 2.0 / pins_id.size(), 1);
                        int temp_net_id = sparse_nets_.size();

                        // Testing
                        // logger_->report("[INFO] Created BundledNet between cluster {} and {}", 
                        //                 cluster_map_[pins_id[i]]->GetName(),
                        //                 cluster_map_[pins_id[j]]->GetName());
                        // std::cout << (2.0 / pins_id.size()) * 100 << std::endl;

                        // Store the net index in the local cluster structure
                        cluster_map_[pins_id[i]]->AddClusterBundledNet(pins_id[j], temp_net_id);
                        // cluster_map_[pins_id[i]]->AddNetBundledNet(net, temp_net_id);
                        cluster_map_[pins_id[j]]->AddClusterBundledNet(pins_id[i], temp_net_id);
                        // cluster_map_[pins_id[j]]->AddNetBundledNet(net, temp_net_id);

                        // Store the net index in the global structure
                        sparse_nets_.push_back(temp_net);
                        net_index_list.push_back(temp_net_id);
                    }
                }

                // Add the net index list to the global storing map
                net_bundled_net_map_.insert(std::make_pair(net, net_index_list));
            }
        }
    }

    // ***************************************************************************
    // * Metrics
    // ***************************************************************************
    /**
     * @brief Traverse Logical Hierarchy. Recursive function to collect the design metrics (num_std_cells, area of std cells)
     *        in the logical hierarchy.
     *        For now we still assume the design has no rams(macros)
    */
    Metrics* HierBlobPlacement::ComputeMetrics(odb::dbModule* module) {
        // Variable definition
        unsigned int num_std_cell = 0;
        float std_cell_area = 0.0;

        for (odb::dbInst* inst: module->getInsts()) {
            const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
            // If the creation failed
            if (liberty_cell == nullptr) continue;

            // Get the master cell of the selected instance
            odb::dbMaster* master = inst->getMaster();

            // Check whether the instance is a pad or cover
            if (master->isPad() || master->isCover()) {
                continue;
            }

            // Get the area, note that if we want to get the area
            // from odb, we have to use the getWidth() and getHeight()
            // no getArea available.
            float inst_area = liberty_cell->area();

            // For now we only support std cells
            num_std_cell += 1;
            std_cell_area += inst_area;
        }

        // Revursively traverse the hierarchical module instances
        for (odb::dbModInst* inst : module->getChildren()) {
            Metrics* metrics = ComputeMetrics(inst->getMaster());
            num_std_cell += metrics->GetNumStdCell();
            std_cell_area += metrics->GetStdCellArea();
        }

        Metrics* metrics
            = new Metrics(num_std_cell, std_cell_area);
        // Add the newly generated metrics to the logical module map
        logical_module_map_[module] = metrics;

        return metrics;
    }

    // Compute the metrics for a cluster
    void HierBlobPlacement::SetClusterMetrics(Cluster* cluster) {
        unsigned int num_std_cell = 0;
        float std_cell_area = 0.0;

        num_std_cell += cluster->GetLeafStdCells().size();
        
        for (odb::dbInst* inst : cluster->GetLeafStdCells()) {
            const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
            std_cell_area += liberty_cell->area();
        }

        Metrics metrics(num_std_cell, std_cell_area);

        for (auto& module : cluster->GetDbModules()) {
            metrics.AddMetrics(*logical_module_map_[module]);
        }

        // Testing
        // std::cout << "[DEBUG] Set Cluster Metrics for cluster: " << cluster->GetName() << std::endl;
        // std::cout << "[DEBUG] Num Std Cell: " << metrics.GetNumStdCell() << std::endl;
        // std::cout << "[DEBUG] Cluster Size: " << metrics.GetStdCellArea() << std::endl;

        // Update metrics
        cluster->SetMetrics(metrics);
    }

    // ***************************************************************************
    // * Simulated Annealing Related Functions
    // ***************************************************************************
    // This function is used to set the cluster id to the corresponding 
    // modules or std cells in the selected level in a BFS manner
    void HierBlobPlacement::ActivateClusters(int level) {
        // Variable Definition
        // Testing
        logger_->report("[STATUS] Set cluster_id to all modules and lead std cells according to level {}", level);

        if (level == 0) {
            SetInstProperty(root_cluster_);
            return;
        }

        // Define queue structure for BFS search
        std::queue<Cluster*> bfs_queue;
        bfs_queue.push(root_cluster_);
        int current_level = 0;

        while(!bfs_queue.empty()) {
            int level_size = bfs_queue.size();

            for (int i = 0; i < level_size; i++) {
                Cluster* sel_cluster = bfs_queue.front();
                bfs_queue.pop();

                if (current_level == level) {
                    // Testing
                    // logger_->report("[INFO] Activate cluster {}", sel_cluster->GetName());

                    // Meet the target level, assign the cluster_id
                    // to corresponding module and leaf std cells
                    SetInstProperty(sel_cluster);
                }

                for (auto child : sel_cluster->GetChildren()) {
                    // Push all child clusters into the queue
                    bfs_queue.push(child);
                }
            }

            if (current_level == level) {
                break; // All clusters in the corresponding level accessed
            }

            ++current_level;
        }

        // Testing
        logger_->report("[STAUS] Finished activating clusters in level {};", level);

    }

    // This function cleans all storing variables for the sliding window structure
    void HierBlobPlacement::CleanSlidingWindowVariable() {
        sparse_dense_map_.clear();
        sparse_nets_.clear();
        dense_nets_.clear();
    }

    // This funciton clears the dbNet to bundled Net index map
    void HierBlobPlacement::ClearBundledNetMap() {
        net_bundled_net_map_.clear();
    }

    // We need to merge the sparse nets into real bundled net connections
    // Traverse all the nets and build the map
    void HierBlobPlacement::CondenseBundledNet() {
        // Variable definition
        std::map<std::pair<int, int>, int> terminal_net_index_map; // Map from termianl pair to the sparse net index
        int counter = 0;

        // Testing
        logger_->report("[STATUS] Start condensing the sparse bundled net");

        // Traverse all sparse net to condense them
        for (auto sel_net : sparse_nets_) {
            // Get the terminal pair
            std::pair<int, int> sel_pair(sel_net.src_cluster_id, sel_net.target_cluster_id);
            std::pair<int, int> reverse_sel_pair(sel_net.target_cluster_id, sel_net.src_cluster_id);

            // Check whether the pair is in the map or not
            if (terminal_net_index_map.find(sel_pair) != terminal_net_index_map.end()) {
                // The mapping exists
                int dense_net_index = terminal_net_index_map[sel_pair];

                dense_nets_[dense_net_index].timing_weight += sel_net.clique_factor * sel_net.timing_weight;

                // Add the mapping information
                sparse_dense_map_[counter] = dense_net_index;
            } else if (terminal_net_index_map.find(reverse_sel_pair) != terminal_net_index_map.end()) {
                // The mapping exists
                int dense_net_index = terminal_net_index_map[reverse_sel_pair];

                dense_nets_[dense_net_index].timing_weight += sel_net.clique_factor * sel_net.timing_weight;

                // Add the mapping information
                sparse_dense_map_[counter] = dense_net_index;
            } else {
                // The selected pair is not in the map
                terminal_net_index_map[sel_pair] = dense_nets_.size();

                // Create the new bundled net
                BundledNet temp_bundled_net(sel_net.src_cluster_id, sel_net.target_cluster_id, 1.0, sel_net.clique_factor * sel_net.timing_weight);
                dense_nets_.push_back(temp_bundled_net);

                // Add the mapping information
                sparse_dense_map_[counter] = dense_nets_.size() - 1;
            }

            counter++;
        }
    }

    // Extract the needed components for the annealer (both vertices and nets)
    void HierBlobPlacement::ExtractBoth(Bbox bbox, int level) {
        // Testing
        logger_->report("[INFO] Extracting moveable Clusters and corresponding nets with in BBOX: {} {} {} {}",
                        bbox.lx,
                        bbox.ly,
                        bbox.ux,
                        bbox.uy);

        // Clean the storing structure
        movable_vertices_.clear();
        selected_dense_nets_.clear();

        // Define mapping structure
        std::map<int, int> cluster_id_vertex_id_map; // Map from cluster id to vertex id

        // Traverse all dense nets
        for (auto sel_net : dense_nets_) {
            // Get the coordinates of both vertices
            int src_cluster_id = sel_net.src_cluster_id;
            int dst_cluster_id = sel_net.target_cluster_id;
            bool src_in_flag = false;
            bool dst_in_flag = false;

            // Get the corresponding coordinates
            float src_cluster_x = cluster_map_[src_cluster_id]->GetX();
            float src_cluster_y = cluster_map_[src_cluster_id]->GetY();
            float src_cluster_z = cluster_map_[src_cluster_id]->GetZ();

            float dst_cluster_x = cluster_map_[dst_cluster_id]->GetX();
            float dst_cluster_y = cluster_map_[dst_cluster_id]->GetY();
            float dst_cluster_z = cluster_map_[dst_cluster_id]->GetZ();

            // Check whether both of vertices are in the chosen Bbox
            if (src_cluster_x > bbox.lx &&
                src_cluster_x < bbox.ux &&
                src_cluster_y > bbox.ly &&
                src_cluster_y < bbox.uy) {
                src_in_flag = true;
            }

            if (level == 1 && !(cluster_map_[src_cluster_id]->GetIOClusterFlag())) src_in_flag = true;

            if (dst_cluster_x > bbox.lx &&
                dst_cluster_x < bbox.ux &&
                dst_cluster_y > bbox.ly &&
                dst_cluster_y < bbox.uy) {
                dst_in_flag = true;
            }

            if (level == 1 && !(cluster_map_[dst_cluster_id]->GetIOClusterFlag())) dst_in_flag = true;

            // Build vertex storing structure
            if (src_in_flag) {
                if (cluster_id_vertex_id_map.find(src_cluster_id) == cluster_id_vertex_id_map.end()) {
                    // Cluster not mapped, create a new one
                    int vertex_id = movable_vertices_.size();
                    Vertex temp_vertex(cluster_map_[src_cluster_id], vertex_id);

                    // Update the cluster_id_vertex_id_map
                    cluster_id_vertex_id_map[src_cluster_id] = vertex_id;

                    movable_vertices_.push_back(temp_vertex);

                    // Add the net to the created vertex
                    movable_vertices_[vertex_id].connected_net_.insert(selected_dense_nets_.size());
                } else {
                    // Cluster already mapped to a vertex
                    int vertex_id = cluster_id_vertex_id_map[src_cluster_id];

                    // Add the net to the created vertex
                    movable_vertices_[vertex_id].connected_net_.insert(selected_dense_nets_.size());
                }
            }

            if (dst_in_flag) {
                if (cluster_id_vertex_id_map.find(dst_cluster_id) == cluster_id_vertex_id_map.end()) {
                    // Cluster not mapped, create a new one
                    int vertex_id = movable_vertices_.size();
                    Vertex temp_vertex(cluster_map_[dst_cluster_id], vertex_id);

                    // Update the cluster_id_vertex_id_map
                    cluster_id_vertex_id_map[dst_cluster_id] = vertex_id;

                    movable_vertices_.push_back(temp_vertex);

                    // Add the net to the created vertex
                    movable_vertices_[vertex_id].connected_net_.insert(selected_dense_nets_.size());
                } else {
                    // Cluster already mapped to a vertex
                    int vertex_id = cluster_id_vertex_id_map[dst_cluster_id];

                    // Add the net to the created vertex
                    movable_vertices_[vertex_id].connected_net_.insert(selected_dense_nets_.size());
                }
            }

            // Create the net storing structure
            if (dst_in_flag && src_in_flag) {
                // Both vertices are movable vertices
                int net_index = selected_dense_nets_.size();
                int src_vertex_index = cluster_id_vertex_id_map[src_cluster_id];
                int dst_vertex_index = cluster_id_vertex_id_map[dst_cluster_id];

                Edge temp_edge(src_vertex_index,
                               dst_vertex_index,
                               net_index,
                               true,
                               -1,
                               -1,
                               -1,
                               sel_net.timing_weight);

                selected_dense_nets_.push_back(temp_edge);
            } else if (src_in_flag) {
                // Src vertex is movable
                int net_index = selected_dense_nets_.size();
                int src_vertex_index = cluster_id_vertex_id_map[src_cluster_id];

                Edge temp_edge(src_vertex_index,
                               -1,
                               net_index,
                               false,
                               dst_cluster_x,
                               dst_cluster_y,
                               dst_cluster_z,
                               sel_net.timing_weight);
                
                selected_dense_nets_.push_back(temp_edge);
            } else if (dst_in_flag) {
                // Dst vertex is movable
                int net_index = selected_dense_nets_.size();
                int dst_vertex_index = cluster_id_vertex_id_map[dst_cluster_id];

                Edge temp_edge(dst_vertex_index,
                               -1,
                               net_index,
                               false,
                               src_cluster_x,
                               src_cluster_y,
                               src_cluster_z,
                               sel_net.timing_weight);
                
                selected_dense_nets_.push_back(temp_edge);
            }

        }

        // Testing
        logger_->report("[INFO] Finished extracting moveable Clusters and corresponding nets with in BBOX: {} {} {} {}",
                        bbox.lx,
                        bbox.ly,
                        bbox.ux,
                        bbox.uy);
    }

    // Extract needed components for the annealer
    // In each round of SA, all clusters outside will be modelded as IO vertices
    // which will be fixed during annealing
    void HierBlobPlacement::ExtractClusters(int level, Bbox bbox) {
        // Testing
        logger_->report("[INFO] Extracting moveable Clusters with in BBOX: {} {} {} {}, level {}",
                        bbox.lx,
                        bbox.ly,
                        bbox.ux,
                        bbox.uy,
                        level);

        // Clear the storing structure
        movable_vertices_.clear();
        
        // Traverse clusters to extract those that meet the specified constraints
        // Define queue structure for BFS search
        std::queue<Cluster*> bfs_queue;
        bfs_queue.push(root_cluster_);
        int current_level = 0;

        // Define mapping structure
        std::map<int, int> cluster_id_vertex_id_map; // Map from cluster id to vertex id

        while(!bfs_queue.empty()) {
            int level_size = bfs_queue.size();

            for (int i = 0; i < level_size; i++) {
                Cluster* sel_cluster = bfs_queue.front();
                bfs_queue.pop();

                if (current_level == level) {
                    // Meet the target level, Check the boundary condition
                    if (level == 1) {
                        
                        // Check whether this is an IO cluster
                        if (sel_cluster->GetIOClusterFlag()) {
                            continue;
                        } else {
                            // Check whether this cluster has already been mapped to a vertex
                            if (cluster_id_vertex_id_map.find(sel_cluster->GetId()) == cluster_id_vertex_id_map.end()) {
                                // Cluster not mapped, create a new one
                                int vertex_id = movable_vertices_.size();
                                Vertex temp_vertex(sel_cluster, vertex_id);

                                // Testing
                                // logger_->report("\t[INFO] Creating vertex {} based on cluster {}", vertex_id, sel_cluster->GetId());

                                // Testing
                                if (sel_cluster->GetId() == 0) {
                                    logger_->report("[ERROR] Root cluster detected, should not happen!");
                                }

                                // Update the cluster_id_vertex_id_map
                                cluster_id_vertex_id_map[sel_cluster->GetId()] = vertex_id;

                                movable_vertices_.push_back(temp_vertex);

                                // Get the neighbor connections
                                std::map<int, std::vector<int>> temp_map; 
                                sel_cluster->GetBundledNetConnections(temp_map);

                                // Add neighbors to the newly created vertex
                                for (auto sel_neighbor_pair : temp_map) {
                                    // Get the cluster id
                                    int neighbor_cluster_id = sel_neighbor_pair.first;
                                    int vertex_id_2 = -1;

                                    // Check whether this is an IO cluster
                                    if (!cluster_map_[neighbor_cluster_id]->GetIOClusterFlag()) {

                                        if (cluster_id_vertex_id_map.find(neighbor_cluster_id) == cluster_id_vertex_id_map.end()) {
                                            vertex_id_2 = movable_vertices_.size();

                                            Vertex temp_vertex_2(cluster_map_[neighbor_cluster_id], vertex_id_2);
                                            movable_vertices_.push_back(temp_vertex_2);

                                            // Testing
                                            // logger_->report("\t\t[INFO] Creating vertex {} based on cluster {} which is a neighbor to vertex {}", vertex_id_2, cluster_map_[neighbor_cluster_id]->GetId(), vertex_id);

                                            // Testing
                                            if (neighbor_cluster_id == 0) {
                                                logger_->report("[ERROR] Root cluster detected, should not happen!");
                                            }

                                            // Update the mapping
                                            cluster_id_vertex_id_map[neighbor_cluster_id] = vertex_id_2;
                                        } else {
                                            vertex_id_2 = cluster_id_vertex_id_map[neighbor_cluster_id];
                                        }

                                        // Testing
                                        // logger_->report("[INFO] Add connection from cluster {} to {}", sel_cluster->GetId(), neighbor_cluster_id);
                                        for (auto net_id : sel_neighbor_pair.second) {
                                            if (sparse_dense_map_[net_id] != sparse_dense_map_[sel_neighbor_pair.second[0]]) {
                                                logger_->report("[ERROR] WRONG sparse to dense bundled net mapping from cluster {} to cluster {}",
                                                                sel_cluster->GetId(),
                                                                neighbor_cluster_id);
                                                logger_->report("\tExpected dense net id {}, got dense net id {}", sparse_dense_map_[sel_neighbor_pair.second[0]], sparse_dense_map_[net_id]);
                                                logger_->report("\tOriginal Sparse net id: {}; Condensed net id: {}\n"
                                                            "\tOriginal Sparse net: src {}, target {}",
                                                            net_id,
                                                            sparse_dense_map_[net_id],
                                                            sparse_nets_[net_id].src_cluster_id,
                                                            sparse_nets_[net_id].target_cluster_id);
                                            }
                                        }

                                        // Testing
                                        // logger_->report("\t\t[INFO] Add neighboring vertex {} (Cluster {})", vertex_id_2, neighbor_cluster_id);

                                        // Get the new dense net id
                                        int dense_net_id = sparse_dense_map_[sel_neighbor_pair.second[0]];
                                        movable_vertices_[vertex_id].moveable_vertex_id_weight_map_[vertex_id_2] = dense_nets_[dense_net_id].timing_weight;
                                    } else {
                                        // This is an IO cluster
                                        int dense_net_id = sparse_dense_map_[sel_neighbor_pair.second[0]];

                                        // Check whether the dense net is already included
                                        if (movable_vertices_[vertex_id].connected_net_.find(dense_net_id) == movable_vertices_[vertex_id].connected_net_.end()) {
                                            // This net has not been added
                                            movable_vertices_[vertex_id].connected_net_.insert(dense_net_id);

                                            // Update the coordinates storing structure inside the corresponding vertex
                                            float second_cluster_x = cluster_map_[neighbor_cluster_id]->GetX();
                                            float second_cluster_y = cluster_map_[neighbor_cluster_id]->GetY();
                                            float second_cluster_z = cluster_map_[neighbor_cluster_id]->GetZ();
                                            std::vector<float> temp_fixed_vertex_info;
                                            temp_fixed_vertex_info.push_back(second_cluster_x);
                                            temp_fixed_vertex_info.push_back(second_cluster_y);
                                            temp_fixed_vertex_info.push_back(second_cluster_z);
                                            temp_fixed_vertex_info.push_back(dense_nets_[dense_net_id].timing_weight);
                                            movable_vertices_[vertex_id].fixed_vertex_connection_map_.push_back(temp_fixed_vertex_info);
                                        }
                                    }
                                }
                                
                            } else {
                                // The corresponding vertex has already been created
                                int vertex_id = cluster_id_vertex_id_map[sel_cluster->GetId()];

                                // Get neighbor connections
                                std::map<int, std::vector<int>> temp_map; 
                                sel_cluster->GetBundledNetConnections(temp_map);

                                // Add neighbors to the existing vertex
                                for (auto sel_neighbor_pair : temp_map) {
                                    // Get the cluster id
                                    int neighbor_cluster_id = sel_neighbor_pair.first;
                                    int vertex_id_2 = -1;

                                    // Check whether this is an IO cluster
                                    if (!cluster_map_[neighbor_cluster_id]->GetIOClusterFlag()) {
                                        if (cluster_id_vertex_id_map.find(neighbor_cluster_id) == cluster_id_vertex_id_map.end()) {
                                            vertex_id_2 = movable_vertices_.size();

                                            Vertex temp_vertex_2(cluster_map_[neighbor_cluster_id], vertex_id_2);
                                            movable_vertices_.push_back(temp_vertex_2);

                                            // Testing
                                            // logger_->report("\t\t[INFO] Creating vertex {} based on cluster {} which is a neighbor to vertex {}", vertex_id_2, cluster_map_[neighbor_cluster_id]->GetId(), vertex_id);

                                            // Testing
                                            if (neighbor_cluster_id == 0) {
                                                logger_->report("[ERROR] Root cluster detected, should not happen!");
                                            }

                                            // Update the mapping
                                            cluster_id_vertex_id_map[neighbor_cluster_id] = vertex_id_2;
                                        } else {
                                            vertex_id_2 = cluster_id_vertex_id_map[neighbor_cluster_id];
                                        }

                                        // Testing
                                        // logger_->report("[INFO] Add connection from cluster {} to {}", sel_cluster->GetId(), neighbor_cluster_id);

                                        // Testing
                                        for (auto net_id : sel_neighbor_pair.second) {
                                            if (sparse_dense_map_[net_id] != sparse_dense_map_[sel_neighbor_pair.second[0]]) {
                                                logger_->report("[ERROR] wrong sparse to dense bundled net mapping from cluster {} to cluster {}",
                                                                sel_cluster->GetId(),
                                                                neighbor_cluster_id);
                                                logger_->report("\tExpected dense net id {}, got dense net id {}", sparse_dense_map_[sel_neighbor_pair.second[0]], sparse_dense_map_[net_id]);
                                                logger_->report("\tOriginal Sparse net id: {}; Condensed net id: {}\n"
                                                            "\tOriginal Sparse net: src {}, target {}",
                                                            net_id,
                                                            sparse_dense_map_[net_id],
                                                            sparse_nets_[net_id].src_cluster_id,
                                                            sparse_nets_[net_id].target_cluster_id);
                                            }
                                        }

                                        // Get the new dense net id
                                        int dense_net_id = sparse_dense_map_[sel_neighbor_pair.second[0]];
                                        movable_vertices_[vertex_id].moveable_vertex_id_weight_map_[vertex_id_2] = dense_nets_[dense_net_id].timing_weight;
                                    } else {
                                        // This is an IO cluster
                                        int dense_net_id = sparse_dense_map_[sel_neighbor_pair.second[0]];

                                        // Check whether the dense net is already included
                                        if (movable_vertices_[vertex_id].connected_net_.find(dense_net_id) == movable_vertices_[vertex_id].connected_net_.end()) {
                                            // This net has not been added
                                            movable_vertices_[vertex_id].connected_net_.insert(dense_net_id);

                                            // Update the coordinates storing structure inside the corresponding vertex
                                            float second_cluster_x = cluster_map_[neighbor_cluster_id]->GetX();
                                            float second_cluster_y = cluster_map_[neighbor_cluster_id]->GetY();
                                            float second_cluster_z = cluster_map_[neighbor_cluster_id]->GetZ();
                                            std::vector<float> temp_fixed_vertex_info;
                                            temp_fixed_vertex_info.push_back(second_cluster_x);
                                            temp_fixed_vertex_info.push_back(second_cluster_y);
                                            temp_fixed_vertex_info.push_back(second_cluster_z);
                                            temp_fixed_vertex_info.push_back(dense_nets_[dense_net_id].timing_weight);
                                            movable_vertices_[vertex_id].fixed_vertex_connection_map_.push_back(temp_fixed_vertex_info);
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        // For levels other than the first level, we need to check the bbox
                        // No IO clusters in these levels, no need to check
                        // Check whether the selected cluster is inside the bbox
                        float cluster_x = sel_cluster->GetX();
                        float cluster_y = sel_cluster->GetY();

                        if (cluster_x > bbox.lx &&
                            cluster_x < bbox.ux &&
                            cluster_y > bbox.ly &&
                            cluster_y < bbox.uy) {
                            // This cluster is inside the selected bbox area
                            if (cluster_id_vertex_id_map.find(sel_cluster->GetId()) == cluster_id_vertex_id_map.end()) {
                                // Not mapped to any existing vertices
                                int vertex_id = movable_vertices_.size();
                                Vertex temp_vertex(sel_cluster, vertex_id);

                                // Testing
                                // logger_->report("\t[INFO] Creating vertex {} based on cluster {}", vertex_id, sel_cluster->GetId());
                                
                                // Update the mapping information
                                cluster_id_vertex_id_map[sel_cluster->GetId()] = vertex_id;

                                movable_vertices_.push_back(temp_vertex);

                                // Get the neighbor connections
                                std::map<int, std::vector<int>> temp_map; 
                                sel_cluster->GetBundledNetConnections(temp_map);

                                // Add neighbors to the newly created vertex
                                for (auto sel_neighbor_pair : temp_map) {
                                    // Get the cluster id
                                    int neighbor_cluster_id = sel_neighbor_pair.first;
                                    float second_cluster_x = cluster_map_[neighbor_cluster_id]->GetX();
                                    float second_cluster_y = cluster_map_[neighbor_cluster_id]->GetY();
                                    int vertex_id_2 = -1;

                                    // Check whether the neighbor is in the region
                                    if (second_cluster_x > bbox.lx && 
                                        second_cluster_x < bbox.ux && 
                                        second_cluster_y > bbox.ly && 
                                        second_cluster_y < bbox.uy &&
                                        !cluster_map_[neighbor_cluster_id]->GetIOClusterFlag()) {
                                        // The neighboring cluster is inside the movable region
                                        if (cluster_id_vertex_id_map.find(neighbor_cluster_id) == cluster_id_vertex_id_map.end()) {
                                            vertex_id_2 = movable_vertices_.size();

                                            Vertex temp_vertex_2(cluster_map_[neighbor_cluster_id], vertex_id_2);
                                            movable_vertices_.push_back(temp_vertex_2);

                                            // Update the mapping
                                            cluster_id_vertex_id_map[neighbor_cluster_id] = vertex_id_2;
                                        } else {
                                            vertex_id_2 = cluster_id_vertex_id_map[neighbor_cluster_id];
                                        }

                                        // Get the new dense net id
                                        int dense_net_id = sparse_dense_map_[sel_neighbor_pair.second[0]];
                                        movable_vertices_[vertex_id].moveable_vertex_id_weight_map_[vertex_id_2] = dense_nets_[dense_net_id].timing_weight;
                                    } else {
                                        // The neighboring cluster is not inside the bbox region
                                        int dense_net_id = sparse_dense_map_[sel_neighbor_pair.second[0]];

                                        // Check whether the dense net is already included
                                        if (movable_vertices_[vertex_id].connected_net_.find(dense_net_id) == movable_vertices_[vertex_id].connected_net_.end()) {
                                            // This net has not been added
                                            movable_vertices_[vertex_id].connected_net_.insert(dense_net_id);

                                            // Update the coordinates storing structure inside the corresponding vertex
                                            float second_cluster_z = cluster_map_[neighbor_cluster_id]->GetZ();
                                            std::vector<float> temp_fixed_vertex_info;
                                            temp_fixed_vertex_info.push_back(second_cluster_x);
                                            temp_fixed_vertex_info.push_back(second_cluster_y);
                                            temp_fixed_vertex_info.push_back(second_cluster_z);
                                            temp_fixed_vertex_info.push_back(dense_nets_[dense_net_id].timing_weight);
                                            movable_vertices_[vertex_id].fixed_vertex_connection_map_.push_back(temp_fixed_vertex_info);
                                        }
                                        
                                    }
                                }

                            } else {
                                // Mapped to an existing vertex
                                int vertex_id = cluster_id_vertex_id_map[sel_cluster->GetId()];

                                // Get neighbor connections
                                std::map<int, std::vector<int>> temp_map; 
                                sel_cluster->GetBundledNetConnections(temp_map);

                                // Add neighbors to the existing vertex
                                for (auto sel_neighbor_pair : temp_map) {
                                    // Get the cluster id
                                    int neighbor_cluster_id = sel_neighbor_pair.first;
                                    float second_cluster_x = cluster_map_[neighbor_cluster_id]->GetX();
                                    float second_cluster_y = cluster_map_[neighbor_cluster_id]->GetY();
                                    int vertex_id_2 = -1;

                                    // Check whether the neighbor is in the region
                                    if (second_cluster_x > bbox.lx && 
                                        second_cluster_x < bbox.ux && 
                                        second_cluster_y > bbox.ly && 
                                        second_cluster_y < bbox.uy &&
                                        !cluster_map_[neighbor_cluster_id]->GetIOClusterFlag()) {
                                        // The neighboring cluster is inside the movable region
                                        if (cluster_id_vertex_id_map.find(neighbor_cluster_id) == cluster_id_vertex_id_map.end()) {
                                            vertex_id_2 = movable_vertices_.size();

                                            Vertex temp_vertex_2(cluster_map_[neighbor_cluster_id], vertex_id_2);
                                            movable_vertices_.push_back(temp_vertex_2);

                                            // Update the mapping information
                                            cluster_id_vertex_id_map[neighbor_cluster_id] = vertex_id_2;
                                        } else {
                                            vertex_id_2 = cluster_id_vertex_id_map[neighbor_cluster_id];
                                        }

                                        // Get the new dense net id
                                        int dense_net_id = sparse_dense_map_[sel_neighbor_pair.second[0]];
                                        movable_vertices_[vertex_id].moveable_vertex_id_weight_map_[vertex_id_2] = dense_nets_[dense_net_id].timing_weight;
                                    } else {
                                        // The neighboring cluster is not inside the bbox region
                                        int dense_net_id = sparse_dense_map_[sel_neighbor_pair.second[0]];

                                        // Check whether the dense net is already included
                                        if (movable_vertices_[vertex_id].connected_net_.find(dense_net_id) == movable_vertices_[vertex_id].connected_net_.end()) {
                                            // This net has not been added
                                            movable_vertices_[vertex_id].connected_net_.insert(dense_net_id);

                                            // Update the coordinates storing structure inside the corresponding vertex
                                            float second_cluster_z = cluster_map_[neighbor_cluster_id]->GetZ();
                                            std::vector<float> temp_fixed_vertex_info;
                                            temp_fixed_vertex_info.push_back(second_cluster_x);
                                            temp_fixed_vertex_info.push_back(second_cluster_y);
                                            temp_fixed_vertex_info.push_back(second_cluster_z);
                                            temp_fixed_vertex_info.push_back(dense_nets_[dense_net_id].timing_weight);
                                            movable_vertices_[vertex_id].fixed_vertex_connection_map_.push_back(temp_fixed_vertex_info);
                                        }
                                    }
                                }
                            }
                        } else {
                            continue; // No need to traverse further
                        }
                    }
                } 

                for (auto child : sel_cluster->GetChildren()) {
                    // Push all child clusters into the queue
                    bfs_queue.push(child);
                }
            }

            if (current_level == level) {
                break; // All clusters in the corresponding level accessed
            }

            ++current_level;
        }

        // Testing
        logger_->report("[INFO] Finished extracting moveable Clusters with in BBOX: {} {} {} {}, level {}",
                        bbox.lx,
                        bbox.ly,
                        bbox.ux,
                        bbox.uy,
                        level);   
    }

    void HierBlobPlacement::SetNumTiers(const int num_tiers) {
        num_tiers_ = num_tiers;
    }

    void HierBlobPlacement::SetTierResolution(const int tier_resolution) {
        tier_resolutions_ = tier_resolution;
    }

    void HierBlobPlacement::SetWireLengthWeight(const float wire_length_weight) {
        wire_length_weight_ = wire_length_weight;
    }

    void HierBlobPlacement::SetVerticalWeight(const float vertical_weight) {
        vertical_connection_weight_ = vertical_weight;
    }

    void HierBlobPlacement::SetAreaUpperBound(const float area_upper_bound) {
        area_upper_bound_ = area_upper_bound;
    }

    void HierBlobPlacement::SetSwapProb(const float swap_prob) {
        swap_prob_ = swap_prob;
    }

    void HierBlobPlacement::SetInitProb(const float init_prob) {
        init_prob_ = init_prob;
    }

    void HierBlobPlacement::SetMaxSteps(const int max_steps) {
        max_steps_ = max_steps;
    }

    void HierBlobPlacement::SetNumPerturbPerStep(const int num_perturb_per_step) {
        num_perturb_per_step_ = num_perturb_per_step;
    }

    void HierBlobPlacement::SetTimingPenalty(const float timing_penalty) {
        timing_penalty_ = timing_penalty;
    }

    void HierBlobPlacement::SetSeed(const unsigned seed) {
        seed_ = seed;
    }
    

    // ***************************************************************************
    // * Helper Functions
    // ***************************************************************************
    // Print Physical Hierarchy Tree
    void HierBlobPlacement::PrintPhysicalHierarchyTree(Cluster* parent, int level) {
        // Variable Definition
        std::string line;

        for (int i = 0; i < level; i++) {
            line += "+---";
        }

        // Format string definition
        line += fmt::format(
            "{}  ({})  num_std_cell :  {}"
            "  std_cell_area :  {}"
            "  location : ({}, {}, {})",
            parent->GetName(),
            parent->GetId(),
            parent->GetNumStdCell(),
            parent->GetStdCellArea(),
            parent->GetX(),
            parent->GetY(),
            parent->GetZ()
        );

        logger_->report("{}\n", line);

        // Iteratively call the print function
        for (auto& cluster : parent->GetChildren()) {
            PrintPhysicalHierarchyTree(cluster, level + 1);
        }
    }

    // Print all group information
    void HierBlobPlacement::PrintGroupInfo(Group sel_group) {
        std::string line = "";
        line += "Group_ID: " + std::to_string(sel_group.group_id) + "; Group_Size: " + std::to_string(sel_group.group_size) + "\n";
        line += "\tIncluded vertices: ";

        // Print all included vertices
        for (int i = 0; i < sel_group.included_vertices.size(); i++) {
            line += std::to_string(sel_group.included_vertices[i]) + " ";
        }

        line += "\n";
        line += "\tNeighboring vertex: ";

        // Print all neighboring vertices
        for (auto neighbor_id : sel_group.neighbor_vetex_id) {
            line += std::to_string(neighbor_id) + " ";
        }

        line += "\n";
        line += "\tConnected Hyperedges: ";

        // Print all connected hyperedges
        for (auto hyperedge_id : sel_group.connected_hyperedge_id) {
            line += std::to_string(hyperedge_id) + " ";
        }

        logger_->report("{}\n", line);

    }

    // Print All the clusters and their statics
    void HierBlobPlacement::PrintClusters() {
        std::string line = "";
        line += "NUM_CLUSTERS  :  " + std::to_string(cluster_map_.size()) + "\n";
        for (auto& [cluster_id, cluster] : cluster_map_) {
            line += cluster->GetName() + "  ";
            line += std::to_string(cluster->GetId()) + "\n";
        }
        logger_->report(line);
    }

    // Print Connection for all the clusters
    void HierBlobPlacement::PrintConnection() {
        std::string line="";
        line += "NUM_CLUSTERS  :  " + std::to_string(cluster_map_.size()) + "\n";

        for (auto& [cluster_id, cluster] : cluster_map_) {
            const std::map<int, float> connections = cluster->GetConnection();
            if (connections.size() == 0) {
                continue;
            }

            line += "cluster " + cluster->GetName() + " : \n";
            for (auto [target, num_nets] : connections) {
                line += "\t\t" + cluster_map_[target]->GetName() + "  ";
                line += std::to_string(static_cast<int>(num_nets)) + "\n";
            }
        }
        logger_->report(line);
    }

    // Print information for the selected timing path
    void HierBlobPlacement::PrintTimingPathInfo(TimingPath sel_timing_path) {
        std::string line = "";
        line += "[INFO] Timing Path Info\n";
        line += "\tPath: ";
        int counter = 0;

        // Print path information
        for (auto sel_path_element : sel_timing_path.path) {
            counter++;
            line += std::to_string(sel_path_element);

            if (counter != sel_timing_path.path.size()) {
                line += " ";
                line += "->";
                line += " ";
            }
        }

        // Print arc information
        line += "\n";
        counter = 0;
        line += "\tArcs: ";

        for (auto sel_arc_element : sel_timing_path.arcs) {
            counter++;
            line += std::to_string(sel_arc_element);

            if (counter != sel_timing_path.arcs.size()) {
                line += " ";
                line += "->";
                line += " ";
            }
        }

        //  Print Cluster id
        line += "\n";
        counter = 0;
        line += "\tClusters: ";

        for (auto sel_cluster_element : sel_timing_path.hyper_arcs) {
            counter++;
            line += std::to_string(sel_cluster_element);

            if (counter != sel_timing_path.hyper_arcs.size()) {
                line += " ";
                line += "->";
                line += " ";
            }
        }

        logger_->report(line);

        // Print slack information
        std::cout << "\tSlack: " << sel_timing_path.slack << std::endl;
    }

    // Writeout the final mapping between cells and the corresponding location
    // This function will write out the final cell mapping file based on the specified level in the physical
    // hierarchy tree.
    void HierBlobPlacement::GenerateFinalOutputFile(std::string output_file_name, int intended_level) {
        // Variable definition
        std::ofstream out_file;

        // Testing
        logger_->report("[STATUS] Write out the final mapping file for clusters at level {}", intended_level);
        
        // Writeout the file - Single output file
        std::string sel_out_file_name = output_file_name + ".txt";
        out_file.open(sel_out_file_name);
        
        // Checkout clusters in the specified level
        std::queue<Cluster*> bfs_queue;
        bfs_queue.push(root_cluster_);
        int current_level = 0;

        while(!bfs_queue.empty()) {
            int level_size = bfs_queue.size();

            for (int i = 0; i < level_size; i++) {
                Cluster* sel_cluster = bfs_queue.front();
                bfs_queue.pop();

                if (current_level == intended_level) {
                    // Write all cell information in the selected cluster
                    WriteClusterCellInfo(out_file, sel_cluster);
                }

                for (auto child : sel_cluster->GetChildren()) {
                    // Push all child clusters into the queue
                    bfs_queue.push(child);
                }
            }

            if (current_level == intended_level) {
                break; // All clusters in the corresponding level accessed
            }

            ++current_level;
        }

        // Close the file 
        out_file.close();

        // Testing
        logger_->report("[STAUS] End of process, finished Generating the output file;");

    }

    // Write out cell information in a cluster
    void HierBlobPlacement::WriteClusterCellInfo(std::ofstream& out_file, Cluster* sel_cluster) {
        // Check all std cell instances
        for (auto& inst : sel_cluster->GetLeafStdCells()) {
            out_file << inst->getName() << " "
                     << sel_cluster->GetZ() + 1 << " "
                     << sel_cluster->GetX() << " "
                     << sel_cluster->GetY() << " r0" << std::endl; 
        }

        // Set properties for all submodules
        for (auto& module : sel_cluster->GetDbModules()) {
            WriteClusterCellInfo(out_file, sel_cluster->GetX(), sel_cluster->GetY(), sel_cluster->GetZ(), module);
        }
    }

    // Write out cell information in a Dbmodule
    void HierBlobPlacement::WriteClusterCellInfo(std::ofstream& out_file, float X, float Y, float Z, odb::dbModule* module) {
        // No macros are included in the design
        for (odb::dbInst* inst : module->getInsts()) {
            out_file << inst->getName() << " "
                     << Z + 1 << " "
                     << X << " "
                     << Y << " r0" << std::endl;
        }

        // Set properties for all submodules
        for (odb::dbModInst* inst : module->getChildren()) {
            WriteClusterCellInfo(out_file, X, Y, Z, inst->getMaster());
        }
    }

    // This function checks the global status of the code
    void HierBlobPlacement::CalGlobalCost() {
        // Variable Definition
        float hpwl_cost = 0;
        float vertical_connstion_cost = 0;

        for (auto sel_net : dense_nets_) {
            float final_timing_weight = sel_net.timing_weight;

            // Get the locations of the selected clusters
            float src_cluster_x_cord = cluster_map_[sel_net.src_cluster_id]->GetX();
            float src_cluster_y_cord = cluster_map_[sel_net.src_cluster_id]->GetY();
            float src_cluster_z_cord = cluster_map_[sel_net.src_cluster_id]->GetZ();
            float dst_cluster_x_cord = cluster_map_[sel_net.target_cluster_id]->GetX();
            float dst_cluster_y_cord = cluster_map_[sel_net.target_cluster_id]->GetY();
            float dst_cluster_z_cord = cluster_map_[sel_net.target_cluster_id]->GetZ();

            // Update timing weight
            if (dst_cluster_z_cord != src_cluster_z_cord) {
                final_timing_weight += timing_penalty_;
            }

            hpwl_cost += final_timing_weight * (std::abs(src_cluster_x_cord - dst_cluster_x_cord) + std::abs(src_cluster_y_cord - dst_cluster_y_cord));
            vertical_connstion_cost += (std::abs(src_cluster_z_cord - dst_cluster_z_cord));
        }

        // Testing
        float weighted_cost = wire_length_weight_ * hpwl_cost + vertical_connection_weight_ * vertical_connstion_cost;
        logger_->report("[GLOBAL_INFO] Overall weighted {}", weighted_cost);
    }

} // namespace mlsa