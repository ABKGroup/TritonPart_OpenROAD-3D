#include "MultilevelSACore.h"

#include <fstream>
#include <iostream>
#include <queue>
#include <unistd.h>
#include <sstream>
#include <iomanip>
#include <cmath>


#include "db_sta/dbNetwork.hh"
#include "object.h"
#include "odb/db.h"
#include "sta/Liberty.hh"
#include "utl/Logger.h"

namespace mlsa {
    SimulatedAnnealingCore::SimulatedAnnealingCore(
        const std::vector<Vertex> vertices,             // Movable vertices
        float area_upper_bound,                         // Ratio of allowed capacities in different partitions
        int num_tiers,                      
        int tier_resolution,
        float bbox_width,
        float bbox_height,
        std::pair<float, float> ll_cord,       // Lower left corner coordinates
        // Weight for different components
        float wire_length_weight,
        float vertical_connection_weight,
        float timing_penalty,
        // Action probability
        float swap_prob,
        // SA parameters
        float init_prob,
        int max_num_steps,
        int num_perturb_per_step,
        unsigned seed
    ) : ll_cord_(ll_cord) {
        /********************* Testing *********************/
        std::cout << "[Param] Total moves " <<  max_num_steps << std::endl;
        std::cout << "[Param] Capacity Ratio " << area_upper_bound << std::endl;
        std::cout << "[Param] Vertical Weight " << vertical_connection_weight << std::endl;
        std::cout << "[Param] Num perturb per step " << num_perturb_per_step << std::endl;
        std::cout << "[Param] Number of tiers " << num_tiers << std::endl;
        std::cout << "[Param] Design fp width: " << bbox_width << std::endl;
        std::cout << "[Param] Design fp height: " << bbox_height << std::endl;
        std::cout << "[Param] Tier resolution " << tier_resolution << std::endl;
        std::cout << "[Param] Random seed " << seed << std::endl;
        /********************* Testing *********************/

        // Floorplan info
        area_upper_bound_ = area_upper_bound;
        num_tiers_ = num_tiers;
        tier_resolution_ = tier_resolution;
        tier_width_ = bbox_width;
        tier_height_ = bbox_height;
        partitions_per_tier_ = tier_resolution_ * tier_resolution_;
        num_partitions_ = num_tiers_ * partitions_per_tier_;

        // Storing structure
        vertices_ = vertices;

        wire_length_weight_ = wire_length_weight;
        vertical_connection_weight_ = vertical_connection_weight;
        timing_penalty_ = timing_penalty;

        init_prob_ = init_prob;
        swap_prob_ = swap_prob;
        move_prob_ = 1 - swap_prob_;

        // SA Hyperparameter
        // TODO: Change the name of the max_num_steps to max_num_moves
        max_steps_ = max_num_steps / num_perturb_per_step;
        num_perturb_per_step_ = num_perturb_per_step;

        // Random number generation
        std::mt19937 rand_gen(seed);
        generator_ = rand_gen;
        std::uniform_real_distribution<float> distribution(0.0, 1.0);
        distribution_ = distribution;

        // Get the overall vertices area
        for (const auto& sel_vertex : vertices_) {
            vertices_area_ += sel_vertex.GetVertexArea();

            // TODO: Fix this
            // Testing, if a vertex is larger than the allowed 

        }

        // Calculate the area for different partitions from input area sizes
        tier_area_ = tier_width_ * tier_height_;
        total_floorplan_area_ = tier_area_ * num_tiers_;

        // Testing, Dynamically change the area upper bound
        if (vertices_area_ / total_floorplan_area_ > area_upper_bound_) {
            if (vertices_area_ / total_floorplan_area_ + 0.03 <= 0.99) {
                area_upper_bound_ = vertices_area_ / total_floorplan_area_ + 0.03;
            } else {
                area_upper_bound_ = 0.99;
            }

            // Doomed run
            doom_run_ = true;
            std::cout << "[ERROR] Doomed run, reuturn false;" << std::endl;
        }

        // Testing
        std::cout << std::endl;
        std::cout << "[Initialization] Lower left corner in the original canvas: (" << ll_cord_.first << ", " << ll_cord.second << ")" << std::endl;
        std::cout << "[Initialization] Number of vertices : " << vertices_.size() << std::endl;
        std::cout << "[Initialization] Total vertex area: " << vertices_area_ << std::endl;
        std::cout << "[Initialization] Tier area: " << tier_area_ << std::endl;
        std::cout << "[Initialization] Total floorplan area: " << total_floorplan_area_ << std::endl;
        std::cout << "[Initialization] Area Upper Bound: " << area_upper_bound_ << std::endl;
        std::cout << "[Initialization] Tier Width: " << tier_width_ << std::endl;
        std::cout << "[Initialization] Tier Height: " << tier_height_ << std::endl;
        std::cout << "[Initialization] Num Tiers: " << num_tiers_ << std::endl;

        // Partition parameters
        const float partition_capacity = tier_area_ / partitions_per_tier_;
        partition_width_ = partition_height_ = std::sqrt(partition_capacity);

        // Initialize all partitions
        // Aspect ratio for all components are assumed to be 1
        // E.g. for tier_resolution_ = 2, one tier has 4 partitions, and are labeled as
        //      2, 3 <- top right corner of the bottom tier
        //      0, 1
        partitions_.clear();

        // Iterate through different tiers
        for (int i = 0; i < num_tiers_; i++) {
            // Iterate through different partitions in one tier 
            for (int j = 0; j < partitions_per_tier_; j++) {
                // Variable definition
                int partition_id = (i * partitions_per_tier_) + j;

                // Calculate the index of different partitions
                int row_index = (j / tier_resolution_);
                int col_index = j % tier_resolution_;

                // Calculate the coordinates of the chosen partition --> center point
                float partition_x_cord = (partition_width_ / 2.0) + (partition_width_ * col_index);
                float partition_y_cord = (partition_height_ / 2.0) + (partition_height_ * row_index);

                // Get the real coordinates
                partition_x_cord += ll_cord_.first;
                partition_y_cord += ll_cord_.second;

                // Create the new partition
                Partition sel_partition(partition_x_cord, partition_y_cord, (float)i, 
                                        partition_width_, partition_height_, partition_capacity, area_upper_bound_ * partition_capacity,
                                        partition_id);
                
                // Testing
                // std::cout << "Create Partition : " << partition_id << ", center location (" << partition_x_cord << ", " << partition_y_cord << ", " << i << ")" << std::endl;

                partitions_.push_back(sel_partition);
            }
        }
    }

    // ***************************************************************************
    // ***************************************************************************
    float SimulatedAnnealingCore::GetWeightCost() const {
        return weighted_cost_;
    }

    int SimulatedAnnealingCore::GetCrossTierNum() const {
        return num_vertical_connections_;
    }

    float SimulatedAnnealingCore::GetVerticalCost() const {
        return vertical_connection_cost_;
    }

    float SimulatedAnnealingCore::GetTierArea() const {
        return tier_area_;
    }

    float SimulatedAnnealingCore::GetTierWidth() const {
        return tier_width_;
    }

    float SimulatedAnnealingCore::GetTierHeight() const {
        return tier_height_;
    }

    float SimulatedAnnealingCore::GetPartitionWidth() const {
        return partition_width_;
    }

    float SimulatedAnnealingCore::GetPartitionHeight() const {
        return partition_height_;
    }

    float SimulatedAnnealingCore::GetWireLength() const {
        return wire_length_;
    }

    void SimulatedAnnealingCore::GetVertices(std::vector<Vertex>& vertices) const {
        vertices = vertices_;
    }

    void SimulatedAnnealingCore::GetPartitions(std::vector<Partition>& partitions) const {
        partitions = partitions_;
    }


    // Check whether this move is valid
    bool SimulatedAnnealingCore::ValidMove(const int sel_vertex, const int partition_id) const {
        // Check whether enough area is available in the selected partition
        float used_partition_area = partitions_[partition_id].GetUsedCapacity();
        float allowed_partition_area = partitions_[partition_id].GetAllowedCapacity();

        // Newly added vertex shouldn't exceed the area limit of the partition
        if (used_partition_area + vertices_[sel_vertex].GetVertexArea() > allowed_partition_area) {
            return false;
        } else {
            return true;
        }
    }

    // Check whether the swap is valid
    bool SimulatedAnnealingCore::ValidSwap(const int delta_capacity, const int influenced_partition_id) const {
        // Check available area
        float used_partition_area = partitions_[influenced_partition_id].GetUsedCapacity();
        float allowed_partition_area = partitions_[influenced_partition_id].GetAllowedCapacity();

        if (used_partition_area + delta_capacity > allowed_partition_area) {
            return false;
        } else {
            return true;
        }
    }

    // Randomly assign the initial solution based on the random seed
    void SimulatedAnnealingCore::InitialAssignmentScratch() {
        // Testing
        std::cout << "Initialize the canvas from scratch" << std::endl;

        // Assign vertices
        for (auto& sel_vertex : vertices_) {
            bool continue_flag = true;

            while (continue_flag) {
                int partition_index = (int) (std::floor(distribution_(generator_) * num_partitions_));

                // Testing
                // std::cout << "\tSelected vertex id: " << sel_vertex.GetVertexId() << "; Vertex Area: " << sel_vertex.GetVertexArea() << std::endl;
                // std::cout << "\tSelected partition id: " << partition_index << "; Partition allowed area before move: " << partitions_[partition_index].GetAllowedCapacity() << "; Used area: " <<  partitions_[partition_index].GetUsedCapacity() << std::endl;
                
                // Check corner cases
                if (partition_index == num_partitions_) partition_index--;

                // Testing
                if (sel_vertex.GetVertexArea() > partitions_[partition_index].GetAllowedCapacity()) {
                    std::cout << "[ERROR] Vertex id: " << sel_vertex.GetVertexId() << "; Vertex Area: " << sel_vertex.GetVertexArea() << " is larger than the allowed area : " << partitions_[partition_index].GetAllowedCapacity() << std::endl;
                    doom_run_ = true;
                    return;
                }

                if (ValidMove(sel_vertex.GetVertexId(), partition_index)) {
                    continue_flag = false;

                    // Testing
                    // std::cout << "\t\tMoved Vertex id: " << sel_vertex.GetVertexId() << "; Vertex Area: " << sel_vertex.GetVertexArea() << std::endl;
                    // std::cout << "\t\tChosen Partition id: " << partition_index << "; Partition allowed area before move: " << partitions_[partition_index].GetAllowedCapacity() << "; Used area: " <<  partitions_[partition_index].GetUsedCapacity() << std::endl;;

                    // Update vertex info
                    sel_vertex.SetPartitionIndex(partition_index);
                    sel_vertex.SetXCord(partitions_[partition_index].GetXCord());
                    sel_vertex.SetYCord(partitions_[partition_index].GetYCord());
                    sel_vertex.SetZCord(partitions_[partition_index].GetZCord());

                    // Testing
                    // std::cout << "Finished updating vertex " << std::endl;

                    // Update partition info
                    partitions_[partition_index].AddCapacity(sel_vertex.GetVertexArea());
                }

            }
        }
    }

    // Follow the original solution from previous SA runs
    void SimulatedAnnealingCore::InitialAssignmentFollow() {
        // We check the closest bin to the selected vertex and set the partition id accordingly
        // Assign locations of different vertices
        for (auto& sel_vertex : vertices_) {
            // Get the flat partition id
            int column_index = (int) (std::floor(sel_vertex.GetXCord() / partition_width_));
            int row_index = (int) (std::floor(sel_vertex.GetYCord() / partition_height_));
            int partition_index = (sel_vertex.GetZCord() * partitions_per_tier_) + (row_index * tier_resolution_) + column_index;

            // Update vertex info
            sel_vertex.SetPartitionIndex(partition_index);
            sel_vertex.SetXCord(partitions_[partition_index].GetXCord());
            sel_vertex.SetYCord(partitions_[partition_index].GetYCord());
            sel_vertex.SetZCord(partitions_[partition_index].GetZCord());


            // Update partition info
            partitions_[partition_index].AddCapacity(sel_vertex.GetVertexArea());
        }
    }

    void SimulatedAnnealingCore::InitializeSA(bool incremental_assignment) {
        // Testing, Doomed runs
        if (doom_run_) {
            std::cout << "[ERROR] Doomed run , return" << std::endl;
            return;
        }
        
        // Give Initial Assignment
        if (incremental_assignment) {
            // Testing
            std::cout << "[INFO] Initialize the SA worker following previous results" << std::endl;

            InitialAssignmentFollow();
        } else {
            // Testing
            std::cout << "[INFO] Initialize the SA worker from scratch" << std::endl;

            InitialAssignmentScratch();
        }

        // Calculate Initial Cost
        CalWeightedCostScratch();

        // Sync-up
        weighted_cost_ = scratch_weighted_cost_;
        wire_length_ = scratch_wire_length_;
        x_wire_length_ = scratch_x_wire_length_;
        y_wire_length_ = scratch_y_wire_length_;
        vertical_connection_cost_ = scratch_vertical_connection_cost_;

        // Variable Definition
        std::vector<float> wire_length_list;
        std::vector<float> vertical_connection_cost_list;

        // Initial runs
        for (int i = 0; i < num_perturb_per_step_; i++) {
            Perturb();

            // Store the current results
            wire_length_list.push_back(wire_length_);
            vertical_connection_cost_list.push_back(vertical_connection_cost_);
        }

        // Calculate the initial temperature
        std::vector<float> overall_cost_list;
        for (int i = 0; i < wire_length_list.size(); i++) {
            wire_length_ = wire_length_list[i];
            vertical_connection_cost_ = vertical_connection_cost_list[i];
            weighted_cost_ = (wire_length_weight_ * wire_length_) + (vertical_connection_weight_ * vertical_connection_cost_);
        
            overall_cost_list.push_back(weighted_cost_);
        }

        float delta_cost = 0.0;
        for (int i = 1; i < overall_cost_list.size(); i++) {
            delta_cost += std::abs(overall_cost_list[i] - overall_cost_list[i - 1]);
        }

        if (overall_cost_list.size() > 1 && delta_cost > 0.0) {
            init_temperature_
                = (-1.0) * (delta_cost / (overall_cost_list.size() - 1)) / log(init_prob_);
        } else {
            init_temperature_ = 1.0;
        }
    }

    void SimulatedAnnealingCore::RunSA() {
        // Check the validity of the input
        if (vertices_.size() == 0) {
            std::cout << "[ERROR] Vertices size is 0! Abort" << std::endl;;
            return;
        }

        if (doom_run_) {
            std::cout << "[ERROR] Doomed run , return" << std::endl;
            return;
        }

        // Variable Definition
        // clock_t start, end;

        // Perturb at the beginning
        Perturb();

        // Record the current status
        float pre_cost = weighted_cost_;
        pre_weighted_cost_ = weighted_cost_;
        float delta_cost = 0.0;

        int step = 1;

        // Update the decreasing factor
        float temperature = init_temperature_;
        const float min_t = 1e-8;
        // temp_decrease_ratio_ = std::exp(std::log(min_t / init_temperature_) / max_steps_);
        move_ = 0;

        // Testing 
        temp_decrease_ratio_ = 0.99;

        // Start timer
        // start = clock();

        /********************* Testing *********************/
        std::cout << "[Initialization] Init_temperature: " << init_temperature_ << std::endl;
        std::cout << "[Initialization] Temp_decrease_rate: " << temp_decrease_ratio_ << std::endl;
        std::cout << "[Initialization] Tier Width: " << tier_width_ << std::endl;
        std::cout << "[Initialization] Tier Height: " << tier_height_ << std::endl;
        std::cout << " Move: " << std::left << std::setw(8) << move_;
        std::cout << " Cost: " << std::left << std::setw(12) << weighted_cost_;
        std::cout << " Vertical_Connection: " << std::left << std::setw(8) << vertical_connection_cost_;
        std::cout << " Vertical Cost: " << std::left << std::setw(8) << vertical_connection_weight_ * vertical_connection_cost_;
        std::cout << " X_HPWL_COST: " << std::left << std::setw(8) << x_wire_length_;
        std::cout << " Y_HPWL_COST: " << std::left << std::setw(8) << y_wire_length_;
        std::cout << " Time: " << 0 << " s" << std::endl;
        /********************* Testing *********************/

        // Run SA
        while (step <= max_steps_) {
            for (unsigned int i = 0; i < num_perturb_per_step_; i++) {
                Perturb();

                // Update cost
                delta_cost = weighted_cost_ - pre_cost;
                const float random_num = distribution_(generator_);
                const float prob
                    = (delta_cost > 0.0) ? exp((-1) * delta_cost / temperature) : 1;
                
                if (random_num < prob) {
                    pre_cost = weighted_cost_;
                } else {
                    // Testing
                    // std::cout << "[STATUS] Action Rejected, restoring to previous status" << std::endl;

                    Restore();
                }

                move_++;

                // Testing
                if (move_ % 25000 == 0) {
                    std::cout << "Move: " << std::left << std::setw(8) << move_;
                    std::cout << " Cost: " << std::left << std::setw(8) << weighted_cost_;
                    std::cout << " Vertical_Connection: " << std::left << std::setw(8) << vertical_connection_cost_;
                    std::cout << " Vertical Cost: " << std::left << std::setw(8) << vertical_connection_weight_ * vertical_connection_cost_;
                    std::cout << " X_HPWL_COST: " << std::left << std::setw(8) << x_wire_length_;
                    std::cout << " Y_HPWL_COST: " << std::left << std::setw(8) << y_wire_length_;
                    std::cout << std::endl;
                }
            }

            // Update Timing
            temperature *= temp_decrease_ratio_;
            step++;
        }
    }

    // This function will calculate the weighted cost from scratch
    void SimulatedAnnealingCore::CalWeightedCostScratch() {
        // Set all costs to zero
        scratch_weighted_cost_ = 0.0;
        scratch_wire_length_ = 0.0;
        scratch_x_wire_length_ = 0.0;
        scratch_y_wire_length_ = 0.0;
        scratch_vertical_connection_cost_ = 0.0;
        num_vertical_connections_ = 0;

        // Check the weight of the two components
        if (wire_length_weight_ <= 0 && vertical_connection_weight_ <=0) {
            return;
        }

        // Iterate through all movable vertices
        for (const auto& sel_vertex : vertices_) {
            // Two parts need to be considered for each vertex
            // 1. connection with fixed clusters
            // 2. Connection with movable clusters
            std::vector<std::vector<float>> fixed_cluster_connections;
            sel_vertex.GetNets(fixed_cluster_connections);
            const std::map<int, float> movable_cluster_map = sel_vertex.GetMap();

            // Testing
            // std::cout << "\t[Scratch] Calculate Cost for vertex (" << sel_vertex.GetVertexId() << ") from scratch." << std::endl;

            // Testing
            // for (const auto sel_vector : fixed_cluster_connections) {
            //     std::cout << "\tX: " << sel_vector[0] << ";Y: " << sel_vector[1] << ";Z: " << sel_vector[2] << "; Timing weight: " << sel_vector[3] << std::endl;
            // }

            // Testing
            // std::cout << "\t[Scratch] Fixed vertices" << std::endl;

            // Step 1: check connection with fixed clusters
            for (const auto sel_fixed_net : fixed_cluster_connections) {
                // Get the corresponding coordinates
                float fixed_x_cord = sel_fixed_net[0];
                float fixed_y_cord = sel_fixed_net[1];
                float fixed_z_cord = sel_fixed_net[2];
                float timing_weight = sel_fixed_net[3];

                // Testing
                // std::cout << "\t\tSelected Fixed vertex cord : (" << fixed_x_cord << ", " << fixed_y_cord << ", " << fixed_z_cord <<")" << std::endl;
                // std::cout << "\t\tTiming weight : " << timing_weight << std::endl;

                // Delta Z
                float delta_z = std::abs(fixed_z_cord - sel_vertex.GetZCord());

                // Update timing weight
                // TODO: Fix the timing penlty calculation part, for now it's not correct. The penalty shall be added to slack instead of the whole timing factor
                if (delta_z > 0) {
                    // Cross tier connections
                    timing_weight += timing_penalty_;
                }

                // Testing
                // std::cout << "\t\tTiming Weight after updating : " << timing_weight << std::endl;

                // HPWL calculation
                scratch_x_wire_length_ += timing_weight * std::abs(fixed_x_cord - sel_vertex.GetXCord());
                scratch_y_wire_length_ += timing_weight * std::abs(fixed_y_cord - sel_vertex.GetYCord());
                scratch_wire_length_ += timing_weight * (std::abs(fixed_x_cord - sel_vertex.GetXCord()) + std::abs(fixed_y_cord - sel_vertex.GetYCord()));

                // Vertical connections
                scratch_vertical_connection_cost_ += delta_z;
                num_vertical_connections_ += delta_z;

                // Testing
                // std::cout << "\t\t[Scratch] scratch wire_length : " << scratch_wire_length_ << std::endl;
                // std::cout << "\t\t[Scratch] scratch x_wire_length : " << scratch_x_wire_length_ << std::endl;
                // std::cout << "\t\t[Scratch] scratch y_wire_length : " << scratch_y_wire_length_ << std::endl;
                // std::cout << "\t\t[Scratch] scratch z_cost : " << scratch_vertical_connection_cost_ << std::endl;
            }

            // Testing
            // std::cout << "\t[Scratch] Movable vertices" << std::endl;

            // Step 2: check all connections with movable clusters
            for (const auto sel_map_pair : sel_vertex.GetMap()) {
                // Get the vertex_id, and HPWL weight
                const int movable_vertex_id = sel_map_pair.first;
                float net_timing_weight = sel_map_pair.second;

                // Get the corresponding coordinates
                float mov_x_cord = vertices_[movable_vertex_id].GetXCord();
                float mov_y_cord = vertices_[movable_vertex_id].GetYCord();
                float mov_z_cord = vertices_[movable_vertex_id].GetZCord();

                // Testing
                // std::cout << "\t\tSelected Movable vertex cord : (" << mov_x_cord << ", " << mov_y_cord << ", " << mov_z_cord <<")" << std::endl;
                // std::cout << "\t\tSelected Movable vertex index : " << movable_vertex_id << std::endl;
                // std::cout << "\t\tTiming weight : " << net_timing_weight << std::endl;

                // Delta Z
                float delta_z = std::abs(mov_z_cord - sel_vertex.GetZCord());

                // Update timing weight
                if (delta_z > 0) {
                    // Cross tier connections
                    net_timing_weight += timing_penalty_;
                }

                // Testing
                // std::cout << "\t\tTiming Weight after updating : " << net_timing_weight << std::endl;

                // Calculate the HPWL
                scratch_x_wire_length_ += net_timing_weight * std::abs(mov_x_cord - sel_vertex.GetXCord());
                scratch_y_wire_length_ += net_timing_weight * std::abs(mov_y_cord - sel_vertex.GetYCord());
                scratch_wire_length_ += net_timing_weight * (std::abs(mov_x_cord - sel_vertex.GetXCord()) + std::abs(mov_y_cord - sel_vertex.GetYCord()));
                
                // Vertical connections
                scratch_vertical_connection_cost_ += delta_z;
                num_vertical_connections_ += delta_z;

                // Testing
                // std::cout << "\t\t[Scratch] scratch wire_length : " << scratch_wire_length_ << std::endl;
                // std::cout << "\t\t[Scratch] scratch x_wire_length : " << scratch_x_wire_length_ << std::endl;
                // std::cout << "\t\t[Scratch] scratch y_wire_length : " << scratch_y_wire_length_ << std::endl;
                // std::cout << "\t\t[Scratch] scratch z_cost : " << scratch_vertical_connection_cost_ << std::endl;
            }
        }

        // Get the overall cost
        scratch_weighted_cost_ = wire_length_weight_ * scratch_wire_length_ + vertical_connection_weight_ * scratch_vertical_connection_cost_;

        // Testing
        // std::cout << "Status After Initialization" << std::endl;
        // std::cout << "\tWeighted Cost: " << scratch_weighted_cost_ << std::endl;
        // std::cout << "\tNumber of Vertical_Connection: " << num_vertical_connections_ << std::endl;
        // std::cout << "\tX_HPWL_COST : " << scratch_x_wire_length_ << std::endl;
        // std::cout << "\tY_HPWL_COST: " << scratch_y_wire_length_ << std::endl;
    }

    // Calculate the incremental cost for the selected vertex
    void SimulatedAnnealingCore::UpdateCostIncre(int sel_vertex, 
        int target_partition, 
        float& delta_cost, 
        float& delta_z_cost, 
        float& delta_x_hpwl, 
        float& delta_y_hpwl, 
        bool swapping_flag, 
        int swap_vertex_id) {
        // Variable definition
        float pre_z_cost = 0.0;
        float pre_cost = 0.0;
        float pre_x_hpwl = 0.0;
        float pre_y_hpwl = 0.0;
        float new_z_cost = 0.0;
        float new_x_hpwl = 0.0;
        float new_y_hpwl = 0.0;
        float new_cost = 0.0;
        const float vertex_x_cord = vertices_[sel_vertex].GetXCord();
        const float vertex_y_cord = vertices_[sel_vertex].GetYCord();
        const float vertex_z_cord = vertices_[sel_vertex].GetZCord();
        const float partition_x_cord = partitions_[target_partition].GetXCord();
        const float partition_y_cord = partitions_[target_partition].GetYCord();
        const float partition_z_cord = partitions_[target_partition].GetZCord();
        std::vector<std::vector<float>> fixed_cluster_connections;
        vertices_[sel_vertex].GetNets(fixed_cluster_connections);
        const std::map<int, float> movable_cluster_map = vertices_[sel_vertex].GetMap();

        // Testing
        std::cout << "\tCalculate incremental cost update for vertex " << sel_vertex << std::endl;
        std::cout << "\tOriginal vertex cord : (" << vertex_x_cord << ", " << vertex_y_cord << ", " << vertex_z_cord <<")" << std::endl;
        std::cout << "\tTarget partition cord : (" << partition_x_cord << ", " << partition_y_cord << ", " << partition_z_cord <<")" << std::endl;
        std::cout << "\t[Incremental-Updating] Fixed vertices" << std::endl;

        // Check the connection with fixed clusters
        for (const auto sel_fixed_net : fixed_cluster_connections) {
            // Get the corresponding coordinates
            float fixed_x_cord = sel_fixed_net[0];
            float fixed_y_cord = sel_fixed_net[1];
            float fixed_z_cord = sel_fixed_net[2];
            float timing_weight = sel_fixed_net[3];
            float new_timing_weight = sel_fixed_net[3];

            // Testing
            std::cout << "\t\tSelected Fixed vertex cord : (" << fixed_x_cord << ", " << fixed_y_cord << ", " << fixed_z_cord <<")" << std::endl;
            std::cout << "\t\tTiming weight : " << timing_weight << std::endl;

            // Delta Z
            float delta_z = std::abs(fixed_z_cord - vertex_z_cord);

            // Update the original timing weight
            if (delta_z > 0) {
                // Cross tier connection
                timing_weight += timing_penalty_;
            }

            // Testing
            std::cout << "\t\tTiming Weight after updating : " << timing_weight << std::endl;

            // Calculate the original costs
            pre_cost = timing_weight * (std::abs(fixed_x_cord - vertex_x_cord) + std::abs(fixed_y_cord - vertex_y_cord));
            pre_x_hpwl = timing_weight * (std::abs(fixed_x_cord - vertex_x_cord));
            pre_y_hpwl = timing_weight * (std::abs(fixed_y_cord - vertex_y_cord));
            pre_z_cost = delta_z;

            // Testing
            std::cout << "\t\t[Original] Pre-cost = " << pre_cost << std::endl;
            std::cout << "\t\t[Original] Pre_x_cost = " << pre_x_hpwl << std::endl;
            std::cout << "\t\t[Original] Pre_y_cost = " << pre_y_hpwl << std::endl;
            std::cout << "\t\t[Original] Pre_z_cost = " << pre_z_cost << std::endl;

            // Calculate the new timing weight
            if (fixed_z_cord != partition_z_cord) {
                // Cross tier connection
                new_timing_weight += timing_penalty_;
            }

            // Calculate the new cost
            new_cost = new_timing_weight * (std::abs(fixed_x_cord - partition_x_cord) + std::abs(fixed_y_cord - partition_y_cord));
            new_x_hpwl = new_timing_weight * (std::abs(fixed_x_cord - partition_x_cord));
            new_y_hpwl = new_timing_weight * (std::abs(fixed_y_cord - partition_y_cord));
            new_z_cost = std::abs(fixed_z_cord - partition_z_cord);

            // Testing
            std::cout << "\t\t[New] New_cost = " << new_cost << std::endl;
            std::cout << "\t\t[New] New_x_cost = " << new_x_hpwl << std::endl;
            std::cout << "\t\t[New] New_y_cost = " << new_y_hpwl << std::endl;
            std::cout << "\t\t[New] New_z_cost = " << new_z_cost << std::endl;

            // Update the delta cost
            delta_cost += 2 * (new_cost - pre_cost);
            delta_x_hpwl += 2 * (new_x_hpwl - pre_x_hpwl);
            delta_y_hpwl += 2 * (new_y_hpwl - pre_y_hpwl);
            delta_z_cost += 2 * (new_z_cost - pre_z_cost);

            // Testing
            std::cout << "\t\t[Delta] delta_cost : " << delta_cost << std::endl;
            std::cout << "\t\t[Delta] delta_x_hpwl : " << delta_x_hpwl << std::endl;
            std::cout << "\t\t[Delta] delta_y_hpwl : " << delta_y_hpwl << std::endl;
            std::cout << "\t\t[Delta] delta_z_cost : " << delta_z_cost << std::endl;
        }

        // Testing
        std::cout << "\t[Incremental-Updating] Movable vertices" << std::endl;

        // Check the connection with movable vertices
        for (const auto sel_map_pair : movable_cluster_map) {
            // if (swapping_flag && sel_map_pair.first == swap_vertex_id) {
            //     continue;
            // }
            // Get the vertex_id, and HPWL weight
            const int movable_vertex_id = sel_map_pair.first;
            float net_timing_weight = sel_map_pair.second;
            float new_net_timing_weight = sel_map_pair.second;

            // Get the connected vertex information
            float mov_x_cord = vertices_[movable_vertex_id].GetXCord();
            float mov_y_cord = vertices_[movable_vertex_id].GetYCord();
            float mov_z_cord = vertices_[movable_vertex_id].GetZCord();

            // Testing
            std::cout << "\t\tSelected Movable vertex cord : (" << mov_x_cord << ", " << mov_y_cord << ", " << mov_z_cord <<")" << std::endl;
            std::cout << "\t\tSelected Movable vertex index : " << movable_vertex_id << std::endl;
            std::cout << "\t\tTiming weight : " << net_timing_weight << std::endl;

            if (swapping_flag && sel_map_pair.first == swap_vertex_id) {
                continue;
            }

            // Delta Z
            float delta_z = std::abs(mov_z_cord - vertex_z_cord);

            // Update the original timing weight
            if (delta_z > 0) {
                // Cross tier connections
                net_timing_weight += timing_penalty_;
            }

            // Testing
            std::cout << "\t\tTiming Weight after updating : " << net_timing_weight << std::endl;

            // Calculate the original cost
            pre_cost = net_timing_weight * (std::abs(mov_x_cord - vertex_x_cord) + std::abs(mov_y_cord - vertex_y_cord));
            pre_x_hpwl = net_timing_weight * (std::abs(mov_x_cord - vertex_x_cord));
            pre_y_hpwl = net_timing_weight * (std::abs(mov_y_cord - vertex_y_cord));
            pre_z_cost = delta_z;

            // Testing
            std::cout << "\t\t[Original] Pre-cost = " << pre_cost << std::endl;
            std::cout << "\t\t[Original] Pre_x_cost = " << pre_x_hpwl << std::endl;
            std::cout << "\t\t[Original] Pre_y_cost = " << pre_y_hpwl << std::endl;
            std::cout << "\t\t[Original] Pre_z_cost = " << pre_z_cost << std::endl;

            // Calculate the new timing weight
            if (mov_z_cord != partition_z_cord) {
                // Cross tier connection
                new_net_timing_weight += timing_penalty_;
            }

            // Calculate the new cost
            new_cost = new_net_timing_weight * (std::abs(mov_x_cord - partition_x_cord) + std::abs(mov_y_cord - partition_y_cord));
            new_x_hpwl = new_net_timing_weight * (std::abs(mov_x_cord - partition_x_cord));
            new_y_hpwl = new_net_timing_weight * (std::abs(mov_y_cord - partition_y_cord));
            new_z_cost = std::abs(mov_z_cord - partition_z_cord);

            // Testing
            std::cout << "\t\t[New] New_cost = " << new_cost << std::endl;
            std::cout << "\t\t[New] New_x_cost = " << new_x_hpwl << std::endl;
            std::cout << "\t\t[New] New_y_cost = " << new_y_hpwl << std::endl;
            std::cout << "\t\t[New] New_z_cost = " << new_z_cost << std::endl;


            // Update the delta cost
            delta_cost += 2 * (new_cost - pre_cost);
            delta_x_hpwl += 2 * (new_x_hpwl - pre_x_hpwl);
            delta_y_hpwl += 2 * (new_y_hpwl - pre_y_hpwl);
            delta_z_cost += 2 * (new_z_cost - pre_z_cost);

            // Testing
            std::cout << "\t\t[Delta] delta_cost : " << delta_cost << std::endl;
            std::cout << "\t\t[Delta] delta_x_hpwl : " << delta_x_hpwl << std::endl;
            std::cout << "\t\t[Delta] delta_y_hpwl : " << delta_y_hpwl << std::endl;
            std::cout << "\t\t[Delta] delta_z_cost : " << delta_z_cost << std::endl;
        }
    }

    float SimulatedAnnealingCore::CalAvg(std::vector<float>& value_list) {
        const auto size = value_list.size();
        if (size == 0) {
            return 0;
        }

        return std::accumulate(value_list.begin(), value_list.end(), 0) / size;
    }

    void SimulatedAnnealingCore::Perturb() {
        // Backup all metrics
        pre_weighted_cost_ = weighted_cost_;
        pre_wire_length_ = wire_length_;
        pre_x_wire_length_ = x_wire_length_;
        pre_y_wire_length_ = y_wire_length_;
        pre_vertical_connection_cost_ = vertical_connection_cost_;
        pre_num_vertical_connections_ = num_vertical_connections_;

        // Generate random number to select the action
        const float op = distribution_(generator_);
        const float action_prob_1 = swap_prob_;
        const float action_prob_2 = action_prob_1 + move_prob_;

        if (op <= action_prob_1) {
            // Conduct the swap operation
            action_id_ = 1;

            // Swap();
            SwapScratch();
        } else {
            // Condut the move operation
            action_id_ = 2;

            // Move();
            MoveScratch();
        }

        // Update the total cost
        weighted_cost_ = wire_length_weight_ * wire_length_ + vertical_connection_weight_ * vertical_connection_cost_;
    }

    void SimulatedAnnealingCore::Swap() {
        // Variable Definition
        bool continue_flag = true;
        float delta_x_hpwl = 0.0;
        float delta_y_hpwl = 0.0;
        float delta_cost = 0.0;
        float delta_z_cost = 0.0;

        // Testing
        std::cout << "\t[INFO] Previous wire_length result : " << pre_wire_length_ << std::endl;
        std::cout << "\t[INFO] Previous vertical_connection_cost result : " << pre_vertical_connection_cost_ << std::endl;
        std::cout << "\t[INFO] Previous x_wire_length result : " << pre_x_wire_length_ << std::endl;
        std::cout << "\t[INFO] Previous y_wire_length result : " << pre_y_wire_length_ << std::endl;

        // Find two valid swap indices
        while (continue_flag) {
            // Generate the swap indices
            int sel_vertex_index_1
                = (int) (std::floor(distribution_(generator_) * vertices_.size()));
            int sel_vertex_index_2
                = (int) (std::floor(distribution_(generator_) * vertices_.size()));

            // Corner case checking
            if (sel_vertex_index_1 == vertices_.size()) sel_vertex_index_1 -= 1;
            if (sel_vertex_index_2 == vertices_.size()) sel_vertex_index_2 -= 1;

            while (sel_vertex_index_1 == sel_vertex_index_2) {
                sel_vertex_index_2 = (int) (std::floor(std::floor(distribution_(generator_) * vertices_.size())));

                // Corner case
                if (sel_vertex_index_2 == vertices_.size()) sel_vertex_index_2 -= 1;
            }

            const int partition_index_1
                = vertices_[sel_vertex_index_1].GetPartitionId();
            const int partition_index_2
                = vertices_[sel_vertex_index_2].GetPartitionId();
            
            // Check whether the swap is needed
            if (partition_index_1 == partition_index_2) {
                continue;
            }

            // Check whether the swap can be made
            const float area_vertex_1 = vertices_[sel_vertex_index_1].GetVertexArea();
            const float area_vertex_2 = vertices_[sel_vertex_index_2].GetVertexArea();

            const bool computation_flag_1 = 
                ValidSwap(area_vertex_1 - area_vertex_2, partition_index_2);
            const bool computation_flag_2 = 
                ValidSwap(area_vertex_2 - area_vertex_1, partition_index_1);

            // Conduct the swap
            if (computation_flag_1 & computation_flag_2) {
                continue_flag = false;

                // Testing
                std::cout << "[DEBUG] Conduct Swapping action" << std::endl;
                std::cout << "\tSelected vertex id 1 : " << sel_vertex_index_1 << "; Selected Vertex area : " << area_vertex_1 << "; Original Partition index : " << partition_index_1 << std::endl;
                std::cout << "\tSelected vertex id 2 : " << sel_vertex_index_2 << "; Selected Vertex area : " << area_vertex_2 << "; Original Partition index : " << partition_index_2 << std::endl;

                // Store the selected partition information
                sel_partition_index_1_ = partition_index_1;
                sel_partition_index_2_ = partition_index_2;

                // Testing
                std::cout << "\t[BEFORE] Partition id : " << partition_index_1 << "; Used Area before swapping : " << partitions_[partition_index_1].GetUsedCapacity() << "; Allowed Area : " << partitions_[partition_index_1].GetAllowedCapacity() << std::endl;
                std::cout << "\t[BEFORE] Partition id : " << partition_index_2 << "; Used Area before swapping : " << partitions_[partition_index_2].GetUsedCapacity() << "; Allowed Area : " << partitions_[partition_index_2].GetAllowedCapacity() << std::endl;

                partitions_[partition_index_1].MinusCapacity(area_vertex_1);
                partitions_[partition_index_1].AddCapacity(area_vertex_2);

                partitions_[partition_index_2].MinusCapacity(area_vertex_2);
                partitions_[partition_index_2].AddCapacity(area_vertex_1);

                // Testing
                std::cout << "\t[AFTER] Partition id : " << partition_index_1 << "; Used Area after swapping : " << partitions_[partition_index_1].GetUsedCapacity() << "; Allowed Area : " << partitions_[partition_index_1].GetAllowedCapacity() << std::endl;
                std::cout << "\t[AFTER] Partition id : " << partition_index_2 << "; Used Area after swapping : " << partitions_[partition_index_2].GetUsedCapacity() << "; Allowed Area : " << partitions_[partition_index_2].GetAllowedCapacity() << std::endl;

                // Update the cost for vertex 1
                UpdateCostIncre(sel_vertex_index_1, 
                                partition_index_2, 
                                delta_cost, 
                                delta_z_cost, 
                                delta_x_hpwl, 
                                delta_y_hpwl,
                                true,
                                sel_vertex_index_2);

                // Update the cost for vertex 2
                UpdateCostIncre(sel_vertex_index_2, 
                                partition_index_1, 
                                delta_cost, 
                                delta_z_cost, 
                                delta_x_hpwl, 
                                delta_y_hpwl,
                                true,
                                sel_vertex_index_1);

                // Update the cost
                wire_length_ += delta_cost;
                vertical_connection_cost_ += delta_z_cost;
                x_wire_length_ += delta_x_hpwl;
                y_wire_length_ += delta_y_hpwl;

                // Testing
                std::cout << "\t[INFO] Final Incremental wire_length result : " << wire_length_ << std::endl;
                std::cout << "\t[INFO] Final Incremental vertical_connection_cost result : " << vertical_connection_cost_ << std::endl;
                std::cout << "\t[INFO] Final Incremental x_wire_length result : " << x_wire_length_ << std::endl;
                std::cout << "\t[INFO] Final Incremental y_wire_length result : " << y_wire_length_ << std::endl;

                // Update the backup info
                ver_1_original_x_ = vertices_[sel_vertex_index_1].GetXCord();
                ver_1_original_y_ = vertices_[sel_vertex_index_1].GetYCord();
                ver_1_original_z_ = vertices_[sel_vertex_index_1].GetZCord();
                ver_2_original_x_ = vertices_[sel_vertex_index_2].GetXCord();
                ver_2_original_y_ = vertices_[sel_vertex_index_2].GetYCord();
                ver_2_original_z_ = vertices_[sel_vertex_index_2].GetZCord();

                // Change the parameters of the corresponding vertices
                vertices_[sel_vertex_index_1].SetXCord(partitions_[partition_index_2].GetXCord());
                vertices_[sel_vertex_index_1].SetYCord(partitions_[partition_index_2].GetYCord());
                vertices_[sel_vertex_index_1].SetZCord(partitions_[partition_index_2].GetZCord());
                vertices_[sel_vertex_index_2].SetXCord(partitions_[partition_index_1].GetXCord());
                vertices_[sel_vertex_index_2].SetYCord(partitions_[partition_index_1].GetYCord());
                vertices_[sel_vertex_index_2].SetZCord(partitions_[partition_index_1].GetZCord());

                vertices_[sel_vertex_index_1].SetPartitionIndex(partition_index_2);
                vertices_[sel_vertex_index_2].SetPartitionIndex(partition_index_1);
                action_vertex_id_1_ = sel_vertex_index_1;
                action_vertex_id_2_ = sel_vertex_index_2;

                // ***************** Testing *****************
                // Calculate cost from scratch
                CalWeightedCostScratch();

                std::cout << "\t[INFO] Final from scratch wire_length result : " << scratch_wire_length_ << std::endl;
                std::cout << "\t[INFO] Final from scratch vertical_connection_cost result : " << scratch_vertical_connection_cost_ << std::endl;
                std::cout << "\t[INFO] Final from scratch x_wire_length result : " << scratch_x_wire_length_ << std::endl;
                std::cout << "\t[INFO] Final from scratch y_wire_length result : " << scratch_y_wire_length_ << std::endl;

                if (wire_length_ != scratch_wire_length_) {
                    std::cout << "[ERROR] From scratch is different from incremental results" << std::endl;
                }
                // ***************** Testing *****************
            }
        } 
    }

    void SimulatedAnnealingCore::SwapScratch() {
        // Variable Definition
        bool continue_flag = true;

        // Find two valid swap indices
        while (continue_flag) {
            // Generate the swap indices
            int sel_vertex_index_1
                = (int) (std::floor(distribution_(generator_) * vertices_.size()));
            int sel_vertex_index_2
                = (int) (std::floor(distribution_(generator_) * vertices_.size()));

            // Corner case checking
            if (sel_vertex_index_1 == vertices_.size()) sel_vertex_index_1 -= 1;
            if (sel_vertex_index_2 == vertices_.size()) sel_vertex_index_2 -= 1;

            while (sel_vertex_index_1 == sel_vertex_index_2) {
                sel_vertex_index_2 = (int) (std::floor(std::floor(distribution_(generator_) * vertices_.size())));

                // Corner case
                if (sel_vertex_index_2 == vertices_.size()) sel_vertex_index_2 -= 1;
            }

            const int partition_index_1
                = vertices_[sel_vertex_index_1].GetPartitionId();
            const int partition_index_2
                = vertices_[sel_vertex_index_2].GetPartitionId();
            
            // Check whether the swap is needed
            if (partition_index_1 == partition_index_2) {
                continue;
            }

            // Check whether the swap can be made
            const float area_vertex_1 = vertices_[sel_vertex_index_1].GetVertexArea();
            const float area_vertex_2 = vertices_[sel_vertex_index_2].GetVertexArea();

            const bool computation_flag_1 = 
                ValidSwap(area_vertex_1 - area_vertex_2, partition_index_2);
            const bool computation_flag_2 = 
                ValidSwap(area_vertex_2 - area_vertex_1, partition_index_1);

            // Conduct the swap
            if (computation_flag_1 & computation_flag_2) {
                continue_flag = false;

                // Store the selected partition information
                sel_partition_index_1_ = partition_index_1;
                sel_partition_index_2_ = partition_index_2;

                partitions_[partition_index_1].MinusCapacity(area_vertex_1);
                partitions_[partition_index_1].AddCapacity(area_vertex_2);

                partitions_[partition_index_2].MinusCapacity(area_vertex_2);
                partitions_[partition_index_2].AddCapacity(area_vertex_1);

                // Update the backup info
                ver_1_original_x_ = vertices_[sel_vertex_index_1].GetXCord();
                ver_1_original_y_ = vertices_[sel_vertex_index_1].GetYCord();
                ver_1_original_z_ = vertices_[sel_vertex_index_1].GetZCord();
                ver_2_original_x_ = vertices_[sel_vertex_index_2].GetXCord();
                ver_2_original_y_ = vertices_[sel_vertex_index_2].GetYCord();
                ver_2_original_z_ = vertices_[sel_vertex_index_2].GetZCord();

                // Change the parameters of the corresponding vertices
                vertices_[sel_vertex_index_1].SetXCord(partitions_[partition_index_2].GetXCord());
                vertices_[sel_vertex_index_1].SetYCord(partitions_[partition_index_2].GetYCord());
                vertices_[sel_vertex_index_1].SetZCord(partitions_[partition_index_2].GetZCord());
                vertices_[sel_vertex_index_2].SetXCord(partitions_[partition_index_1].GetXCord());
                vertices_[sel_vertex_index_2].SetYCord(partitions_[partition_index_1].GetYCord());
                vertices_[sel_vertex_index_2].SetZCord(partitions_[partition_index_1].GetZCord());

                vertices_[sel_vertex_index_1].SetPartitionIndex(partition_index_2);
                vertices_[sel_vertex_index_2].SetPartitionIndex(partition_index_1);
                action_vertex_id_1_ = sel_vertex_index_1;
                action_vertex_id_2_ = sel_vertex_index_2;

                // Calculate the cost from scratch
                CalWeightedCostScratch();

                // Sync-up
                weighted_cost_ = scratch_weighted_cost_;
                wire_length_ = scratch_wire_length_;
                x_wire_length_ = scratch_x_wire_length_;
                y_wire_length_ = scratch_y_wire_length_;
                vertical_connection_cost_ = scratch_vertical_connection_cost_;
            }
        }
    }

    void SimulatedAnnealingCore::MoveScratch() {
        // Variable definition
        bool continue_flag = true;

        // Find a valid move action
        while (continue_flag) {
            int sel_vertex_index
                = (int) (std::floor(distribution_(generator_) * vertices_.size()));
            
            int sel_partition_index
                = (int) (std::floor(distribution_(generator_) * partitions_.size()));

            // Corner case checking
            if (sel_vertex_index == vertices_.size()) sel_vertex_index -= 1;
            if (sel_partition_index == partitions_.size()) sel_partition_index -= 1;

            const int original_partition_index
                = vertices_[sel_vertex_index].GetPartitionId();

            if (sel_partition_index == original_partition_index) continue;

            if (ValidMove(sel_vertex_index, sel_partition_index)) {
                continue_flag = false;

                // Change the partition info
                sel_partition_index_1_ = original_partition_index;
                sel_partition_index_2_ = sel_partition_index;

                partitions_[original_partition_index].MinusCapacity(vertices_[sel_vertex_index].GetVertexArea());
                partitions_[sel_partition_index].AddCapacity(vertices_[sel_vertex_index].GetVertexArea());

                // Update backup info
                ver_1_original_x_ = vertices_[sel_vertex_index].GetXCord();
                ver_1_original_y_ = vertices_[sel_vertex_index].GetYCord();
                ver_1_original_z_ = vertices_[sel_vertex_index].GetZCord();

                // Change the parameters of the corresponding vertex
                vertices_[sel_vertex_index].SetXCord(partitions_[sel_partition_index].GetXCord());
                vertices_[sel_vertex_index].SetYCord(partitions_[sel_partition_index].GetYCord());
                vertices_[sel_vertex_index].SetZCord(partitions_[sel_partition_index].GetZCord());

                // Move the vertex
                vertices_[sel_vertex_index].SetPartitionIndex(sel_partition_index);
                action_vertex_id_1_ = sel_vertex_index;

                // Calculate the cost from scratch
                CalWeightedCostScratch();

                // Sync-up
                weighted_cost_ = scratch_weighted_cost_;
                wire_length_ = scratch_wire_length_;
                x_wire_length_ = scratch_x_wire_length_;
                y_wire_length_ = scratch_y_wire_length_;
                vertical_connection_cost_ = scratch_vertical_connection_cost_;
            }
        }
    }

    void SimulatedAnnealingCore::Move() {
        // Variable definition
        bool continue_flag = true;
        float delta_cost = 0.0;
        float delta_z_cost = 0.0;
        float delta_x_hpwl = 0.0;
        float delta_y_hpwl = 0.0;

        // Testing
        std::cout << "\t[INFO] Previous wire_length result : " << pre_wire_length_ << std::endl;
        std::cout << "\t[INFO] Previous vertical_connection_cost result : " << pre_vertical_connection_cost_ << std::endl;
        std::cout << "\t[INFO] Previous x_wire_length result : " << pre_x_wire_length_ << std::endl;
        std::cout << "\t[INFO] Previous y_wire_length result : " << pre_y_wire_length_ << std::endl;


        // Find a valid move action
        while (continue_flag) {
            int sel_vertex_index
                = (int) (std::floor(distribution_(generator_) * vertices_.size()));
            
            int sel_partition_index
                = (int) (std::floor(distribution_(generator_) * partitions_.size()));

            // Corner case checking
            if (sel_vertex_index == vertices_.size()) sel_vertex_index -= 1;
            if (sel_partition_index == partitions_.size()) sel_partition_index -= 1;

            const int original_partition_index
                = vertices_[sel_vertex_index].GetPartitionId();

            if (sel_partition_index == original_partition_index) continue;

            if (ValidMove(sel_vertex_index, sel_partition_index)) {
                continue_flag = false;

                // Testing
                std::cout << "[DEBUG] Conduct Move action" << std::endl;
                std::cout << "\tSelected vertex id : " << sel_vertex_index << "; Selected Vertex area : " << vertices_[sel_vertex_index].GetVertexArea() << "; Original Partition index : " << sel_partition_index << std::endl;

                // Change the partition info
                sel_partition_index_1_ = original_partition_index;
                sel_partition_index_2_ = sel_partition_index;

                // Testing
                std::cout << "\t[BEFORE] Partition id : " << sel_partition_index_1_ << "; Used Area before swapping : " << partitions_[sel_partition_index_1_].GetUsedCapacity() << "; Allowed Area : " << partitions_[sel_partition_index_1_].GetAllowedCapacity() << std::endl;
                std::cout << "\t[BEFORE] Partition id : " << sel_partition_index_2_ << "; Used Area before swapping : " << partitions_[sel_partition_index_2_].GetUsedCapacity() << "; Allowed Area : " << partitions_[sel_partition_index_2_].GetAllowedCapacity() << std::endl;

                partitions_[original_partition_index].MinusCapacity(vertices_[sel_vertex_index].GetVertexArea());
                partitions_[sel_partition_index].AddCapacity(vertices_[sel_vertex_index].GetVertexArea());

                // Testing
                std::cout << "\t[BEFORE] Partition id : " << sel_partition_index_1_ << "; Used Area before swapping : " << partitions_[sel_partition_index_1_].GetUsedCapacity() << "; Allowed Area : " << partitions_[sel_partition_index_1_].GetAllowedCapacity() << std::endl;
                std::cout << "\t[BEFORE] Partition id : " << sel_partition_index_2_ << "; Used Area before swapping : " << partitions_[sel_partition_index_2_].GetUsedCapacity() << "; Allowed Area : " << partitions_[sel_partition_index_2_].GetAllowedCapacity() << std::endl;

                // Incremental cost update
                UpdateCostIncre(sel_vertex_index, 
                                sel_partition_index, 
                                delta_cost, 
                                delta_z_cost, 
                                delta_x_hpwl, 
                                delta_y_hpwl,
                                false,
                                0);

                // Update backup info
                ver_1_original_x_ = vertices_[sel_vertex_index].GetXCord();
                ver_1_original_y_ = vertices_[sel_vertex_index].GetYCord();
                ver_1_original_z_ = vertices_[sel_vertex_index].GetZCord();

                // Update cost
                wire_length_ += delta_cost;
                vertical_connection_cost_ += delta_z_cost;
                x_wire_length_ += delta_x_hpwl;
                y_wire_length_ += delta_y_hpwl;

                // Testing
                std::cout << "\t[INFO] Final Incremental wire_length result : " << wire_length_ << std::endl;
                std::cout << "\t[INFO] Final Incremental vertical_connection_cost result : " << vertical_connection_cost_ << std::endl;
                std::cout << "\t[INFO] Final Incremental x_wire_length result : " << x_wire_length_ << std::endl;
                std::cout << "\t[INFO] Final Incremental y_wire_length result : " << y_wire_length_ << std::endl;

                // Change the parameters of the corresponding vertex
                vertices_[sel_vertex_index].SetXCord(partitions_[sel_partition_index].GetXCord());
                vertices_[sel_vertex_index].SetYCord(partitions_[sel_partition_index].GetYCord());
                vertices_[sel_vertex_index].SetZCord(partitions_[sel_partition_index].GetZCord());

                // Move the vertex
                vertices_[sel_vertex_index].SetPartitionIndex(sel_partition_index);
                action_vertex_id_1_ = sel_vertex_index;

                // ***************** Testing *****************
                // Calculate cost from scratch
                CalWeightedCostScratch();

                std::cout << "\t[INFO] Final from scratch wire_length result : " << scratch_wire_length_ << std::endl;
                std::cout << "\t[INFO] Final from scratch vertical_connection_cost result : " << scratch_vertical_connection_cost_ << std::endl;
                std::cout << "\t[INFO] Final from scratch x_wire_length result : " << scratch_x_wire_length_ << std::endl;
                std::cout << "\t[INFO] Final from scratch y_wire_length result : " << scratch_y_wire_length_ << std::endl;

                if (wire_length_ != scratch_wire_length_) {
                    std::cout << "[ERROR] From scratch is different from incremental results" << std::endl;
                }
                // ***************** Testing *****************
            }
        }
    }

    void SimulatedAnnealingCore::Restore() {
        // Check the conducted action
        if (action_id_ == 1) {
            // Restore Swap operation
            vertices_[action_vertex_id_1_].SetXCord(ver_1_original_x_);
            vertices_[action_vertex_id_1_].SetYCord(ver_1_original_y_);
            vertices_[action_vertex_id_1_].SetZCord(ver_1_original_z_);
            vertices_[action_vertex_id_2_].SetXCord(ver_2_original_x_);
            vertices_[action_vertex_id_2_].SetYCord(ver_2_original_y_);
            vertices_[action_vertex_id_2_].SetZCord(ver_2_original_z_);
            vertices_[action_vertex_id_1_].SetPartitionIndex(sel_partition_index_1_);
            vertices_[action_vertex_id_2_].SetPartitionIndex(sel_partition_index_2_);
            // Restore partition information
            partitions_[sel_partition_index_1_].AddCapacity(vertices_[action_vertex_id_1_].GetVertexArea());
            partitions_[sel_partition_index_1_].MinusCapacity(vertices_[action_vertex_id_2_].GetVertexArea());
            partitions_[sel_partition_index_2_].AddCapacity(vertices_[action_vertex_id_2_].GetVertexArea());
            partitions_[sel_partition_index_2_].MinusCapacity(vertices_[action_vertex_id_1_].GetVertexArea());
        } else {
            // Restore Move operation
            vertices_[action_vertex_id_1_].SetXCord(ver_1_original_x_);
            vertices_[action_vertex_id_1_].SetYCord(ver_1_original_y_);
            vertices_[action_vertex_id_1_].SetZCord(ver_1_original_z_);
            vertices_[action_vertex_id_1_].SetPartitionIndex(sel_partition_index_1_);
            // Restore Partition information
            partitions_[sel_partition_index_1_].AddCapacity(vertices_[action_vertex_id_1_].GetVertexArea());
            partitions_[sel_partition_index_2_].MinusCapacity(vertices_[action_vertex_id_1_].GetVertexArea());
        }

        // Restore all costs
        weighted_cost_ = pre_weighted_cost_;
        wire_length_ = pre_wire_length_;
        x_wire_length_ = pre_x_wire_length_;
        y_wire_length_ = pre_y_wire_length_;
        vertical_connection_cost_ = pre_vertical_connection_cost_;
        num_vertical_connections_ = pre_num_vertical_connections_;

    }
}