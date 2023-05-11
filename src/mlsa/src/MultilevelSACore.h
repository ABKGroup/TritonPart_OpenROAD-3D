#pragma once

#include <random>
#include <vector>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <iostream>

#include "odb/dbTypes.h"
#include "object.h"

namespace odb {
    class dbBlock;
    class dbBTerm;
    class dbDatabase;
    class dbITerm;
    class dbInst;
    class dbModule;
} // namespace odb

namespace sta {
    class dbNetwork;
    class dbSta;
} // namespace sta

namespace utl {
    class Logger;
} // namespace utl

namespace mlsa {

    /**
     * @class  SimulatedAnnealingCore
     * @author Jiantao Liu
     * @date   04/12/2023
     * @brief  Implement the simulated annealing framework
    */
    class SimulatedAnnealingCore {
        public:
            SimulatedAnnealingCore(
                const std::vector<Vertex> vertices,             // Movable vertices
                float area_upper_bound,                         // Ratio of allowed capacities in different partitions
                int num_tiers,                      
                int tier_resolution,
                float bbox_width,
                float bbox_height,
                std::pair<float, float> ll_cord,                // Lower left corner coordinates
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
            );

            // Get info
            float GetWeightCost() const;
            int GetCrossTierNum() const;
            float GetVerticalCost() const;
            float GetTierArea() const;
            float GetTierWidth() const;
            float GetTierHeight() const;
            float GetPartitionWidth() const;
            float GetPartitionHeight() const;
            float GetWireLength() const;
            void GetVertices(std::vector<Vertex>& vertices) const;
            void GetPartitions(std::vector<Partition>& partitions) const;
            bool ValidMove(const int sel_vertex, const int partition_id) const;
            bool ValidSwap(const int delta_capacity, const int influenced_partition_id) const;

            // Initialize placement
            void InitialAssignmentScratch();
            void InitialAssignmentFollow();

            // Initialize SA
            void InitializeSA(bool incremental_assignment);

            // Run the SA worker
            void RunSA();

            // Update From Scratch Cost
            void CalWeightedCostScratch();

            // Incremental cost Update
            void UpdateCostIncre(int sel_vertex, 
                                 int target_partition, 
                                 float& delta_cost, 
                                 float& delta_z_cost, 
                                 float& delta_x_hpwl, 
                                 float& delta_y_hpwl, 
                                 bool swapping_flag, 
                                 int swap_vertex_id);

            // Helper function
            float CalAvg(std::vector<float>& value_list);

            // Actions
            void Perturb();
            void Swap();
            void Move();
            void Restore();
            void MoveScratch();
            void SwapScratch();

        private:
            // Floorplan info
            float total_floorplan_area_ = 0.0;
            float vertices_area_ = 0.0;
            float tier_area_ = 0.0;
            float tier_width_ = 0.0;
            float tier_height_ = 0.0;
            float partition_width_ = 0.0;
            float partition_height_ = 0.0;
            float area_upper_bound_ = 0.0;
            std::pair<float, float> ll_cord_; // Global coordinates for the lower left corner

            // Storing structure
            const std::vector<BundledNet> nets_;               // Vector from the outer physical hierarchy tree
            std::vector<Vertex> vertices_;
            std::vector<Partition> partitions_;
            
            // Parameters
            int num_tiers_ = 0;
            int num_partitions_ = 0;
            int tier_resolution_ = 0;
            int partitions_per_tier_ = 0;

            // Weights
            float wire_length_weight_ = 0.0;
            float vertical_connection_weight_ = 0.0;
            float timing_penalty_ = 0.0;

            // SA parameters
            float init_temperature_ = 0.0;
            float temp_decrease_ratio_ = 0.0;
            int max_steps_ = 0;
            int num_perturb_per_step_ = 0;

            // Seed for reproduciability
            std::mt19937 generator_;
            std::uniform_real_distribution<float> distribution_;

            // Last state 
            float pre_weighted_cost_ = 0.0;
            float pre_wire_length_ = 0.0;
            float pre_x_wire_length_ = 0.0;
            float pre_y_wire_length_ = 0.0;
            float pre_vertical_connection_cost_ = 0.0;
            int pre_num_vertical_connections_ = 0;
            int sel_partition_index_1_ = 0;
            int sel_partition_index_2_ = 0;
            float ver_1_original_x_ = 0.0;
            float ver_1_original_y_ = 0.0;
            float ver_1_original_z_ = 0.0;
            float ver_2_original_x_ = 0.0;
            float ver_2_original_y_ = 0.0;
            float ver_2_original_z_ = 0.0;

            // Book keeping
            int action_vertex_id_1_ = -1;
            int action_vertex_id_2_ = -1;
            int action_id_ = -1;

            // Metrics
            float scratch_weighted_cost_ = 0.0;
            float scratch_wire_length_ = 0.0;
            float scratch_x_wire_length_ = 0.0;
            float scratch_y_wire_length_ = 0.0;
            float scratch_vertical_connection_cost_ = 0.0;
            float weighted_cost_ = 0.0;
            float x_wire_length_ = 0.0;
            float y_wire_length_ = 0.0;
            float wire_length_ = 0.0;
            float vertical_connection_cost_ = 0.0;
            int num_vertical_connections_ = 0;

            // Constraints
            float balanced_area_upper_ = 0.0;
            float balanced_area_lower_ = 0.0;

            // Prob of different actions
            float init_prob_ = 0.0;
            float swap_prob_ = 0.0;
            float move_prob_ = 0.0;

            // Debug parameters
            int move_ = 0;

            // Debug - Fault flag
            bool doom_run_ = false;


    }; // End of class SimulatedAnnealingCore


} // End of namespace mlsa