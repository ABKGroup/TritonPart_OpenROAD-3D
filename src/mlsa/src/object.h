#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <queue>

#include "odb/dbTypes.h"

namespace odb {
    class dbInst;
    class dbModule;
    class dbNet;
}   // namespace odb

namespace utl {
    class Logger;
} // namespace utl

namespace mlsa {
    // ***************************************************************************
    // * Basic Functions
    // ***************************************************************************
    // Conversion between dbu and microns
    float DbuToMicron(int metric, float dbu);
    int MicronToDbu(float metric, float dbu);

    // Define the type for clusters, for now we only support StdCellClusters, but it
    // can be extended to include macros. A StdCellCluster only has std cells or dbModule_
    // instances
    enum ClusterType {
        StdCellCluster
    };

    // ***************************************************************************
    // * Metrics Class
    // ***************************************************************************
    // Metrics class for logical modules and clusters
    // For now we assue there's no macros in the design
    class Metrics {
        public:
            Metrics() = default;
            Metrics(unsigned int num_std_cells,
                    float std_cell_area);

            void AddMetrics(const Metrics& metrics);
            void InflateStdCellArea(float std_cell_util);

            // Set info
            void SetNumStdCell(unsigned int num_std_cell);
            void SetStdCellArea(float std_cell_area);
            void SetInflateStdCellArea(float percentage);

            // Get info
            unsigned int GetNumStdCell() const;
            float GetStdCellArea() const;
            float GetInflateStdCellArea() const;

        private:
            // During the autoclustering stage, we cluster the 
            // logical modules or clusters based on num_std_cell.
            unsigned int num_std_cell_ = 0;

            // For conducting the operations in SA, we need the 
            // area for std cells and the inflated area
            float std_cell_area_ = 0.0;
            float inflate_std_cell_area_ = 0.0;
    }; // Metrics Class

    // ***************************************************************************
    // * Cluster Class
    // ***************************************************************************
    /**
     * @class Cluster
     * @author Zhiang Wang, Jiantao Liu
     * @date 03/17/2023
     * @brief Basic units that will be shared by both the autoclustering and Multilevel-SA
     *        We convert the original gate-level netlist into a cluster-level netlist. All 
     *        physical information will be stored in Clusters, Partitions in SA will just store
     *        high-level information about the annealing
     */
    class Cluster {
        public:
            // Constructor
            Cluster(int cluster_id, utl::Logger* logger); //cluster name can be updated
            Cluster(int cluster_id, const std::string& cluster_name, utl::Logger* logger);

            // Cluster id is unique for each cluster and can not be changed once assigned
            int GetId() const;
            // Cluster names
            const std::string GetName() const;
            void SetName(const std::string& name);
            // Set cluster type, for now we only support standard cell clusters.
            void SetClusterType(const ClusterType& cluster_type);
            const ClusterType GetClusterType() const;

            // Instances in the cluster, two types
            // 1) dbInst
            // 2) dbModule
            void AddDbModule(odb::dbModule* db_module);
            void AddLeafStdCell(odb::dbInst* leaf_std_cell);
            const std::vector<odb::dbModule*> GetDbModules() const;
            const std::vector<odb::dbInst*> GetLeafStdCells() const;
            void ClearDbModules();
            void ClearLeafStdCells();
            void CopyInstances(const Cluster& cluster); // Only based on cluster type
            void GetBundledNetConnections(std::map<int, std::vector<int>>& map) { map = bundled_net_map_; }
            void InitBundledNetConnection();

            // IO cluster
            // If the cluster is specified as an IO cluster,
            // the exact location of the cluster must be given
            void SetIOClusterFlag(const float x_cord,
                                  const float y_cord,
                                  const float z_cord,
                                  const float width,
                                  const float height);

            bool GetIOClusterFlag() const;

            // Metric supports
            void SetMetrics(const Metrics& metrics);
            void SetInflatedArea(float percentage);
            const Metrics GetMetrics() const;
            int GetNumStdCell() const;
            float GetStdCellArea() const;
            float GetInflatedArea() const;

            // Physical Location Information
            // Aspect ratio is assumed to be 1
            float GetWidth() const;
            float GetHeight() const;
            float GetX() const;
            float GetY() const;
            float GetZ() const;
            int GetPartitionId() const;
            void SetX(float x_cord);
            void SetY(float y_cord);
            void SetZ(float z_cord);
            void SetPartitionId(int partition_id);
            // [0]: x_cord; [1]: y_cord, [2]: z_cord
            void GetLocation(std::vector<float>& location) const;

            // Hierarchy Support
            void SetParent(Cluster* parent);
            void AddChild(Cluster* child);
            void RemoveChild(const Cluster* child);
            void AddChildren(const std::vector<Cluster*>& children);
            void RemoveChildren();
            Cluster* GetParent() const;
            std::vector<Cluster*> GetChildren() const;

            // If the cluster is a leaf cluster
            bool IsLeaf() const;
            bool MergeCluster(Cluster& cluster,
                              bool& delete_flag); // return true if succeed

            // Connection signature support
            void InitConnection();
            void AddConnection(int cluster_id, float weight);
            const std::map<int, float> GetConnection() const;
            bool IsSameConnectionSignature(const Cluster& cluster,
                                           float net_threshold);
            // Get closely-connected cluster if such cluster exists
            // For example, if a small cluster A is closely connected to a
            // well-formed cluster B, (there are also other well-formed clusters
            // C, D), A is only connected to B and A has no connection with C, D
            int GetCloseCluster(const std::vector<int>& candidate_cluster,
                                float net_threshold);
            
            // Virtual Connections
            const std::vector<std::pair<int, int>> GetVirtualConnections() const;
            void AddVirtualConnection(int src, int target);

            // SA related functions
            void AddClusterBundledNet(int target_cluster_id, int bundled_net_id);
            void AddNetBundledNet(odb::dbNet* sel_net, int bundled_net_id);

            // Print Basic Information
            void PrintBasicInformation(utl::Logger* logger) const;

        private:
            // Private variables
            int id_ = -1; // Cluster id (a valid cluster id should be nonnegative)
            std::string name_ = "";             // cluster name
            ClusterType type_ = StdCellCluster; // For now, we only support standard cell 
            
            // Two types of instances will be in the cluster
            // 1) Sub logical modules
            // 2) Leaf standard cells
            // dbModule is an object representing logical module in the OpenDB
            std::vector<odb::dbModule*> db_modules_;
            // Standard cell instances in the cluster
            std::vector<odb::dbInst*> leaf_std_cells_;

            // The boundled IOs are modeled as a cluster with no area
            // The position is set to the center of thsoe IOs
            bool io_cluster_flag_ = false;

            // Each clusters uses metrics to store its statistics --> Internal information
            Metrics metrics_;

            // Each cluster is a node in the physical hierarchy tree
            Cluster* parent_ = nullptr;         // Parent of the current cluster
            std::vector<Cluster*> children_;    // Children of the current cluster

            // Connection_map_ is used for grouping small clusters during auto-clustering
            std::map<int, float> connection_map_; // cluster_id, number of connections

            // Mapping structure for the bundled net building process
            std::map<int, std::vector<int>> bundled_net_map_; // cluster_id, vector of bundled net index
            // std::map<odb::dbNet*, std::vector<int>> net_to_bundled_net_map_; // dbNet , vector of bundled_net index map

            // Store the virtual connections between children.
            // The virtual connections are used to embed the timing information between different
            // clusters.
            std::vector<std::pair<int, int>> virtual_connections_;

            utl::Logger* logger_;

            // For now we just have soft macros (standard cells), we don't differentiate between the 
            // physical entity and the logical abstraction.
            // Physical information
            float x_ = 0.0;
            float y_ = 0.0;
            float z_ = 0.0;
            int partition_id_ = -1; // Valid partition_id shall always be positive
            float height_ = 0.0;    // Optional, for now our constraint is area instead of the actual shape
            float width_ = 0.0;     

            // Status flag
            bool fixed_ = false; // Do not touch this cluster
    }; // Cluster Class

    // ***************************************************************************
    // * BundledNet Struct
    // ***************************************************************************
    // We use a clique model to represent all multi pin net in the netlist
    // Thus, we just have two-pin nets
    struct BundledNet {
        BundledNet(int src, int dst, float clique_factor, float timing_weight) {
            this->src_cluster_id = src;
            this->target_cluster_id = dst;
            this->timing_weight = timing_weight;
            this->clique_factor = clique_factor;
            this->hpwl = 0;
        }

        bool operator==(const BundledNet& net) {
            return (src_cluster_id == net.src_cluster_id) &&
                    (target_cluster_id == net.target_cluster_id);
        }

        int GetOpposite(int sel_cluster_id) {
            if (sel_cluster_id == src_cluster_id) {
                return target_cluster_id;
            } else {
                return src_cluster_id;
            }
        }

        // Number of bundled connections (can incorporate the weight of timing info)
        float clique_factor = 1.0;
        float timing_weight;
        // SA calculation support
        float hpwl; // HPWL of the net

        int src_cluster_id = -1;
        int target_cluster_id = -1;
    };

    // ***************************************************************************
    // * Timing Related Struct
    // ***************************************************************************
    // Class used for sta timing extraction
    // Three types of information can be used to represent a timing path
    // A timing path can be
    // 1) A sequence of vertices(insts): a -> b -> c;
    // 2) A sequence of hyperedges: e1 -> e2 -> e3;
    // 3) A sequence of cluster-ids : c1 -> c2 -> c3;
    // Here we separate vertices and clusters just to ease the search and weighting process 
    struct TimingPath {
        std::vector<int> path;       // A list of vertex ids -> path-based method
        std::vector<int> arcs;       // A list of hyperedge ids -> net-based methods
        std::vector<int> hyper_arcs; // A list of cluster ids -> net-based methods
        // Slack for the corresponding critical timing paths
        // Normalized to the specified clock period
        float slack = 0.0;
    };

    // ***************************************************************************
    // * Bottom-up grouing 
    // ***************************************************************************
    // Hypergraph representation 
    struct Hypergraph {
        // Define data structures for the intended hypergraph
        std::map<int, int> cluster_vertex_id_map;               // Fixed vertices map form cluster_id to vertex_id
        std::map<odb::dbInst*, int> inst_vertex_id_map;         // Map from instance to vertex_id
        std::map<odb::dbNet*, int> net_hyperedge_id_map;        // Map from net to hyperedge_id
        std::map<int, odb::dbInst*> vertex_id_instance_map;     // Map from vertex id to the corresponding instance
        std::vector<std::vector<int>> hyperedges;               // Structure for hyperedge representation
        std::vector<float> vertex_weight;                       // Vertex weight will be the sum of areas of the vertex
        std::vector<float> hyperedge_weight;                    // Hyperedge weight is a combination of connectivity and path slack
        std::map<int, float> hyperedge_id_net_slack_weight_map; // Map from hyperedge_id to net slack weight
        std::map<int, int> vertex_id_group_id_map;              // Map from movable vertex_id to group id
    };

    // Group representations for a cluster of vertices in the hypergraph
    struct Group {
        // Define the data structures for the intended group
        std::vector<int> included_vertices;  // Movable vertices that belong to this, always the first level
        std::set<int> neighbor_vetex_id;     // Vertex_id or group_id that's are connected with this group
        std::set<int> connected_hyperedge_id;// Hyperedges that are connected with this group;

        // Indexing parameters
        int group_id = -1;
        float group_size = 0;
    };

    // Clustering score class that will be used in the priority queue
    // Three elements will be included in the Priority Queue Element
    // 1) pair<int, int> : which stores the group_id pairs. (pair(u, v))
    // 2) d(u, v) : which stores the clustering score.
    // 3) valid check: which stores whether the element is valid or not
    class PQElement {
        public:
            PQElement(std::pair<int, int> group_pair,
                    float clustering_score,
                    bool valid_flag)
                    : group_pair(group_pair), valid(valid_flag) {
                cluster_score = clustering_score;
            }
            
            bool ValidEntry() { return valid; }

            std::pair<int, int> group_pair;
            float cluster_score;
            bool valid;
    };

    // Custom comparator functor for handling PQElement class
    struct PQComparator {
        bool operator()(const PQElement& a, const PQElement& b) {
            if (a.cluster_score == b.cluster_score) {
                return true;
            }
            return a.cluster_score < b.cluster_score;
        }
    };

    // Define the custom priority queue
    template<typename T, typename Container = std::vector<T>, typename Compare = PQComparator>
    class AccessiblePriorityQueue: public std::priority_queue<T, Container, Compare> {
        public:
            // Add a method to access the underlying container
            Container& get_container() {
                return this->c;
            }
    };


    // ***************************************************************************
    // * Simulated Annealing Related Representation Definition
    // ***************************************************************************
    // Struct used for representing the bounding box
    struct Bbox {
        // Define the bottom left and top right coordinates
        float lx = 0.0;
        float ly = 0.0;
        float ux = 0.0;
        float uy = 0.0;
    };

    /**
     * @class Vertex
     * @author Jiantao Liu
     * @date 04/10/2023
     * @brief This class is used to represent the original cluster class in a much lighter way
     *        in Simulated Annealing
    */
    class Vertex {
        public:
            // Constructor
            Vertex(Cluster* cluster, int vertex_id) {
                cluster_id_ = cluster->GetId();
                vertex_id_ = vertex_id;
                std_area_ = cluster->GetStdCellArea();
                x_ = cluster->GetX();
                y_ = cluster->GetY();
                z_ = cluster->GetZ();
                io_cluster_flag_ = cluster->GetIOClusterFlag();
            } // Create the Vertex from a cluster

            // Helper function
            void PrintVertexInfo() {
                std::cout << "Vertex Id : " << vertex_id_ << "; Cluster Id: " << cluster_id_ << std::endl;
                std::cout << "\tStd Cell Area: " << std_area_ << std::endl;
                PrintConnections();
            }

            // Set inner parameters
            void SetVertexId(const int vertex_id) {
                vertex_id_ = vertex_id;
            }

            void SetXCord(const float x_cord) {
                x_ = x_cord;
            }

            void SetYCord(const float y_cord) {
                y_ = y_cord;
            }

            void SetZCord(const float z_cord) {
                z_ = z_cord;
            }

            void SetPartitionIndex(const int new_partition_index) {
                partition_id_ = new_partition_index;
            }

            void SetIOFlag(bool IO_flag) {
                io_cluster_flag_ = IO_flag;
            }

            void AddNet(const int bundled_net_index) {
                connected_net_.insert(bundled_net_index);
            }

            // Get info
            float GetVertexArea() const {
                return std_area_;
            }

            int GetVertexClusterId() const {
                return cluster_id_;
            }

            float GetVertexId() const {
                return vertex_id_;
            }

            int GetPartitionId() const {
                return partition_id_;
            }

            float GetXCord() const {
                return x_;
            }

            float GetYCord() const {
                return y_;
            }

            float GetZCord() const {
                return z_;
            }

            void GetNets(std::vector<std::vector<float>>& influenced_nets) const {
                influenced_nets = fixed_vertex_connection_map_;
            }

            const std::map<int, float>& GetMap() const {
                return moveable_vertex_id_weight_map_;
            }

            void PrintConnections() {
                std::cout << "\t[Net with fixed cluster] " << ": ";

                for (const auto sel_net : connected_net_) {
                    std::cout << sel_net;
                    std::cout << " ";
                }

                std::cout << std::endl;

                std::cout << "\t[Net with movable vertices] : ";

                for (const auto sel_net : moveable_vertex_id_weight_map_) {
                    std::cout << sel_net.first;
                    std::cout << " ";
                }

                std::cout << std::endl;
            }

            // Connected bundled net
            std::set<int> connected_net_; // Index of the connected bundled net - fixed
            std::vector<std::vector<float>> fixed_vertex_connection_map_; // vector elements: x_cord, y_cord, z_cord, timing_weight;
            std::map<int, float> moveable_vertex_id_weight_map_; // Map from movable vertex id to timing weight.

        private:
            int cluster_id_ = -1;
            int vertex_id_ = -1;
            float std_area_ = 0.0;
            float x_ = 0.0;
            float y_ = 0.0;
            float z_ = 0.0;
            int partition_id_ = -1;

            // IO vertex flag, which will be fixed during the annealing process
            bool io_cluster_flag_ = false;
    };

    class Edge {
        public:
            Edge(int src_vertex_id, 
                 int dst_vertex_id,
                 int edge_index, 
                 bool full_movable, 
                 float fix_x,
                 float fix_y,
                 float fix_z,
                 float timing_weight) {
                    full_movable_ = full_movable;
                    edge_index_ = edge_index;

                    // All movable vertex will be set to src first
                    src_vertex_id_ = src_vertex_id;
                    dst_vertex_id_ = dst_vertex_id;

                    fix_x_cord_ = fix_x;
                    fix_y_cord_ = fix_y;
                    fix_z_cord_ = fix_z;

                    timing_weight_ = timing_weight;
                }

            // Set info
            void SetSrcVertex(int src_vertex_id) {
                src_vertex_id_ = src_vertex_id;
            }

            void SetDstVertex(int dst_vertex_id) {
                dst_vertex_id_ = dst_vertex_id;
            }

            void SetMovableFlag(bool movable_flag) {
                full_movable_ = movable_flag;
            }

            void SetEdgeIndex(int edge_index) {
                edge_index_ = edge_index;
            }

            void SetXCord(float x_cord) {
                fix_x_cord_ = x_cord;
            }

            void SetYCord(float y_cord) {
                fix_y_cord_ = y_cord;
            }

            void SetZCord(float z_cord) {
                fix_z_cord_ = z_cord;
            }

            void SetTimingWeight (float timing_weight) {
                timing_weight_ = timing_weight;
            }

            // Get Info
            bool GetFullMovFlag() {
                return full_movable_;
            }

            int GetEdgeIndex() {
                return edge_index_;
            }

            int GetSrcIndex() {
                return src_vertex_id_;
            }

            int GetDstIndex() {
                return dst_vertex_id_;
            }

            float GetXCord() {
                return fix_x_cord_;
            }

            float GetYCord() {
                return fix_y_cord_;
            }

            float GetZCord() {
                return fix_z_cord_;
            }

            float GetTimingWeight() {
                return timing_weight_;
            }
    
        private:
            bool full_movable_ = false;
            int edge_index_ = 0;

            int src_vertex_id_ = -1;
            int dst_vertex_id_ = -1;
            float fix_x_cord_ = 0.0;
            float fix_y_cord_ = 0.0;
            float fix_z_cord_ = 0.0;
            float timing_weight_ = 0.0;
    };

    class Partition {
        public:
            Partition(
                float x_center_cord,
                float y_center_cord,
                float z_center_cord,
                float width,
                float height,
                float overall_capacity,
                float allowed_capacity,
                int partition_id
            ) {
                x_center_cord_ = x_center_cord;
                y_center_cord_ = y_center_cord;
                z_center_cord_ = z_center_cord;
                width_ = width;
                height_ = height;
                overall_capacity_ = overall_capacity;
                allowed_capacity_ = allowed_capacity;
                partition_id_ = partition_id;
            }

            // Set inner variables
            void SetXCord(float new_x_cord) {
                x_center_cord_ = new_x_cord;
            }

            void SetYCord(float new_y_cord) {
                y_center_cord_ = new_y_cord;
            }

            void SetZCord(float new_z_cord) {
                z_center_cord_ = new_z_cord;
            }

            void SetWidth(float new_width) {
                width_ = new_width;
            }

            void SetHeight(float new_height) {
                height_ = new_height;
            }

            void SetPartitionID(int partition_id) {
                partition_id_ = partition_id;
            }

            void SetOverallCapacity(float new_capacity) {
                overall_capacity_ = new_capacity;
            }

            void SetAllowedCapacity(float allowed_capacity) {
                allowed_capacity_ = allowed_capacity;
            }

            void AddCapacity(const float added_capacity) {
                used_capacity_ += added_capacity;
            }

            void MinusCapacity(const float capacity) {
                used_capacity_ -= capacity;
            }

            // Get inner info
            float GetXCord() const {
                return x_center_cord_;
            }

            float GetYCord() const {
                return y_center_cord_;
            }

            float GetZCord() const {
                return z_center_cord_;
            }

            float GetWidth() const {
                return width_;
            }

            float GetHeight() const {
                return height_;
            }

            float GetOverallCapacity() const {
                return overall_capacity_;
            }

            float GetUsedCapacity() const {
                return used_capacity_;
            }

            float GetAllowedCapacity() const {
                return allowed_capacity_;
            }

            int GetPartitionID() const {
                return partition_id_;
            }

        private:
            float x_center_cord_ = 0.0;
            float y_center_cord_ = 0.0;
            float z_center_cord_ = 0.0;
            float width_ = 0.0;
            float height_ = 0.0;
            float overall_capacity_ = 0.0;
            float allowed_capacity_ = 0.0;
            float used_capacity_ = 0.0;
            int partition_id_ = 0;
    };

    
} // namespace mlsa