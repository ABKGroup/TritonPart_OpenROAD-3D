#include "object.h"

#include <iostream>

#include "odb/db.h"
#include "utl/Logger.h"


namespace mlsa {
    // ***************************************************************************
    // * Basic Functions
    // ***************************************************************************

    // Conversion between dbu and microns
    float DbuToMicron(int metric, float dbu) {
        return metric/dbu;
    }

    int MicronToDbu(float metric, float dbu) {
        return std::round(metric * dbu);
    }

    // ***************************************************************************
    // * Metrics Class
    // ***************************************************************************
    Metrics::Metrics(
        unsigned int num_std_cell,
        float std_cell_area
    ) {
        num_std_cell_ = num_std_cell;
        std_cell_area_ = std_cell_area;
    }

    void Metrics::AddMetrics(const Metrics& metrics) {
        num_std_cell_ += metrics.GetNumStdCell();
        std_cell_area_ += metrics.GetStdCellArea();
        inflate_std_cell_area_ += metrics.GetInflateStdCellArea();
    }

    void Metrics::InflateStdCellArea(float std_cell_util) {
        if ((std_cell_util > 0.0) && (std_cell_util < 1.0)) {
            inflate_std_cell_area_ /= std_cell_util;
        }
    }

    void Metrics::SetNumStdCell(unsigned int num_std_cell) {
        num_std_cell_ = num_std_cell;
    }

    void Metrics::SetStdCellArea(float std_cell_area) {
        std_cell_area_ = std_cell_area;
    }

    void Metrics::SetInflateStdCellArea(float percentage) {
        inflate_std_cell_area_ = std_cell_area_ * ( 1 + percentage);
    }

    unsigned int Metrics::GetNumStdCell() const {
        return num_std_cell_;
    }

    float Metrics::GetStdCellArea() const {
        return std_cell_area_;
    }

    float Metrics::GetInflateStdCellArea() const {
        return inflate_std_cell_area_;
    }

    // ***************************************************************************
    // * Cluster Class
    // ***************************************************************************
    // Constructors and Destructors
    Cluster::Cluster(int cluster_id, utl::Logger* logger) {
        id_ = cluster_id;
        logger_ = logger;
    }

    Cluster::Cluster(int cluster_id,
                     const std::string& cluster_name,
                     utl::Logger* logger) {
        id_ = cluster_id;
        name_ = cluster_name;
        logger_ = logger;
    }

    // Cluster id
    int Cluster::GetId() const {
        return id_;
    }

    const std::string Cluster::GetName() const {
        return name_;
    }

    void Cluster::SetName(const std::string& name) {
        name_ = name;
    }

    // Cluster Type
    // For now we only support StdCellCluster
    void Cluster::SetClusterType(const ClusterType& cluster_type) {
        type_ = cluster_type;
    }

    const ClusterType Cluster::GetClusterType() const {
        return type_;
    }

    // Instances in the Cluster
    // StdCellCluster: dbModule & dbInst
    void Cluster::AddDbModule(odb::dbModule* db_module) {
        db_modules_.push_back(db_module);
    }

    void Cluster::AddLeafStdCell(odb::dbInst* leaf_std_cell) {
        leaf_std_cells_.push_back(leaf_std_cell);
    }

    const std::vector<odb::dbModule*> Cluster::GetDbModules() const {
        return db_modules_;
    }

    const std::vector<odb::dbInst*> Cluster::GetLeafStdCells() const {
        return leaf_std_cells_;
    }

    void Cluster::ClearDbModules() {
        db_modules_.clear();
    }

    void Cluster::ClearLeafStdCells() {
        leaf_std_cells_.clear();
    }

    // Copy all instances based on the given cluster
    void Cluster::CopyInstances(const Cluster& cluster) {
        // Clear all
        db_modules_.clear();
        leaf_std_cells_.clear();

        // Insert new elements
        // For now we only support StdCellCluster
        // if (type_ == StdCellCluster) {
            leaf_std_cells_.insert(leaf_std_cells_.end(),
                                   cluster.leaf_std_cells_.begin(),
                                   cluster.leaf_std_cells_.end());
            db_modules_.insert(db_modules_.end(),
                               cluster.db_modules_.begin(),
                               cluster.db_modules_.end());
        // }
    }

    // Supports for bundled IO clusters
    void Cluster::SetIOClusterFlag(const float x_cord,
                          const float y_cord,
                          const float z_cord,
                          const float width,
                          const float height) {
        io_cluster_flag_ = true;

        // Set the location of the IO cluster
        // Bundled IO cluster will be set to the center point of the actual
        // bundled IO. For now we assume there is no area related to the bundled IO
        x_ = x_cord;
        y_ = y_cord;
        z_ = z_cord;
        width_ = width;
        height_ = height;
    }

    bool Cluster::GetIOClusterFlag() const {
        return io_cluster_flag_;
    }

    // Functions related to Metrics
    void Cluster::SetMetrics(const Metrics& metrics) {
        metrics_ = metrics;
    }

    void Cluster::SetInflatedArea(float percentage) {
        metrics_.SetInflateStdCellArea(percentage);
    }

    const Metrics Cluster::GetMetrics() const {
        return metrics_;
    }

    int Cluster::GetNumStdCell() const {
        // For now we only support StdCellClusters
        return metrics_.GetNumStdCell();
    }

    float Cluster::GetStdCellArea() const {
        return metrics_.GetStdCellArea();
    }

    float Cluster::GetInflatedArea() const {
        return metrics_.GetInflateStdCellArea();
    }

    // Functions related to physical locations
    float Cluster::GetWidth() const {
        return width_;
    }

    float Cluster::GetHeight() const {
        return height_;
    }

    float Cluster::GetX() const {
        return x_;
    }

    float Cluster::GetY() const {
        return y_;
    }

    float Cluster::GetZ() const {
        return z_;
    }

    int Cluster::GetPartitionId() const {
        return partition_id_;
    }

    void Cluster::SetX(float x_cord) {
        x_ = x_cord;

        // Set X location for all it's childs, if any
        for (auto sel_child : GetChildren()) {
            sel_child->SetX(x_cord);
        }
    }

    void Cluster::SetY(float y_cord) {
        y_ = y_cord;

        // Set Y location for all it's childs, if any
        for (auto sel_child : GetChildren()) {
            sel_child->SetY(y_cord);
        }
    }

    void Cluster::SetZ(float z_cord) {
        z_ = z_cord;

        // Set Z location for all it's childs, if any
        for (auto sel_child : GetChildren()) {
            sel_child->SetZ(z_cord);
        }
    }

    void Cluster::SetPartitionId(int partition_id) {
        partition_id_ = partition_id;
    }

    // Get the X, Y, Z location of the chosen cluster
    // [0]: x, [1]: y, [2]: z
    void Cluster::GetLocation(std::vector<float>& location) const {
        location.push_back(x_);
        location.push_back(y_);
        location.push_back(z_);
    }

    // Functions related to the hierarchy structure
    void Cluster::SetParent(Cluster* parent) {
        parent_ = parent;
    }

    void Cluster::AddChild(Cluster* child) {
        children_.push_back(child);
    }

    void Cluster::RemoveChild(const Cluster* child) {
        children_.erase(std::find(children_.begin(), children_.end(), child));
    }

    void Cluster::AddChildren(const std::vector<Cluster*>& children) {
        std::copy(children.begin(), children.end(), std::back_inserter(children_));
    }

    void Cluster::RemoveChildren() {
        children_.clear();
    }

    Cluster* Cluster::GetParent() const {
        return parent_;
    }

    std::vector<Cluster*> Cluster::GetChildren() const {
        return children_;
    }

    // Check whether the cluster is a leaf cluster
    bool Cluster::IsLeaf() const {
        return (children_.size() == 0);
    }

    bool Cluster::MergeCluster(Cluster& cluster, bool& delete_flag) {
        // We only merge clusters that belong to the same parent cluster
        if (parent_ != cluster.parent_) {
            return false;
        }

        parent_->RemoveChild(&cluster);
        metrics_.AddMetrics(cluster.metrics_);

        // Change the name of the cluster
        name_ += "||" + cluster.name_;
        // First handle the components in these two clusters
        leaf_std_cells_.insert(leaf_std_cells_.end(),
                               cluster.leaf_std_cells_.begin(),
                               cluster.leaf_std_cells_.end());
        db_modules_.insert(db_modules_.end(),
                           cluster.db_modules_.begin(),
                           cluster.db_modules_.end());
        delete_flag = true;

        // Check whether current cluster is a leaf cluster or not.
        // TODO: Think about other ways to handle this merging.
        // If current cluster is not a leaf cluster
        if (children_.size() > 0) {
            children_.push_back(&cluster);
            cluster.SetParent(this);
            delete_flag = false;
        }

        return true;
    }

    // Functions related to connection signature
    void Cluster::InitConnection() {
        connection_map_.clear();
    }

    void Cluster::InitBundledNetConnection() {
        bundled_net_map_.clear();
    }

    void Cluster::AddConnection(int cluster_id, float weight) {
        if (connection_map_.find(cluster_id) == connection_map_.end()) {
            connection_map_[cluster_id] = weight;
        } else {
            connection_map_[cluster_id] += weight;
        }
    }

    // Add new bundled net id to the inner target cluster to bundled net index mapping structure
    void Cluster::AddClusterBundledNet(int target_cluster_id, int bundled_net_id) {
        if (bundled_net_map_.find(target_cluster_id) == bundled_net_map_.end()) {
            bundled_net_map_[target_cluster_id] = {bundled_net_id};
        } else {
            bundled_net_map_[target_cluster_id].push_back(bundled_net_id);
        }
    }

    const std::map<int, float> Cluster::GetConnection() const {
        return connection_map_;
    }

    // The connection signature is based on connection topology
    // if the number of connections between two clusters is larger than the
    // net_threshold we think the two clusters are connected, otherwise they
    // are disconnected
    bool Cluster::IsSameConnectionSignature(const Cluster& cluster,
                                            float net_threshold) {
        // Neighbors of the current cluster
        std::vector<int> neighbors;
        // Neighbors of the input cluster
        std::vector<int> cluster_neighbors;

        // Build connection signature for the current cluster
        for (auto& [cluster_id, weight] : connection_map_) {
            if ((cluster_id != id_) && cluster_id != cluster.id_
                 && (weight >= net_threshold)) {
                    neighbors.push_back(cluster_id);
                }
        }

        // Build connection signature for the input cluster
        for (auto& [cluster_id, weight] : cluster.connection_map_) {
            if ((cluster_id != id_) && (cluster_id != cluster.id_)
                && (weight >= net_threshold)) {
                    cluster_neighbors.push_back(cluster_id);
                }
        }

        // Check the connection signature
        if (neighbors.size() != cluster_neighbors.size()) {
            return false;
        } else {
            std::sort(neighbors.begin(), neighbors.end());
            std::sort(cluster_neighbors.begin(), cluster_neighbors.end());
            for (int i = 0; i < neighbors.size(); i++) {
                if (neighbors[i] != cluster_neighbors[i]) {
                    return false;
                }
            }
        }

        return true;
    }

    //
    // Get closely-connected cluster if such cluster exists
    // For example, if a small cluster A is closely connected to a
    // well-formed cluster B, (there are also other well-formed clusters
    // C, D), A is only connected to B and A has no connection with C, D
    //
    // candidate_clusters are small clusters that need to be merged,
    // any cluster not in candidate_clusters is a well-formed cluster
    //
    int Cluster::GetCloseCluster(const std::vector<int>& candidate_cluster,
                                float net_threshold) {
        int close_cluster = -1;
        int num_close_clusters = 0;

        for (auto& [cluster_id, num_nets] : connection_map_) {
            if (num_nets > net_threshold &&
                std::find(candidate_cluster.begin(), candidate_cluster.end(), cluster_id)
                        == candidate_cluster.end()) {
                num_close_clusters++;
                close_cluster = cluster_id;                
            }
        }

        if (num_close_clusters == 1) {
            return close_cluster;
        } else {
            return -1;
        }
    }

    // Functions related to virtual connections
    const std::vector<std::pair<int, int>> Cluster::GetVirtualConnections() const {
        return virtual_connections_;
    }

    void Cluster::AddVirtualConnection(int src, int target) {
        virtual_connections_.push_back(std::pair<int, int>(src, target));
    }

    // Print all basic information
    void Cluster::PrintBasicInformation(utl::Logger* logger) const {
        std::string line = "\n";

        line += std::string(80, '*') + "\n";
        line += "[INFO] cluster_name :  " + name_ + "   ";
        line += "cluster_id : " + std::to_string(id_) + " \n";
        line += "num_std_cell : " + std::to_string(GetNumStdCell()) + "  ";
        line += "area : " + std::to_string(GetStdCellArea()) + "  ";
        line += "location : (" + std::to_string(GetX()) + " , ";
        line += std::to_string(GetY()) + " , " + std::to_string(GetZ()) + ")\n";

        // Call Logger
        logger->report(line); 
    }

} // Namespace mlsa