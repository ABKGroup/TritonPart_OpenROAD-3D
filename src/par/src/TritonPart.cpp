///////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2022, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////
// High-level description
// This is the interfaces for TritonPart
// It works in two ways:
// 1) default mode :  read a verilog netlist and extract the hypergraph based on
// netlist 2) classical mode : read the hypergraph file in hmetis format
///////////////////////////////////////////////////////////////////////////////
#include "TritonPart.h"

#include <iostream>
#include <set>
#include <string>

#include "TPCoarsener.h"
#include "TPHypergraph.h"
#include "TPMultilevel.h"
#include "TPPartitioner.h"
#include "TPRefiner.h"
#include "TritonPart.h"
#include "Utilities.h"
#include "odb/db.h"
#include "sta/ArcDelayCalc.hh"
#include "sta/Bfs.hh"
#include "sta/Corner.hh"
#include "sta/DcalcAnalysisPt.hh"
#include "sta/ExceptionPath.hh"
#include "sta/FuncExpr.hh"
#include "sta/Graph.hh"
#include "sta/GraphDelayCalc.hh"
#include "sta/Liberty.hh"
#include "sta/Network.hh"
#include "sta/PathAnalysisPt.hh"
#include "sta/PathEnd.hh"
#include "sta/PathExpanded.hh"
#include "sta/PathRef.hh"
#include "sta/PatternMatch.hh"
#include "sta/PortDirection.hh"
#include "sta/Sdc.hh"
#include "sta/Search.hh"
#include "sta/SearchPred.hh"
#include "sta/Sequential.hh"
#include "sta/Sta.hh"
#include "sta/Units.hh"
#include "utl/Logger.h"

using utl::PAR;

namespace par {

// -----------------------------------------------------------------------------------
// Public functions
// -----------------------------------------------------------------------------------

// The function for partitioning a hypergraph
// This is used for replacing hMETIS
// Key supports:
// (1) fixed vertices constraint in fixed_file
// (2) community attributes in community_file (This can be used to guide the
// partitioning process) (3) stay together attributes in group_file. (4)
// placement information is specified in placement file The format is that each
// line cooresponds to a group fixed vertices, community and placement
// attributes both follows the hMETIS format
void TritonPart::PartitionHypergraph(unsigned int num_parts_arg,
                                     float balance_constraint_arg,
                                     unsigned int seed_arg,
                                     int vertex_dimension_arg,
                                     int hyperedge_dimension_arg,
                                     int placement_dimension_arg,
                                     const char* hypergraph_file_arg,
                                     const char* fixed_file_arg,
                                     const char* community_file_arg,
                                     const char* group_file_arg,
                                     const char* placement_file_arg)
{
  logger_->report("========================================");
  logger_->report("[STATUS] Starting TritonPart Partitioner");
  logger_->report("========================================");
  logger_->report("[INFO] Partitioning parameters**** ");
  // Parameters
  num_parts_ = num_parts_arg;
  ub_factor_ = balance_constraint_arg;
  seed_ = seed_arg;
  vertex_dimensions_ = vertex_dimension_arg;
  hyperedge_dimensions_ = hyperedge_dimension_arg;
  placement_dimensions_ = placement_dimension_arg;
  // local parameters
  std::string hypergraph_file = hypergraph_file_arg;
  std::string fixed_file = fixed_file_arg;
  std::string community_file = community_file_arg;
  std::string group_file = group_file_arg;
  std::string placement_file = placement_file_arg;
  // solution file
  std::string solution_file
      = hypergraph_file + std::string(".part.") + std::to_string(num_parts_);
  logger_->report("[PARAM] Number of partitions = {}", num_parts_);
  logger_->report("[PARAM] UBfactor = {}", ub_factor_);
  logger_->report("[PARAM] Seed = {}", seed_);
  logger_->report("[PARAM] Vertex dimensions = {}", vertex_dimensions_);
  logger_->report("[PARAM] Hyperedge dimensions = {}", hyperedge_dimensions_);
  logger_->report("[PARAM] Placement dimensions = {}", placement_dimensions_);
  logger_->report("[PARAM] Hypergraph file = {}", hypergraph_file);
  logger_->report("[PARAM] Solution file = {}", solution_file);
  logger_->report("[PARAM] Global net threshold = {}", global_net_threshold_);
  if (!fixed_file.empty()) {
    logger_->report("[PARAM] Fixed file  = {}", fixed_file);
  }
  if (!community_file.empty()) {
    logger_->report("[PARAM] Community file = {}", community_file);
  }
  if (!group_file.empty()) {
    logger_->report("[PARAM] Group file = {}", group_file);
  }
  if (!placement_file.empty()) {
    logger_->report("[PARAM] Placement file = {}", placement_file);
  }

  // set the random seed
  srand(seed_);  // set the random seed

  timing_aware_flag_ = false;
  logger_->report(
      "[WARNING] Reset the timing_aware_flag to false. Timing-driven mode is "
      "not supported");

  // build hypergraph: read the basic hypergraph information and other
  // constraints
  ReadHypergraph(
      hypergraph_file, fixed_file, community_file, group_file, placement_file);

  // call the multilevel partitioner to partition hypergraph_
  // but the evaluation is the original_hypergraph_
  MultiLevelPartition();

  // write the solution in hmetis format
  std::ofstream solution_file_output;
  solution_file_output.open(solution_file);
  for (auto part_id : solution_) {
    solution_file_output << part_id << std::endl;
  }
  solution_file_output.close();

  // finish hypergraph partitioning
  logger_->report("===============================================");
  logger_->report("Exiting TritonPart");
}

// Top level interface
// The function for partitioning a hypergraph
// This is the main API for TritonPart
// Key supports:
// (1) fixed vertices constraint in fixed_file
// (2) community attributes in community_file (This can be used to guide the
// partitioning process) (3) stay together attributes in group_file. (4)
// timing-driven partitioning (5) fence-aware partitioning (6) placement-aware
// partitioning, placement information is extracted from OpenDB
void TritonPart::PartitionDesign(unsigned int num_parts_arg,
                                 float balance_constraint_arg,
                                 unsigned int seed_arg,
                                 bool timing_aware_flag_arg,
                                 int top_n_arg,
                                 bool placement_flag_arg,
                                 bool fence_flag_arg,
                                 float fence_lx_arg,
                                 float fence_ly_arg,
                                 float fence_ux_arg,
                                 float fence_uy_arg,
                                 const char* fixed_file_arg,
                                 const char* community_file_arg,
                                 const char* group_file_arg,
                                 const char* solution_filename_arg)
{
  logger_->report("========================================");
  logger_->report("[STATUS] Starting TritonPart Partitioner");
  logger_->report("========================================");
  logger_->report("[INFO] Partitioning parameters**** ");
  const int dbu = db_->getTech()->getDbUnitsPerMicron();
  block_ = db_->getChip()->getBlock();
  // Parameters
  num_parts_ = num_parts_arg;
  ub_factor_ = balance_constraint_arg;
  seed_ = seed_arg;
  vertex_dimensions_ = 1;  // for design partitioning, vertex weight is the area
                           // of the instance
  hyperedge_dimensions_
      = 1;  // for design partitioning, hyperedge weight is the connectivity
  timing_aware_flag_ = timing_aware_flag_arg;
  if (timing_aware_flag_ == false) {
    top_n_ = 0;  // timing driven flow is disabled
  } else {
    top_n_ = top_n_arg;  // extract the top_n critical timing paths
  }
  placement_flag_ = placement_flag_arg;
  if (placement_flag_ == false) {
    placement_dimensions_ = 0;  // no placement information
  } else {
    placement_dimensions_ = 2;  // 2D canvas
  }
  fence_flag_ = fence_flag_arg;
  fence_ = Rect(fence_lx_arg * dbu,
                fence_ly_arg * dbu,
                fence_ux_arg * dbu,
                fence_uy_arg * dbu);
  if (fence_flag_ == false || fence_.IsValid() == false) {
    fence_.Reset();
  }
  // local parameters
  std::string fixed_file = fixed_file_arg;
  std::string community_file = community_file_arg;
  std::string group_file = group_file_arg;
  std::string solution_file = solution_filename_arg;
  logger_->report("[PARAM] Number of partitions = {}", num_parts_);
  logger_->report("[PARAM] UBfactor = {}", ub_factor_);
  logger_->report("[PARAM] Seed = {}", seed_);
  logger_->report("[PARAM] Vertex dimensions = {}", vertex_dimensions_);
  logger_->report("[PARAM] Hyperedge dimensions = {}", hyperedge_dimensions_);
  logger_->report("[PARAM] Placement dimensions = {}", placement_dimensions_);
  logger_->report("[PARAM] Timing aware flag = {}", timing_aware_flag_);
  logger_->report("[PARAM] Global net threshold = {}", global_net_threshold_);
  logger_->report("[PARAM] Top {} critical timing paths are extracted.",
                  top_n_);
  logger_->report("[PARAM] Fence aware flag = {}", fence_flag_);
  if (fence_flag_ == true) {
    logger_->report(
        "[PARAM] fence_lx = {}, fence_ly = {}, fence_ux = {}, fence_uy = {}",
        fence_.lx / dbu,
        fence_.ly / dbu,
        fence_.ux / dbu,
        fence_.uy / dbu);
  }
  if (!fixed_file.empty()) {
    logger_->report("[PARAM] Fixed file  = {}", fixed_file);
  }
  if (!community_file.empty()) {
    logger_->report("[PARAM] Community file = {}", community_file);
  }
  if (!group_file.empty()) {
    logger_->report("[PARAM] Group file = {}", group_file);
  }
  if (!solution_file.empty()) {
    logger_->report("[PARAM] Solution file = {}", solution_file);
  }

  // set the random seed
  srand(seed_);  // set the random seed

  // read the netlist from OpenDB
  // for IO port and insts (std cells and macros),
  // there is an attribute for vertex_id
  logger_->report("========================================");
  logger_->report("[STATUS] Reading netlist**** ");
  // if the fence_flag_ is true, only consider the instances within the fence
  ReadNetlist(fixed_file, community_file, group_file);
  logger_->report("[STATUS] Finish reading netlist****");

  // call the multilevel partitioner to partition hypergraph_
  // but the evaluation is the original_hypergraph_
  MultiLevelPartition();

  // Write out the solution.
  // Format 1: write the clustered netlist in verilog directly
  for (auto term : block_->getBTerms()) {
    auto vertex_id_property = odb::dbIntProperty::find(term, "vertex_id");
    const int vertex_id = vertex_id_property->getValue();
    if (vertex_id == -1) {
      continue;  // This instance is not used
    }
    const int partition_id = solution_[vertex_id];
    if (auto property = odb::dbIntProperty::find(term, "partition_id")) {
      property->setValue(partition_id);
    } else {
      odb::dbIntProperty::create(term, "partition_id", partition_id);
    }
  }

  for (auto inst : block_->getInsts()) {
    auto vertex_id_property = odb::dbIntProperty::find(inst, "vertex_id");
    const int vertex_id = vertex_id_property->getValue();
    if (vertex_id == -1) {
      continue;  // This instance is not used
    }
    const int partition_id = solution_[vertex_id];
    if (auto property = odb::dbIntProperty::find(inst, "partition_id")) {
      property->setValue(partition_id);
    } else {
      odb::dbIntProperty::create(inst, "partition_id", partition_id);
    }
  }

  // Format 2: write the explicit solution
  // each line :  instance_name  partition_id
  if (!solution_file.empty()) {
    std::string solution_file_name = solution_file;
    if (fence_flag_ == true) {
      // if the fence_flag_ is set to true, we need to update the solution file
      // to reflect the fence
      std::stringstream str_ss;
      str_ss.setf(std::ios::fixed);
      str_ss.precision(3);
      str_ss << ".lx_" << fence_.lx / dbu;
      str_ss << ".ly_" << fence_.ly / dbu;
      str_ss << ".ux_" << fence_.ux / dbu;
      str_ss << ".uy_" << fence_.uy / dbu;
      solution_file_name = solution_file_name + str_ss.str();
    }
    logger_->report("[INFO] Updated solution file name = {}",
                    solution_file_name);
    std::ofstream file_output;
    file_output.open(solution_file_name);

    for (auto term : block_->getBTerms()) {
      if (auto property = odb::dbIntProperty::find(term, "partition_id")) {
        file_output << term->getName() << "  ";
        file_output << property->getValue() << "  ";
        file_output << std::endl;
      }
    }

    for (auto inst : block_->getInsts()) {
      if (auto property = odb::dbIntProperty::find(inst, "partition_id")) {
        file_output << inst->getName() << "  ";
        file_output << property->getValue() << "  ";
        file_output << std::endl;
      }
    }
    file_output.close();
  }

  logger_->report("===============================================");
  logger_->report("Exiting TritonPart");
}

void TritonPart::EvaluateHypergraphSolution(unsigned int num_parts_arg,
                                            float balance_constraint_arg,
                                            int vertex_dimension_arg,
                                            int hyperedge_dimension_arg,
                                            const char* hypergraph_file_arg,
                                            const char* fixed_file_arg,
                                            const char* group_file_arg,
                                            const char* solution_file_arg)

{
  logger_->report("========================================");
  logger_->report("[STATUS] Starting Evaluating Hypergraph Solution");
  logger_->report("========================================");
  logger_->report("[INFO] Partitioning parameters**** ");
  // Parameters
  num_parts_ = num_parts_arg;
  ub_factor_ = balance_constraint_arg;
  seed_ = 0;  // use the default random seed (no meaning in this function)
  vertex_dimensions_ = vertex_dimension_arg;
  hyperedge_dimensions_ = hyperedge_dimension_arg;
  placement_dimensions_
      = 0;  // use the default value (no meaning in this function)
  // local parameters
  std::string hypergraph_file = hypergraph_file_arg;
  std::string fixed_file = fixed_file_arg;
  std::string group_file = group_file_arg;
  // solution file
  std::string solution_file = solution_file_arg;
  std::string community_file(
      "");  // no community file is used (no meaning in this function)
  std::string placement_file(
      "");  // no placement file is used (no meaning in this function)
  logger_->report("[PARAM] Number of partitions = {}", num_parts_);
  logger_->report("[PARAM] UBfactor = {}", ub_factor_);
  logger_->report("[PARAM] Seed = {}", seed_);
  logger_->report("[PARAM] Vertex dimensions = {}", vertex_dimensions_);
  logger_->report("[PARAM] Hyperedge dimensions = {}", hyperedge_dimensions_);
  logger_->report("[PARAM] Placement dimensions = {}", placement_dimensions_);
  logger_->report("[PARAM] Hypergraph file = {}", hypergraph_file);
  logger_->report("[PARAM] Solution file = {}", solution_file);
  if (!fixed_file.empty()) {
    logger_->report("[PARAM] Fixed file  = {}", fixed_file);
  }
  if (!community_file.empty()) {
    logger_->report("[PARAM] Community file = {}", community_file);
  }
  if (!group_file.empty()) {
    logger_->report("[PARAM] Group file = {}", group_file);
  }
  if (!placement_file.empty()) {
    logger_->report("[PARAM] Placement file = {}", placement_file);
  }

  int part_id = -1;
  std::ifstream solution_file_input(solution_file);
  if (!solution_file_input.is_open()) {
    logger_->error(
        PAR, 2511, "Can not open the solution file : {}", solution_file);
  }
  while (solution_file_input >> part_id) {
    solution_.push_back(part_id);
  }
  solution_file_input.close();

  // set the random seed
  srand(seed_);  // set the random seed

  timing_aware_flag_ = false;
  logger_->report(
      "[WARNING] Reset the timing_aware_flag to false. Timing-driven mode is "
      "not supported");

  // build hypergraph: read the basic hypergraph information and other
  // constraints
  ReadHypergraph(
      hypergraph_file, fixed_file, community_file, group_file, placement_file);

  // check the weighting scheme
  if (static_cast<int>(e_wt_factors_.size()) != hyperedge_dimensions_) {
    logger_->report(
        "[WARNING] no hyperedge weighting is specified. Use default value of "
        "1.");
    e_wt_factors_.clear();
    e_wt_factors_.resize(hyperedge_dimensions_);
    std::fill(e_wt_factors_.begin(), e_wt_factors_.end(), 1.0);
  }
  logger_->report("[PARAM] hyperedge weight factor : [ {} ]",
                  GetVectorString(e_wt_factors_));

  if (static_cast<int>(v_wt_factors_.size()) != vertex_dimensions_) {
    logger_->report(
        "[WARNING] no vertex weighting is specified. Use default value of 1.");
    v_wt_factors_.clear();
    v_wt_factors_.resize(vertex_dimensions_);
    std::fill(v_wt_factors_.begin(), v_wt_factors_.end(), 1.0);
  }
  logger_->report("[PARAM] vertex weight factor : [ {} ]",
                  GetVectorString(v_wt_factors_));

  if (static_cast<int>(placement_wt_factors_.size()) != placement_dimensions_) {
    if (placement_dimensions_ <= 0) {
      placement_wt_factors_.clear();
    } else {
      logger_->report(
          "[WARNING] no placement weighting is specified. Use default value of "
          "1.");
      placement_wt_factors_.clear();
      placement_wt_factors_.resize(placement_dimensions_);
      std::fill(
          placement_wt_factors_.begin(), placement_wt_factors_.end(), 1.0f);
    }
  }
  logger_->report("[PARAM] placement weight factor : [ {} ]",
                  GetVectorString(placement_wt_factors_));

  // following parameters are not used
  net_timing_factor_ = 0.0;
  path_timing_factor_ = 0.0;
  path_snaking_factor_ = 0.0;
  timing_exp_factor_ = 0.0;
  extra_delay_ = 0.0;

  // print all the weighting parameters
  logger_->report("[PARAM] net_timing_factor : {}", net_timing_factor_);
  logger_->report("[PARAM] path_timing_factor : {}", path_timing_factor_);
  logger_->report("[PARAM] path_snaking_factor : {}", path_snaking_factor_);
  logger_->report("[PARAM] timing_exp_factor : {}", timing_exp_factor_);
  logger_->report("[PARAM] extra_delay : {}", extra_delay_);

  // create the evaluator class
  TP_evaluator_ptr evaluator
      = std::make_shared<GoldenEvaluator>(num_parts_,
                                          // weight vectors
                                          e_wt_factors_,
                                          v_wt_factors_,
                                          placement_wt_factors_,
                                          // timing related weight
                                          net_timing_factor_,
                                          path_timing_factor_,
                                          path_snaking_factor_,
                                          timing_exp_factor_,
                                          extra_delay_,
                                          original_hypergraph_,
                                          logger_);

  evaluator->ConstraintAndCutEvaluator(
      original_hypergraph_, solution_, ub_factor_, group_attr_, true);

  logger_->report("===============================================");
  logger_->report("Exiting Evaluating Hypergraph Solution");
}

// Function to evaluate the hypergraph partitioning solution
// This can be used to write the timing-weighted hypergraph
// and evaluate the solution.
// If the solution file is empty, then this function is to write the
// solution. If the solution file is not empty, then this function is to
// evaluate the solution without writing the hypergraph again This function is
// only used for testing
void TritonPart::EvaluatePartDesignSolution(
    unsigned int num_parts_arg,
    float balance_constraint_arg,
    bool timing_aware_flag_arg,
    int top_n_arg,
    bool fence_flag_arg,
    float fence_lx_arg,
    float fence_ly_arg,
    float fence_ux_arg,
    float fence_uy_arg,
    const char* fixed_file_arg,
    const char* community_file_arg,
    const char* group_file_arg,
    const char* hypergraph_file_arg,
    const char* hypergraph_int_weight_file_arg,
    const char* solution_filename_arg)
{
  logger_->report("========================================");
  logger_->report("[STATUS] Starting TritonPart Partitioner");
  logger_->report("========================================");
  logger_->report("[INFO] Partitioning parameters**** ");
  const int dbu = db_->getTech()->getDbUnitsPerMicron();
  block_ = db_->getChip()->getBlock();
  // Parameters
  num_parts_ = num_parts_arg;
  ub_factor_ = balance_constraint_arg;
  seed_ = 0;               // This parameter is not used.  just a random value
  vertex_dimensions_ = 1;  // for design partitioning, vertex weight is the area
                           // of the instance
  hyperedge_dimensions_
      = 1;  // for design partitioning, hyperedge weight is the connectivity
  timing_aware_flag_ = timing_aware_flag_arg;
  if (timing_aware_flag_ == false) {
    top_n_ = 0;  // timing driven flow is disabled
  } else {
    top_n_ = top_n_arg;  // extract the top_n critical timing paths
  }
  placement_flag_ = false;  // We do not need this parameter here
  if (placement_flag_ == false) {
    placement_dimensions_ = 0;  // no placement information
  } else {
    placement_dimensions_ = 2;  // 2D canvas
  }
  fence_flag_ = fence_flag_arg;
  fence_ = Rect(fence_lx_arg * dbu,
                fence_ly_arg * dbu,
                fence_ux_arg * dbu,
                fence_uy_arg * dbu);
  if (fence_flag_ == false || fence_.IsValid() == false) {
    fence_.Reset();
  }
  // local parameters
  std::string fixed_file = fixed_file_arg;
  std::string community_file = community_file_arg;
  std::string group_file = group_file_arg;
  std::string solution_file = solution_filename_arg;
  std::string hypergraph_file = hypergraph_file_arg;
  std::string hypergraph_int_weight_file = hypergraph_int_weight_file_arg;
  logger_->report("[PARAM] Number of partitions = {}", num_parts_);
  logger_->report("[PARAM] UBfactor = {}", ub_factor_);
  logger_->report("[PARAM] Seed = {}", seed_);
  logger_->report("[PARAM] Vertex dimensions = {}", vertex_dimensions_);
  logger_->report("[PARAM] Hyperedge dimensions = {}", hyperedge_dimensions_);
  logger_->report("[PARAM] Placement dimensions = {}", placement_dimensions_);
  logger_->report("[PARAM] Timing aware flag = {}", timing_aware_flag_);
  logger_->report("[PARAM] Global net threshold = {}", global_net_threshold_);
  logger_->report("[PARAM] Top {} critical timing paths are extracted.",
                  top_n_);
  logger_->report("[PARAM] Fence aware flag = {}", fence_flag_);
  if (fence_flag_ == true) {
    logger_->report(
        "[PARAM] fence_lx = {}, fence_ly = {}, fence_ux = {}, fence_uy = {}",
        fence_.lx / dbu,
        fence_.ly / dbu,
        fence_.ux / dbu,
        fence_.uy / dbu);
  }
  if (!fixed_file.empty()) {
    logger_->report("[PARAM] Fixed file  = {}", fixed_file);
  }
  if (!community_file.empty()) {
    logger_->report("[PARAM] Community file = {}", community_file);
  }
  if (!group_file.empty()) {
    logger_->report("[PARAM] Group file = {}", group_file);
  }
  if (!hypergraph_file.empty()) {
    logger_->report("[PARAM] Hypergraph file = {}", hypergraph_file);
  }
  if (!hypergraph_int_weight_file.empty()) {
    logger_->report("[PARAM] Hypergraph_int_weight_file = {}",
                    hypergraph_int_weight_file);
  }
  if (!solution_file.empty()) {
    logger_->report("[PARAM] Solution file = {}", solution_file);
  }

  // set the random seed
  srand(seed_);  // set the random seed

  // read the netlist from OpenDB
  // for IO port and insts (std cells and macros),
  // there is an attribute for vertex_id
  logger_->report("========================================");
  logger_->report("[STATUS] Reading netlist**** ");
  // if the fence_flag_ is true, only consider the instances within the fence
  ReadNetlist(fixed_file, community_file, group_file);
  logger_->report("[STATUS] Finish reading netlist****");

  // check the weighting scheme
  if (static_cast<int>(e_wt_factors_.size()) != hyperedge_dimensions_) {
    logger_->report(
        "[WARNING] no hyperedge weighting is specified. Use default value of "
        "1.");
    e_wt_factors_.clear();
    e_wt_factors_.resize(hyperedge_dimensions_);
    std::fill(e_wt_factors_.begin(), e_wt_factors_.end(), 1.0);
  }
  logger_->report("[PARAM] hyperedge weight factor : [ {} ]",
                  GetVectorString(e_wt_factors_));

  if (static_cast<int>(v_wt_factors_.size()) != vertex_dimensions_) {
    logger_->report(
        "[WARNING] no vertex weighting is specified. Use default value of 1.");
    v_wt_factors_.clear();
    v_wt_factors_.resize(vertex_dimensions_);
    std::fill(v_wt_factors_.begin(), v_wt_factors_.end(), 1.0);
  }
  logger_->report("[PARAM] vertex weight factor : [ {} ]",
                  GetVectorString(v_wt_factors_));

  if (static_cast<int>(placement_wt_factors_.size()) != placement_dimensions_) {
    if (placement_dimensions_ <= 0) {
      placement_wt_factors_.clear();
    } else {
      logger_->report(
          "[WARNING] no placement weighting is specified. Use default value of "
          "1.");
      placement_wt_factors_.clear();
      placement_wt_factors_.resize(placement_dimensions_);
      std::fill(
          placement_wt_factors_.begin(), placement_wt_factors_.end(), 1.0f);
    }
  }
  logger_->report("[PARAM] placement weight factor : [ {} ]",
                  GetVectorString(placement_wt_factors_));

  // print all the weighting parameters
  logger_->report("[PARAM] net_timing_factor : {}", net_timing_factor_);
  logger_->report("[PARAM] path_timing_factor : {}", path_timing_factor_);
  logger_->report("[PARAM] path_snaking_factor : {}", path_snaking_factor_);
  logger_->report("[PARAM] timing_exp_factor : {}", timing_exp_factor_);
  logger_->report("[PARAM] extra_delay : {}", extra_delay_);

  // create the evaluator class
  TP_evaluator_ptr evaluator
      = std::make_shared<GoldenEvaluator>(num_parts_,
                                          // weight vectors
                                          e_wt_factors_,
                                          v_wt_factors_,
                                          placement_wt_factors_,
                                          // timing related weight
                                          net_timing_factor_,
                                          path_timing_factor_,
                                          path_snaking_factor_,
                                          timing_exp_factor_,
                                          extra_delay_,
                                          original_hypergraph_,
                                          logger_);

  evaluator->InitializeTiming(original_hypergraph_);

  if (hypergraph_file.empty() == false) {
    evaluator->WriteWeightedHypergraph(original_hypergraph_, hypergraph_file);
    logger_->report("Finish writing hypergraph");
  }

  // This is for hMETIS. hMETIS only accept integer weight
  if (hypergraph_int_weight_file.empty() == false) {
    evaluator->WriteIntWeightHypergraph(original_hypergraph_,
                                        hypergraph_int_weight_file);
    logger_->report("Finish writing integer weight hypergraph");
  }

  if (solution_file.empty() == false) {
    int part_id = -1;
    std::ifstream solution_file_input(solution_file);
    if (!solution_file_input.is_open()) {
      logger_->error(
          PAR, 2514, "Can not open the solution file : {}", solution_file);
    }
    while (solution_file_input >> part_id) {
      solution_.push_back(part_id);
    }
    solution_file_input.close();
    evaluator->ConstraintAndCutEvaluator(
        original_hypergraph_, solution_, ub_factor_, group_attr_, true);

    // generate the timing report
    if (timing_aware_flag_ == true) {
      logger_->report("[STATUS] Displaying timing path cuts statistics");
      PathStats path_stats = evaluator->GetTimingCuts(original_hypergraph_, solution_);
      evaluator->PrintPathStats(path_stats);
    }

    logger_->report("===============================================");
    logger_->report("Exiting TritonPart");
    return;
  }
}

// k-way partitioning used by Hier-RTLMP
std::vector<int> TritonPart::PartitionKWaySimpleMode(
    unsigned int num_parts_arg,
    float balance_constraint_arg,
    unsigned int seed_arg,
    const std::vector<std::vector<int>>& hyperedges,
    const std::vector<float>& vertex_weights,
    const std::vector<float>& hyperedge_weights)
{
  num_parts_ = num_parts_arg;
  ub_factor_ = balance_constraint_arg;
  seed_ = seed_arg;
  vertex_dimensions_ = 1;  // for design partitioning, vertex weight is the area
                           // of the instance
  hyperedge_dimensions_
      = 1;  // for design partitioning, hyperedge weight is the connectivity
  timing_aware_flag_ = false;
  placement_flag_ = false;
  placement_dimensions_ = 0;
  fence_flag_ = false;
  hyperedges_ = hyperedges;
  fixed_attr_.clear();
  community_attr_.clear();
  group_attr_.clear();

  // convert vertex and hyperedge weights
  for (const auto& weight : vertex_weights) {
    std::vector<float> v_wt{weight};
    vertex_weights_.push_back(v_wt);
  }

  for (const auto& weight : hyperedge_weights) {
    std::vector<float> e_wt{weight};
    hyperedge_weights_.push_back(e_wt);
  }

  // Build the original hypergraph first
  original_hypergraph_ = std::make_shared<TPHypergraph>(vertex_dimensions_,
                                                        hyperedge_dimensions_,
                                                        placement_dimensions_,
                                                        hyperedges_,
                                                        vertex_weights_,
                                                        hyperedge_weights_,
                                                        fixed_attr_,
                                                        community_attr_,
                                                        placement_attr_,
                                                        logger_);

  // call the multilevel partitioner to partition hypergraph_
  // but the evaluation is the original_hypergraph_
  MultiLevelPartition();

  return solution_;
}

// --------------------------------------------------------------------------------------
// Private functions
// --------------------------------------------------------------------------------------

// for hypergraph partitioning
// Read hypergraph from input files and related constraint files
void TritonPart::ReadHypergraph(std::string hypergraph_file,
                                std::string fixed_file,
                                std::string community_file,
                                std::string group_file,
                                std::string placement_file)
{
  // read hypergraph file
  std::ifstream hypergraph_file_input(hypergraph_file);
  if (!hypergraph_file_input.is_open()) {
    logger_->error(PAR,
                   2500,
                   "Can not open the input hypergraph file : {}",
                   hypergraph_file);
  }
  // Check the number of vertices, number of hyperedges, weight flag
  std::string cur_line;
  std::getline(hypergraph_file_input, cur_line);
  std::istringstream cur_line_buf(cur_line);
  std::vector<int> stats{std::istream_iterator<int>(cur_line_buf),
                         std::istream_iterator<int>()};
  num_hyperedges_ = stats[0];
  num_vertices_ = stats[1];
  bool hyperedge_weight_flag = false;
  bool vertex_weight_flag = false;
  if (stats.size() == 3) {
    if ((stats[2] % 10) == 1) {
      hyperedge_weight_flag = true;
    }
    if (stats[2] >= 10) {
      vertex_weight_flag = true;
    }
  }

  // clear the related vectors
  hyperedges_.clear();
  hyperedge_weights_.clear();
  vertex_weights_.clear();
  hyperedges_.reserve(num_hyperedges_);
  hyperedge_weights_.reserve(num_hyperedges_);
  vertex_weights_.reserve(num_vertices_);

  // Read hyperedge information
  for (int i = 0; i < num_hyperedges_; i++) {
    std::getline(hypergraph_file_input, cur_line);
    if (hyperedge_weight_flag == true) {
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> hvec{std::istream_iterator<float>(cur_line_buf),
                              std::istream_iterator<float>()};
      std::vector<float>::iterator breakpoint{hvec.begin()
                                              + hyperedge_dimensions_};
      // read first hyperedge_dimensions_ elements as hyperege weights
      std::vector<float> hwts(hvec.begin(), breakpoint);
      // read remaining elements as hyperedge
      std::vector<int> hyperedge(breakpoint, hvec.end());
      for (auto& value : hyperedge) {
        value--;  // the vertex id starts from 1 in the hypergraph file
      }
      hyperedge_weights_.push_back(hwts);
      hyperedges_.push_back(hyperedge);
    } else {
      std::istringstream cur_line_buf(cur_line);
      std::vector<int> hyperedge{std::istream_iterator<int>(cur_line_buf),
                                 std::istream_iterator<int>()};
      for (auto& value : hyperedge) {
        value--;  // the vertex id starts from 1 in the hypergraph file
      }
      std::vector<float> hwts(hyperedge_dimensions_,
                              1.0);  // each dimension has the same weight
      hyperedge_weights_.push_back(hwts);
      hyperedges_.push_back(hyperedge);
    }
  }

  // Read weight for vertices
  for (int i = 0; i < num_vertices_; i++) {
    if (vertex_weight_flag == true) {
      std::getline(hypergraph_file_input, cur_line);
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> vwts{std::istream_iterator<float>(cur_line_buf),
                              std::istream_iterator<float>()};
      vertex_weights_.push_back(vwts);
    } else {
      std::vector<float> vwts(vertex_dimensions_, 1.0);
      vertex_weights_.push_back(vwts);
    }
  }

  // Read fixed vertices
  if (fixed_file.size() > 0) {
    int part_id = -1;
    std::ifstream fixed_file_input(fixed_file);
    if (!fixed_file_input.is_open()) {
      logger_->error(PAR, 2501, "Can not open the fixed file : {}", fixed_file);
    }
    while (fixed_file_input >> part_id) {
      fixed_attr_.push_back(part_id);
    }
    fixed_file_input.close();
    if (static_cast<int>(fixed_attr_.size()) != num_vertices_) {
      logger_->report("[WARNING] Reset the fixed attributes to NONE !");
      fixed_attr_.clear();
    }
  }

  // Read community file
  if (community_file.size() > 0) {
    int part_id = -1;
    std::ifstream community_file_input(community_file);
    if (!community_file_input.is_open()) {
      logger_->error(
          PAR, 2502, "Can not open the community file : {}", community_file);
    }
    while (community_file_input >> part_id) {
      community_attr_.push_back(part_id);
    }
    community_file_input.close();
    if (static_cast<int>(community_attr_.size()) != num_vertices_) {
      logger_->report("[WARNING] Reset the community attributes to NONE !");
      community_attr_.clear();
    }
  }

  // read group file
  if (group_file.size() > 0) {
    std::ifstream group_file_input(group_file);
    if (!group_file_input.is_open()) {
      logger_->error(PAR, 2503, "Can not open the group file : {}", group_file);
    }
    group_attr_.clear();
    std::string cur_line;
    while (std::getline(group_file_input, cur_line)) {
      std::istringstream cur_line_buf(cur_line);
      std::vector<int> group_info{std::istream_iterator<int>(cur_line_buf),
                                  std::istream_iterator<int>()};
      // reduce by 1 because definition of hMETIS
      for (auto& value : group_info) {
        value--;
      }
      if (group_info.size() > 1) {
        group_attr_.push_back(group_info);
      }
    }
    group_file_input.close();
  }

  // Read placement file
  if (placement_file.size() > 0) {
    std::ifstream placement_file_input(placement_file);
    if (!placement_file_input.is_open()) {
      logger_->error(
          PAR, 2504, "Can not open the placement file : {}", placement_file);
    }
    std::string cur_line;
    // We assume the embedding has been normalized
    const float max_placement_value
        = 1.0;  // we assume the embedding has been normalized,
                // so we assume max_value is 1.0
    std::vector<std::vector<float>> temp_placement_attr;
    while (std::getline(placement_file_input, cur_line)) {
      std::vector<std::string> elements = SplitLine(
          cur_line);  // split the line based on deliminator empty space, ','
      std::vector<float> vertex_placement;
      for (auto& ele : elements) {
        if (ele == "NaN" || ele == "nan" || ele == "NAN") {
          vertex_placement.push_back(max_placement_value);
        } else {
          vertex_placement.push_back(std::stof(ele));
        }
      }
      temp_placement_attr.push_back(vertex_placement);
    }
    placement_file_input.close();
    // Here comes the very important part for placement-driven clustering
    // Since we have so many vertices, the embedding value for each single
    // vertex usually very small, around e-5 - e-7 So we need to normalize the
    // placement embedding based on average distance again Here we randomly
    // sample num_vertices of pairs to compute the average norm
    float avg_distance = 0.0;
    int num_random_pair = 0;
    for (int i = 0; i < num_vertices_; i++) {
      std::mt19937 gen;
      std::uniform_real_distribution<> dist(0.0, 1.0);
      gen.seed(seed_);
      const int u = num_vertices_ * dist(gen);
      const int v = num_vertices_ * dist(gen);
      if (u != v) {
        avg_distance += norm2(temp_placement_attr[u] - temp_placement_attr[v]);
        num_random_pair++;
      }
    }
    avg_distance = avg_distance / num_random_pair;
    logger_->report(
        "[INFO] Normalize the placement embedding with the value of {} !",
        avg_distance);

    // perform normalization
    for (auto& ele : temp_placement_attr) {
      placement_attr_.push_back(DivideFactor(ele, avg_distance));
    }

    /*
    while (std::getline(placement_file_input, cur_line)) {
      std::istringstream cur_line_buf(cur_line);
      std::vector<float> vertex_placement{
          std::istream_iterator<float>(cur_line_buf),
          std::istream_iterator<float>()};
      if (static_cast<int>(vertex_placement.size()) == placement_dimensions_) {
        placement_attr_.push_back(vertex_placement);
      }
    }
    */
    if (static_cast<int>(placement_attr_.size()) != num_vertices_) {
      logger_->report("[WARNING] Reset the placement attributes to NONE !");
      placement_attr_.clear();
    }
  }

  // Build the original hypergraph first
  original_hypergraph_ = std::make_shared<TPHypergraph>(vertex_dimensions_,
                                                        hyperedge_dimensions_,
                                                        placement_dimensions_,
                                                        hyperedges_,
                                                        vertex_weights_,
                                                        hyperedge_weights_,
                                                        fixed_attr_,
                                                        community_attr_,
                                                        placement_attr_,
                                                        logger_);

  // show the status of hypergraph
  logger_->report("[INFO] Hypergraph Information**");
  logger_->report("[INFO] Vertices = {}", original_hypergraph_->num_vertices_);
  logger_->report("[INFO] Hyperedges = {}",
                  original_hypergraph_->num_hyperedges_);
}

// for design partitioning
// Convert the netlist into hypergraphs
// read fixed_file, community_file and group_file
// read placement information
// read timing information
void TritonPart::ReadNetlist(std::string fixed_file,
                             std::string community_file,
                             std::string group_file)
{
  // assign vertex_id property of each instance and each IO port
  // the vertex_id property will be removed after the partitioning
  vertex_weights_.clear();
  vertex_types_.clear();
  fixed_attr_.clear();
  community_attr_.clear();
  group_attr_.clear();
  placement_attr_.clear();
  // traverse all the instances
  int vertex_id = 0;
  // check if the fence constraint is specified
  if (fence_flag_ == true) {
    // check IO ports
    for (auto term : block_->getBTerms()) {
      // -1 means that the instance is not used by the partitioner
      odb::dbIntProperty::create(term, "vertex_id", -1);
      odb::Rect box = term->getBBox();
      if (box.xMin() >= fence_.lx && box.xMax() <= fence_.ux
          && box.yMin() >= fence_.ly && box.yMax() <= fence_.uy) {
        odb::dbIntProperty::create(term, "vertex_id", vertex_id++);
        std::vector<float> vwts(vertex_dimensions_,
                                0.0);  // IO port has no area
        vertex_weights_.emplace_back(vwts);
        vertex_types_.emplace_back(PORT);
        odb::dbIntProperty::find(term, "vertex_id")->setValue(vertex_id++);
        if (placement_flag_ == true) {
          std::vector<float> loc{(box.xMin() + box.xMax()) / 2.0,
                                 (box.yMin() + box.yMax()) / 2.0};
          placement_attr_.emplace_back(loc);
        }
      }
    }
    // check instances
    for (auto inst : block_->getInsts()) {
      // -1 means that the instance is not used by the partitioner
      odb::dbIntProperty::create(inst, "vertex_id", -1);
      const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
      if (liberty_cell == nullptr) {
        continue;  // ignore the instance with no liberty
      }
      odb::dbMaster* master = inst->getMaster();
      // check if the instance is a pad or a cover macro
      if (master->isPad() || master->isCover()) {
        continue;
      }
      odb::dbBox* box = inst->getBBox();
      // check if the inst is within the fence
      if (box->xMin() >= fence_.lx && box->xMax() <= fence_.ux
          && box->yMin() >= fence_.ly && box->yMax() <= fence_.uy) {
        const float area = liberty_cell->area();
        std::vector<float> vwts(vertex_dimensions_, area);
        vertex_weights_.emplace_back(vwts);
        if (master->isBlock()) {
          vertex_types_.emplace_back(MACRO);
        } else if (liberty_cell->hasSequentials()) {
          vertex_types_.emplace_back(SEQ_STD_CELL);
        } else {
          vertex_types_.emplace_back(COMB_STD_CELL);
        }
        if (placement_flag_ == true) {
          std::vector<float> loc{(box->xMin() + box->xMax()) / 2.0,
                                 (box->yMin() + box->yMax()) / 2.0};
          placement_attr_.emplace_back(loc);
        }
        odb::dbIntProperty::find(inst, "vertex_id")->setValue(vertex_id++);
      }
    }
  } else {
    for (auto term : block_->getBTerms()) {
      odb::dbIntProperty::create(term, "vertex_id", vertex_id++);
      vertex_types_.emplace_back(PORT);
      std::vector<float> vwts(vertex_dimensions_, 0.0);
      vertex_weights_.push_back(vwts);
      if (placement_flag_ == true) {
        odb::Rect box = term->getBBox();
        std::vector<float> loc{(box.xMin() + box.xMax()) / 2.0,
                               (box.yMin() + box.yMax()) / 2.0};
        placement_attr_.emplace_back(loc);
      }
    }

    for (auto inst : block_->getInsts()) {
      // -1 means that the instance is not used by the partitioner
      odb::dbIntProperty::create(inst, "vertex_id", -1);
      const sta::LibertyCell* liberty_cell = network_->libertyCell(inst);
      if (liberty_cell == nullptr) {
        continue;  // ignore the instance with no liberty
      }
      odb::dbMaster* master = inst->getMaster();
      // check if the instance is a pad or a cover macro
      if (master->isPad() || master->isCover()) {
        continue;
      }
      const float area = liberty_cell->area();
      std::vector<float> vwts(vertex_dimensions_, area);
      vertex_weights_.emplace_back(vwts);
      if (master->isBlock()) {
        vertex_types_.emplace_back(MACRO);
      } else if (liberty_cell->hasSequentials()) {
        vertex_types_.emplace_back(SEQ_STD_CELL);
      } else {
        vertex_types_.emplace_back(COMB_STD_CELL);
      }
      odb::dbIntProperty::find(inst, "vertex_id")->setValue(vertex_id++);
      if (placement_flag_ == true) {
        odb::dbBox* box = inst->getBBox();
        std::vector<float> loc{(box->xMin() + box->xMax()) / 2.0,
                               (box->yMin() + box->yMax()) / 2.0};
        placement_attr_.emplace_back(loc);
      }
    }
  }

  num_vertices_ = vertex_id;

  // read fixed instance file
  if (fixed_file.empty() == false) {
    std::ifstream file_input(fixed_file);
    if (!file_input.is_open()) {
      logger_->report("[WARNING] Cannot open the fixed instance file : {}",
                      fixed_file);
    } else {
      fixed_attr_.resize(num_vertices_);
      std::fill(fixed_attr_.begin(), fixed_attr_.end(), -1);
      std::string cur_line;
      while (std::getline(file_input, cur_line)) {
        std::stringstream ss(cur_line);
        std::string inst_name;
        int partition_id = -1;
        ss >> inst_name;
        ss >> partition_id;
        auto db_inst = block_->findInst(inst_name.c_str());
        const int vertex_id
            = odb::dbIntProperty::find(db_inst, "vertex_id")->getValue();
        if (vertex_id > -1) {
          fixed_attr_[vertex_id] = partition_id;
        }
      }
    }
    file_input.close();
  }

  // read community attribute file
  if (community_file.empty() == false) {
    std::ifstream file_input(community_file);
    if (!file_input.is_open()) {
      logger_->report("[WARNING] Cannot open the community file : {}",
                      community_file);
    } else {
      community_attr_.resize(num_vertices_);
      std::fill(community_attr_.begin(), community_attr_.end(), -1);
      std::string cur_line;
      while (std::getline(file_input, cur_line)) {
        std::stringstream ss(cur_line);
        std::string inst_name;
        int partition_id = -1;
        ss >> inst_name;
        ss >> partition_id;
        auto db_inst = block_->findInst(inst_name.c_str());
        const int vertex_id
            = odb::dbIntProperty::find(db_inst, "vertex_id")->getValue();
        if (vertex_id > -1) {
          community_attr_[vertex_id] = partition_id;
        }
      }
    }
    file_input.close();
  }

  // read the group file
  if (group_file.empty() == false) {
    std::ifstream file_input(group_file);
    if (!file_input.is_open()) {
      logger_->report("[WARNING] Cannot open the group file : {}", group_file);
    } else {
      group_attr_.clear();
      std::string cur_line;
      while (std::getline(file_input, cur_line)) {
        std::stringstream ss(cur_line);
        std::string inst_name;
        std::vector<int> inst_group;
        while (ss >> inst_name) {
          auto db_inst = block_->findInst(inst_name.c_str());
          const int vertex_id
              = odb::dbIntProperty::find(db_inst, "vertex_id")->getValue();
          if (vertex_id > -1) {
            inst_group.push_back(vertex_id);
          }
        }
        if (inst_group.size() > 1) {
          group_attr_.push_back(inst_group);
        }
      }
    }
    file_input.close();
  }

  // Check all the hyperedges,
  // we do not check the parallel hyperedges
  // because we need to consider timing graph
  hyperedges_.clear();
  hyperedge_weights_.clear();
  // Each net correponds to an hyperedge
  // Traverse the hyperedge and assign hyperedge_id to each net
  // the hyperedge_id property will be removed after partitioning
  int hyperedge_id = 0;
  for (auto net : block_->getNets()) {
    odb::dbIntProperty::create(net, "hyperedge_id", -1);
    // ignore all the power net
    if (net->getSigType().isSupply())
      continue;
    // check the hyperedge
    int driver_id = -1;      // vertex id of the driver instance
    std::set<int> loads_id;  // vertex id of sink instances
    // check the connected instances
    for (odb::dbITerm* iterm : net->getITerms()) {
      odb::dbInst* inst = iterm->getInst();
      const int vertex_id
          = odb::dbIntProperty::find(inst, "vertex_id")->getValue();
      if (vertex_id == -1) {
        continue;  // the current instance is not used
      }
      if (iterm->getIoType() == odb::dbIoType::OUTPUT) {
        driver_id = vertex_id;
      } else {
        loads_id.insert(vertex_id);
      }
    }
    // check the connected IO pins
    for (odb::dbBTerm* bterm : net->getBTerms()) {
      const int vertex_id
          = odb::dbIntProperty::find(bterm, "vertex_id")->getValue();
      if (vertex_id == -1) {
        continue;  // the current bterm is not used
      }
      if (bterm->getIoType() == odb::dbIoType::INPUT) {
        driver_id = vertex_id;
      } else {
        loads_id.insert(vertex_id);
      }
    }
    // check the hyperedges
    std::vector<int> hyperedge;
    if (driver_id != -1 && loads_id.size() > 0) {
      hyperedge.push_back(driver_id);
      for (auto& load_id : loads_id) {
        if (load_id != driver_id) {
          hyperedge.push_back(load_id);
        }
      }
    }
    // Ignore all the single-vertex hyperedge and large global netthreshold
    // if (hyperedge.size() > 1 && hyperedge.size() <= global_net_threshold_) {
    if (hyperedge.size() > 1) {
      hyperedges_.push_back(hyperedge);
      hyperedge_weights_.push_back(
          std::vector<float>(hyperedge_dimensions_, 1.0));
      odb::dbIntProperty::find(net, "hyperedge_id")->setValue(hyperedge_id++);
    }
  }  // finish hyperedge
  num_hyperedges_ = static_cast<int>(hyperedges_.size());

  // add timing features
  if (timing_aware_flag_ == true) {
    logger_->report("[STATUS] Extracting timing paths**** ");
    BuildTimingPaths();  // create timing paths
  }

  if (num_vertices_ == 0 || num_hyperedges_ == 0) {
    logger_->error(utl::PAR, 2677, "There is no vertices and hyperedges");
  }

  // build the timing graph
  // map each net to the timing arc in the timing graph
  std::vector<std::set<int>> hyperedges_arc_set;
  for (int e = 0; e < num_hyperedges_; e++) {
    const std::set<int> arc_set{e};
    hyperedges_arc_set.push_back(arc_set);
  }

  original_hypergraph_ = std::make_shared<TPHypergraph>(vertex_dimensions_,
                                                        hyperedge_dimensions_,
                                                        placement_dimensions_,
                                                        hyperedges_,
                                                        vertex_weights_,
                                                        hyperedge_weights_,
                                                        fixed_attr_,
                                                        community_attr_,
                                                        placement_attr_,
                                                        vertex_types_,
                                                        hyperedge_slacks_,
                                                        hyperedges_arc_set,
                                                        timing_paths_,
                                                        logger_);
  // show the status of hypergraph
  logger_->report("[INFO] Netlist Information**");
  logger_->report("[INFO] Vertices = {}", original_hypergraph_->num_vertices_);
  logger_->report("[INFO] Hyperedges = {}",
                  original_hypergraph_->num_hyperedges_);
  logger_->report("[INFO] Number of timing paths = {}", timing_paths_.size());
}

// Find all the critical timing paths
// The codes below similar to gui/src/staGui.cpp
// Please refer to sta/Search/ReportPath.cc for how to check the timing path
// Currently we can only consider single-clock design
// TODO:  how to handle multi-clock design
void TritonPart::BuildTimingPaths()
{
  if (timing_aware_flag_ == false || top_n_ <= 0) {
    logger_->report("[WARNING] Timing driven partitioning is disabled");
    return;
  }
  sta_->ensureGraph();     // Ensure that the timing graph has been built
  sta_->searchPreamble();  // Make graph and find delays
  sta_->ensureLevelized();
  // Step 1:  find the top_n critical timing paths
  sta::ExceptionFrom* e_from = nullptr;
  sta::ExceptionThruSeq* e_thrus = nullptr;
  sta::ExceptionTo* e_to = nullptr;
  bool include_unconstrained = false;
  bool get_max = true;  // max for setup check, min for hold check
  // Timing paths are grouped into path groups according to the clock
  // associated with the endpoint of the path, for example, path group for clk
  int group_count = top_n_;
  int endpoint_count = 1;  // The number of paths to report for each endpoint.
  // Definition for findPathEnds function in Search.hh
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
  sta::PathEndSeq path_ends = sta_->search()->findPathEnds(  // from, thrus, to,
                                                             // unconstrained
      e_from,   // return paths from a list of clocks/instances/ports/register
                // clock pins or latch data pins
      e_thrus,  // return paths through a list of instances/ports/nets
      e_to,     // return paths to a list of clocks/instances/ports or pins
      include_unconstrained,  // return unconstrained paths
      // corner, min_max,
      sta_->cmdCorner(),  // return paths for a process corner
      get_max ? sta::MinMaxAll::max()
              : sta::MinMaxAll::min(),  // return max/min paths checks
      // group_count, endpoint_count, unique_pins
      group_count,     // number of paths in total
      endpoint_count,  // number of paths for each endpoint
      true,
      -sta::INF,
      sta::INF,  // slack_min, slack_max,
      true,      // sort_by_slack
      nullptr,   // group_names
      // setup, hold, recovery, removal,
      get_max,
      !get_max,
      false,
      false,
      // clk_gating_setup, clk_gating_hold
      false,
      false);
  // check all the timing paths
  for (auto& path_end : path_ends) {
    // Printing timing paths to logger
    // sta_->reportPathEnd(path_end);
    auto* path = path_end->path();
    TimingPath timing_path;                     // create the timing path
    const float slack = path_end->slack(sta_);  // slack information
    // TODO: to be deleted.  We should not
    // normalize the slack according to the clock period for multi-clock design
    const float clock_period = path_end->targetClk(sta_)->period();
    maximum_clock_period_ = std::max(maximum_clock_period_, clock_period);
    timing_path.slack = slack;
    // logger_->report("clock_period = {}, slack = {}", maximum_clock_period_,
    // timing_path.slack);
    sta::PathExpanded expand(path, sta_);
    expand.path(expand.size() - 1);
    for (size_t i = 0; i < expand.size(); i++) {
      // PathRef is reference to a path vertex
      sta::PathRef* ref = expand.path(i);
      sta::Pin* pin = ref->vertex(sta_)->pin();
      // Nets connect pins at a level of the hierarchy
      auto net = network_->net(pin);  // sta::Net*
      // Check if the pin is connected to a net
      if (net == nullptr) {
        continue;  // check if the net exists
      }
      if (network_->isTopLevelPort(pin) == true) {
        auto bterm = block_->findBTerm(network_->pathName(pin));
        const int vertex_id
            = odb::dbIntProperty::find(bterm, "vertex_id")->getValue();
        if (vertex_id == -1) {
          continue;
        }
        if (timing_path.path.empty() == true
            || timing_path.path.back() != vertex_id) {
          timing_path.path.push_back(vertex_id);
        }
      } else {
        auto inst = network_->instance(pin);
        auto db_inst = block_->findInst(network_->pathName(inst));
        const int vertex_id
            = odb::dbIntProperty::find(db_inst, "vertex_id")->getValue();
        if (vertex_id == -1) {
          continue;
        }
        if (timing_path.path.empty() == true
            || timing_path.path.back() != vertex_id) {
          timing_path.path.push_back(vertex_id);
        }
      }
      auto db_net = block_->findNet(
          network_->pathName(net));  // convert sta::Net* to dbNet*
      const int hyperedge_id
          = odb::dbIntProperty::find(db_net, "hyperedge_id")->getValue();
      if (hyperedge_id == -1) {
        continue;
      }
      timing_path.arcs.push_back(hyperedge_id);
    }
    // add timing path
    if (timing_path.arcs.size() > 0) {
      timing_paths_.push_back(timing_path);
    }
  }

  // normalize all the slack
  logger_->report("[INFO] maximum_clock_period : {} second",
                  maximum_clock_period_);
  extra_delay_ = extra_delay_ / maximum_clock_period_;
  logger_->report("[INFO] normalized extra delay : {}", extra_delay_);
  for (auto& timing_path : timing_paths_) {
    timing_path.slack = timing_path.slack / maximum_clock_period_;
  }
  logger_->report(
      "[INFO] We normalized the slack of each path based on maximum clock "
      "period");
  // resize the hyperedge_slacks_
  hyperedge_slacks_.clear();
  hyperedge_slacks_.resize(num_hyperedges_);
  std::fill(hyperedge_slacks_.begin(),
            hyperedge_slacks_.end(),
            maximum_clock_period_);
  logger_->report(
      "[INFO] We normalized the slack of each net based on maximum clock "
      "period");
  int num_unconstrained_hyperedges = 0;
  // check the slack on each net
  for (auto db_net : block_->getNets()) {
    const int hyperedge_id
        = odb::dbIntProperty::find(db_net, "hyperedge_id")->getValue();
    if (hyperedge_id == -1) {
      continue;  // this net is not used
    }
    sta::Net* net = network_->dbToSta(db_net);
    const float slack = sta_->netSlack(net, sta::MinMax::max());
    // set the slack of unconstrained net to max_clock_period_
    if (slack > maximum_clock_period_) {
      num_unconstrained_hyperedges++;
      hyperedge_slacks_[hyperedge_id] = 1.0;
    } else {
      hyperedge_slacks_[hyperedge_id] = slack / maximum_clock_period_;
    }
  }
  logger_->report("[STATUS] Finish traversing timing graph");
  logger_->report("[WARNING] {} unconstrained hyperedges !",
                  num_unconstrained_hyperedges);
  logger_->report(
      "[WARNING] Reset the slack of all unconstrained hyperedges to {} seconds",
      maximum_clock_period_);
}

// Partition the hypergraph_ with the multilevel methodology
// the return value is the partitioning solution
void TritonPart::MultiLevelPartition()
{
  auto start_time_stamp_global = std::chrono::high_resolution_clock::now();

  // check the weighting scheme
  if (static_cast<int>(e_wt_factors_.size()) != hyperedge_dimensions_) {
    logger_->report(
        "[WARNING] no hyperedge weighting is specified. Use default value of "
        "1.");
    e_wt_factors_.clear();
    e_wt_factors_.resize(hyperedge_dimensions_);
    std::fill(e_wt_factors_.begin(), e_wt_factors_.end(), 1.0);
  }
  logger_->report("[PARAM] hyperedge weight factor : [ {} ]",
                  GetVectorString(e_wt_factors_));

  if (static_cast<int>(v_wt_factors_.size()) != vertex_dimensions_) {
    logger_->report(
        "[WARNING] no vertex weighting is specified. Use default value of 1.");
    v_wt_factors_.clear();
    v_wt_factors_.resize(vertex_dimensions_);
    std::fill(v_wt_factors_.begin(), v_wt_factors_.end(), 1.0);
  }
  logger_->report("[PARAM] vertex weight factor : [ {} ]",
                  GetVectorString(v_wt_factors_));

  if (static_cast<int>(placement_wt_factors_.size()) != placement_dimensions_) {
    if (placement_dimensions_ <= 0) {
      placement_wt_factors_.clear();
    } else {
      logger_->report(
          "[WARNING] no placement weighting is specified. Use default value of "
          "1.");
      placement_wt_factors_.clear();
      placement_wt_factors_.resize(placement_dimensions_);
      std::fill(
          placement_wt_factors_.begin(), placement_wt_factors_.end(), 1.0f);
    }
  }
  logger_->report("[PARAM] placement weight factor : [ {} ]",
                  GetVectorString(placement_wt_factors_));
  // print all the weighting parameters
  logger_->report("[PARAM] net_timing_factor : {}", net_timing_factor_);
  logger_->report("[PARAM] path_timing_factor : {}", path_timing_factor_);
  logger_->report("[PARAM] path_snaking_factor : {}", path_snaking_factor_);
  logger_->report("[PARAM] timing_exp_factor : {}", timing_exp_factor_);
  // coarsening related parameters
  logger_->report("[PARAM] coarsen order : {}", ToString(coarsen_order_));
  logger_->report("[PARAM] thr_coarsen_hyperedge_size_skip : {}",
                  thr_coarsen_hyperedge_size_skip_);
  logger_->report("[PARAM] thr_coarsen_vertices : {}", thr_coarsen_vertices_);
  logger_->report("[PARAM] thr_coarsen_hyperedges : {}",
                  thr_coarsen_hyperedges_);
  logger_->report("[PARAM] coarsening_ratio : {}", coarsening_ratio_);
  logger_->report("[PARAM] max_coarsen_iters : {}", max_coarsen_iters_);
  logger_->report("[PARAM] adj_diff_ratio : {}", adj_diff_ratio_);
  logger_->report("[PARAM] min_num_vertcies_each_part : {}",
                  min_num_vertices_each_part_);
  // initial partitioning parameter
  logger_->report("[PARAM] num_initial_solutions : {}", num_initial_solutions_);
  logger_->report("[PARAM] num_best_initial_solutions : {}",
                  num_best_initial_solutions_);
  // refinement related parameters
  logger_->report("[PARAM] refine_iters : {}", refiner_iters_);
  logger_->report("[PARAM] max_moves (FM or greedy refinement) : {}",
                  max_moves_);
  logger_->report("[PARAM] early_stop_ratio : {}", early_stop_ratio_);
  logger_->report("[PARAM] total_corking_passes : {}", total_corking_passes_);
  logger_->report("[PARAM] v_cycle_flag : {}", v_cycle_flag_);
  logger_->report("[PARAM] max_num_vcycle : {}", max_num_vcycle_);
  logger_->report("[PARAM] num_coarsen_solutions : {}", num_coarsen_solutions_);
  logger_->report("[PARAM] num_vertices_threshold_ilp : {}",
                  num_vertices_threshold_ilp_);

  // create the evaluator class
  TP_evaluator_ptr tritonpart_evaluator
      = std::make_shared<GoldenEvaluator>(num_parts_,
                                          // weight vectors
                                          e_wt_factors_,
                                          v_wt_factors_,
                                          placement_wt_factors_,
                                          // timing related weight
                                          net_timing_factor_,
                                          path_timing_factor_,
                                          path_snaking_factor_,
                                          timing_exp_factor_,
                                          extra_delay_,
                                          original_hypergraph_,
                                          logger_);

  // create the balance constraint
  MATRIX<float> upper_block_balance
      = original_hypergraph_->GetUpperVertexBalance(num_parts_, ub_factor_);
  MATRIX<float> lower_block_balance
      = original_hypergraph_->GetLowerVertexBalance(num_parts_, ub_factor_);

  // Step 1 : create all the coarsening, partitionig and refinement class
  const std::vector<float> thr_cluster_weight
      = DivideFactor(original_hypergraph_->GetTotalVertexWeights(),
                     min_num_vertices_each_part_ * num_parts_);

  // create the coarsener cluster
  TP_coarsening_ptr tritonpart_coarsener
      = std::make_shared<TPcoarsener>(num_parts_,
                                      thr_coarsen_hyperedge_size_skip_,
                                      thr_coarsen_vertices_,
                                      thr_coarsen_hyperedges_,
                                      coarsening_ratio_,
                                      max_coarsen_iters_,
                                      adj_diff_ratio_,
                                      thr_cluster_weight,
                                      seed_,
                                      coarsen_order_,
                                      tritonpart_evaluator,
                                      logger_);

  // create the initial partitioning class
  TP_partitioning_ptr tritonpart_partitioner = std::make_shared<TPpartitioner>(
      num_parts_, seed_, tritonpart_evaluator, logger_);

  // create the refinement classes
  // We have four types of refiner
  // (1) greedy refinement. try to one entire hyperedge each time
  TP_greedy_refiner_ptr greedy_refiner
      = std::make_shared<TPgreedyRefine>(num_parts_,
                                         refiner_iters_,
                                         path_timing_factor_,
                                         path_snaking_factor_,
                                         max_moves_,
                                         tritonpart_evaluator,
                                         logger_);

  // (2) ILP-based partitioning (only for two-way since k-way ILP partitioning
  // is too timing-consuming)
  TP_ilp_refiner_ptr ilp_refiner
      = std::make_shared<TPilpRefine>(num_parts_,
                                      refiner_iters_,
                                      path_timing_factor_,
                                      path_snaking_factor_,
                                      max_moves_,
                                      tritonpart_evaluator,
                                      logger_);

  // (3) direct k-way FM
  TP_k_way_fm_refiner_ptr k_way_fm_refiner
      = std::make_shared<TPkWayFMRefine>(num_parts_,
                                         refiner_iters_,
                                         path_timing_factor_,
                                         path_snaking_factor_,
                                         max_moves_,
                                         total_corking_passes_,
                                         tritonpart_evaluator,
                                         logger_);

  // (4) k-way pair-wise FM
  TP_k_way_pm_refiner_ptr k_way_pm_refiner
      = std::make_shared<TPkWayPMRefine>(num_parts_,
                                         refiner_iters_,
                                         path_timing_factor_,
                                         path_snaking_factor_,
                                         max_moves_,
                                         total_corking_passes_,
                                         tritonpart_evaluator,
                                         logger_);

  // create the multi-level class
  TP_multi_level_partitioner tritonpart_mlevel_partitioner
      = std::make_shared<TPmultilevelPartitioner>(num_parts_,
                                                  ub_factor_,
                                                  v_cycle_flag_,
                                                  num_initial_solutions_,
                                                  num_best_initial_solutions_,
                                                  num_vertices_threshold_ilp_,
                                                  max_num_vcycle_,
                                                  num_coarsen_solutions_,
                                                  seed_,
                                                  timing_aware_flag_,
                                                  tritonpart_coarsener,
                                                  tritonpart_partitioner,
                                                  k_way_fm_refiner,
                                                  k_way_pm_refiner,
                                                  greedy_refiner,
                                                  ilp_refiner,
                                                  tritonpart_evaluator,
                                                  logger_);

  if (timing_aware_flag_ == true) {
    // Initialize the timing on original_hypergraph_
    tritonpart_evaluator->InitializeTiming(original_hypergraph_);
  }
  // Use coarsening to do preprocessing step
  // Build the hypergraph used to call multi-level partitioner
  // (1) remove single-vertex hyperedge
  // (2) remove lager hyperedge
  // (3) detect parallel hyperedges
  // (4) handle group information
  // (5) group fixed vertices based on each block
  // group vertices based on group_attr_
  // the original_hypergraph_ will be modified by the Group Vertices command
  // We will store the mapping relationship of vertices between
  // original_hypergraph_ and hypergraph_
  tritonpart_coarsener->SetThrCoarsenHyperedgeSizeSkip(global_net_threshold_);
  hypergraph_
      = tritonpart_coarsener->GroupVertices(original_hypergraph_, group_attr_);
  tritonpart_coarsener->SetThrCoarsenHyperedgeSizeSkip(
      thr_coarsen_hyperedge_size_skip_);

  // partition on the processed hypergraph
  std::vector<int> solution = tritonpart_mlevel_partitioner->Partition(
      hypergraph_, upper_block_balance, lower_block_balance);

  // Translate the solution of hypergraph to original_hypergraph_
  // solution to solution_
  solution_.clear();
  solution_.resize(original_hypergraph_->num_vertices_);
  std::fill(solution_.begin(), solution_.end(), -1);
  for (int cluster_id = 0; cluster_id < hypergraph_->num_vertices_;
       cluster_id++) {
    const int part_id = solution[cluster_id];
    for (const auto& v : hypergraph_->vertex_c_attr_[cluster_id]) {
      solution_[v] = part_id;
    }
  }

  // Perform the last-minute refinement
  tritonpart_coarsener->SetThrCoarsenHyperedgeSizeSkip(global_net_threshold_);
  tritonpart_mlevel_partitioner->VcycleRefinement(
      hypergraph_, upper_block_balance, lower_block_balance, solution_);

  // evaluate on the original hypergraph
  // tritonpart_evaluator->CutEvaluator(original_hypergraph_, solution_, true);
  tritonpart_evaluator->ConstraintAndCutEvaluator(
      original_hypergraph_, solution_, ub_factor_, group_attr_, true);

  // generate the timing report
  if (timing_aware_flag_ == true) {
    logger_->report("[STATUS] Displaying timing path cuts statistics");
    PathStats path_stats = tritonpart_evaluator->GetTimingCuts(original_hypergraph_, solution_);
    tritonpart_evaluator->PrintPathStats(path_stats);
  }
  
  // print the runtime
  auto end_timestamp_global = std::chrono::high_resolution_clock::now();
  double total_global_time
      = std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_timestamp_global - start_time_stamp_global)
            .count();
  total_global_time *= 1e-9;
  logger_->report("[INFO] The runtime of multi-level partitioner : {} seconds",
                  total_global_time);
}

}  // namespace par
