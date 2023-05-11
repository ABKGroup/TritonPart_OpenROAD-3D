sta::define_cmd_args "rtl_blob_placer" { -tolerance tolerance \
                                          -max_num_level max_num_level \
                                          -coarsening_ratio coarsening_ratio \
                                          -num_bundled_ios num_bundled_ios \
                                          -large_net_threshould large_net_threshould \
                                          -signature_net_threshold signature_net_threshold \
                                          -intended_num_clusters intended_num_clusters \
                                          -num_tiers num_tiers \
                                          -tier_resolution tier_resolution \
                                          -wire_length_weight wire_length_weight \
                                          -vertical_weight vertical_weight \
                                          -area_upper_bound area_upper_bound \
                                          -swap_prob swap_prob \
                                          -init_prob init_prob \
                                          -max_steps max_steps \
                                          -num_perturb_per_step num_perturb_per_step \
                                          -timing_penalty timing_penalty \
                                          -seed seed \
                                        }
proc rtl_blob_placer { args } {
    sta::parse_key_args "rtl_blob_placer" args keys {
        -tolerance -max_num_level \
        -coarsening_ratio -num_bundled_ios -large_net_threshould \
        -signature_net_threshold -intended_num_clusters \
        -num_tiers -tier_resolution -wire_length_weight -vertical_weight \
        -area_upper_bound -swap_prob -init_prob -max_steps -num_perturb_per_step \
        -timing_penalty \
        -seed \
    } flag { }

# Set the default parameters for the blob_placer
# Set auto defaults for min/max std cells based on design
    set tolerance 0.1
    set max_num_level 2
    set coarsening_ratio 10.0
    set num_bundled_ios 4
    set large_net_threshould 50
    set signature_net_threshold 50
    set intended_num_clusters 300
    set num_tiers 2
    set tier_resolution 4
    set wire_length_weight 1.0
    set vertical_weight 1.0
    set area_upper_bound 0.8
    set swap_prob 0.5
    set init_prob 0.8
    set max_steps 1000000
    set num_perturb_per_step 10000
    set timing_penalty 10.0
    set seed 10

    if { [info exists keys(-tolerance)] } {
        set tolerance $keys(-tolerance)
    }
    if { [info exists keys(-max_num_level)] } {
        set max_num_level $keys(-max_num_level)
    }
    if { [info exists keys(-coarsening_ratio)] } {
        set coarsening_ratio $keys(-coarsening_ratio)
    }
    if { [info exists keys(-num_bundled_ios)] } {
        set num_bundled_ios $keys(-num_bundled_ios)
    }
    if { [info exists keys(-large_net_threshould)] } {
        set large_net_threshould $keys(-large_net_threshould)
    }
    if { [info exists keys(-signature_net_threshold)] } {
        set signature_net_threshold $keys(-signature_net_threshold)
    }
    if { [info exists keys(-intended_num_clusters)] } {
        set intended_num_clusters $keys(-intended_num_clusters)
    }
    if { [info exists keys(-num_tiers)] } {
        set num_tiers $keys(-num_tiers)
    }
    if { [info exists keys(-tier_resolution)] } {
        set tier_resolution $keys(-tier_resolution)
    }
    if { [info exists keys(-wire_length_weight)] } {
        set wire_length_weight $keys(-wire_length_weight)
    }
    if { [info exists keys(-vertical_weight)] } {
        set vertical_weight $keys(-vertical_weight)
    }
    if { [info exists keys(-area_upper_bound)] } {
        set area_upper_bound $keys(-area_upper_bound)
    }
    if { [info exists keys(-swap_prob)] } {
        set swap_prob $keys(-swap_prob)
    }
    if { [info exists keys(-init_prob)] } {
        set init_prob $keys(-init_prob)
    }
    if { [info exists keys(-max_steps)] } {
        set max_steps $keys(-max_steps)
    }
    if { [info exists keys(-num_perturb_per_step)] } {
        set num_perturb_per_step $keys(-num_perturb_per_step)
    }
    if { [info exists keys(-timing_penalty)] } {
        set timing_penalty $keys(-timing_penalty)
    }
    if { [info exists keys(-seed)] } {
        set seed $keys(-seed)
    }

    # TODO: Add the report directory creation part
    if {![mlsa::rtl_blob_placer_cmd $tolerance \
                                    $max_num_level \
                                    $coarsening_ratio \
                                    $num_bundled_ios \
                                    $large_net_threshould \
                                    $signature_net_threshold \
                                    $intended_num_clusters \
                                    $num_tiers \
                                    $tier_resolution \
                                    $wire_length_weight \
                                    $vertical_weight \
                                    $area_upper_bound \
                                    $swap_prob \
                                    $init_prob \
                                    $max_steps \
                                    $num_perturb_per_step \
                                    $timing_penalty \
                                    $seed \
                                    ]} {
        return false
    }

    return true
}