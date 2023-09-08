using Revise

using Oiler

# options
geol_slip_rate_weight = 2.0
save_results = false

# load data

# EARS
ear_block_file = "../block_data/ears_blocks.geojson"
ear_fault_file = "../block_data/ears_faults.geojson"  
gsrm_af_file = "../strain_data/gsrm_gps_af.geojson"

# C Asia
cea_block_file = "../../c_asia_blocks/block_data/c_asia_blocks.geojson"
cea_fault_file = "../../c_asia_blocks/block_data/c_asia_faults.geojson"
cea_slip_rate_file = "../../c_asia_blocks/block_data/c_asia_geol_slip_rates.geojson"
cea_tris_file = "../block_data/c_asia_sub_tris.geojson" # not using yet

# global
glo_block_file = "../../global_scale_plates/global_scale_plates.geojson"
glo_fault_file = "../../global_scale_plates/global_scale_faults.geojson"
glo_slip_rates_file = "../../global_scale_plates/global_scale_slip_rates.geojson"

# Anatolia
ana_block_file = "../../anatolia_blocks/block_data/anatolia_blocks.geojson"
ana_fault_file = "../../anatolia_blocks/block_data/anatolia_faults.geojson"
weiss_vel_field_file = "../../anatolia_blocks/geod_data/weiss_et_al_2020_vels_down_100.geojson" # not using


# boundary
ears_boundary_file = "../block_data/ears_bounds.geojson"

@info "joining blocks"
ear_block_df = Oiler.IO.gis_vec_file_to_df(ear_block_file)
glo_block_df = Oiler.IO.gis_vec_file_to_df(glo_block_file; fid_drop=["ant"])
cea_block_df = Oiler.IO.gis_vec_file_to_df(cea_block_file)
ana_block_df = Oiler.IO.gis_vec_file_to_df(ana_block_file)
block_df = vcat(ear_block_df,
                glo_block_df,
                cea_block_df,
                ana_block_df,
                cols=:union)

@info "culling blocks"
println("n blocks before ", size(block_df, 1))
bound_df = Oiler.IO.gis_vec_file_to_df(ears_boundary_file)
block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df)
println("n blocks after ", size(block_df, 1))

@info "doing faults"
fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        ear_fault_file,
                                                        glo_fault_file,
                                                        cea_fault_file,
                                                        ana_fault_file;
                                                        block_df=block_df,
                                                        subset_in_bounds=true,
                                                        check_blocks=true,
                                                        )
fault_df[:,:fid] = string.(fault_df[:,:fid])
println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))

@info "doing non-fault block boundaries"
@time non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bound_vels = vcat(map(b->Oiler.Boundaries.boundary_to_vels(b, ee=2.0, en=2.0), 
                      non_fault_bounds)...)
println("n non-fault-bound vels: ", length(bound_vels))

@info "doing GNSS"
gsrm_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_af_file)

@info "doing GSRM vels"
@time gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

gnss_vels = gsrm_vels

println("n gnss vels: ", length(gnss_vels))

vels = vcat(fault_vels,
            bound_vels,
            gnss_vels, 
            )

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                                weighted=true,
                                                elastic_floor=1e-4,
                                                check_nans=true,
                                                predict_vels=true,
                                                check_closures=true,
                                                pred_se=true,
                                                )
                                                
map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults)

Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 ref_pole="1111", directory="../web_viewer")

show()
