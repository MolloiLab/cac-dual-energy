### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ d6965c35-ace0-4ff6-92f8-435853937946
using Pkg; Pkg.instantiate()

# ╔═╡ 024de588-307b-452d-812d-3a58b9265089
using DrWatson

# ╔═╡ 3c5e1804-23dd-4629-baa5-5106aacd74a6
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ 0c4ba5cc-8e40-40c6-b596-785973bbc74e
using PlutoUI, CairoMakie, Statistics, CSV, DataFrames, DICOM, CSVFiles, Printf

# ╔═╡ a9c366d0-e90a-406c-b970-3d5a76bc7973
using StatsBase: quantile!, rmsd

# ╔═╡ 11835515-0d70-4a4b-b7f0-742153e4c8e2
using CalciumScoring

# ╔═╡ 65212870-6665-4ce8-b1f6-f80f9304bf04
using MaterialDecomposition

# ╔═╡ d7e045fc-ac9c-4f53-836e-aedcadad1fbe
using GLM, MLJBase

# ╔═╡ f1708909-5df5-477b-9657-82029aab0a50
include(srcdir("masks.jl")); include(srcdir("dicom_utils.jl"));

# ╔═╡ 95509262-a9a1-41b6-972d-e7dcf7bb2371
include(srcdir("helper_functions.jl"));

# ╔═╡ b98db832-06bd-492f-8ca4-38e298086c8a
include(srcdir("plot_utils.jl"));

# ╔═╡ 937d6c6d-80c7-4bfc-9c61-9d1830c84b66
TableOfContents()

# ╔═╡ c69df360-2e83-4d92-89a5-beaa1a20a1c5
md"""
# Low Density
"""

# ╔═╡ 2a8133c5-d7c0-4f85-8091-9a6c9d2aceb4
begin
	sizes_cal = [30]
	energies = [80, 120, 135]
end;

# ╔═╡ 0a82a5ee-d54d-41fe-9cd3-336d8d9ca8e9
md"""
## Volume Fraction & Agatston Calibration
"""

# ╔═╡ ebec3a75-49ff-4ab0-b738-59cac17a001f
path = joinpath(datadir("dcms", "cal", "100", "30"), string(energies[2]))

# ╔═╡ 3f9c0eb2-b727-4b9a-91eb-f9c482d6809f
dcm_viz = load_dcm_array(dcmdir_parse(path));

# ╔═╡ d3e09d9a-187f-4db0-a8c5-5fc0dc408d25
heatmap(dcm_viz[:, :, 2], colormap = :grays)

# ╔═╡ f1575684-58ac-4606-92ef-993377c8309d
begin
	dcm_cal = dcmdir_parse(path)
	dcm_array_cal = load_dcm_array(dcm_cal)
	offset = 5

	center_insert1, center_insert2 = 187, 318
	calibration_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array_cal, 3))
	
	for z in axes(dcm_array_cal, 3)
		rows, cols, depth = size(dcm_array_cal)
		half_row, half_col = center_insert1, center_insert2
		row_range = half_row-offset:half_row+offset
		col_range = half_col-offset:half_col+offset	
		calibration_rod[:, :, z] .= dcm_array_cal[row_range, col_range, z];
	end

	

	hu_calcium_100 = mean(calibration_rod)
	ρ_calcium_100 = 0.100 # mg/cm^3
end

# ╔═╡ 82712845-c2ba-49f1-a8ad-9127cbc321b2
md"""
## Calibration
"""

# ╔═╡ a1289bbb-21b8-471e-9116-86f04f53c085
#second best option i found:
#densities_cal_all = [3, 6,10,50,100,200,300,400,500,600]

#Worse option:
#densities_cal_all = [3,25,50,75,100,200,300,400,500,600]

#Original points work well now.
densities_cal_all = [0, 50, 100, 200, 300, 400, 500, 600]

# ╔═╡ a5935d6a-c825-4f8e-b559-c22104bb151e
optimal_cal_low = densities_cal_all

# ╔═╡ ea5853a5-5cd8-4017-824c-c1a9c945b357
begin
	densities_cal_low = optimal_cal_low
	
	# Low Energy
	means_80_low = Dict(:density => densities_cal_low, :means => zeros(length(densities_cal_low)))
	for (i, density) in enumerate(densities_cal_low)
		
		
		if i ==1 #for background
			path = joinpath(datadir("dcms", "cal", "100", "30"), string(energies[1]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
			offset = 5
			bkg_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array, 3))
			center_insert_bkg, center_insert_bkg2 = 250, 318

			for z in axes(dcm_array, 3)
					rows, cols, depth = size(dcm_array)
					half_row, half_col = center_insert_bkg, center_insert_bkg2
					row_range = half_row-offset:half_row+offset
					col_range = half_col-offset:half_col+offset	
					calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
			
			means_80_low[:means][i] = mean(calibration_rod)

		else
			for _size in sizes_cal
				path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[1]))
				dcm = dcmdir_parse(path)
				dcm_array = load_dcm_array(dcm)
			
				center_insert1, center_insert2 = 187, 318
				offset = 5
				calibration_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array, 3))
				
				for z in axes(dcm_array, 3)
					rows, cols, depth = size(dcm_array)
					half_row, half_col = center_insert1, center_insert2
					row_range = half_row-offset:half_row+offset
					col_range = half_col-offset:half_col+offset	
					calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
				end
			
				means_80_low[:means][i] = mean(calibration_rod)
			end
		end

			
			
	end

	# High Energy
	means_135_low = Dict(:density => densities_cal_low, :means => zeros(length(densities_cal_low)))
	for (i, density) in enumerate(densities_cal_low)


		if i ==1 #for background
			path = joinpath(datadir("dcms", "cal", "50", "30"), string(energies[3]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
			offset = 5
			bkg_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array, 3))
			center_insert_bkg, center_insert_bkg2 = 250, 318

			for z in axes(dcm_array, 3)
					rows, cols, depth = size(dcm_array)
					half_row, half_col = center_insert_bkg, center_insert_bkg2
					row_range = half_row-offset:half_row+offset
					col_range = half_col-offset:half_col+offset	
					calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
			
				means_135_low[:means][i] = mean(calibration_rod)

		else 
			for _size in sizes_cal
				path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[3]))
				dcm = dcmdir_parse(path)
				dcm_array = load_dcm_array(dcm)
			
				center_insert1, center_insert2 = 187, 318
				
				offset = 5
				calibration_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array, 3))
				
				for z in axes(dcm_array, 3)
					rows, cols, depth = size(dcm_array)
					half_row, half_col = center_insert1, center_insert2
					row_range = half_row-offset:half_row+offset
					col_range = half_col-offset:half_col+offset	
					calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
				end
			
				means_135_low[:means][i] = mean(calibration_rod)
			end
		end 
	end

	# Fit Parameters
	calculated_intensities_low = hcat(means_80_low[:means], means_135_low[:means])
	ps_low = fit_calibration(calculated_intensities_low, densities_cal_low)
end

# ╔═╡ 66d07988-d98a-4bab-8243-33d6f4485cf0
md"""
### Check Results
"""

# ╔═╡ c93775d8-3422-4ad1-9f40-6d57477d8425
begin
	predicted_densities_low = []
	
	for i in 1:length(densities_cal_low)
		append!(
			predicted_densities_low, 
			quantify(calculated_intensities_low[i, 1], calculated_intensities_low[i, 2], ps_low
			)
		)
	end
end

# ╔═╡ 636abeda-6c10-4d21-ae08-10e1274b77a1
calculated_intensities_low

# ╔═╡ d0b1613e-890e-417d-b7eb-57b0bc9fafe2
df_low = DataFrame(
	densities = densities_cal_low,
	predicted_densities = predicted_densities_low
)

# ╔═╡ 93679b54-9364-4579-8767-53e36499b6d8
md"""
## Validation
"""

# ╔═╡ 853d895d-91d4-487a-a927-e74424bba9ed
densities_val_all = [
	"15_18_22"
	"26_29_36"
	"52_59_73"
	"110_210_310"
	"410_610_780"
]

# ╔═╡ d6a65460-51e0-4b18-83e5-1844f2566cc6
begin
	densities_val_low = densities_val_all[1:3]
	sizes_val = ["small", "medium", "large"]
end;

# ╔═╡ e45adb9a-0519-4063-b8ae-2cb7e1005397
energies_val = [80, 120, 135]

# ╔═╡ 8a774374-ecce-40d7-8b64-609124242d1f
md"""
### Load masks
"""

# ╔═╡ 38bf3544-c7a4-46d6-9df2-6fe269b24cb8
begin
	# Load masks
	masks_large = Dict(
		:mask_L_HD => Array(CSV.read(datadir("julia_arrays", "large", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_M_HD => Array(CSV.read(datadir("julia_arrays", "large", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_S_HD => Array(CSV.read(datadir("julia_arrays", "large", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_L_MD => Array(CSV.read(datadir("julia_arrays", "large", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_M_MD => Array(CSV.read(datadir("julia_arrays", "large", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_S_MD => Array(CSV.read(datadir("julia_arrays", "large", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_M_LD => Array(CSV.read(datadir("julia_arrays", "large", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_L_LD => Array(CSV.read(datadir("julia_arrays", "large", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_S_LD => Array(CSV.read(datadir("julia_arrays", "large", "mask_S_LD.csv"), DataFrame; header=false)),
	)
	
	masks_medium = Dict(
		:mask_L_HD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_M_HD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_S_HD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_L_MD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_M_MD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_S_MD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_M_LD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_L_LD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_S_LD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_S_LD.csv"), DataFrame; header=false)),
	)

	masks_small = Dict(
		:mask_L_HD => Array(CSV.read(datadir("julia_arrays", "small", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_M_HD => Array(CSV.read(datadir("julia_arrays", "small", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_S_HD => Array(CSV.read(datadir("julia_arrays", "small", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_L_MD => Array(CSV.read(datadir("julia_arrays", "small", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_M_MD => Array(CSV.read(datadir("julia_arrays", "small", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_S_MD => Array(CSV.read(datadir("julia_arrays", "small", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_M_LD => Array(CSV.read(datadir("julia_arrays", "small", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_L_LD => Array(CSV.read(datadir("julia_arrays", "small", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_S_LD => Array(CSV.read(datadir("julia_arrays", "small", "mask_S_LD.csv"), DataFrame; header=false)),
	)
end;

# ╔═╡ 5588545b-f070-4e02-ac12-3e1b5a150d88
md"""
### Predict
"""

# ╔═╡ 20a32f24-2c4b-4308-a6e1-802a72322fa7
begin
	dfs_low_md = []
	dfs_low_vf = []
	dfs_low_a = []
	for density in densities_val_low
		for _size in sizes_val
			if _size == "small"
				center_insert = [175, 320]
				mask_L_HD = masks_small[:mask_L_HD]
				mask_M_HD = masks_small[:mask_M_HD]
				mask_S_HD = masks_small[:mask_S_HD]
				mask_L_MD = masks_small[:mask_L_MD]
				mask_M_MD = masks_small[:mask_M_MD]
				mask_S_MD = masks_small[:mask_S_MD]
				mask_L_LD = masks_small[:mask_L_LD]
				mask_M_LD = masks_small[:mask_M_LD]
				mask_S_LD = masks_small[:mask_S_LD]
			elseif _size == "medium"
				center_insert = [230, 370]
				mask_L_HD = masks_medium[:mask_L_HD]
				mask_M_HD = masks_medium[:mask_M_HD]
				mask_S_HD = masks_medium[:mask_S_HD]
				mask_L_MD = masks_medium[:mask_L_MD]
				mask_M_MD = masks_medium[:mask_M_MD]
				mask_S_MD = masks_medium[:mask_S_MD]
				mask_L_LD = masks_medium[:mask_L_LD]
				mask_M_LD = masks_medium[:mask_M_LD]
				mask_S_LD = masks_medium[:mask_S_LD]
			else
				center_insert = [285, 420]
				mask_L_HD = masks_large[:mask_L_HD]
				mask_M_HD = masks_large[:mask_M_HD]
				mask_S_HD = masks_large[:mask_S_HD]
				mask_L_MD = masks_large[:mask_L_MD]
				mask_M_MD = masks_large[:mask_M_MD]
				mask_S_MD = masks_large[:mask_S_MD]
				mask_L_LD = masks_large[:mask_L_LD]
				mask_M_LD = masks_large[:mask_M_LD]
				mask_S_LD = masks_large[:mask_S_LD]
			end

			# Large Insert Masks
			dilated_mask_L_HD = dilate_mask_large(mask_L_HD)
			dilated_mask_L_HD_3D = cat(dilated_mask_L_HD, dilated_mask_L_HD, dilated_mask_L_HD, dims=3)

			ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)
			ring_mask_L_HD_3D = cat(ring_mask_L_HD, ring_mask_L_HD, ring_mask_L_HD, dims=3)

			dilated_mask_L_MD = dilate_mask_large(mask_L_MD)
			dilated_mask_L_MD_3D = cat(dilated_mask_L_MD, dilated_mask_L_MD, dilated_mask_L_MD, dims=3)

			ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)
			ring_mask_L_MD_3D = cat(ring_mask_L_MD, ring_mask_L_MD, ring_mask_L_MD, dims=3)

			dilated_mask_L_LD = dilate_mask_large(mask_L_LD)
			dilated_mask_L_LD_3D = cat(dilated_mask_L_LD, dilated_mask_L_LD, dilated_mask_L_LD, dims=3)

			ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)
			ring_mask_L_LD_3D = cat(ring_mask_L_LD, ring_mask_L_LD, ring_mask_L_LD, dims=3)
			

			# Medium Insert Masks
			dilated_mask_M_HD = dilate_mask_medium(mask_M_HD)
			dilated_mask_M_HD_3D = cat(dilated_mask_M_HD, dilated_mask_M_HD, dilated_mask_M_HD, dims=3)

			ring_mask_M_HD = ring_mask_large(dilated_mask_M_HD)
			ring_mask_M_HD_3D = cat(ring_mask_M_HD, ring_mask_M_HD, ring_mask_M_HD, dims=3)

			dilated_mask_M_MD = dilate_mask_medium(mask_M_MD)
			dilated_mask_M_MD_3D = cat(dilated_mask_M_MD, dilated_mask_M_MD, dilated_mask_M_MD, dims=3)

			ring_mask_M_MD = ring_mask_large(dilated_mask_M_MD)
			ring_mask_M_MD_3D = cat(ring_mask_M_MD, ring_mask_M_MD, ring_mask_M_MD, dims=3)

			dilated_mask_M_LD = dilate_mask_medium(mask_M_LD)
			dilated_mask_M_LD_3D = cat(dilated_mask_M_LD, dilated_mask_M_LD, dilated_mask_M_LD, dims=3)

			ring_mask_M_LD = ring_mask_large(dilated_mask_M_LD)
			ring_mask_M_LD_3D = cat(ring_mask_M_LD, ring_mask_M_LD, ring_mask_M_LD, dims=3)

			
			# Small Insert Masks
			dilated_mask_S_HD = dilate_mask_small(mask_S_HD)
			dilated_mask_S_HD_3D = cat(dilated_mask_S_HD, dilated_mask_S_HD, dilated_mask_S_HD, dims=3)

			ring_mask_S_HD = ring_mask_large(dilated_mask_S_HD)
			ring_mask_S_HD_3D = cat(ring_mask_S_HD, ring_mask_S_HD, ring_mask_S_HD, dims=3)

			dilated_mask_S_MD = dilate_mask_small(mask_S_MD)
			dilated_mask_S_MD_3D = cat(dilated_mask_S_MD, dilated_mask_S_MD, dilated_mask_S_MD, dims=3)

			ring_mask_S_MD = ring_mask_large(dilated_mask_S_MD)
			ring_mask_S_MD_3D = cat(ring_mask_S_MD, ring_mask_S_MD, ring_mask_S_MD, dims=3)

			dilated_mask_S_LD = dilate_mask_small(mask_S_LD)
			dilated_mask_S_LD_3D = cat(dilated_mask_S_LD, dilated_mask_S_LD, dilated_mask_S_LD, dims=3)

			ring_mask_S_LD = ring_mask_large(dilated_mask_S_LD)
			ring_mask_S_LD_3D = cat(ring_mask_S_LD, ring_mask_S_LD, ring_mask_S_LD, dims=3)

			#------- Material Decomposition -------#
			
			# Low Energy
			path_80 = datadir("dcms", "val", density, _size, string(energies[1]))
			dcm_80 = dcmdir_parse(path_80)
			dcm_array_80 = load_dcm_array(dcm_80)
			pixel_size = get_pixel_size(dcm_80[1].meta)
			
			means_80 = [
				mean(dcm_array_80[dilated_mask_L_HD_3D]), mean(dcm_array_80[dilated_mask_L_MD_3D]), mean(dcm_array_80[dilated_mask_L_LD_3D]),
				mean(dcm_array_80[dilated_mask_M_HD_3D]), mean(dcm_array_80[dilated_mask_M_MD_3D]), mean(dcm_array_80[dilated_mask_M_LD_3D]),
				mean(dcm_array_80[dilated_mask_S_HD_3D]), mean(dcm_array_80[dilated_mask_S_MD_3D]), mean(dcm_array_80[dilated_mask_S_LD_3D])
			]

			# High Energy
			path_135 = datadir("dcms", "val", density, _size, string(energies[3]))
			dcm_135 = dcmdir_parse(path_135)
			dcm_array_135 = load_dcm_array(dcm_135)
			
			means_135 = [
				mean(dcm_array_135[dilated_mask_L_HD_3D]), mean(dcm_array_135[dilated_mask_L_MD_3D]), mean(dcm_array_135[dilated_mask_L_LD_3D]),
				mean(dcm_array_135[dilated_mask_M_HD_3D]), mean(dcm_array_135[dilated_mask_M_MD_3D]), mean(dcm_array_135[dilated_mask_M_LD_3D]),
				mean(dcm_array_135[dilated_mask_S_HD_3D]), mean(dcm_array_135[dilated_mask_S_MD_3D]), mean(dcm_array_135[dilated_mask_S_LD_3D])
			]

			# Background masks
			background_mask = zeros(size(dcm_array_80))
			background_mask[
				(center_insert[1]-5):(center_insert[1]+5),
				(center_insert[2]-5):(center_insert[2]+5),
				1:3,
			] .= 1

			dilated_mask_L_bkg = Bool.(dilate_mask_large_bkg(background_mask))
			ring_mask_L_bkg = Bool.(ring_mask_large(dilated_mask_L_bkg))

			dilated_mask_M_bkg = Bool.(dilate_mask_medium_bkg(background_mask))
			ring_mask_M_bkg = Bool.(ring_mask_medium(dilated_mask_M_bkg))

			dilated_mask_S_bkg = Bool.(dilate_mask_small_bkg(background_mask))
			ring_mask_S_bkg = Bool.(ring_mask_small(dilated_mask_S_bkg))

			# Background means
			means_80_bkg = [
				mean(dcm_array_80[dilated_mask_L_bkg]),
				mean(dcm_array_80[dilated_mask_M_bkg]),
				mean(dcm_array_80[dilated_mask_S_bkg])
			]

			means_135_bkg = [
				mean(dcm_array_135[dilated_mask_L_bkg]),
				mean(dcm_array_135[dilated_mask_M_bkg]),
				mean(dcm_array_135[dilated_mask_S_bkg])
			]
			
			
			# Calculate ground truth mass 
			# π * radius_mm^2 * slice_thickness_mm * number of slices
			vol_small_insert_gt = π * (1/2)^2 * (pixel_size[3] * 3) # mm^3
			vol_medium_insert_gt = π * (3/2)^2 * (pixel_size[3] * 3) # mm^3
			vol_large_insert_gt = π * (5/2)^2 * (pixel_size[3] * 3) # mm^3
			
			volumes_inserts_mm3 = [vol_small_insert_gt, vol_medium_insert_gt, 
			vol_large_insert_gt] # mm^3
			volumes_inserts_cm3 = volumes_inserts_mm3 * 1e-3 # cm^3
			
			gt_density = parse.(Int, split(density, "_"))
			gt_mass_large_inserts = gt_density .* volumes_inserts_cm3[3]
			gt_mass_medium_inserts = gt_density .* volumes_inserts_cm3[2]
			gt_mass_small_inserts = gt_density .* volumes_inserts_cm3[1]

			# Calculate predicted mass
			calculated_intensities = hcat(means_80, means_135)
			predicted_densities = zeros(size(means_135))
			for i in eachindex(predicted_densities)
				predicted_densities[i] = quantify(means_80[i], means_80[i], ps_low)
			end

			# Calculate background mass
			calculated_intensities_bkg = hcat(means_80_bkg, means_135_bkg)
			predicted_densities_bkg = zeros(size(means_80_bkg))
			for i in eachindex(predicted_densities_bkg)
				predicted_densities_bkg[i] = predict_concentration(means_80_bkg[i], means_135_bkg[i], ps_low)
			end
			
			voxel_size_mm3 = (pixel_size[1] * pixel_size[2] * pixel_size[3]) # mm^3
			voxel_size_cm3 = voxel_size_mm3 * 1e-3 # cm^3
			vol_large = count(dilated_mask_L_HD_3D) * voxel_size_cm3
			vol_medium = count(dilated_mask_M_HD_3D) * voxel_size_cm3
			vol_small = count(dilated_mask_S_HD_3D) * voxel_size_cm3
			
			predicted_mass_large_inserts = predicted_densities[1:3] .* vol_large
			predicted_mass_medium_inserts = predicted_densities[4:6] .* vol_medium
			predicted_mass_small_inserts = predicted_densities[7:9] .* vol_small

			vol_large_bkg = count(dilated_mask_L_HD_3D) * voxel_size_cm3
			vol_medium_bkg = count(dilated_mask_M_HD_3D) * voxel_size_cm3
			vol_small_bkg = count(dilated_mask_S_HD_3D) * voxel_size_cm3

			predicted_masses_bkg = [
				predicted_densities_bkg[1] * vol_large_bkg,
				predicted_densities_bkg[2] * vol_medium_bkg,
				predicted_densities_bkg[3] * vol_small_bkg,
			]

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				insert_densities = [:low_density, :medium_density, :high_density],
				gt_mass_large_inserts = gt_mass_large_inserts,
				predicted_mass_large_inserts = predicted_mass_large_inserts,
				gt_mass_medium_inserts = gt_mass_medium_inserts,
				predicted_mass_medium_inserts = predicted_mass_medium_inserts,
				gt_mass_small_inserts = gt_mass_small_inserts,
				predicted_mass_small_inserts = predicted_mass_small_inserts,
				bkg_mass = predicted_masses_bkg
			)
			push!(dfs_low_md, df_results)

			#------- Volume Fraction -------#
			
			path = datadir("dcms", "val", density, _size, string(energies[2]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
			voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

			# Score Inserts
			hu_heart_tissue_large_hd = mean(dcm_array[ring_mask_L_HD_3D])
			mass_large_hd = score(dcm_array[dilated_mask_L_HD_3D], hu_calcium_100, hu_heart_tissue_large_hd, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_large_md = mean(dcm_array[ring_mask_L_MD_3D])
			mass_large_md = score(dcm_array[dilated_mask_L_MD_3D], hu_calcium_100, hu_heart_tissue_large_md, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_large_ld = mean(dcm_array[ring_mask_L_LD_3D])
			mass_large_ld = score(dcm_array[dilated_mask_L_LD_3D], hu_calcium_100, hu_heart_tissue_large_ld, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_medium_hd = mean(dcm_array[ring_mask_M_HD_3D])
			mass_medium_hd = score(dcm_array[dilated_mask_M_HD_3D], hu_calcium_100, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_medium_md = mean(dcm_array[ring_mask_M_MD_3D])
			mass_medium_md = score(dcm_array[dilated_mask_M_MD_3D], hu_calcium_100, hu_heart_tissue_medium_md, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_medium_ld = mean(dcm_array[ring_mask_M_LD_3D])
			mass_medium_ld = score(dcm_array[dilated_mask_M_LD_3D], hu_calcium_100, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_small_hd = mean(dcm_array[ring_mask_S_HD_3D])
			mass_small_hd = score(dcm_array[dilated_mask_S_HD_3D], hu_calcium_100, hu_heart_tissue_large_hd, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_small_md = mean(dcm_array[ring_mask_S_MD_3D])
			mass_small_md = score(dcm_array[dilated_mask_S_MD_3D], hu_calcium_100, hu_heart_tissue_large_md, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_small_ld = mean(dcm_array[ring_mask_S_LD_3D])
			mass_small_ld = score(dcm_array[dilated_mask_S_LD_3D], hu_calcium_100, hu_heart_tissue_large_ld, voxel_size, ρ_calcium_100, VolumeFraction())

			predicted_mass_large_inserts = [mass_large_hd, mass_large_md, mass_large_ld]
			predicted_mass_medium_inserts = [mass_medium_hd, mass_medium_md, mass_medium_ld]
			predicted_mass_small_inserts = [mass_small_hd, mass_small_md, mass_small_ld]

			# Background
			hu_heart_tissue_large_bkg = mean(dcm_array[ring_mask_L_bkg])
			mass_large_bkg = score(dcm_array[dilated_mask_L_bkg], hu_calcium_100, hu_heart_tissue_large_bkg, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_medium_bkg = mean(dcm_array[ring_mask_M_bkg])
			mass_medium_bkg = score(dcm_array[dilated_mask_M_bkg], hu_calcium_100, hu_heart_tissue_medium_bkg, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_small_bkg = mean(dcm_array[ring_mask_S_bkg])
			mass_small_bkg = score(dcm_array[dilated_mask_S_bkg], hu_calcium_100, hu_heart_tissue_small_bkg, voxel_size, ρ_calcium_100, VolumeFraction())

			mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				insert_densities = [:low_density, :medium_density, :high_density],
				gt_mass_large_inserts = gt_mass_large_inserts,
				predicted_mass_large_inserts = predicted_mass_large_inserts,
				gt_mass_medium_inserts = gt_mass_medium_inserts,
				predicted_mass_medium_inserts = predicted_mass_medium_inserts,
				gt_mass_small_inserts = gt_mass_small_inserts,
				predicted_mass_small_inserts = predicted_mass_small_inserts,
				bkg_mass = mass_bkg
			)
			push!(dfs_low_vf, df_results)

			#------- Agatston -------#
			
			alg = Agatston()
			mass_cal_factor = ρ_calcium_100 / hu_calcium_100
			kV = 120

			overlayed_bkg_mask_L = create_mask(dcm_array, dilated_mask_L_bkg)
			overlayed_bkg_mask_M = create_mask(dcm_array, dilated_mask_M_bkg)
			overlayed_bkg_mask_S = create_mask(dcm_array, dilated_mask_S_bkg)
			
			_, _, mass_bkg_large = score(
				overlayed_bkg_mask_L,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			_, _, mass_bkg_medium = score(
				overlayed_bkg_mask_M,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			_, _, mass_bkg_small = score(
				overlayed_bkg_mask_S,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			mass_bkg = [mass_bkg_large, mass_bkg_medium, mass_bkg_small]
				
			overlayed_mask_l_hd = create_mask(dcm_array, dilated_mask_L_HD_3D)
			agat_l_hd, vol_l_hd, mass_l_hd = score(
				overlayed_mask_l_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_l_md = create_mask(dcm_array, dilated_mask_L_MD_3D)
			agat_l_md, vol_l_md, mass_l_md = score(
				overlayed_mask_l_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_l_ld = create_mask(dcm_array, dilated_mask_L_LD_3D)
			agat_l_ld, vol_l_ld, mass_l_ld = score(
				overlayed_mask_l_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			
			overlayed_mask_m_hd = create_mask(dcm_array, dilated_mask_M_HD_3D)
			agat_m_hd, vol_m_hd, mass_m_hd = score(
				overlayed_mask_m_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_m_md = create_mask(dcm_array, dilated_mask_M_MD_3D)
			agat_m_md, vol_m_md, mass_m_md = score(
				overlayed_mask_m_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_m_ld = create_mask(dcm_array, dilated_mask_M_LD_3D)
			agat_m_ld, vol_m_ld, mass_m_ld = score(
				overlayed_mask_m_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_s_hd = create_mask(dcm_array, dilated_mask_S_HD_3D)
			agat_s_hd, vol_s_hd, mass_s_hd = score(
				overlayed_mask_s_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_s_md = create_mask(dcm_array, dilated_mask_S_MD_3D)
			agat_s_md, vol_s_md, mass_s_md = score(
				overlayed_mask_s_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_s_ld = create_mask(dcm_array, dilated_mask_S_LD_3D)
			agat_s_ld, vol_s_ld, mass_s_ld = score(
				overlayed_mask_s_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			predicted_mass_large_inserts = [mass_l_hd, mass_l_md, mass_l_ld]
			predicted_mass_medium_inserts = [mass_m_hd, mass_m_md, mass_m_ld]
			predicted_mass_small_inserts = [mass_s_hd, mass_s_md, mass_s_ld]

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				insert_densities = [:low_density, :medium_density, :high_density],
				gt_mass_large_inserts = gt_mass_large_inserts,
				predicted_mass_large_inserts = predicted_mass_large_inserts,
				gt_mass_medium_inserts = gt_mass_medium_inserts,
				predicted_mass_medium_inserts = predicted_mass_medium_inserts,
				gt_mass_small_inserts = gt_mass_small_inserts,
				predicted_mass_small_inserts = predicted_mass_small_inserts,
				bkg_mass = mass_bkg
			)
			push!(dfs_low_a, df_results)	
		end
	end
end

# ╔═╡ c2db3887-3be1-499a-ba8c-d0a9f3f3f645
md"""
# High Density
"""

# ╔═╡ 7cfd0261-208d-41b6-b356-63d9e02ca13a
md"""
## Calibration
"""

# ╔═╡ cac42faf-836c-42a2-bed2-f0bba2a13be0
# optimal_cal_high = [10, 25, 45, 75, 150, 250, 350, 450, 550, 650]
optimal_cal_high = densities_cal_all

# ╔═╡ bd61d48e-886b-40c6-a54a-76468c52cd4e
begin
	densities_cal_high = optimal_cal_high

	# Low Energy
	means_80_high = Dict(:density => densities_cal_high, :means => zeros(length(densities_cal_high)))
	for (i, density) in enumerate(densities_cal_high)
		if i ==1 #for background
			path = joinpath(datadir("dcms", "cal", "100", "30"), string(energies[1]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
			offset = 5
			bkg_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array, 3))
			center_insert_bkg, center_insert_bkg2 = 250, 318

			for z in axes(dcm_array, 3)
					rows, cols, depth = size(dcm_array)
					half_row, half_col = center_insert_bkg, center_insert_bkg2
					row_range = half_row-offset:half_row+offset
					col_range = half_col-offset:half_col+offset	
					calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
			
			means_80_high[:means][i] = mean(calibration_rod)

		else
			for _size in sizes_cal
				path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[1]))
				dcm = dcmdir_parse(path)
				dcm_array = load_dcm_array(dcm)
			
				center_insert1, center_insert2 = 187, 318
				offset = 5
				calibration_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array, 3))
				
				for z in axes(dcm_array, 3)
					rows, cols, depth = size(dcm_array)
					half_row, half_col = center_insert1, center_insert2
					row_range = half_row-offset:half_row+offset
					col_range = half_col-offset:half_col+offset	
					calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
				end
			
				means_80_high[:means][i] = mean(calibration_rod)
			end
		end 
	end

	# High energy	
	means_135_high = Dict(:density => densities_cal_high, :means => zeros(length(densities_cal_high)))
	for (i, density) in enumerate(densities_cal_high)
		
		if i ==1 #for background
			path = joinpath(datadir("dcms", "cal", "200", "30"), string(energies[3]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
			offset = 5
			bkg_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array, 3))
			center_insert_bkg, center_insert_bkg2 = 250, 318

			for z in axes(dcm_array, 3)
					rows, cols, depth = size(dcm_array)
					half_row, half_col = center_insert_bkg, center_insert_bkg2
					row_range = half_row-offset:half_row+offset
					col_range = half_col-offset:half_col+offset	
					calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
			
			means_135_high[:means][i] = mean(calibration_rod)

		else
			for _size in sizes_cal
				path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[3]))
				dcm = dcmdir_parse(path)
				dcm_array = load_dcm_array(dcm)
			
				center_insert1, center_insert2 = 187, 318
				
				offset = 5
				calibration_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array, 3))
				
				for z in axes(dcm_array, 3)
					rows, cols, depth = size(dcm_array)
					half_row, half_col = center_insert1, center_insert2
					row_range = half_row-offset:half_row+offset
					col_range = half_col-offset:half_col+offset	
					calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
				end
			
				means_135_high[:means][i] = mean(calibration_rod)
			end
		end 
	end

	# Fit Parameters
	calculated_intensities_high = hcat(means_80_high[:means], means_135_high[:means]) # low energy, high energy
		
	ps_high = fit_calibration(calculated_intensities_high, densities_cal_high)
end

# ╔═╡ 6c01fb2c-c8cf-4484-8509-3fc540c038a0
md"""
### Check Results
"""

# ╔═╡ b410ca06-2fe0-4671-925e-bbdaf9786849
begin
	predicted_densities_high = []
	
	for i in 1:length(densities_cal_high)
		append!(
			predicted_densities_high, 
			quantify(calculated_intensities_high[i, 1], calculated_intensities_high[i, 2], ps_high
			)
		)
	end
end

# ╔═╡ e7e5385b-4016-4bb3-8fd8-c203f8e5f8d8
df_high = DataFrame(
	densities = densities_cal_high,
	predicted_densities = predicted_densities_high,
)

# ╔═╡ d01812f8-e0ec-4b2d-9992-e7db2610d4fe
md"""
## Validation
"""

# ╔═╡ 397e2b5a-8d37-4398-b92b-678b4ba080ee
densities_val_high = densities_val_all[length(densities_val_low):end];

# ╔═╡ e0c01bd6-0cec-408d-b534-cf36cb5bd78a
md"""
### Predict
"""

# ╔═╡ 6952fc87-c733-4e67-b1aa-ec9dd46ebe3d
begin
	dfs_high_md = []
	dfs_high_vf = []
	dfs_high_a = []
	for density in densities_val_high
		for _size in sizes_val
			if _size == "small"
				center_insert = [175, 320]
				mask_L_HD = masks_small[:mask_L_HD]
				mask_M_HD = masks_small[:mask_M_HD]
				mask_S_HD = masks_small[:mask_S_HD]
				mask_L_MD = masks_small[:mask_L_MD]
				mask_M_MD = masks_small[:mask_M_MD]
				mask_S_MD = masks_small[:mask_S_MD]
				mask_L_LD = masks_small[:mask_L_LD]
				mask_M_LD = masks_small[:mask_M_LD]
				mask_S_LD = masks_small[:mask_S_LD]
			elseif _size == "medium"
				center_insert = [230, 370]
				mask_L_HD = masks_medium[:mask_L_HD]
				mask_M_HD = masks_medium[:mask_M_HD]
				mask_S_HD = masks_medium[:mask_S_HD]
				mask_L_MD = masks_medium[:mask_L_MD]
				mask_M_MD = masks_medium[:mask_M_MD]
				mask_S_MD = masks_medium[:mask_S_MD]
				mask_L_LD = masks_medium[:mask_L_LD]
				mask_M_LD = masks_medium[:mask_M_LD]
				mask_S_LD = masks_medium[:mask_S_LD]
			else
				center_insert = [285, 420]
				mask_L_HD = masks_large[:mask_L_HD]
				mask_M_HD = masks_large[:mask_M_HD]
				mask_S_HD = masks_large[:mask_S_HD]
				mask_L_MD = masks_large[:mask_L_MD]
				mask_M_MD = masks_large[:mask_M_MD]
				mask_S_MD = masks_large[:mask_S_MD]
				mask_L_LD = masks_large[:mask_L_LD]
				mask_M_LD = masks_large[:mask_M_LD]
				mask_S_LD = masks_large[:mask_S_LD]
			end

			# Large Insert Masks
			dilated_mask_L_HD = dilate_mask_large(mask_L_HD)
			dilated_mask_L_HD_3D = cat(dilated_mask_L_HD, dilated_mask_L_HD, dilated_mask_L_HD, dims=3)

			ring_mask_L_HD = ring_mask_large(dilated_mask_L_HD)
			ring_mask_L_HD_3D = cat(ring_mask_L_HD, ring_mask_L_HD, ring_mask_L_HD, dims=3)

			dilated_mask_L_MD = dilate_mask_large(mask_L_MD)
			dilated_mask_L_MD_3D = cat(dilated_mask_L_MD, dilated_mask_L_MD, dilated_mask_L_MD, dims=3)

			ring_mask_L_MD = ring_mask_large(dilated_mask_L_MD)
			ring_mask_L_MD_3D = cat(ring_mask_L_MD, ring_mask_L_MD, ring_mask_L_MD, dims=3)

			dilated_mask_L_LD = dilate_mask_large(mask_L_LD)
			dilated_mask_L_LD_3D = cat(dilated_mask_L_LD, dilated_mask_L_LD, dilated_mask_L_LD, dims=3)

			ring_mask_L_LD = ring_mask_large(dilated_mask_L_LD)
			ring_mask_L_LD_3D = cat(ring_mask_L_LD, ring_mask_L_LD, ring_mask_L_LD, dims=3)
			

			# Medium Insert Masks
			dilated_mask_M_HD = dilate_mask_medium(mask_M_HD)
			dilated_mask_M_HD_3D = cat(dilated_mask_M_HD, dilated_mask_M_HD, dilated_mask_M_HD, dims=3)

			ring_mask_M_HD = ring_mask_large(dilated_mask_M_HD)
			ring_mask_M_HD_3D = cat(ring_mask_M_HD, ring_mask_M_HD, ring_mask_M_HD, dims=3)

			dilated_mask_M_MD = dilate_mask_medium(mask_M_MD)
			dilated_mask_M_MD_3D = cat(dilated_mask_M_MD, dilated_mask_M_MD, dilated_mask_M_MD, dims=3)

			ring_mask_M_MD = ring_mask_large(dilated_mask_M_MD)
			ring_mask_M_MD_3D = cat(ring_mask_M_MD, ring_mask_M_MD, ring_mask_M_MD, dims=3)

			dilated_mask_M_LD = dilate_mask_medium(mask_M_LD)
			dilated_mask_M_LD_3D = cat(dilated_mask_M_LD, dilated_mask_M_LD, dilated_mask_M_LD, dims=3)

			ring_mask_M_LD = ring_mask_large(dilated_mask_M_LD)
			ring_mask_M_LD_3D = cat(ring_mask_M_LD, ring_mask_M_LD, ring_mask_M_LD, dims=3)

			
			# Small Insert Masks
			dilated_mask_S_HD = dilate_mask_small(mask_S_HD)
			dilated_mask_S_HD_3D = cat(dilated_mask_S_HD, dilated_mask_S_HD, dilated_mask_S_HD, dims=3)

			ring_mask_S_HD = ring_mask_large(dilated_mask_S_HD)
			ring_mask_S_HD_3D = cat(ring_mask_S_HD, ring_mask_S_HD, ring_mask_S_HD, dims=3)

			dilated_mask_S_MD = dilate_mask_small(mask_S_MD)
			dilated_mask_S_MD_3D = cat(dilated_mask_S_MD, dilated_mask_S_MD, dilated_mask_S_MD, dims=3)

			ring_mask_S_MD = ring_mask_large(dilated_mask_S_MD)
			ring_mask_S_MD_3D = cat(ring_mask_S_MD, ring_mask_S_MD, ring_mask_S_MD, dims=3)

			dilated_mask_S_LD = dilate_mask_small(mask_S_LD)
			dilated_mask_S_LD_3D = cat(dilated_mask_S_LD, dilated_mask_S_LD, dilated_mask_S_LD, dims=3)

			ring_mask_S_LD = ring_mask_large(dilated_mask_S_LD)
			ring_mask_S_LD_3D = cat(ring_mask_S_LD, ring_mask_S_LD, ring_mask_S_LD, dims=3)

			#------- Material Decomposition -------#
			
			# Low Energy
			path_80 = datadir("dcms", "val", density, _size, string(energies[1]))
			dcm_80 = dcmdir_parse(path_80)
			dcm_array_80 = load_dcm_array(dcm_80)
			pixel_size = get_pixel_size(dcm_80[1].meta)
			
			means_80 = [
				mean(dcm_array_80[dilated_mask_L_HD_3D]), mean(dcm_array_80[dilated_mask_L_MD_3D]), mean(dcm_array_80[dilated_mask_L_LD_3D]),
				mean(dcm_array_80[dilated_mask_M_HD_3D]), mean(dcm_array_80[dilated_mask_M_MD_3D]), mean(dcm_array_80[dilated_mask_M_LD_3D]),
				mean(dcm_array_80[dilated_mask_S_HD_3D]), mean(dcm_array_80[dilated_mask_S_MD_3D]), mean(dcm_array_80[dilated_mask_S_LD_3D])
			]

			# High Energy
			path_135 = datadir("dcms", "val", density, _size, string(energies[3]))
			dcm_135 = dcmdir_parse(path_135)
			dcm_array_135 = load_dcm_array(dcm_135)
			
			means_135 = [
				mean(dcm_array_135[dilated_mask_L_HD_3D]), mean(dcm_array_135[dilated_mask_L_MD_3D]), mean(dcm_array_135[dilated_mask_L_LD_3D]),
				mean(dcm_array_135[dilated_mask_M_HD_3D]), mean(dcm_array_135[dilated_mask_M_MD_3D]), mean(dcm_array_135[dilated_mask_M_LD_3D]),
				mean(dcm_array_135[dilated_mask_S_HD_3D]), mean(dcm_array_135[dilated_mask_S_MD_3D]), mean(dcm_array_135[dilated_mask_S_LD_3D])
			]

			# Background masks
			background_mask = zeros(size(dcm_array_80))
			background_mask[
				(center_insert[1]-5):(center_insert[1]+5),
				(center_insert[2]-5):(center_insert[2]+5),
				1:3,
			] .= 1

			dilated_mask_L_bkg = Bool.(dilate_mask_large_bkg(background_mask))
			ring_mask_L_bkg = Bool.(ring_mask_large(dilated_mask_L_bkg))

			dilated_mask_M_bkg = Bool.(dilate_mask_medium_bkg(background_mask))
			ring_mask_M_bkg = Bool.(ring_mask_medium(dilated_mask_M_bkg))

			dilated_mask_S_bkg = Bool.(dilate_mask_small_bkg(background_mask))
			ring_mask_S_bkg = Bool.(ring_mask_small(dilated_mask_S_bkg))

			# Background means
			means_80_bkg = [
				mean(dcm_array_80[dilated_mask_L_bkg]),
				mean(dcm_array_80[dilated_mask_M_bkg]),
				mean(dcm_array_80[dilated_mask_S_bkg])
			]

			means_135_bkg = [
				mean(dcm_array_135[dilated_mask_L_bkg]),
				mean(dcm_array_135[dilated_mask_M_bkg]),
				mean(dcm_array_135[dilated_mask_S_bkg])
			]
			
			
			# Calculate ground truth mass 
			# π * radius_mm^2 * slice_thickness_mm * number of slices
			vol_small_insert_gt = π * (1/2)^2 * (pixel_size[3] * 3) # mm^3
			vol_medium_insert_gt = π * (3/2)^2 * (pixel_size[3] * 3) # mm^3
			vol_large_insert_gt = π * (5/2)^2 * (pixel_size[3] * 3) # mm^3
			
			volumes_inserts_mm3 = [vol_small_insert_gt, vol_medium_insert_gt, 
			vol_large_insert_gt] # mm^3
			volumes_inserts_cm3 = volumes_inserts_mm3 * 1e-3 # cm^3
			
			gt_density = parse.(Int, split(density, "_"))
			gt_mass_large_inserts = gt_density .* volumes_inserts_cm3[3]
			gt_mass_medium_inserts = gt_density .* volumes_inserts_cm3[2]
			gt_mass_small_inserts = gt_density .* volumes_inserts_cm3[1]

			# Calculate predicted mass
			calculated_intensities = hcat(means_80, means_135)
			predicted_densities = zeros(size(means_135))
			for i in eachindex(predicted_densities)
				predicted_densities[i] = quantify(means_80[i], means_80[i], ps_high)
			end

			# Calculate background mass
			calculated_intensities_bkg = hcat(means_80_bkg, means_135_bkg)
			predicted_densities_bkg = zeros(size(means_80_bkg))
			for i in eachindex(predicted_densities_bkg)
				predicted_densities_bkg[i] = predict_concentration(means_80_bkg[i], means_135_bkg[i], ps_high)
			end
			
			voxel_size_mm3 = (pixel_size[1] * pixel_size[2] * pixel_size[3]) # mm^3
			voxel_size_cm3 = voxel_size_mm3 * 1e-3 # cm^3
			vol_large = count(dilated_mask_L_HD_3D) * voxel_size_cm3
			vol_medium = count(dilated_mask_M_HD_3D) * voxel_size_cm3
			vol_small = count(dilated_mask_S_HD_3D) * voxel_size_cm3
			
			predicted_mass_large_inserts = predicted_densities[1:3] .* vol_large
			predicted_mass_medium_inserts = predicted_densities[4:6] .* vol_medium
			predicted_mass_small_inserts = predicted_densities[7:9] .* vol_small

			vol_large_bkg = count(dilated_mask_L_HD_3D) * voxel_size_cm3
			vol_medium_bkg = count(dilated_mask_M_HD_3D) * voxel_size_cm3
			vol_small_bkg = count(dilated_mask_S_HD_3D) * voxel_size_cm3

			predicted_masses_bkg = [
				predicted_densities_bkg[1] * vol_large_bkg,
				predicted_densities_bkg[2] * vol_medium_bkg,
				predicted_densities_bkg[3] * vol_small_bkg,
			]

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				insert_densities = [:low_density, :medium_density, :high_density],
				gt_mass_large_inserts = gt_mass_large_inserts,
				predicted_mass_large_inserts = predicted_mass_large_inserts,
				gt_mass_medium_inserts = gt_mass_medium_inserts,
				predicted_mass_medium_inserts = predicted_mass_medium_inserts,
				gt_mass_small_inserts = gt_mass_small_inserts,
				predicted_mass_small_inserts = predicted_mass_small_inserts,
				bkg_mass = predicted_masses_bkg
			)
			push!(dfs_high_md, df_results)

			#------- Volume Fraction -------#
			
			path = datadir("dcms", "val", density, _size, string(energies[2]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
			voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

			# Score Inserts
			hu_heart_tissue_large_hd = mean(dcm_array[ring_mask_L_HD_3D])
			mass_large_hd = score(dcm_array[dilated_mask_L_HD_3D], hu_calcium_100, hu_heart_tissue_large_hd, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_large_md = mean(dcm_array[ring_mask_L_MD_3D])
			mass_large_md = score(dcm_array[dilated_mask_L_MD_3D], hu_calcium_100, hu_heart_tissue_large_md, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_large_ld = mean(dcm_array[ring_mask_L_LD_3D])
			mass_large_ld = score(dcm_array[dilated_mask_L_LD_3D], hu_calcium_100, hu_heart_tissue_large_ld, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_medium_hd = mean(dcm_array[ring_mask_M_HD_3D])
			mass_medium_hd = score(dcm_array[dilated_mask_M_HD_3D], hu_calcium_100, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_medium_md = mean(dcm_array[ring_mask_M_MD_3D])
			mass_medium_md = score(dcm_array[dilated_mask_M_MD_3D], hu_calcium_100, hu_heart_tissue_medium_md, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_medium_ld = mean(dcm_array[ring_mask_M_LD_3D])
			mass_medium_ld = score(dcm_array[dilated_mask_M_LD_3D], hu_calcium_100, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_small_hd = mean(dcm_array[ring_mask_S_HD_3D])
			mass_small_hd = score(dcm_array[dilated_mask_S_HD_3D], hu_calcium_100, hu_heart_tissue_large_hd, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_small_md = mean(dcm_array[ring_mask_S_MD_3D])
			mass_small_md = score(dcm_array[dilated_mask_S_MD_3D], hu_calcium_100, hu_heart_tissue_large_md, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_small_ld = mean(dcm_array[ring_mask_S_LD_3D])
			mass_small_ld = score(dcm_array[dilated_mask_S_LD_3D], hu_calcium_100, hu_heart_tissue_large_ld, voxel_size, ρ_calcium_100, VolumeFraction())

			predicted_mass_large_inserts = [mass_large_hd, mass_large_md, mass_large_ld]
			predicted_mass_medium_inserts = [mass_medium_hd, mass_medium_md, mass_medium_ld]
			predicted_mass_small_inserts = [mass_small_hd, mass_small_md, mass_small_ld]

			# Background
			hu_heart_tissue_large_bkg = mean(dcm_array[ring_mask_L_bkg])
			mass_large_bkg = score(dcm_array[dilated_mask_L_bkg], hu_calcium_100, hu_heart_tissue_large_bkg, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_medium_bkg = mean(dcm_array[ring_mask_M_bkg])
			mass_medium_bkg = score(dcm_array[dilated_mask_M_bkg], hu_calcium_100, hu_heart_tissue_medium_bkg, voxel_size, ρ_calcium_100, VolumeFraction())

			hu_heart_tissue_small_bkg = mean(dcm_array[ring_mask_S_bkg])
			mass_small_bkg = score(dcm_array[dilated_mask_S_bkg], hu_calcium_100, hu_heart_tissue_small_bkg, voxel_size, ρ_calcium_100, VolumeFraction())

			mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				insert_densities = [:low_density, :medium_density, :high_density],
				gt_mass_large_inserts = gt_mass_large_inserts,
				predicted_mass_large_inserts = predicted_mass_large_inserts,
				gt_mass_medium_inserts = gt_mass_medium_inserts,
				predicted_mass_medium_inserts = predicted_mass_medium_inserts,
				gt_mass_small_inserts = gt_mass_small_inserts,
				predicted_mass_small_inserts = predicted_mass_small_inserts,
				bkg_mass = mass_bkg
			)
			push!(dfs_high_vf, df_results)

			#------- Agatston -------#
			
			alg = Agatston()
			mass_cal_factor = ρ_calcium_100 / hu_calcium_100
			kV = 120

			overlayed_bkg_mask_L = create_mask(dcm_array, dilated_mask_L_bkg)
			overlayed_bkg_mask_M = create_mask(dcm_array, dilated_mask_M_bkg)
			overlayed_bkg_mask_S = create_mask(dcm_array, dilated_mask_S_bkg)
			
			_, _, mass_bkg_large = score(
				overlayed_bkg_mask_L,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			_, _, mass_bkg_medium = score(
				overlayed_bkg_mask_M,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			_, _, mass_bkg_small = score(
				overlayed_bkg_mask_S,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			mass_bkg = [mass_bkg_large, mass_bkg_medium, mass_bkg_small]
				
			overlayed_mask_l_hd = create_mask(dcm_array, dilated_mask_L_HD_3D)
			agat_l_hd, vol_l_hd, mass_l_hd = score(
				overlayed_mask_l_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_l_md = create_mask(dcm_array, dilated_mask_L_MD_3D)
			agat_l_md, vol_l_md, mass_l_md = score(
				overlayed_mask_l_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_l_ld = create_mask(dcm_array, dilated_mask_L_LD_3D)
			agat_l_ld, vol_l_ld, mass_l_ld = score(
				overlayed_mask_l_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)
			
			overlayed_mask_m_hd = create_mask(dcm_array, dilated_mask_M_HD_3D)
			agat_m_hd, vol_m_hd, mass_m_hd = score(
				overlayed_mask_m_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_m_md = create_mask(dcm_array, dilated_mask_M_MD_3D)
			agat_m_md, vol_m_md, mass_m_md = score(
				overlayed_mask_m_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_m_ld = create_mask(dcm_array, dilated_mask_M_LD_3D)
			agat_m_ld, vol_m_ld, mass_m_ld = score(
				overlayed_mask_m_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_s_hd = create_mask(dcm_array, dilated_mask_S_HD_3D)
			agat_s_hd, vol_s_hd, mass_s_hd = score(
				overlayed_mask_s_hd,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_s_md = create_mask(dcm_array, dilated_mask_S_MD_3D)
			agat_s_md, vol_s_md, mass_s_md = score(
				overlayed_mask_s_md,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			overlayed_mask_s_ld = create_mask(dcm_array, dilated_mask_S_LD_3D)
			agat_s_ld, vol_s_ld, mass_s_ld = score(
				overlayed_mask_s_ld,
				pixel_size,
				mass_cal_factor,
				alg;
				kV=kV
			)

			predicted_mass_large_inserts = [mass_l_hd, mass_l_md, mass_l_ld]
			predicted_mass_medium_inserts = [mass_m_hd, mass_m_md, mass_m_ld]
			predicted_mass_small_inserts = [mass_s_hd, mass_s_md, mass_s_ld]

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				insert_densities = [:low_density, :medium_density, :high_density],
				gt_mass_large_inserts = gt_mass_large_inserts,
				predicted_mass_large_inserts = predicted_mass_large_inserts,
				gt_mass_medium_inserts = gt_mass_medium_inserts,
				predicted_mass_medium_inserts = predicted_mass_medium_inserts,
				gt_mass_small_inserts = gt_mass_small_inserts,
				predicted_mass_small_inserts = predicted_mass_small_inserts,
				bkg_mass = mass_bkg
			)
			push!(dfs_high_a, df_results)
		end
	end
end

# ╔═╡ 1a70a6e7-108d-477e-b1a5-7e755cbc0619
md"""
# Analyze
"""

# ╔═╡ 3baf66b0-5e25-43d2-9b53-6b2d35f268a7
medphys_theme = Theme(
    Axis = (
        backgroundcolor = :white,
		xgridcolor = :gray,
		xgridwidth = 0.1,
		xlabelsize = 15,
		xticklabelsize = 15,
		ygridcolor = :gray,
		ygridwidth = 0.1,
		ylabelsize = 15,
		yticklabelsize = 15,
		bottomsplinecolor = :black,
		leftspinecolor = :black,
		titlesize = 25
	)
);

# ╔═╡ c3919e4d-c474-4aeb-bdc8-76e139e1d8c2
md"""
## Low Density
"""

# ╔═╡ fef76b74-ddfe-4578-9ccc-a8ff68b60137
md"""
### Accuracy
"""

# ╔═╡ 6d84cb36-80cf-43af-b100-caab847b347a
new_df_low_md = vcat(dfs_low_md[1:length(dfs_low_md)]...);

# ╔═╡ fd222d84-8b90-418e-b405-929557aaff8a
new_df_low_vf = vcat(dfs_low_vf[1:length(dfs_low_vf)]...);

# ╔═╡ a32f1676-cb1c-4b57-99a0-b6f7b542dbaf
new_df_low_a = vcat(dfs_low_a[1:length(dfs_low_a)]...);

# ╔═╡ c9f166c4-1269-42c5-85fb-21e0786a2ec5
co_1_low_md, r_squared_1_low_md, rms_values_1_low_md, pred_1_low_md = calculate_coefficients(new_df_low_md);

# ╔═╡ 2f742443-0dd4-4d25-90ff-4c7941cba0a2
co_1_low_vf, r_squared_1_low_vf, rms_values_1_low_vf, pred_1_low_vf = calculate_coefficients(new_df_low_vf);

# ╔═╡ 7b6def43-c729-4698-a752-9d37f182e816
co_1_low_a, r_squared_1_low_a, rms_values_1_low_a, pred_1_low_a = calculate_coefficients(new_df_low_a);

# ╔═╡ 1cb3b222-e86d-41a5-98f0-d4faa9363e87
function accuracy_low()
	f = Figure()

	##-- A --##
	ax = Axis(
		f[1, 1],
		xticks = collect(-5:10:25),
		yticks = collect(-5:10:25),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Material Decomposition (Low Density)",
	)
	
	df = new_df_low_md
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_low_md, linestyle=:dashdot)
	create_textbox(f[1, 1], co_1_low_md, r_squared_1_low_md, rms_values_1_low_md)
	
	xlims!(ax, low=-5, high=25)
	ylims!(ax, low=-5, high=25)

	ax = Axis(
		f[2, 1],
		xticks = collect(-5:10:25),
		yticks = collect(-5:10:25),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Volume Fraction (Low Density)",
	)
	
	df = new_df_low_vf
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_low_vf, linestyle=:dashdot)
	create_textbox(f[2, 1], co_1_low_vf, r_squared_1_low_vf, rms_values_1_low_vf)
	
	xlims!(ax, low=-5, high=25)
	ylims!(ax, low=-5, high=25)


	ax = Axis(
		f[3, 1],
		xticks = collect(-5:10:25),
		yticks = collect(-5:10:25),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Agatston (Low Density)",
	)
	
	df = new_df_low_a
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_low_a, linestyle=:dashdot)
	create_textbox(f[3, 1], co_1_low_a, r_squared_1_low_a, rms_values_1_low_a)
	
	xlims!(ax, low=-5, high=25)
	ylims!(ax, low=-5, high=25)

	#-- LABELS --##
	f[2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip([["A"], ["B"], ["C"]], [f[1,1], f[2, 1], f[3, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end
	f
end

# ╔═╡ c23d6948-d28a-4c66-ae2d-78686082e5f3
with_theme(accuracy_low, medphys_theme)

# ╔═╡ 6fb4896c-a086-4ef1-bbb7-1f1603338b3c
md"""
### Sensitivity & Specificity
"""

# ╔═╡ 7d218b4c-d467-4317-b239-a96de01a9d81
std_level = 1.5

# ╔═╡ d1d90206-1f7d-46b7-9063-8d6cd530a6e1
md"""
#### False-Negatives
"""

# ╔═╡ 82baca25-0d76-40a9-b19f-4838e6e84d72
begin
	# Agatston
	array_a = hcat(new_df_low_a[!, :predicted_mass_large_inserts], new_df_low_a[!, :predicted_mass_medium_inserts], new_df_low_a[!, :predicted_mass_small_inserts])
	total_cac = length(array_a)
	total_zero_a = length(findall(x -> x <= 0, array_a))

	# Material Decomposition
	total_zero_md = 0
	local df = new_df_low_md
	for i in range(start = 1, stop = nrow(df), step = 3)
		j = i:i+2
		μ, σ = mean(df[j, :bkg_mass]), std(df[j, :bkg_mass]) * std_level 
		array = hcat(
			df[j, :predicted_mass_large_inserts],
			df[j, :predicted_mass_medium_inserts],
			df[j, :predicted_mass_small_inserts]
		)
		negatives = length(findall(x -> x <= μ + σ, array))
		total_zero_md += negatives
	end

	# Volume Fraction
	total_zero_vf = 0
	local df = new_df_low_vf
	for i in range(start = 1, stop = nrow(df), step = 3)
		j = i:i+2
		μ, σ = mean(df[j, :bkg_mass]), std(df[j, :bkg_mass]) * std_level 
		array = hcat(
			df[j, :predicted_mass_large_inserts],
			df[j, :predicted_mass_medium_inserts],
			df[j, :predicted_mass_small_inserts]
		)
		negatives = length(findall(x -> x <= μ + σ, array))
		total_zero_vf += negatives
	end
end

# ╔═╡ a3b18554-6f32-4502-8fd4-d8550ddbbf07
total_zero_md, total_zero_vf, total_zero_a

# ╔═╡ f35612ba-8130-4c01-b91c-68886bccf358
md"""
#### False-Positives
"""

# ╔═╡ 7873647b-d7f8-4d8c-89d7-f15c8d64653a
begin
	# Agatston
	array_a_pos = new_df_low_a[!, :bkg_mass]
	total_cac_pos = length(array_a_pos)
	total_pos_a = length(findall(x -> x > 0, array_a_pos))

	# Material Decomposition
	total_pos_md = 0
	local df = new_df_low_md
	for i in range(start = 1, stop = nrow(df), step = 3)
		j = i:i+2
		μ, σ = mean(df[j, :bkg_mass]), std(df[j, :bkg_mass]) * std_level 
		array = df[j, :bkg_mass]
		positives = length(findall(x -> x > μ + σ, array))
		total_pos_md += positives
	end

	# Material Decomposition
	total_pos_vf = 0
	local df = new_df_low_vf
	for i in range(start = 1, stop = nrow(df), step = 3)
		j = i:i+2
		μ, σ = mean(df[j, :bkg_mass]), std(df[j, :bkg_mass]) * std_level 
		array = df[j, :bkg_mass]
		positives = length(findall(x -> x > μ + σ, array))
		total_pos_vf += positives
	end
end

# ╔═╡ 5ca5e077-3900-4948-9136-8e1b907758b1
total_pos_md, total_pos_vf, total_pos_a

# ╔═╡ ca5d95fe-ae81-41e5-9a54-ddc81962f2a3
function sensitivity_specificity()
    f = Figure()
    colors = Makie.wong_colors()

	labels = ["Material Decomposition", "Volume Fraction", "Agatston"]

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (eachindex(labels), labels),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_zero_md / total_cac) * 100
	h2 = (total_zero_vf / total_cac) * 100
	h3 = (total_zero_a / total_cac) * 100 
	heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[eachindex(labels)], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)

	##-- B --##
    ax = Axis(
		f[2, 1]; 
		xticks = (eachindex(labels), labels),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_pos_md / total_cac_pos) * 100
	h2 = (total_pos_vf / total_cac_pos) * 100
	h3 = (total_pos_a / total_cac_pos) * 100
    heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[eachindex(labels)], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	Label(
		f[1:2, 0],
		"Low Density",
		fontsize = 20,
		rotation = pi/2
	)

	# save(plotsdir("sensitivity_specificity_low.png"), f)

    return f
end

# ╔═╡ f1a8465b-aea5-46a5-ab95-ef181096a682
with_theme(sensitivity_specificity, medphys_theme)

# ╔═╡ 2c0270a0-36e9-4b9d-855b-834aa443a808
md"""
## High Density
"""

# ╔═╡ d9e33929-635a-450b-aa08-7ce0a6ee1170
md"""
### Accuracy
"""

# ╔═╡ 8a315d66-ff37-4e7d-9439-10567aabcd90
new_df_high_md = vcat(dfs_high_md[1:length(dfs_high_md)]...);

# ╔═╡ 603ae0fc-4fc5-4a2f-8dd7-e3b868bba417
new_df_high_vf = vcat(dfs_high_vf[1:length(dfs_high_vf)]...);

# ╔═╡ fe9b5fcb-1221-4898-976b-051b4896e470
new_df_high_a = vcat(dfs_high_a[1:length(dfs_high_a)]...);

# ╔═╡ 1d176c7b-802f-462d-94ca-97849df1e529
co_1_high_md, r_squared_1_high_md, rms_values_1_high_md, pred_1_high_md = calculate_coefficients(new_df_high_md);

# ╔═╡ ed9ec9d9-7262-4be2-99de-32c9fea8fc33
co_1_high_vf, r_squared_1_high_vf, rms_values_1_high_vf, pred_1_high_vf = calculate_coefficients(new_df_high_vf);

# ╔═╡ 9f3b22d1-e52a-4064-98e7-77856b5923bd
co_1_high_a, r_squared_1_high_a, rms_values_1_high_a, pred_1_high_a = calculate_coefficients(new_df_high_a);

# ╔═╡ 708ea6a0-2129-42b9-b05e-b2766bcaefd7
function accuracy_high()
	f = Figure()

	ax = Axis(
		f[1, 1],
		xticks = collect(-50:100:250),
		yticks = collect(-50:100:250),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Material Decomposition (High Density)",
	)
	
	df = new_df_high_md
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_high_md, linestyle=:dashdot)
	create_textbox(f[1, 1], co_1_high_md, r_squared_1_high_md, rms_values_1_high_md)
	
	xlims!(ax, low=-50, high=250)
	ylims!(ax, low=-50, high=250)

	ax = Axis(
		f[2, 1],
		xticks = collect(-50:100:250),
		yticks = collect(-50:100:250),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Volume Fraction (High Density)",
	)
	
	df = new_df_high_vf
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_high_vf, linestyle=:dashdot)
	create_textbox(f[2, 1], co_1_high_vf, r_squared_1_high_vf, rms_values_1_high_vf)
	
	xlims!(ax, low=-50, high=250)
	ylims!(ax, low=-50, high=250)

	ax = Axis(
		f[3, 1],
		xticks = collect(-50:100:250),
		yticks = collect(-50:100:250),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Agatston (High Density)",
	)
	
	df = new_df_high_a
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_high_a, linestyle=:dashdot)
	create_textbox(f[3, 1], co_1_high_a, r_squared_1_high_a, rms_values_1_high_a)
	
	xlims!(ax, low=-50, high=250)
	ylims!(ax, low=-50, high=250)

	#-- LABELS --##
	f[2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip([["A"], ["B"], ["C"]], [f[1,1], f[2, 1], f[3, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end
	f
end

# ╔═╡ 40854076-fcc6-481a-b1c8-12ded408f242
with_theme(accuracy_high, medphys_theme)

# ╔═╡ 2e0c8139-fe10-4b4e-b00d-c713208d8555
md"""
### Sensitivity & Specificity
"""

# ╔═╡ eaab5224-3302-4dc4-9d09-c2728b547f23
md"""
#### False-Negatives
"""

# ╔═╡ 986e1b6a-5052-4a1d-b3f8-65d05e8bf3bf
begin
	# Agatston
	array_a_high = hcat(new_df_high_a[!, :predicted_mass_large_inserts], new_df_high_a[!, :predicted_mass_medium_inserts], new_df_high_a[!, :predicted_mass_small_inserts])
	total_cac_high = length(array_a_high)
	total_zero_a_high = length(findall(x -> x <= 0, array_a_high))

	# Material Decomposition
	total_zero_md_high = 0
	local df = new_df_high_md
	for i in range(start = 1, stop = nrow(df), step = 3)
		j = i:i+2
		μ, σ = mean(df[j, :bkg_mass]), std(df[j, :bkg_mass]) * std_level 
		array = hcat(
			df[j, :predicted_mass_large_inserts],
			df[j, :predicted_mass_medium_inserts],
			df[j, :predicted_mass_small_inserts]
		)
		negatives = length(findall(x -> x <= μ + σ, array))
		total_zero_md_high += negatives
	end

	# Volume Fraction
	total_zero_vf_high = 0
	local df = new_df_high_vf
	for i in range(start = 1, stop = nrow(df), step = 3)
		j = i:i+2
		μ, σ = mean(df[j, :bkg_mass]), std(df[j, :bkg_mass]) * std_level 
		array = hcat(
			df[j, :predicted_mass_large_inserts],
			df[j, :predicted_mass_medium_inserts],
			df[j, :predicted_mass_small_inserts]
		)
		negatives = length(findall(x -> x <= μ + σ, array))
		total_zero_vf_high += negatives
	end
end

# ╔═╡ 6da99d21-dde5-4e21-aae9-f798e6d54a2b
total_zero_md_high, total_zero_vf_high, total_zero_a_high

# ╔═╡ 812e634d-ff38-43b9-8a8b-8a21a80654b5
md"""
#### False-Positives
"""

# ╔═╡ 7b0777fc-d104-480f-abd6-2f192ec539da
begin
	# Agatston
	array_a_pos_high = new_df_high_a[!, :bkg_mass]
	total_cac_pos_high = length(array_a_pos_high)
	total_pos_a_high = length(findall(x -> x > 0, array_a_pos_high))

	# Material Decomposition
	total_pos_md_high = 0
	local df = new_df_high_md
	for i in range(start = 1, stop = nrow(df), step = 3)
		j = i:i+2
		μ, σ = mean(df[j, :bkg_mass]), std(df[j, :bkg_mass]) * std_level 
		array = df[j, :bkg_mass]
		positives = length(findall(x -> x > μ + σ, array))
		total_pos_md_high += positives
	end

	# Material Decomposition
	total_pos_vf_high = 0
	local df = new_df_high_vf
	for i in range(start = 1, stop = nrow(df), step = 3)
		j = i:i+2
		μ, σ = mean(df[j, :bkg_mass]), std(df[j, :bkg_mass]) * std_level 
		array = df[j, :bkg_mass]
		positives = length(findall(x -> x > μ + σ, array))
		total_pos_vf_high += positives
	end
end

# ╔═╡ 77cf61eb-939e-42a7-917d-8fe4811522d5
total_pos_md, total_pos_vf, total_pos_a

# ╔═╡ f41d1a38-08f6-4d3e-8932-ef4d6783e065
function sensitivity_specificity_high()
    f = Figure()
    colors = Makie.wong_colors()

	labels = ["Material Decomposition", "Volume Fraction", "Agatston"]

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		xticks = (eachindex(labels), labels),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_zero_md_high / total_cac_high) * 100
	h2 = (total_zero_vf_high / total_cac_high) * 100
	h3 = (total_zero_a_high / total_cac_high) * 100 
	heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[eachindex(labels)], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)

	##-- B --##
    ax = Axis(
		f[2, 1]; 
		xticks = (eachindex(labels), labels),
		title = "False-Negative (CAC=0)",
		ylabel = "False-Negative (%)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_pos_md_high / total_cac_pos_high) * 100
	h2 = (total_pos_vf_high / total_cac_pos_high) * 100
	h3 = (total_pos_a_high / total_cac_pos_high) * 100
    heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[eachindex(labels)], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)

	
	for (label, layout) in zip(["A", "B"], [f[1,1], f[2,1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end

	Label(
		f[1:2, 0],
		"High Density",
		fontsize = 20,
		rotation = pi/2
	)

	# save(plotsdir("sensitivity_specificity_high.png"), f)

    return f
end

# ╔═╡ cb0a3ad1-de3c-475d-a727-1dcbe2ffefb4
with_theme(sensitivity_specificity_high, medphys_theme)

# ╔═╡ 39e42f15-de5c-445e-bd07-67745da578d7
md"""
## Combined
"""

# ╔═╡ 503fadca-83b4-444c-a3c9-a48fa6a67694
md"""
### Accuracy
"""

# ╔═╡ 0b0094a2-e267-4820-89e8-1cf95145975e
function accuracy_combined()
	f = Figure()

	#------ Low Density ------#
	ax = Axis(
		f[1, 1],
		xticks = collect(-5:10:25),
		yticks = collect(-5:10:25),
		title = "Material Decomposition (Low)",
	)
	
	df = new_df_low_md
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_low_md, linestyle=:dashdot)
	create_textbox(f[1, 1], co_1_low_md, r_squared_1_low_md, rms_values_1_low_md)
	
	xlims!(ax, low=-5, high=25)
	ylims!(ax, low=-5, high=25)

	ax = Axis(
		f[2, 1],
		xticks = collect(-5:10:25),
		yticks = collect(-5:10:25),
		title = "Volume Fraction (Low)",
	)
	
	df = new_df_low_vf
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_low_vf, linestyle=:dashdot)
	create_textbox(f[2, 1], co_1_low_vf, r_squared_1_low_vf, rms_values_1_low_vf)
	
	xlims!(ax, low=-5, high=25)
	ylims!(ax, low=-5, high=25)


	ax = Axis(
		f[3, 1],
		xticks = collect(-5:10:25),
		yticks = collect(-5:10:25),
		title = "Agatston (Low)",
	)
	
	df = new_df_low_a
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_low_a, linestyle=:dashdot)
	create_textbox(f[3, 1], co_1_low_a, r_squared_1_low_a, rms_values_1_low_a)
	
	xlims!(ax, low=-5, high=25)
	ylims!(ax, low=-5, high=25)


	#------ High Density ------#
	ax = Axis(
		f[1, 2],
		xticks = collect(-50:100:250),
		yticks = collect(-50:100:250),
		title = "Material Decomposition (High)",
	)
	
	df = new_df_high_md
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_high_md, linestyle=:dashdot)
	create_textbox(f[1, 2], co_1_high_md, r_squared_1_high_md, rms_values_1_high_md)
	
	xlims!(ax, low=-50, high=250)
	ylims!(ax, low=-50, high=250)

	ax = Axis(
		f[2, 2],
		xticks = collect(-50:100:250),
		yticks = collect(-50:100:250),
		title = "Volume Fraction (High)",
	)
	
	df = new_df_high_vf
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_high_vf, linestyle=:dashdot)
	create_textbox(f[2, 2], co_1_high_vf, r_squared_1_high_vf, rms_values_1_high_vf)
	
	xlims!(ax, low=-50, high=250)
	ylims!(ax, low=-50, high=250)

	ax = Axis(
		f[3, 2],
		xticks = collect(-50:100:250),
		yticks = collect(-50:100:250),
		title = "Agatston (High)",
	)
	
	df = new_df_high_a
	sc1 = scatter!(
		df[!, :gt_mass_large_inserts], df[!, :predicted_mass_large_inserts]
	)
	sc2 = scatter!(
		df[!, :gt_mass_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	sc3 = scatter!(
		df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	)
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_high_a, linestyle=:dashdot)
	create_textbox(f[3, 2], co_1_high_a, r_squared_1_high_a, rms_values_1_high_a)
	
	xlims!(ax, low=-50, high=250)
	ylims!(ax, low=-50, high=250)

	#-- LABELS --##
	f[2, 3] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	Label(
		f[2, 0],
		"Calculated Mass (mg)",
		rotation = pi/2,
		fontsize = 12
	)
	Label(
		f[end+1, 1:2],
		"Known Mass (mg)",
		fontsize = 12
	)
	
	save(plotsdir("accuracy_combined.png"), f)
	
	f
end

# ╔═╡ c27811d2-d651-4198-bf90-246c4dafd8ad
accuracy_combined()

# ╔═╡ 8060a585-765a-47ac-bc66-04808bd401f6
md"""
### Sensitivity & Specificity
"""

# ╔═╡ 7d7afb9b-2dbd-4373-b745-c245f86e48f7
function sensitivity_specificity_combined()
    f = Figure()
    colors = Makie.wong_colors()

	labels = ["Material \nDecomposition", "Volume \nFraction", "Agatston"]

    ##-- A --##
    ax = Axis(
		f[1, 1]; 
		title = "False-Negative (Low)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_zero_md / total_cac) * 100
	h2 = (total_zero_vf / total_cac) * 100
	h3 = (total_zero_a / total_cac) * 100 
	heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[eachindex(labels)], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)
	hidexdecorations!(ax)

	##-- B --##
    ax = Axis(
		f[2, 1]; 
		title = "False-Positive (Low)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_pos_md / total_cac_pos) * 100
	h2 = (total_pos_vf / total_cac_pos) * 100
	h3 = (total_pos_a / total_cac_pos) * 100
    heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(
		table, heights1; 
		color=colors[eachindex(labels)], bar_labels=[l1, l2, l3]
	)

    ylims!(ax; low=0, high=100)
	hidexdecorations!(ax)

	##-- A --##
    ax = Axis(
		f[1, 2]; 
		title = "False-Negative (High)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_zero_md_high / total_cac_high) * 100
	h2 = (total_zero_vf_high / total_cac_high) * 100
	h3 = (total_zero_a_high / total_cac_high) * 100 
	heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[eachindex(labels)], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)
	hidexdecorations!(ax)

	##-- B --##
    ax = Axis(
		f[2, 2]; 
		title = "False-Positive (High)",
		yticks = [0, 25, 50, 75, 100]
	)

    table = [1, 2, 3]
	h1 = (total_pos_md_high / total_cac_pos_high) * 100
	h2 = (total_pos_vf_high / total_cac_pos_high) * 100
	h3 = (total_pos_a_high / total_cac_pos_high) * 100
    heights1 = [h1, h2, h3]
	l1 = @sprintf "%.2f" h1
	l2 = @sprintf "%.2f" h2
	l3 = @sprintf "%.2f" h3
    barplot!(table, heights1; color=colors[eachindex(labels)], bar_labels=[l1, l2, l3])

    ylims!(ax; low=0, high=100)
	hidexdecorations!(ax)

	group_color = [PolyElement(color = color, strokecolor = :transparent)
    for color in colors[eachindex(labels)]]
	
	f[1:2, 3] = Legend(
		f,
		group_color,
		labels,
		framevisible = false
	)

	Label(
		f[1:2, 0],
		"Percentage (%)",
		fontsize = 12,
		rotation = pi/2
	)

	save(plotsdir("sensitivity_specificity_combined.png"), f)

    return f
end

# ╔═╡ 6dd9b880-0e55-46d9-a52a-66e16f1bc7d9
sensitivity_specificity_combined()

# ╔═╡ Cell order:
# ╠═024de588-307b-452d-812d-3a58b9265089
# ╠═3c5e1804-23dd-4629-baa5-5106aacd74a6
# ╠═d6965c35-ace0-4ff6-92f8-435853937946
# ╠═0c4ba5cc-8e40-40c6-b596-785973bbc74e
# ╠═a9c366d0-e90a-406c-b970-3d5a76bc7973
# ╠═11835515-0d70-4a4b-b7f0-742153e4c8e2
# ╠═65212870-6665-4ce8-b1f6-f80f9304bf04
# ╠═f1708909-5df5-477b-9657-82029aab0a50
# ╠═937d6c6d-80c7-4bfc-9c61-9d1830c84b66
# ╟─c69df360-2e83-4d92-89a5-beaa1a20a1c5
# ╠═2a8133c5-d7c0-4f85-8091-9a6c9d2aceb4
# ╠═a5935d6a-c825-4f8e-b559-c22104bb151e
# ╟─0a82a5ee-d54d-41fe-9cd3-336d8d9ca8e9
# ╠═ebec3a75-49ff-4ab0-b738-59cac17a001f
# ╠═3f9c0eb2-b727-4b9a-91eb-f9c482d6809f
# ╠═d3e09d9a-187f-4db0-a8c5-5fc0dc408d25
# ╠═f1575684-58ac-4606-92ef-993377c8309d
# ╟─82712845-c2ba-49f1-a8ad-9127cbc321b2
# ╠═a1289bbb-21b8-471e-9116-86f04f53c085
# ╠═ea5853a5-5cd8-4017-824c-c1a9c945b357
# ╠═66d07988-d98a-4bab-8243-33d6f4485cf0
# ╠═c93775d8-3422-4ad1-9f40-6d57477d8425
# ╠═636abeda-6c10-4d21-ae08-10e1274b77a1
# ╠═d0b1613e-890e-417d-b7eb-57b0bc9fafe2
# ╟─93679b54-9364-4579-8767-53e36499b6d8
# ╠═853d895d-91d4-487a-a927-e74424bba9ed
# ╠═d6a65460-51e0-4b18-83e5-1844f2566cc6
# ╠═e45adb9a-0519-4063-b8ae-2cb7e1005397
# ╟─8a774374-ecce-40d7-8b64-609124242d1f
# ╠═38bf3544-c7a4-46d6-9df2-6fe269b24cb8
# ╟─5588545b-f070-4e02-ac12-3e1b5a150d88
# ╠═20a32f24-2c4b-4308-a6e1-802a72322fa7
# ╟─c2db3887-3be1-499a-ba8c-d0a9f3f3f645
# ╟─7cfd0261-208d-41b6-b356-63d9e02ca13a
# ╠═cac42faf-836c-42a2-bed2-f0bba2a13be0
# ╠═bd61d48e-886b-40c6-a54a-76468c52cd4e
# ╟─6c01fb2c-c8cf-4484-8509-3fc540c038a0
# ╠═b410ca06-2fe0-4671-925e-bbdaf9786849
# ╠═e7e5385b-4016-4bb3-8fd8-c203f8e5f8d8
# ╟─d01812f8-e0ec-4b2d-9992-e7db2610d4fe
# ╠═397e2b5a-8d37-4398-b92b-678b4ba080ee
# ╟─e0c01bd6-0cec-408d-b534-cf36cb5bd78a
# ╠═6952fc87-c733-4e67-b1aa-ec9dd46ebe3d
# ╟─1a70a6e7-108d-477e-b1a5-7e755cbc0619
# ╠═d7e045fc-ac9c-4f53-836e-aedcadad1fbe
# ╠═95509262-a9a1-41b6-972d-e7dcf7bb2371
# ╠═b98db832-06bd-492f-8ca4-38e298086c8a
# ╠═3baf66b0-5e25-43d2-9b53-6b2d35f268a7
# ╟─c3919e4d-c474-4aeb-bdc8-76e139e1d8c2
# ╟─fef76b74-ddfe-4578-9ccc-a8ff68b60137
# ╠═6d84cb36-80cf-43af-b100-caab847b347a
# ╠═fd222d84-8b90-418e-b405-929557aaff8a
# ╠═a32f1676-cb1c-4b57-99a0-b6f7b542dbaf
# ╠═c9f166c4-1269-42c5-85fb-21e0786a2ec5
# ╠═2f742443-0dd4-4d25-90ff-4c7941cba0a2
# ╠═7b6def43-c729-4698-a752-9d37f182e816
# ╟─1cb3b222-e86d-41a5-98f0-d4faa9363e87
# ╟─c23d6948-d28a-4c66-ae2d-78686082e5f3
# ╟─6fb4896c-a086-4ef1-bbb7-1f1603338b3c
# ╠═7d218b4c-d467-4317-b239-a96de01a9d81
# ╟─d1d90206-1f7d-46b7-9063-8d6cd530a6e1
# ╠═82baca25-0d76-40a9-b19f-4838e6e84d72
# ╠═a3b18554-6f32-4502-8fd4-d8550ddbbf07
# ╟─f35612ba-8130-4c01-b91c-68886bccf358
# ╠═7873647b-d7f8-4d8c-89d7-f15c8d64653a
# ╠═5ca5e077-3900-4948-9136-8e1b907758b1
# ╟─ca5d95fe-ae81-41e5-9a54-ddc81962f2a3
# ╠═f1a8465b-aea5-46a5-ab95-ef181096a682
# ╟─2c0270a0-36e9-4b9d-855b-834aa443a808
# ╟─d9e33929-635a-450b-aa08-7ce0a6ee1170
# ╠═8a315d66-ff37-4e7d-9439-10567aabcd90
# ╠═603ae0fc-4fc5-4a2f-8dd7-e3b868bba417
# ╠═fe9b5fcb-1221-4898-976b-051b4896e470
# ╠═1d176c7b-802f-462d-94ca-97849df1e529
# ╠═ed9ec9d9-7262-4be2-99de-32c9fea8fc33
# ╠═9f3b22d1-e52a-4064-98e7-77856b5923bd
# ╟─708ea6a0-2129-42b9-b05e-b2766bcaefd7
# ╟─40854076-fcc6-481a-b1c8-12ded408f242
# ╟─2e0c8139-fe10-4b4e-b00d-c713208d8555
# ╟─eaab5224-3302-4dc4-9d09-c2728b547f23
# ╠═986e1b6a-5052-4a1d-b3f8-65d05e8bf3bf
# ╠═6da99d21-dde5-4e21-aae9-f798e6d54a2b
# ╟─812e634d-ff38-43b9-8a8b-8a21a80654b5
# ╠═7b0777fc-d104-480f-abd6-2f192ec539da
# ╠═77cf61eb-939e-42a7-917d-8fe4811522d5
# ╟─f41d1a38-08f6-4d3e-8932-ef4d6783e065
# ╟─cb0a3ad1-de3c-475d-a727-1dcbe2ffefb4
# ╟─39e42f15-de5c-445e-bd07-67745da578d7
# ╟─503fadca-83b4-444c-a3c9-a48fa6a67694
# ╟─0b0094a2-e267-4820-89e8-1cf95145975e
# ╟─c27811d2-d651-4198-bf90-246c4dafd8ad
# ╟─8060a585-765a-47ac-bc66-04808bd401f6
# ╟─7d7afb9b-2dbd-4373-b745-c245f86e48f7
# ╟─6dd9b880-0e55-46d9-a52a-66e16f1bc7d9
