### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ e5248e4b-76bc-4a45-b38d-2c7f55f8d8f7
using DrWatson

# ╔═╡ 67c7c5e9-e181-4560-b77a-50e9d7977348
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ 155ddbb8-f3a2-4db9-b256-7fd771785d66
using PlutoUI, CairoMakie, Statistics, CSV, DataFrames, DICOM, CSVFiles, Printf

# ╔═╡ 90830c78-6f24-41a3-a2eb-1fce3756bcaf
using StatsBase: quantile!, rmsd

# ╔═╡ 68ec4904-d908-4c1e-98e1-00ec7ba06f59
# ╠═╡ show_logs = false
using CalciumScoring

# ╔═╡ 00c7bf5b-58b1-4184-8ddc-b97a81eca736
# ╠═╡ show_logs = false
using MaterialDecomposition

# ╔═╡ d394450d-825d-4071-9102-239c1dd79b98
using GLM, MLJBase

# ╔═╡ 29f9cf29-9437-42a6-8dab-89105273c187
include(srcdir("masks.jl")); include(srcdir("dicom_utils.jl"));

# ╔═╡ a08b9637-f6d3-4671-9cfa-2a67b27f6585
include(srcdir("helper_functions.jl"));

# ╔═╡ 4b9aa903-625f-4d4e-89de-aa175b33cbf5
include(srcdir("plot_utils.jl"));

# ╔═╡ 1b568480-5467-4ffd-9099-de81066e407e
TableOfContents()

# ╔═╡ d1502f07-edc5-4baf-b7e9-808d37aeb6b3
md"""
# Low Density
"""

# ╔═╡ 668d1999-ce4e-4510-9dd4-84a26ab26dcf
begin
	sizes_cal = [30]
	energies = [80, 120, 135]
end;

# ╔═╡ b38de269-be3f-4e5d-adfb-43b1db438937
optimal_cal_low = [16, 25, 31, 45, 50, 75, 100, 150, 200, 250, 300, 350]

# ╔═╡ 66025838-3122-4316-ba66-35b89e9d510c
md"""
## Volume Fraction & Agatston Calibration
"""

# ╔═╡ 00e313f8-1f2f-4e16-ae14-85c1f082c843
path = joinpath(datadir("dcms", "cal", "100", "30"), string(energies[2]))

# ╔═╡ 17ced35c-de57-4633-be0f-754ccfbfe180
begin
	dcm_cal = dcmdir_parse(path)
	dcm_array_cal = load_dcm_array(dcm_cal)

	center_insert1, center_insert2 = 187, 318
	offset = 5
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

# ╔═╡ c02cc808-ac3b-479a-b5a1-9abb36b93a03
md"""
## Calibration
"""

# ╔═╡ b9753033-46ef-4509-8102-d4d294171257
densities_cal_all = [3, 6, 9, 10, 16, 25, 31, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800] # calcium densities

# ╔═╡ c30d5213-abcc-482d-887a-a63b54eed243
begin
	densities_cal_low = optimal_cal_low
	
	# Low Energy
	means_80_low = Dict(:density => densities_cal_low, :means => zeros(length(densities_cal_low)))
	for (i, density) in enumerate(densities_cal_low)
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

	# High Energy
	means_135_low = Dict(:density => densities_cal_low, :means => zeros(length(densities_cal_low)))
	for (i, density) in enumerate(densities_cal_low)
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

	# Fit Parameters
	calculated_intensities_low = hcat(means_80_low[:means], means_135_low[:means])
	ps_low = fit_calibration(calculated_intensities_low, densities_cal_low)
end

# ╔═╡ 905ccae5-0ee5-4242-b58a-25ef4bcbb8b9
md"""
### Check Results
"""

# ╔═╡ 2402b78f-580b-47eb-962f-707c90d6e74c
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

# ╔═╡ 896312ad-3284-487f-b549-bd34ac8b1f67
calculated_intensities_low

# ╔═╡ c0095d39-14f2-464d-8f55-ea1739a7b313
df_low = DataFrame(
	densities = densities_cal_low,
	predicted_densities = predicted_densities_low,
)

# ╔═╡ cb0fc64e-8e36-4250-9ea6-427b5dccf27d
md"""
## Validation
"""

# ╔═╡ 0c05f4a7-45ad-4dc3-8d1c-cfd4a48fce95
densities_val_all = [
	"15_18_22"
	"26_29_36"
	"52_59_73"
	"110_210_310"
	"410_610_780"
] 

# ╔═╡ af661a05-e9ec-4b0b-8d7e-0ce39af016ec
begin
	densities_val_low = densities_val_all[1:3]
	sizes_val = ["small", "medium", "large"]
end;

# ╔═╡ 06780785-7a4d-4e8f-bc3e-7163ac6374f4
energies_val = [80, 120, 135]

# ╔═╡ 4c2bf021-6a1d-4c8c-a3f4-7aee6c39c361
md"""
### Load masks
"""

# ╔═╡ 0029b07b-6c4f-49f3-841b-3e75e2f2a34f
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

# ╔═╡ e35e963d-1b3f-4d6e-8598-ea6c104935c2
md"""
### Predict
"""

# ╔═╡ 42fd7e43-96fe-4c85-9807-9239cc66f10b
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

# ╔═╡ 0666cb65-8c79-4766-a4f4-4d57797b98f3
md"""
# High Density
"""

# ╔═╡ feddd76b-bb09-4f05-894f-a89b1bcfbcff
md"""
## Calibration
"""

# ╔═╡ 2e882999-ded4-43f5-aa45-66f6438dc641
optimal_cal_high = [10, 25, 45, 75, 150, 250, 350, 450, 550, 650]

# ╔═╡ b63a9eee-7d14-4b4a-92b5-6b960aad90c1
begin
	densities_cal_high = optimal_cal_high

	# Low Energy
	means_80_high = Dict(:density => densities_cal_high, :means => zeros(length(densities_cal_high)))
	for (i, density) in enumerate(densities_cal_high)
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

	# High energy	
	means_135_high = Dict(:density => densities_cal_high, :means => zeros(length(densities_cal_high)))
	for (i, density) in enumerate(densities_cal_high)
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

	# Fit Parameters
	calculated_intensities_high = hcat(means_80_high[:means], means_135_high[:means]) # low energy, high energy
		
	ps_high = fit_calibration(calculated_intensities_high, densities_cal_high)
end

# ╔═╡ 776b150d-ddbe-494a-9e53-632a33bb502b
md"""
### Check Results
"""

# ╔═╡ 3d45ea72-bc18-4651-8a0d-286294d1725e
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

# ╔═╡ 7f0eb143-a2bf-4623-9d3e-9698f03f5311
df_high = DataFrame(
	densities = densities_cal_high,
	predicted_densities = predicted_densities_high,
)

# ╔═╡ 397c56dc-79c0-4a18-8a04-fded2676c658
md"""
## Validation
"""

# ╔═╡ ccce7595-bc08-4855-819b-9dd3e9e88198
densities_val_high = densities_val_all[length(densities_val_low):end];

# ╔═╡ ca68557b-b343-4f1c-96cb-5cccb10636ad
md"""
### Predict
"""

# ╔═╡ 044e548f-f81e-4cf8-a872-8e5440e74ddb
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

# ╔═╡ 5ffcf3a8-e3d5-4093-8794-7eb133252d04
md"""
# Analyze
"""

# ╔═╡ c3f49d87-b42f-4e72-9f8f-02c8aa72af4f
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

# ╔═╡ 1b6ce948-c725-46e8-a414-c64a7ca49239
md"""
## Low Density
"""

# ╔═╡ 4355944f-2562-48d1-bddc-82cb6e23472c
md"""
### Accuracy
"""

# ╔═╡ ca81d1ea-252a-4d43-aab8-4a08f95176dc
new_df_low_md = vcat(dfs_low_md[1:length(dfs_low_md)]...);

# ╔═╡ b0af4257-b1a9-448d-83ec-fd5494595aa2
new_df_low_vf = vcat(dfs_low_vf[1:length(dfs_low_vf)]...);

# ╔═╡ 9efebaa2-458a-4d5a-a44d-f01f1dddac18
new_df_low_a = vcat(dfs_low_a[1:length(dfs_low_a)]...);

# ╔═╡ ed6a9d9d-ce6a-4eab-955b-c6e676103603
co_1_low_md, r_squared_1_low_md, rms_values_1_low_md, pred_1_low_md = calculate_coefficients(new_df_low_md);

# ╔═╡ 073282e1-f842-41a0-98dc-c1b4bf8959b0
co_1_low_vf, r_squared_1_low_vf, rms_values_1_low_vf, pred_1_low_vf = calculate_coefficients(new_df_low_vf);

# ╔═╡ d3459d52-d9a0-492b-b540-342ea732486b
co_1_low_a, r_squared_1_low_a, rms_values_1_low_a, pred_1_low_a = calculate_coefficients(new_df_low_a);

# ╔═╡ 5c002c0a-10cf-4eb3-bac1-08ac25aab4cf
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

# ╔═╡ 41387ee8-07fa-4a54-bb1d-cbe78d56e08b
with_theme(accuracy_low, medphys_theme)

# ╔═╡ 3387dc7d-4c94-4bf2-a917-fa9a576b1661
md"""
### Sensitivity & Specificity
"""

# ╔═╡ d8f08ec3-bddc-4846-bd55-381373a51301
std_level = 1.5

# ╔═╡ a2bf8116-6ddb-456e-b0a6-2a1145ac3cc9
md"""
#### False-Negatives
"""

# ╔═╡ 3e424011-f270-45d5-832d-117c2bd2f47c
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

# ╔═╡ a1fb39ec-e333-4c35-9c72-29281becf9f4
total_zero_md, total_zero_vf, total_zero_a

# ╔═╡ d2be0a3d-7d2a-469f-8beb-9806db8d2a92
md"""
#### False-Positives
"""

# ╔═╡ 00669225-9347-46c8-9e03-56c49205d0d4
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

# ╔═╡ 34a76785-5316-4d9b-8734-6b9f8ff429bc
total_pos_md, total_pos_vf, total_pos_a

# ╔═╡ 120a6073-5447-4377-816f-cff88d153e58
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

# ╔═╡ fe51bc0a-78f9-4508-8171-0a056db887a4
with_theme(sensitivity_specificity, medphys_theme)

# ╔═╡ a2261f23-10a8-4bc1-b5cb-c96da65e4b25
md"""
## High Density
"""

# ╔═╡ b78bace6-e275-4950-8608-35888c339e12
md"""
### Accuracy
"""

# ╔═╡ 066d0e67-a4db-4acc-b1c8-9fcb6fd76e88
new_df_high_md = vcat(dfs_high_md[1:length(dfs_high_md)]...);

# ╔═╡ 0263f7e5-e98d-44e1-8dd5-5d1d4380e296
new_df_high_vf = vcat(dfs_high_vf[1:length(dfs_high_vf)]...);

# ╔═╡ 4a40ce3f-40dd-4d9d-af16-2bebebe9cb68
new_df_high_a = vcat(dfs_high_a[1:length(dfs_high_a)]...);

# ╔═╡ 4cadf68d-5cce-4331-8e96-8757a42667ef
co_1_high_md, r_squared_1_high_md, rms_values_1_high_md, pred_1_high_md = calculate_coefficients(new_df_high_md);

# ╔═╡ 3a799ff0-593c-4754-8180-c79c8a9c40cf
co_1_high_vf, r_squared_1_high_vf, rms_values_1_high_vf, pred_1_high_vf = calculate_coefficients(new_df_high_vf);

# ╔═╡ 6a4c5316-459d-4791-ae18-c9c756c9145a
co_1_high_a, r_squared_1_high_a, rms_values_1_high_a, pred_1_high_a = calculate_coefficients(new_df_high_a);

# ╔═╡ 2eccc861-3b08-48a4-942a-e51feadb8e89
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

# ╔═╡ 323e72a2-6bd9-4398-9fe3-4d8cb4347979
with_theme(accuracy_high, medphys_theme)

# ╔═╡ cd48d27b-3d24-499f-baec-7a5bfd3338e7
md"""
### Sensitivity & Specificity
"""

# ╔═╡ 30a6f69f-ca39-4a49-a37d-2f504e9f9d14
md"""
#### False-Negatives
"""

# ╔═╡ b8c2e2cc-10c8-4096-abf7-bd352d3e4cf0
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

# ╔═╡ 0869a394-8065-4620-b337-49d01fa3202d
total_zero_md_high, total_zero_vf_high, total_zero_a_high

# ╔═╡ e90a1cb1-d104-4a21-8837-47a9764a96a8
md"""
#### False-Positives
"""

# ╔═╡ 78103305-4811-4178-a70c-1ceb7101d8d0
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

# ╔═╡ 4f784cc2-7dd9-4a09-9eb0-e3ae2a012f88
total_pos_md, total_pos_vf, total_pos_a

# ╔═╡ 476630f6-6630-4cf3-9a9e-a415423eefd7
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

# ╔═╡ 9c0dd631-6575-4bec-bc01-66d0579ef1a9
with_theme(sensitivity_specificity_high, medphys_theme)

# ╔═╡ 0c842b80-8caf-4251-baa1-f99da406dc50
md"""
## Combined
"""

# ╔═╡ 7f9b5200-487b-4788-b6fe-e2ef810dc61a
md"""
### Accuracy
"""

# ╔═╡ 1c660d27-f650-4cef-871f-de28927721d0
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

# ╔═╡ 306e9179-23e3-4d06-87d0-71c91bf33367
accuracy_combined()

# ╔═╡ 3f672521-eda4-42f7-8afe-3e6ebd85cdd1
md"""
### Sensitivity & Specificity
"""

# ╔═╡ c7ab1789-59d0-4307-b09d-44a2f33cfb43
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

# ╔═╡ 85c7f5d9-fac4-4c36-97a6-3fa9d7f7ca0d
sensitivity_specificity_combined()

# ╔═╡ Cell order:
# ╠═e5248e4b-76bc-4a45-b38d-2c7f55f8d8f7
# ╠═67c7c5e9-e181-4560-b77a-50e9d7977348
# ╠═155ddbb8-f3a2-4db9-b256-7fd771785d66
# ╠═90830c78-6f24-41a3-a2eb-1fce3756bcaf
# ╠═68ec4904-d908-4c1e-98e1-00ec7ba06f59
# ╠═00c7bf5b-58b1-4184-8ddc-b97a81eca736
# ╠═29f9cf29-9437-42a6-8dab-89105273c187
# ╠═1b568480-5467-4ffd-9099-de81066e407e
# ╟─d1502f07-edc5-4baf-b7e9-808d37aeb6b3
# ╠═668d1999-ce4e-4510-9dd4-84a26ab26dcf
# ╠═b38de269-be3f-4e5d-adfb-43b1db438937
# ╟─66025838-3122-4316-ba66-35b89e9d510c
# ╠═00e313f8-1f2f-4e16-ae14-85c1f082c843
# ╠═17ced35c-de57-4633-be0f-754ccfbfe180
# ╟─c02cc808-ac3b-479a-b5a1-9abb36b93a03
# ╠═b9753033-46ef-4509-8102-d4d294171257
# ╠═c30d5213-abcc-482d-887a-a63b54eed243
# ╟─905ccae5-0ee5-4242-b58a-25ef4bcbb8b9
# ╠═2402b78f-580b-47eb-962f-707c90d6e74c
# ╠═896312ad-3284-487f-b549-bd34ac8b1f67
# ╠═c0095d39-14f2-464d-8f55-ea1739a7b313
# ╟─cb0fc64e-8e36-4250-9ea6-427b5dccf27d
# ╠═0c05f4a7-45ad-4dc3-8d1c-cfd4a48fce95
# ╠═af661a05-e9ec-4b0b-8d7e-0ce39af016ec
# ╠═06780785-7a4d-4e8f-bc3e-7163ac6374f4
# ╟─4c2bf021-6a1d-4c8c-a3f4-7aee6c39c361
# ╠═0029b07b-6c4f-49f3-841b-3e75e2f2a34f
# ╟─e35e963d-1b3f-4d6e-8598-ea6c104935c2
# ╠═42fd7e43-96fe-4c85-9807-9239cc66f10b
# ╟─0666cb65-8c79-4766-a4f4-4d57797b98f3
# ╟─feddd76b-bb09-4f05-894f-a89b1bcfbcff
# ╠═2e882999-ded4-43f5-aa45-66f6438dc641
# ╠═b63a9eee-7d14-4b4a-92b5-6b960aad90c1
# ╟─776b150d-ddbe-494a-9e53-632a33bb502b
# ╠═3d45ea72-bc18-4651-8a0d-286294d1725e
# ╠═7f0eb143-a2bf-4623-9d3e-9698f03f5311
# ╟─397c56dc-79c0-4a18-8a04-fded2676c658
# ╠═ccce7595-bc08-4855-819b-9dd3e9e88198
# ╟─ca68557b-b343-4f1c-96cb-5cccb10636ad
# ╠═044e548f-f81e-4cf8-a872-8e5440e74ddb
# ╟─5ffcf3a8-e3d5-4093-8794-7eb133252d04
# ╠═d394450d-825d-4071-9102-239c1dd79b98
# ╠═a08b9637-f6d3-4671-9cfa-2a67b27f6585
# ╠═4b9aa903-625f-4d4e-89de-aa175b33cbf5
# ╠═c3f49d87-b42f-4e72-9f8f-02c8aa72af4f
# ╟─1b6ce948-c725-46e8-a414-c64a7ca49239
# ╟─4355944f-2562-48d1-bddc-82cb6e23472c
# ╠═ca81d1ea-252a-4d43-aab8-4a08f95176dc
# ╠═b0af4257-b1a9-448d-83ec-fd5494595aa2
# ╠═9efebaa2-458a-4d5a-a44d-f01f1dddac18
# ╠═ed6a9d9d-ce6a-4eab-955b-c6e676103603
# ╠═073282e1-f842-41a0-98dc-c1b4bf8959b0
# ╠═d3459d52-d9a0-492b-b540-342ea732486b
# ╟─5c002c0a-10cf-4eb3-bac1-08ac25aab4cf
# ╟─41387ee8-07fa-4a54-bb1d-cbe78d56e08b
# ╟─3387dc7d-4c94-4bf2-a917-fa9a576b1661
# ╠═d8f08ec3-bddc-4846-bd55-381373a51301
# ╟─a2bf8116-6ddb-456e-b0a6-2a1145ac3cc9
# ╠═3e424011-f270-45d5-832d-117c2bd2f47c
# ╠═a1fb39ec-e333-4c35-9c72-29281becf9f4
# ╟─d2be0a3d-7d2a-469f-8beb-9806db8d2a92
# ╠═00669225-9347-46c8-9e03-56c49205d0d4
# ╠═34a76785-5316-4d9b-8734-6b9f8ff429bc
# ╟─120a6073-5447-4377-816f-cff88d153e58
# ╟─fe51bc0a-78f9-4508-8171-0a056db887a4
# ╟─a2261f23-10a8-4bc1-b5cb-c96da65e4b25
# ╟─b78bace6-e275-4950-8608-35888c339e12
# ╠═066d0e67-a4db-4acc-b1c8-9fcb6fd76e88
# ╠═0263f7e5-e98d-44e1-8dd5-5d1d4380e296
# ╠═4a40ce3f-40dd-4d9d-af16-2bebebe9cb68
# ╠═4cadf68d-5cce-4331-8e96-8757a42667ef
# ╠═3a799ff0-593c-4754-8180-c79c8a9c40cf
# ╠═6a4c5316-459d-4791-ae18-c9c756c9145a
# ╟─2eccc861-3b08-48a4-942a-e51feadb8e89
# ╟─323e72a2-6bd9-4398-9fe3-4d8cb4347979
# ╟─cd48d27b-3d24-499f-baec-7a5bfd3338e7
# ╟─30a6f69f-ca39-4a49-a37d-2f504e9f9d14
# ╠═b8c2e2cc-10c8-4096-abf7-bd352d3e4cf0
# ╠═0869a394-8065-4620-b337-49d01fa3202d
# ╟─e90a1cb1-d104-4a21-8837-47a9764a96a8
# ╠═78103305-4811-4178-a70c-1ceb7101d8d0
# ╠═4f784cc2-7dd9-4a09-9eb0-e3ae2a012f88
# ╟─476630f6-6630-4cf3-9a9e-a415423eefd7
# ╟─9c0dd631-6575-4bec-bc01-66d0579ef1a9
# ╟─0c842b80-8caf-4251-baa1-f99da406dc50
# ╟─7f9b5200-487b-4788-b6fe-e2ef810dc61a
# ╟─1c660d27-f650-4cef-871f-de28927721d0
# ╟─306e9179-23e3-4d06-87d0-71c91bf33367
# ╟─3f672521-eda4-42f7-8afe-3e6ebd85cdd1
# ╟─c7ab1789-59d0-4307-b09d-44a2f33cfb43
# ╟─85c7f5d9-fac4-4c36-97a6-3fa9d7f7ca0d
