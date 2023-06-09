### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e5248e4b-76bc-4a45-b38d-2c7f55f8d8f7
using DrWatson

# ╔═╡ 67c7c5e9-e181-4560-b77a-50e9d7977348
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ 155ddbb8-f3a2-4db9-b256-7fd771785d66
using PlutoUI, CairoMakie, Statistics, CSV, DataFrames, DICOM, CSVFiles

# ╔═╡ 90830c78-6f24-41a3-a2eb-1fce3756bcaf
using StatsBase: quantile!, rmsd

# ╔═╡ 68ec4904-d908-4c1e-98e1-00ec7ba06f59
# ╠═╡ show_logs = false
using CalciumScoring

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

# ╔═╡ c02cc808-ac3b-479a-b5a1-9abb36b93a03
md"""
## Calibration
"""

# ╔═╡ b9753033-46ef-4509-8102-d4d294171257
densities_cal_all = [3, 6, 9, 10, 16, 25, 31, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800] # calcium densities

# ╔═╡ aefb07e8-0023-45d0-ae74-c33a665c40a2
densities_cal_low = densities_cal_all[1:12]

# ╔═╡ a38fee20-1c9c-41ce-991a-a78a5354474c
densities_cal_high = densities_cal_all[length(densities_cal_low)-5:end-10]

# ╔═╡ ee95b544-2e02-4de3-bd1c-019d8fc3cdd6
md"""
### Low Energy
"""

# ╔═╡ 668d1999-ce4e-4510-9dd4-84a26ab26dcf
begin
	sizes_cal = [30]
	energies = [80, 135]
end;

# ╔═╡ c30d5213-abcc-482d-887a-a63b54eed243
begin
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
end

# ╔═╡ 37181932-1398-4c45-9ff0-fb8d488ab562
md"""
### High Energy
"""

# ╔═╡ 84ab3b1d-cd05-43b1-b889-fdcf94529023
begin
	means_135_low = Dict(:density => densities_cal_low, :means => zeros(length(densities_cal_low)))
	for (i, density) in enumerate(densities_cal_low)
		for _size in sizes_cal
			path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[2]))
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

# ╔═╡ 6760a907-7817-4692-abd3-2c1b66dc3a50
md"""
### Fit Parameters
"""

# ╔═╡ 727de9ea-d78c-4c5f-9223-351a2fac45e3
calculated_intensities_low = hcat(means_80_low[:means], means_135_low[:means]) # low energy, high energy

# ╔═╡ 38910d55-20d0-445e-98ea-c2f36aaa5255
ps_low = fit_calibration(calculated_intensities_low, densities_cal_low)

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
			score(calculated_intensities_low[i, 1], calculated_intensities_low[i, 2], ps_low, MaterialDecomposition()
			)
		)
	end
end

# ╔═╡ c0095d39-14f2-464d-8f55-ea1739a7b313
df = DataFrame(
	densities = densities_cal_low,
	predicted_densities = predicted_densities_low,
	mean_intensities_low = means_80_low[:means],
	mean_intensities_high = means_135_low[:means],
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

# ╔═╡ 42fd7e43-96fe-4c85-9807-9239cc66f10b
begin
	dfs_low = []
	for density in densities_val_low
		for _size in sizes_val
			if _size == "small"
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

			dilated_mask_L_HD = dilate_recursively(mask_L_HD, 2)
			dilated_mask_L_HD_3D = cat(dilated_mask_L_HD, dilated_mask_L_HD, dilated_mask_L_HD, dims=3)

			dilated_mask_L_MD = dilate_recursively(mask_L_MD, 2)
			dilated_mask_L_MD_3D = cat(dilated_mask_L_MD, dilated_mask_L_MD, dilated_mask_L_MD, dims=3)

			dilated_mask_L_LD = dilate_recursively(mask_L_LD, 2)
			dilated_mask_L_LD_3D = cat(dilated_mask_L_LD, dilated_mask_L_LD, dilated_mask_L_LD, dims=3)
			

			dilated_mask_M_HD = dilate_recursively(mask_M_HD, 2)
			dilated_mask_M_HD_3D = cat(dilated_mask_M_HD, dilated_mask_M_HD, dilated_mask_M_HD, dims=3)

			dilated_mask_M_MD = dilate_recursively(mask_M_MD, 2)
			dilated_mask_M_MD_3D = cat(dilated_mask_M_MD, dilated_mask_M_MD, dilated_mask_M_MD, dims=3)

			dilated_mask_M_LD = dilate_recursively(mask_M_LD, 2)
			dilated_mask_M_LD_3D = cat(dilated_mask_M_LD, dilated_mask_M_LD, dilated_mask_M_LD, dims=3)
			

			dilated_mask_S_HD = dilate_recursively(mask_S_HD, 2)
			dilated_mask_S_HD_3D = cat(dilated_mask_S_HD, dilated_mask_S_HD, dilated_mask_S_HD, dims=3)

			dilated_mask_S_MD = dilate_recursively(mask_S_MD, 2)
			dilated_mask_S_MD_3D = cat(dilated_mask_S_MD, dilated_mask_S_MD, dilated_mask_S_MD, dims=3)

			dilated_mask_S_LD = dilate_recursively(mask_S_LD, 2)
			dilated_mask_S_LD_3D = cat(dilated_mask_S_LD, dilated_mask_S_LD, dilated_mask_S_LD, dims=3)

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
			path_135 = datadir("dcms", "val", density, _size, string(energies[2]))
			dcm_135 = dcmdir_parse(path_135)
			dcm_array_135 = load_dcm_array(dcm_135)
			
			means_135 = [
				mean(dcm_array_135[dilated_mask_L_HD_3D]), mean(dcm_array_135[dilated_mask_L_MD_3D]), mean(dcm_array_135[dilated_mask_L_LD_3D]),
				mean(dcm_array_135[dilated_mask_M_HD_3D]), mean(dcm_array_135[dilated_mask_M_MD_3D]), mean(dcm_array_135[dilated_mask_M_LD_3D]),
				mean(dcm_array_135[dilated_mask_S_HD_3D]), mean(dcm_array_135[dilated_mask_S_MD_3D]), mean(dcm_array_135[dilated_mask_S_LD_3D])
			]

			calculated_intensities = hcat(means_80, means_135)
			predicted_densities = zeros(size(means_135))
			for i in eachindex(predicted_densities)
				predicted_densities[i] = score(means_80[i], means_80[i], ps_low, MaterialDecomposition())
			end
			
			## Calculate predicted mass
			voxel_size_mm3 = (pixel_size[1] * pixel_size[2] * pixel_size[3]) # mm^3
			voxel_size_cm3 = voxel_size_mm3 * 1e-3 # cm^3
			vol_large = count(dilated_mask_L_HD_3D) * voxel_size_cm3
			vol_medium = count(dilated_mask_M_HD_3D) * voxel_size_cm3
			vol_small = count(dilated_mask_S_HD_3D) * voxel_size_cm3
			
			predicted_mass_large_inserts = predicted_densities[1:3] .* vol_large
			predicted_mass_medium_inserts = predicted_densities[4:6] .* vol_medium
			predicted_mass_small_inserts = predicted_densities[7:9] .* vol_small

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

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				gt_mass_large_inserts = gt_mass_large_inserts,
				predicted_mass_large_inserts = predicted_mass_large_inserts,
				gt_mass_medium_inserts = gt_mass_medium_inserts,
				predicted_mass_medium_inserts = predicted_mass_medium_inserts,
				gt_mass_small_inserts = gt_mass_small_inserts,
				predicted_mass_small_inserts = predicted_mass_small_inserts,
			)
			push!(dfs_low, df_results)
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

# ╔═╡ 124b7573-113c-4cdf-9a5c-e0cc9c78dfaa
md"""
### Low Energy
"""

# ╔═╡ b63a9eee-7d14-4b4a-92b5-6b960aad90c1
begin
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
end

# ╔═╡ 6bc5781a-d025-46d8-bd70-19cb374d3c4e
md"""
### High Energy
"""

# ╔═╡ 390c73cc-8253-4547-afb3-c225b58fb2c4
begin
	means_135_high = Dict(:density => densities_cal_high, :means => zeros(length(densities_cal_high)))
	for (i, density) in enumerate(densities_cal_high)
		for _size in sizes_cal
			path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[2]))
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

# ╔═╡ 07e1b431-84c6-4cfc-ad2b-edc54fdc59ca
md"""
### Fit Parameters
"""

# ╔═╡ 565eade7-574f-473a-a9a4-6aeca2a812e6
calculated_intensities_high = hcat(means_80_high[:means], means_135_high[:means]) # low energy, high energy

# ╔═╡ da5d5ca7-5e6d-425c-823b-881c59ff8a3f
ps_high = fit_calibration(calculated_intensities_high, densities_cal_high)

# ╔═╡ 776b150d-ddbe-494a-9e53-632a33bb502b
md"""
### Check Results
"""

# ╔═╡ 7accb307-22f1-4035-951b-14bd04bcf5b2
begin
	predicted_densities_high = []
	
	for i in 1:length(densities_cal_high)
		append!(
			predicted_densities_high, 
			score(calculated_intensities_high[i, 1], calculated_intensities_high[i, 2], ps_high, MaterialDecomposition()
			)
		)
	end
end

# ╔═╡ 3226a396-db02-4f50-ac67-9f8c57e7566b
df_high = DataFrame(
	densities = densities_cal_high,
	predicted_densities = predicted_densities_high,
	mean_intensities_low = means_80_high[:means],
	mean_intensities_high = means_135_high[:means],
)

# ╔═╡ 397c56dc-79c0-4a18-8a04-fded2676c658
md"""
## Validation
"""

# ╔═╡ ccce7595-bc08-4855-819b-9dd3e9e88198
densities_val_high = densities_val_all[length(densities_val_low):end];

# ╔═╡ 044e548f-f81e-4cf8-a872-8e5440e74ddb
begin
	dfs_high = []
	for density in densities_val_high
		for _size in sizes_val
			if _size == "small"
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

			dilated_mask_L_HD = dilate_recursively(mask_L_HD, 2)
			dilated_mask_L_HD_3D = cat(dilated_mask_L_HD, dilated_mask_L_HD, dilated_mask_L_HD, dims=3)

			dilated_mask_L_MD = dilate_recursively(mask_L_MD, 2)
			dilated_mask_L_MD_3D = cat(dilated_mask_L_MD, dilated_mask_L_MD, dilated_mask_L_MD, dims=3)

			dilated_mask_L_LD = dilate_recursively(mask_L_LD, 2)
			dilated_mask_L_LD_3D = cat(dilated_mask_L_LD, dilated_mask_L_LD, dilated_mask_L_LD, dims=3)
			

			dilated_mask_M_HD = dilate_recursively(mask_M_HD, 2)
			dilated_mask_M_HD_3D = cat(dilated_mask_M_HD, dilated_mask_M_HD, dilated_mask_M_HD, dims=3)

			dilated_mask_M_MD = dilate_recursively(mask_M_MD, 2)
			dilated_mask_M_MD_3D = cat(dilated_mask_M_MD, dilated_mask_M_MD, dilated_mask_M_MD, dims=3)

			dilated_mask_M_LD = dilate_recursively(mask_M_LD, 2)
			dilated_mask_M_LD_3D = cat(dilated_mask_M_LD, dilated_mask_M_LD, dilated_mask_M_LD, dims=3)
			

			dilated_mask_S_HD = dilate_recursively(mask_S_HD, 2)
			dilated_mask_S_HD_3D = cat(dilated_mask_S_HD, dilated_mask_S_HD, dilated_mask_S_HD, dims=3)

			dilated_mask_S_MD = dilate_recursively(mask_S_MD, 2)
			dilated_mask_S_MD_3D = cat(dilated_mask_S_MD, dilated_mask_S_MD, dilated_mask_S_MD, dims=3)

			dilated_mask_S_LD = dilate_recursively(mask_S_LD, 2)
			dilated_mask_S_LD_3D = cat(dilated_mask_S_LD, dilated_mask_S_LD, dilated_mask_S_LD, dims=3)

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
			path_135 = datadir("dcms", "val", density, _size, string(energies[2]))
			dcm_135 = dcmdir_parse(path_135)
			dcm_array_135 = load_dcm_array(dcm_135)
			
			means_135 = [
				mean(dcm_array_135[dilated_mask_L_HD_3D]), mean(dcm_array_135[dilated_mask_L_MD_3D]), mean(dcm_array_135[dilated_mask_L_LD_3D]),
				mean(dcm_array_135[dilated_mask_M_HD_3D]), mean(dcm_array_135[dilated_mask_M_MD_3D]), mean(dcm_array_135[dilated_mask_M_LD_3D]),
				mean(dcm_array_135[dilated_mask_S_HD_3D]), mean(dcm_array_135[dilated_mask_S_MD_3D]), mean(dcm_array_135[dilated_mask_S_LD_3D])
			]

			calculated_intensities = hcat(means_80, means_135)
			predicted_densities = zeros(size(means_135))
			for i in eachindex(predicted_densities)
				predicted_densities[i] = score(means_80[i], means_80[i], ps_high, MaterialDecomposition())
			end
			
			## Calculate predicted mass
			voxel_size_mm3 = (pixel_size[1] * pixel_size[2] * pixel_size[3]) # mm^3
			voxel_size_cm3 = voxel_size_mm3 * 1e-3 # cm^3
			vol_large = count(dilated_mask_L_HD_3D) * voxel_size_cm3
			vol_medium = count(dilated_mask_M_HD_3D) * voxel_size_cm3
			vol_small = count(dilated_mask_S_HD_3D) * voxel_size_cm3
			
			predicted_mass_large_inserts = predicted_densities[1:3] .* vol_large
			predicted_mass_medium_inserts = predicted_densities[4:6] .* vol_medium
			predicted_mass_small_inserts = predicted_densities[7:9] .* vol_small

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

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				gt_mass_large_inserts = gt_mass_large_inserts,
				predicted_mass_large_inserts = predicted_mass_large_inserts,
				gt_mass_medium_inserts = gt_mass_medium_inserts,
				predicted_mass_medium_inserts = predicted_mass_medium_inserts,
				gt_mass_small_inserts = gt_mass_small_inserts,
				predicted_mass_small_inserts = predicted_mass_small_inserts,
			)
			push!(dfs_high, df_results)
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

# ╔═╡ ca81d1ea-252a-4d43-aab8-4a08f95176dc
new_df_low = vcat(dfs_low[1:length(dfs_low)]...);

# ╔═╡ ed6a9d9d-ce6a-4eab-955b-c6e676103603
co_1_low, r_squared_1_low, rms_values_1_low, pred_1_low = calculate_coefficients(new_df_low);

# ╔═╡ 066d0e67-a4db-4acc-b1c8-9fcb6fd76e88
new_df_high = vcat(dfs_high[1:length(dfs_high)]...);

# ╔═╡ 4cadf68d-5cce-4331-8e96-8757a42667ef
co_1_high, r_squared_1_high, rms_values_1_high, pred_1_high = calculate_coefficients(new_df_high);

# ╔═╡ 5c002c0a-10cf-4eb3-bac1-08ac25aab4cf
function accuracy()
	f = Figure()

	##-- A --##
	ax = Axis(
		f[1, 1],
		xticks = collect(-5:5:25),
		yticks = collect(-5:5:25),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Material Decomposition (Low Density)",
	)
	
	df = new_df_low
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
	ln2 = lines!(collect(1:1000), pred_1_low, linestyle=:dashdot)
	create_textbox(f[1, 1], co_1_low, r_squared_1_low, rms_values_1_low)
	
	xlims!(ax, low=-5, high=25)
	ylims!(ax, low=-5, high=25)

	ax = Axis(
		f[2, 1],
		xticks = collect(-50:50:250),
		yticks = collect(-50:50:250),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Material Decomposition (High Density)",
	)
	
	df = new_df_high
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
	ln2 = lines!(collect(1:1000), pred_1_high, linestyle=:dashdot)
	create_textbox(f[2, 1], co_1_high, r_squared_1_high, rms_values_1_high)
	
	xlims!(ax, low=-50, high=250)
	ylims!(ax, low=-50, high=250)

	#-- LABELS --##
	f[1:2, 2] = Legend(f, [sc1, sc2, sc3, ln1, ln2], ["Large Inserts", "Medium Inserts", "Small Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip([["A"], ["B"]], [f[1,1], f[2, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end
	f
end

# ╔═╡ 41387ee8-07fa-4a54-bb1d-cbe78d56e08b
with_theme(accuracy, medphys_theme)

# ╔═╡ Cell order:
# ╠═e5248e4b-76bc-4a45-b38d-2c7f55f8d8f7
# ╠═67c7c5e9-e181-4560-b77a-50e9d7977348
# ╠═155ddbb8-f3a2-4db9-b256-7fd771785d66
# ╠═90830c78-6f24-41a3-a2eb-1fce3756bcaf
# ╠═68ec4904-d908-4c1e-98e1-00ec7ba06f59
# ╠═29f9cf29-9437-42a6-8dab-89105273c187
# ╠═1b568480-5467-4ffd-9099-de81066e407e
# ╟─d1502f07-edc5-4baf-b7e9-808d37aeb6b3
# ╟─c02cc808-ac3b-479a-b5a1-9abb36b93a03
# ╠═b9753033-46ef-4509-8102-d4d294171257
# ╠═aefb07e8-0023-45d0-ae74-c33a665c40a2
# ╠═a38fee20-1c9c-41ce-991a-a78a5354474c
# ╟─ee95b544-2e02-4de3-bd1c-019d8fc3cdd6
# ╠═668d1999-ce4e-4510-9dd4-84a26ab26dcf
# ╠═c30d5213-abcc-482d-887a-a63b54eed243
# ╟─37181932-1398-4c45-9ff0-fb8d488ab562
# ╠═84ab3b1d-cd05-43b1-b889-fdcf94529023
# ╟─6760a907-7817-4692-abd3-2c1b66dc3a50
# ╠═727de9ea-d78c-4c5f-9223-351a2fac45e3
# ╠═38910d55-20d0-445e-98ea-c2f36aaa5255
# ╟─905ccae5-0ee5-4242-b58a-25ef4bcbb8b9
# ╠═2402b78f-580b-47eb-962f-707c90d6e74c
# ╠═c0095d39-14f2-464d-8f55-ea1739a7b313
# ╟─cb0fc64e-8e36-4250-9ea6-427b5dccf27d
# ╠═0c05f4a7-45ad-4dc3-8d1c-cfd4a48fce95
# ╠═af661a05-e9ec-4b0b-8d7e-0ce39af016ec
# ╟─4c2bf021-6a1d-4c8c-a3f4-7aee6c39c361
# ╠═0029b07b-6c4f-49f3-841b-3e75e2f2a34f
# ╠═42fd7e43-96fe-4c85-9807-9239cc66f10b
# ╟─0666cb65-8c79-4766-a4f4-4d57797b98f3
# ╟─feddd76b-bb09-4f05-894f-a89b1bcfbcff
# ╟─124b7573-113c-4cdf-9a5c-e0cc9c78dfaa
# ╠═b63a9eee-7d14-4b4a-92b5-6b960aad90c1
# ╟─6bc5781a-d025-46d8-bd70-19cb374d3c4e
# ╠═390c73cc-8253-4547-afb3-c225b58fb2c4
# ╟─07e1b431-84c6-4cfc-ad2b-edc54fdc59ca
# ╠═565eade7-574f-473a-a9a4-6aeca2a812e6
# ╠═da5d5ca7-5e6d-425c-823b-881c59ff8a3f
# ╟─776b150d-ddbe-494a-9e53-632a33bb502b
# ╠═7accb307-22f1-4035-951b-14bd04bcf5b2
# ╠═3226a396-db02-4f50-ac67-9f8c57e7566b
# ╟─397c56dc-79c0-4a18-8a04-fded2676c658
# ╠═ccce7595-bc08-4855-819b-9dd3e9e88198
# ╠═044e548f-f81e-4cf8-a872-8e5440e74ddb
# ╟─5ffcf3a8-e3d5-4093-8794-7eb133252d04
# ╠═d394450d-825d-4071-9102-239c1dd79b98
# ╠═a08b9637-f6d3-4671-9cfa-2a67b27f6585
# ╠═4b9aa903-625f-4d4e-89de-aa175b33cbf5
# ╠═c3f49d87-b42f-4e72-9f8f-02c8aa72af4f
# ╠═ca81d1ea-252a-4d43-aab8-4a08f95176dc
# ╠═ed6a9d9d-ce6a-4eab-955b-c6e676103603
# ╠═066d0e67-a4db-4acc-b1c8-9fcb6fd76e88
# ╠═4cadf68d-5cce-4331-8e96-8757a42667ef
# ╟─5c002c0a-10cf-4eb3-bac1-08ac25aab4cf
# ╟─41387ee8-07fa-4a54-bb1d-cbe78d56e08b
