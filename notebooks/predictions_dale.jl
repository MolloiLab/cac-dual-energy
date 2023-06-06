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

# ╔═╡ 29f9cf29-9437-42a6-8dab-89105273c187
include(srcdir("masks.jl")); include(srcdir("dicom_utils.jl"));

# ╔═╡ 1b568480-5467-4ffd-9099-de81066e407e
TableOfContents()

# ╔═╡ c02cc808-ac3b-479a-b5a1-9abb36b93a03
md"""
# Calibration
"""

# ╔═╡ 668d1999-ce4e-4510-9dd4-84a26ab26dcf
begin
	densities_cal = [10, 16, 25, 45, 50, 75, 100, 200, 300, 400, 500, 600, 800]
	sizes_cal = [30]
	energies = [80, 135]
end

# ╔═╡ ee95b544-2e02-4de3-bd1c-019d8fc3cdd6
md"""
## Low Energy
"""

# ╔═╡ 73ab1438-d7b2-483a-a69b-2cb2d942a2fa
begin
	means_80 = Dict(:density => densities_cal, :means => zeros(length(densities_cal)))
	for (i, density) in enumerate(densities_cal)
		for _size in sizes_cal
			path = joinpath(datadir("data_new","dcms", "cal"), string(density), string(_size), string(energies[1]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
		
			center_insert1, center_insert2 = 187, 318
			
			calibration_rod = zeros(25, 25, size(dcm_array, 3))
			
			for z in axes(dcm_array, 3)
				rows, cols, depth = size(dcm_array)
				half_row, half_col = center_insert1, center_insert2
				offset = 12
				row_range = half_row-offset:half_row+offset
				col_range = half_col-offset:half_col+offset	
				calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
		
			means_80[:means][i] = mean(calibration_rod)
		end
	end
end

# ╔═╡ 37181932-1398-4c45-9ff0-fb8d488ab562
md"""
## High Energy
"""

# ╔═╡ 4a8cfe9e-6b9a-4645-9f6d-13afa87d8213
begin
	means_135 = Dict(:density => densities_cal, :means => zeros(length(densities_cal)))
	for (i, density) in enumerate(densities_cal)
		for _size in sizes_cal
			path = joinpath(datadir("data_new","dcms", "cal"), string(density), string(_size), string(energies[2]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
		
			center_insert1, center_insert2 = 187, 318
			
			calibration_rod = zeros(25, 25, size(dcm_array, 3))
			
			for z in axes(dcm_array, 3)
				rows, cols, depth = size(dcm_array)
				half_row, half_col = center_insert1, center_insert2
				offset = 12
				row_range = half_row-offset:half_row+offset
				col_range = half_col-offset:half_col+offset	
				calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
		
			means_135[:means][i] = mean(calibration_rod)
		end
	end
end

# ╔═╡ 6760a907-7817-4692-abd3-2c1b66dc3a50
md"""
## Fit Parameters
"""

# ╔═╡ 727de9ea-d78c-4c5f-9223-351a2fac45e3
calculated_intensities = hcat(means_80[:means], means_135[:means]) # low energy, high energy

# ╔═╡ 38910d55-20d0-445e-98ea-c2f36aaa5255
ps = fit_calibration(calculated_intensities, densities_cal)

# ╔═╡ 905ccae5-0ee5-4242-b58a-25ef4bcbb8b9
md"""
## Check Results
"""

# ╔═╡ 2402b78f-580b-47eb-962f-707c90d6e74c
begin
	predicted_densities = []
	
	for i in 1:length(densities_cal)
		append!(
			predicted_densities, 
			score(calculated_intensities[i, 1], calculated_intensities[i, 2], ps, MaterialDecomposition()
			)
		)
	end
end

# ╔═╡ c0095d39-14f2-464d-8f55-ea1739a7b313
df = DataFrame(
	densities = densities_cal,
	predicted_densities = predicted_densities,
	mean_intensities_low = means_80[:means],
	mean_intensities_high = means_135[:means],
)

# ╔═╡ cb0fc64e-8e36-4250-9ea6-427b5dccf27d
md"""
# Validation
"""

# ╔═╡ af661a05-e9ec-4b0b-8d7e-0ce39af016ec
begin
	densities_val = [
		"15_18_22"
		"26_29_36"
		"52_59_73"
		"110_210_310"
		"410_610_780"
	] 

	sizes_val = ["small", "medium", "large"]
end;

# ╔═╡ 4c2bf021-6a1d-4c8c-a3f4-7aee6c39c361
md"""
## Load masks
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
	dfs = []
	for density in densities_val
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
			path_80 = datadir("data_new","dcms", "val", density, _size, string(energies[1]))
			dcm_80 = dcmdir_parse(path_80)
			dcm_array_80 = load_dcm_array(dcm_80)
			pixel_size = get_pixel_size(dcm_80[1].meta)
			
			means_80 = [
				mean(dcm_array_80[dilated_mask_L_HD_3D]), mean(dcm_array_80[dilated_mask_L_MD_3D]), mean(dcm_array_80[dilated_mask_L_LD_3D]),
				mean(dcm_array_80[dilated_mask_M_HD_3D]), mean(dcm_array_80[dilated_mask_M_MD_3D]), mean(dcm_array_80[dilated_mask_M_LD_3D]),
				mean(dcm_array_80[dilated_mask_S_HD_3D]), mean(dcm_array_80[dilated_mask_S_MD_3D]), mean(dcm_array_80[dilated_mask_S_LD_3D])
			]

			# High Energy
			path_135 = datadir("data_new","dcms", "val", density, _size, string(energies[2]))
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
				predicted_densities[i] = score(means_80[i], means_80[i], ps, MaterialDecomposition())
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
			push!(dfs, df_results)
		end
	end
end

# ╔═╡ 5e6aa434-51c9-460c-9675-ac33d061379d
dfs

# ╔═╡ Cell order:
# ╠═e5248e4b-76bc-4a45-b38d-2c7f55f8d8f7
# ╠═67c7c5e9-e181-4560-b77a-50e9d7977348
# ╠═155ddbb8-f3a2-4db9-b256-7fd771785d66
# ╠═90830c78-6f24-41a3-a2eb-1fce3756bcaf
# ╠═68ec4904-d908-4c1e-98e1-00ec7ba06f59
# ╠═29f9cf29-9437-42a6-8dab-89105273c187
# ╠═1b568480-5467-4ffd-9099-de81066e407e
# ╟─c02cc808-ac3b-479a-b5a1-9abb36b93a03
# ╠═668d1999-ce4e-4510-9dd4-84a26ab26dcf
# ╟─ee95b544-2e02-4de3-bd1c-019d8fc3cdd6
# ╠═73ab1438-d7b2-483a-a69b-2cb2d942a2fa
# ╟─37181932-1398-4c45-9ff0-fb8d488ab562
# ╠═4a8cfe9e-6b9a-4645-9f6d-13afa87d8213
# ╟─6760a907-7817-4692-abd3-2c1b66dc3a50
# ╠═727de9ea-d78c-4c5f-9223-351a2fac45e3
# ╠═38910d55-20d0-445e-98ea-c2f36aaa5255
# ╟─905ccae5-0ee5-4242-b58a-25ef4bcbb8b9
# ╠═2402b78f-580b-47eb-962f-707c90d6e74c
# ╠═c0095d39-14f2-464d-8f55-ea1739a7b313
# ╟─cb0fc64e-8e36-4250-9ea6-427b5dccf27d
# ╠═af661a05-e9ec-4b0b-8d7e-0ce39af016ec
# ╟─4c2bf021-6a1d-4c8c-a3f4-7aee6c39c361
# ╠═0029b07b-6c4f-49f3-841b-3e75e2f2a34f
# ╠═42fd7e43-96fe-4c85-9807-9239cc66f10b
# ╠═5e6aa434-51c9-460c-9675-ac33d061379d
