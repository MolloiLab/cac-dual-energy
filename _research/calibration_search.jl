### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 47c381ac-37d0-4566-8602-66d707ec8a91
using DrWatson

# ╔═╡ 4fa4a2f9-b837-42bd-80e0-ded9b38c5623
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ b1eb1bd7-9e57-4228-b5e2-d45d2858bf56
using PlutoUI, CairoMakie, Statistics, CSV, DataFrames, DICOM, CSVFiles

# ╔═╡ e06ebcc3-5cd7-49f5-80ea-dc202b39e16f
using StatsBase: quantile!, rmsd

# ╔═╡ 33b118ab-c37a-4875-bc83-2336d4061b7c
using CalciumScoring

# ╔═╡ c952efc2-42a2-4d6e-9632-0584a91c9fd6
using GLM, MLJBase

# ╔═╡ eb662ecd-1b18-4f3e-8690-e62f10090148
include(srcdir("helper_functions.jl"));

# ╔═╡ 4a64afb2-b424-4519-a4b5-cc0ea4f9932a
include(srcdir("plot_utils.jl"));

# ╔═╡ 3994b2c3-0a09-45e2-bc81-51b8a48049ed
using OrderedCollections

# ╔═╡ 5f9cca1a-0478-409e-ba14-5af1a32483fc
include(srcdir("masks.jl")); include(srcdir("dicom_utils.jl"));

# ╔═╡ 9e42e006-e13d-477f-a6fa-8fc93f181160
TableOfContents()

# ╔═╡ 2e71d6be-5dc0-4a76-8d15-6e425e35277a
md"""
# Low Density
"""

# ╔═╡ f0c9c98c-0dc5-45cd-9ff4-9ad7bda22049
begin
	sizes_cal = [30]
	energies = [80, 120, 135]
end;

# ╔═╡ a565d7c9-7f6a-4598-a986-cde2e5afa983
md"""
## Volume Fraction & Agatston Calibration
"""

# ╔═╡ 0a2622c8-0094-4143-a9e8-512d171e0577
path = joinpath(datadir("dcms", "cal", "100", "30"), string(energies[2]))

# ╔═╡ bd5558ab-3632-46b4-8059-ff490667348e
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

# ╔═╡ 86f21a96-d2c7-4726-a1e0-3823d7056ed2
md"""
## Calibration
"""

# ╔═╡ 54b47e68-ef05-493d-b29e-e7f72b7c87bb
densities_cal_all = [3, 6, 9, 10, 16, 25, 31, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800] # calcium densities

# ╔═╡ faf5c214-5a5e-420b-93ec-e084cff6966e
begin
	_a1 = 9
	_density_ranges1 = [
		range(i; length = _a1, step = 1) for i in 1:(25 - (_a1 - 1))
	]
	_density_ranges2 = [
		range(i; length = _a1, step = 2) for i in 1:(25 - 2 * (_a1 - 1))
	]
	_density_ranges3 = [
		range(i; length = _a1, step = 3) for i in 1:(25 - 3 * (_a1 - 1))
	]

	_a1 = 10
	_density_ranges4 = [
		range(i; length = _a1, step = 1) for i in 1:(25 - (_a1 - 1))
	]
	_density_ranges5 = [
		range(i; length = _a1, step = 2) for i in 1:(25 - 2 * (_a1 - 1))
	]
	_density_ranges6 = [
		range(i; length = _a1, step = 3) for i in 1:(25 - 3 * (_a1 - 1))
	]

	_a1 = 12
	_density_ranges7 = [
		range(i; length = _a1, step = 1) for i in 1:(25 - (_a1 - 1))
	]
	_density_ranges8 = [
		range(i; length = _a1, step = 2) for i in 1:(25 - 2 * (_a1 - 1))
	]
	_density_ranges9 = [
		range(i; length = _a1, step = 3) for i in 1:(25 - 3 * (_a1 - 1))
	]

	_a1 = 14
	_density_ranges10 = [
		range(i; length = _a1, step = 1) for i in 1:(25 - (_a1 - 1))
	]
	_density_ranges11 = [
		range(i; length = _a1, step = 2) for i in 1:(25 - 2 * (_a1 - 1))
	]
	_density_ranges12 = [
		range(i; length = _a1, step = 3) for i in 1:(25 - 3 * (_a1 - 1))
	]
	
	density_ranges = vcat(
		_density_ranges1, 
		_density_ranges2, 
		_density_ranges3,
		_density_ranges4, 
		_density_ranges5, 
		_density_ranges6,
		_density_ranges7, 
		_density_ranges8, 
		_density_ranges9,
		_density_ranges10, 
		_density_ranges11, 
		_density_ranges12,
	)
end

# ╔═╡ ba9ac2da-c13d-43d2-b38d-6c23819d03f8
begin
	calculated_intensities_low_tuple = []
	ps_low_tuple = []
	for d_range in density_ranges
		densities_cal_low = densities_cal_all[d_range]
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
		calculated_intensities_low = hcat(means_80_low[:means], means_135_low[:means]) # low energy, high energy
		push!(calculated_intensities_low_tuple, calculated_intensities_low)
		
		ps_low = fit_calibration(calculated_intensities_low, densities_cal_low)
		push!(ps_low_tuple, ps_low)
	end
end

# ╔═╡ 054e438f-b47f-46a2-9105-f0068b699e8e
md"""
### Check Results
"""

# ╔═╡ afae797e-deb4-4411-9aea-115df49326b1
md"""
## Validation
"""

# ╔═╡ 003049a8-a27f-4d2e-b079-c3641eb02354
densities_val_all = [
	"15_18_22"
	"26_29_36"
	"52_59_73"
	"110_210_310"
	"410_610_780"
]

# ╔═╡ d7a5af04-87f8-4399-9595-3fc1b1baf014
begin
	densities_val_low = densities_val_all[1:3]
	sizes_val = ["small", "medium", "large"]
end;

# ╔═╡ c58a3410-9341-425d-a698-2adffb91b8bb
energies_val = [80, 120, 135]

# ╔═╡ c9ff2870-7c17-4fa5-a0ab-ea4301abcc41
md"""
### Load masks
"""

# ╔═╡ b18a49a4-eee0-4556-ac5c-594f3ef39d8e
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

# ╔═╡ c212cfb4-ccb2-4926-9b04-c2f8992ca184
md"""
### Predict
"""

# ╔═╡ d761d60c-af91-4265-ba4d-ce432b9f8947
begin
	dfs_low_md = []
	dfs_low_vf = []
	dfs_low_a = []
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
			
			for r in eachindex(density_ranges)
				calculated_intensities = hcat(means_80, means_135)
				predicted_densities = zeros(size(means_135))
				for i in eachindex(predicted_densities)
					predicted_densities[i] = score(means_80[i], means_80[i], ps_low_tuple[r], MaterialDecomposition())
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
	
				df_results = DataFrame(
					phantom_size = _size,
					density = density,
					insert_densities = [:low_density, :medium_density, :high_density],
					calibration_points = [densities_cal_all[density_ranges[r]], densities_cal_all[density_ranges[r]], densities_cal_all[density_ranges[r]]],
					gt_mass_large_inserts = gt_mass_large_inserts,
					predicted_mass_large_inserts = predicted_mass_large_inserts,
					gt_mass_medium_inserts = gt_mass_medium_inserts,
					predicted_mass_medium_inserts = predicted_mass_medium_inserts,
					gt_mass_small_inserts = gt_mass_small_inserts,
					predicted_mass_small_inserts = predicted_mass_small_inserts,
				)
				push!(dfs_low_md, df_results)
			end

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
			)
			push!(dfs_low_vf, df_results)

			#------- Agatston -------#
			
			alg = Agatston()
			mass_cal_factor = ρ_calcium_100 / hu_calcium_100
			kV = 120
				
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
			)
			push!(dfs_low_a, df_results)
			
		end
	end
end

# ╔═╡ 55e8841e-d73c-4cb9-b498-79b887033cdc
md"""
# High Density
"""

# ╔═╡ c766a18c-109d-479a-8e75-e55ea23bdf4b
md"""
## Calibration
"""

# ╔═╡ f2aaa90a-803a-4979-aa55-5366779b84d1
begin
	calculated_intensities_high_tuple = []
	ps_high_tuple = []
	for d_range in density_ranges
		densities_cal_high = densities_cal_all[d_range]
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
		push!(calculated_intensities_high_tuple, calculated_intensities_high)
			
		ps_high = fit_calibration(calculated_intensities_high, densities_cal_high)
		push!(ps_high_tuple, ps_high)
	end
end

# ╔═╡ ea132b81-c2ba-45ae-8722-5ebb2fed33b1
md"""
### Check Results
"""

# ╔═╡ 4f180915-a6a1-4f92-a00a-6e8c58433431
md"""
## Validation
"""

# ╔═╡ b818c0a3-fc94-4fed-a8c5-ddf4448c03bf
densities_val_high = densities_val_all[length(densities_val_low):end];

# ╔═╡ 7b51d8aa-7479-4839-9a10-83521f71bc3d
md"""
### Predict
"""

# ╔═╡ 93b5b7a7-1712-42d7-a381-4efc0aed495d
begin
	dfs_high_md = []
	dfs_high_vf = []
	dfs_high_a = []
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
			
			for r in eachindex(density_ranges)
				calculated_intensities = hcat(means_80, means_135)
				predicted_densities = zeros(size(means_135))
				for i in eachindex(predicted_densities)
					predicted_densities[i] = score(means_80[i], means_80[i], ps_high_tuple[r], MaterialDecomposition())
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
	
				df_results = DataFrame(
					phantom_size = _size,
					density = density,
					insert_densities = [:low_density, :medium_density, :high_density],
					calibration_points = [densities_cal_all[density_ranges[r]], densities_cal_all[density_ranges[r]], densities_cal_all[density_ranges[r]]],
					gt_mass_large_inserts = gt_mass_large_inserts,
					predicted_mass_large_inserts = predicted_mass_large_inserts,
					gt_mass_medium_inserts = gt_mass_medium_inserts,
					predicted_mass_medium_inserts = predicted_mass_medium_inserts,
					gt_mass_small_inserts = gt_mass_small_inserts,
					predicted_mass_small_inserts = predicted_mass_small_inserts,
				)
				push!(dfs_high_md, df_results)
			end

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
			)
			push!(dfs_high_vf, df_results)

			#------- Agatston -------#
			
			alg = Agatston()
			mass_cal_factor = ρ_calcium_100 / hu_calcium_100
			kV = 120
				
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
			)
			push!(dfs_high_a, df_results)
			
		end
	end
end

# ╔═╡ fecb08f9-46f6-4a96-abc2-ba48b7d31aae
md"""
# Analyze
"""

# ╔═╡ 7b7d3805-4ed2-4a5e-8fd2-9ac7d5eb71ee
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

# ╔═╡ 1d1475c5-6d2d-4ea6-bb79-be9fbab57ced
md"""
## Low Density
"""

# ╔═╡ 889317b0-afaa-485f-8f3e-3096c2f8a118
md"""
### Accuracy
"""

# ╔═╡ 6c30e840-0ccc-4822-a3a5-4c52b226bbe3
new_df_low_md = vcat(dfs_low_md[1:length(dfs_low_md)]...);

# ╔═╡ e42e347c-7734-4b45-835a-74c9535a4912
new_df_low_vf = vcat(dfs_low_vf[1:length(dfs_low_vf)]...);

# ╔═╡ 18ae2522-e304-47c4-8dc1-3281ce93d647
new_df_low_a = vcat(dfs_low_a[1:length(dfs_low_a)]...);

# ╔═╡ 4080cdbc-3ded-41eb-980f-e1f1a19947e3
co_1_low_md, r_squared_1_low_md, rms_values_1_low_md, pred_1_low_md = calculate_coefficients(new_df_low_md);

# ╔═╡ 5ca129f6-dcbd-4fdf-bf3e-5e053bd287b2
co_1_low_vf, r_squared_1_low_vf, rms_values_1_low_vf, pred_1_low_vf = calculate_coefficients(new_df_low_vf);

# ╔═╡ 922a6d18-0576-4381-85bf-e90e5291fd7b
co_1_low_a, r_squared_1_low_a, rms_values_1_low_a, pred_1_low_a = calculate_coefficients(new_df_low_a);

# ╔═╡ a119fa83-553e-444f-ab22-56279fba7fb4
function calculate_coefficients_calibration(df;
    label1=:gt_mass_large_inserts,
    label2=:gt_mass_medium_inserts,
    label3=:gt_mass_small_inserts,
    label4=:predicted_mass_large_inserts,
    label5=:predicted_mass_medium_inserts,
    label6=:predicted_mass_small_inserts
)
    gt_array = vec(hcat(df[!, label1], df[!, label2], df[!, label3]))
    calc_array = vec(hcat(df[!, label4], df[!, label5], df[!, label6]))
    data = DataFrame(X=gt_array, Y=calc_array)
    model = lm(@formula(Y ~ X), data)
    rms_values = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model))
    ]

    return rms_values
end

# ╔═╡ 9915bd09-7b40-443d-8cb7-01752d221507
begin
	cal_low_md_dict = OrderedDict()
	for i in range(start = 1; stop = size(new_df_low_md, 1), step = 3)
		j = i:i+2
		key = new_df_low_md[i, :calibration_points]
		if !haskey(cal_low_md_dict, :calibration_points)
			cal_low_md_dict[key] = OrderedDict{Symbol, Float64}()
		end
		rms_values = calculate_coefficients_calibration(new_df_low_md[j, :])

		cal_low_md_dict[key][:rmse] = rms_values[1]
		cal_low_md_dict[key][:rmsd] = rms_values[2]
	end
end

# ╔═╡ 0ad3ff78-db47-4905-b1b1-1f65e23e2aa0
@bind lim_ld PlutoUI.Slider(1.5:5, show_value = true, default = 2)

# ╔═╡ 29e7e252-4619-48b8-bfe6-9d2911e7a090
let
	f = Figure()
	labels = [string(key) for (key, value) in cal_low_md_dict]
	ax = CairoMakie.Axis(
		f[1, 1],
		title = "RMSE vs Calibration Points (Low Density)",
		# xticks = (1:length(labels), labels),
		# xticks = 1:length(labels),
		xticklabelrotation = π/2
	)

	i = 0
	for (key, value) in cal_low_md_dict
		i += 1
		barplot!(i, Float32.(value[:rmse]), label = labels[i])
	end

	ylims!(ax; low = 0, high = lim_ld)

	f
end

# ╔═╡ e9d081d6-b724-467a-a880-840cbbbc4dd3
labels = [string(key) for (key, value) in cal_low_md_dict]

# ╔═╡ 8616d09d-db8e-4f9e-a961-f3671fb30908
@bind lim_hd PlutoUI.Slider(1:5:100, show_value = true, default = 5)

# ╔═╡ 428a9351-3510-41a9-a320-bb61bf047f24
@bind r_low PlutoUI.Slider(1:length(density_ranges), show_value = true, default = 6)

# ╔═╡ 248416a4-a38d-4b05-9716-25b3b738b32d
begin
	predicted_densities_low = []
	
	for i in 1:length(densities_cal_all[density_ranges[r_low]])
		append!(
			predicted_densities_low, 
			score(calculated_intensities_low_tuple[r_low][i, 1], calculated_intensities_low_tuple[r_low][i, 2], ps_low_tuple[r_low], MaterialDecomposition()
			)
		)
	end
end

# ╔═╡ 0f109831-ff67-4f37-9d14-1507633a6d2d
df_low = DataFrame(
	densities = densities_cal_all[density_ranges[r_low]],
	predicted_densities = predicted_densities_low,
)

# ╔═╡ cbf9d4db-6fa5-428f-9e59-dc7c6ef2b9ca
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

# ╔═╡ 7eb7cedf-32ae-43a2-8f1b-1487aa4fdf55
with_theme(accuracy_low, medphys_theme)

# ╔═╡ 345553f6-007f-43b5-ab7c-6e1e2752234a
md"""
## High Density
"""

# ╔═╡ 8e4a9c1b-1cdb-4909-b855-02f36e0fea7f
md"""
### Accuracy
"""

# ╔═╡ a3783f08-f870-42a0-9e75-bd6567c91f96
new_df_high_md = vcat(dfs_high_md[1:length(dfs_high_md)]...);

# ╔═╡ c72067e2-fb3d-460d-b208-92ec8e6d76f6
new_df_high_vf = vcat(dfs_high_vf[1:length(dfs_high_vf)]...);

# ╔═╡ 5eca5057-5dc1-4e10-8a37-4dd21c4f147f
new_df_high_a = vcat(dfs_high_a[1:length(dfs_high_a)]...);

# ╔═╡ 89265fb5-367f-4278-9675-035bad46e97d
co_1_high_md, r_squared_1_high_md, rms_values_1_high_md, pred_1_high_md = calculate_coefficients(new_df_high_md);

# ╔═╡ 29a9f8c1-8b08-4f5e-8410-6c5311180b4f
co_1_high_vf, r_squared_1_high_vf, rms_values_1_high_vf, pred_1_high_vf = calculate_coefficients(new_df_high_vf);

# ╔═╡ 8253ed0a-a613-4875-a676-9349d7dad5be
co_1_high_a, r_squared_1_high_a, rms_values_1_high_a, pred_1_high_a = calculate_coefficients(new_df_high_a);

# ╔═╡ f50b9ccf-e31f-4b25-a4df-020f47ee38ad
begin
	cal_high_md_dict = OrderedDict()
	for i in range(start = 1; stop = size(new_df_high_md, 1), step = 3)
		j = i:i+2
		key = new_df_high_md[i, :calibration_points]
		if !haskey(cal_high_md_dict, :calibration_points)
			cal_high_md_dict[key] = OrderedDict{Symbol, Float64}()
		end
		rms_values = calculate_coefficients_calibration(new_df_high_md[j, :])

		cal_high_md_dict[key][:rmse] = rms_values[1]
		cal_high_md_dict[key][:rmsd] = rms_values[2]
	end
end

# ╔═╡ 7be78672-d78a-4566-8635-1cacbf77828d
let
	f = Figure()
	labels = [string(key) for (key, value) in cal_high_md_dict]
	ax = CairoMakie.Axis(
		f[1, 1],
		title = "RMSE vs Calibration Points (High Density)",
		# xticks = (1:length(labels), labels),
		xticklabelrotation = π/2
	)

	i = 0
	for (key, value) in cal_high_md_dict
		i += 1
		barplot!(i, Float32.(value[:rmse]), label = labels[i])
	end

	ylims!(ax; low = 0, high = lim_hd)

	f
end

# ╔═╡ d923ce17-dbbd-4505-9861-d562bc8d4ee1
cal_high_md_dict

# ╔═╡ ffa20df9-b847-4c3d-8228-33d9aea05ddf
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

# ╔═╡ 0813fc04-868b-4982-88e0-d90f201cf3c3
@bind r_high PlutoUI.Slider(1:length(density_ranges), show_value = true, default = 6)

# ╔═╡ d342309b-e7d3-4f2d-b90c-c88d4347a04f
begin
	predicted_densities_high = []
	
	for i in 1:length(densities_cal_all[density_ranges[r_high]])
		append!(
			predicted_densities_high, 
			score(calculated_intensities_high_tuple[r_high][i, 1], calculated_intensities_high_tuple[r_high][i, 2], ps_high_tuple[r_high], MaterialDecomposition()
			)
		)
	end
end

# ╔═╡ 0c6e4e21-42c7-4490-afff-db8ffce5217e
df_high = DataFrame(
	densities = densities_cal_all[density_ranges[r_high]],
	predicted_densities = predicted_densities_high,
)

# ╔═╡ aa6c6c88-f645-44a5-aa4e-4e6e3494be8d
with_theme(accuracy_high, medphys_theme)

# ╔═╡ Cell order:
# ╠═47c381ac-37d0-4566-8602-66d707ec8a91
# ╠═4fa4a2f9-b837-42bd-80e0-ded9b38c5623
# ╠═b1eb1bd7-9e57-4228-b5e2-d45d2858bf56
# ╠═3994b2c3-0a09-45e2-bc81-51b8a48049ed
# ╠═e06ebcc3-5cd7-49f5-80ea-dc202b39e16f
# ╠═33b118ab-c37a-4875-bc83-2336d4061b7c
# ╠═5f9cca1a-0478-409e-ba14-5af1a32483fc
# ╠═9e42e006-e13d-477f-a6fa-8fc93f181160
# ╟─2e71d6be-5dc0-4a76-8d15-6e425e35277a
# ╠═f0c9c98c-0dc5-45cd-9ff4-9ad7bda22049
# ╠═a565d7c9-7f6a-4598-a986-cde2e5afa983
# ╠═0a2622c8-0094-4143-a9e8-512d171e0577
# ╠═bd5558ab-3632-46b4-8059-ff490667348e
# ╠═86f21a96-d2c7-4726-a1e0-3823d7056ed2
# ╠═54b47e68-ef05-493d-b29e-e7f72b7c87bb
# ╠═faf5c214-5a5e-420b-93ec-e084cff6966e
# ╠═ba9ac2da-c13d-43d2-b38d-6c23819d03f8
# ╠═054e438f-b47f-46a2-9105-f0068b699e8e
# ╠═248416a4-a38d-4b05-9716-25b3b738b32d
# ╠═0f109831-ff67-4f37-9d14-1507633a6d2d
# ╠═afae797e-deb4-4411-9aea-115df49326b1
# ╠═003049a8-a27f-4d2e-b079-c3641eb02354
# ╠═d7a5af04-87f8-4399-9595-3fc1b1baf014
# ╠═c58a3410-9341-425d-a698-2adffb91b8bb
# ╠═c9ff2870-7c17-4fa5-a0ab-ea4301abcc41
# ╠═b18a49a4-eee0-4556-ac5c-594f3ef39d8e
# ╠═c212cfb4-ccb2-4926-9b04-c2f8992ca184
# ╠═d761d60c-af91-4265-ba4d-ce432b9f8947
# ╠═55e8841e-d73c-4cb9-b498-79b887033cdc
# ╠═c766a18c-109d-479a-8e75-e55ea23bdf4b
# ╠═f2aaa90a-803a-4979-aa55-5366779b84d1
# ╠═ea132b81-c2ba-45ae-8722-5ebb2fed33b1
# ╠═d342309b-e7d3-4f2d-b90c-c88d4347a04f
# ╠═0c6e4e21-42c7-4490-afff-db8ffce5217e
# ╠═4f180915-a6a1-4f92-a00a-6e8c58433431
# ╠═b818c0a3-fc94-4fed-a8c5-ddf4448c03bf
# ╠═7b51d8aa-7479-4839-9a10-83521f71bc3d
# ╠═93b5b7a7-1712-42d7-a381-4efc0aed495d
# ╠═fecb08f9-46f6-4a96-abc2-ba48b7d31aae
# ╠═c952efc2-42a2-4d6e-9632-0584a91c9fd6
# ╠═eb662ecd-1b18-4f3e-8690-e62f10090148
# ╠═4a64afb2-b424-4519-a4b5-cc0ea4f9932a
# ╠═7b7d3805-4ed2-4a5e-8fd2-9ac7d5eb71ee
# ╠═1d1475c5-6d2d-4ea6-bb79-be9fbab57ced
# ╠═889317b0-afaa-485f-8f3e-3096c2f8a118
# ╠═6c30e840-0ccc-4822-a3a5-4c52b226bbe3
# ╠═e42e347c-7734-4b45-835a-74c9535a4912
# ╠═18ae2522-e304-47c4-8dc1-3281ce93d647
# ╠═4080cdbc-3ded-41eb-980f-e1f1a19947e3
# ╠═5ca129f6-dcbd-4fdf-bf3e-5e053bd287b2
# ╠═922a6d18-0576-4381-85bf-e90e5291fd7b
# ╠═a119fa83-553e-444f-ab22-56279fba7fb4
# ╠═9915bd09-7b40-443d-8cb7-01752d221507
# ╠═0ad3ff78-db47-4905-b1b1-1f65e23e2aa0
# ╠═29e7e252-4619-48b8-bfe6-9d2911e7a090
# ╠═e9d081d6-b724-467a-a880-840cbbbc4dd3
# ╠═8616d09d-db8e-4f9e-a961-f3671fb30908
# ╠═7be78672-d78a-4566-8635-1cacbf77828d
# ╠═428a9351-3510-41a9-a320-bb61bf047f24
# ╠═cbf9d4db-6fa5-428f-9e59-dc7c6ef2b9ca
# ╠═7eb7cedf-32ae-43a2-8f1b-1487aa4fdf55
# ╠═345553f6-007f-43b5-ab7c-6e1e2752234a
# ╠═8e4a9c1b-1cdb-4909-b855-02f36e0fea7f
# ╠═a3783f08-f870-42a0-9e75-bd6567c91f96
# ╠═c72067e2-fb3d-460d-b208-92ec8e6d76f6
# ╠═5eca5057-5dc1-4e10-8a37-4dd21c4f147f
# ╠═89265fb5-367f-4278-9675-035bad46e97d
# ╠═29a9f8c1-8b08-4f5e-8410-6c5311180b4f
# ╠═8253ed0a-a613-4875-a676-9349d7dad5be
# ╠═f50b9ccf-e31f-4b25-a4df-020f47ee38ad
# ╠═d923ce17-dbbd-4505-9861-d562bc8d4ee1
# ╠═ffa20df9-b847-4c3d-8228-33d9aea05ddf
# ╠═0813fc04-868b-4982-88e0-d90f201cf3c3
# ╠═aa6c6c88-f645-44a5-aa4e-4e6e3494be8d
