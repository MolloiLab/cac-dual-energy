### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 6c99eb9f-99ee-40d9-b43b-e220d0c991ab
# ╠═╡ show_logs = false
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ efdd5526-1c01-438f-a4d2-47adb669289c
begin
    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!, rmsd
end

# ╔═╡ 6b14198c-41bd-4f8f-8e2c-f9d31082f17c
include(srcdir("helper_functions.jl")); include(srcdir("masks.jl"));

# ╔═╡ e716418b-ddee-4bc2-987a-d4e4b314d93a
begin	
	ENERGIES = [80, 135]
	SIZES = ["Small", "Medium", "Large", "Small1","Medium1","Large1"]
	DENSITIES = ["Density1","Density2","Density3"]
end

# ╔═╡ 40521a9a-6e20-43e9-abb0-f0cc174bc06b
begin
	small_param = CSV.read(datadir("calibration_params", "Small.csv"), DataFrame)
	med_param = CSV.read(datadir("calibration_params", "Medium.csv"), DataFrame)
	large_param = CSV.read(datadir("calibration_params", "Large.csv"), DataFrame)
end;

# ╔═╡ 57f0f760-2925-4aaf-a1c4-7ca882cb06d1
ca_dens = [
	[733, 733, 733, 411, 411, 411, 151, 151, 151],
	[669, 669, 669, 370, 370, 370, 90, 90, 90],
	[552, 552, 552, 222, 222, 222, 52, 52, 52],
	[797, 797, 797, 101, 101, 101, 37, 37, 37], 
	[403, 403, 403, 48, 48, 48, 32, 32, 32],
	[199, 199, 199, 41, 41, 41, 27, 27, 27]
]

# ╔═╡ f9232844-0fdf-4f86-93d9-e57560962453
begin
	dfs_m = []
	dfs_a = []
	dfs_v = []
	for _size in SIZES 
		for density in DENSITIES
			@info _size, density

			#---------------- Reusable Pieces ----------------#
			#Load masks 
			if (_size == SIZES[1] || _size == SIZES[4])
				_SIZE = "small"
				_SIZE2 = "Small"
			elseif (_size == SIZES[2] || _size == SIZES[5])
				_SIZE = "medium"
				_SIZE2 = "Medium"
			elseif (_size == SIZES[3] || _size == SIZES[6])
				_SIZE = "large"
				_SIZE2 = "Large"
			end
				
			root_new = datadir("julia_arrays", _SIZE)
	
			mask_L_HD = Array(CSV.read(joinpath(root_new, "mask_L_HD.csv"), 			DataFrame; header=false))
			mask_M_HD = Array(CSV.read(joinpath(root_new, "mask_M_HD.csv"), 
			DataFrame; header=false))
			mask_S_HD = Array(CSV.read(joinpath(root_new, "mask_S_HD.csv"), 			DataFrame; header=false))
			mask_L_MD = Array(CSV.read(joinpath(root_new, "mask_L_MD.csv"), DataFrame; header=false))
			mask_M_MD = Array(CSV.read(joinpath(root_new, "mask_M_MD.csv"), DataFrame; header=false))
			mask_S_MD = Array(CSV.read(joinpath(root_new, "mask_S_MD.csv"), DataFrame; header=false))
			mask_L_LD = Array(CSV.read(joinpath(root_new, "mask_L_LD.csv"), 			DataFrame; header=false))
			mask_M_LD = Array(CSV.read(joinpath(root_new, "mask_M_LD.csv"), 			DataFrame; header=false))
			mask_S_LD = Array(CSV.read(joinpath(root_new, "mask_S_LD.csv"), 			DataFrame; header=false))

			masks = mask_L_HD+mask_L_MD+mask_L_LD+mask_M_HD+mask_M_MD+mask_M_LD+mask_S_HD+mask_S_MD+mask_S_LD;

			## energy1
			pth = datadir("dcms_measurement_combined", _size, density, string(ENERGIES[1]))		
			dcm = dcmdir_parse(pth)
			dcm_array = load_dcm_array(dcm)
			
			
			dilate_mask_S_HD = dilate(dilate(mask_S_HD))
			dilate_mask_S_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_S_HD_3D[:, :, z] = dilate_mask_S_HD
			end
		
			dilate_mask_M_HD = dilate(dilate(mask_M_HD))
			dilate_mask_M_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_M_HD_3D[:, :, z] = dilate_mask_M_HD
			end
		
			dilate_mask_L_HD = dilate(dilate(mask_L_HD))
			dilate_mask_L_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_L_HD_3D[:, :, z] = dilate_mask_L_HD
			end
			
			dilate_mask_S_MD = dilate(dilate(mask_S_MD))
			dilate_mask_S_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_S_MD_3D[:, :, z] = dilate_mask_S_MD
			end
		
			dilate_mask_M_MD = dilate(dilate(mask_M_MD))
			dilate_mask_M_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_M_MD_3D[:, :, z] = dilate_mask_M_MD
			end
			
			dilate_mask_L_MD = dilate(dilate(mask_L_MD))
			dilate_mask_L_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_L_MD_3D[:, :, z] = dilate_mask_L_MD
			end
			
			dilate_mask_S_LD = dilate(dilate(mask_S_LD))
			dilate_mask_S_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_S_LD_3D[:, :, z] = dilate_mask_S_LD
			end
		
			dilate_mask_M_LD = dilate(dilate(mask_M_LD))
			dilate_mask_M_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_M_LD_3D[:, :, z] = dilate_mask_M_LD
			end
			
			dilate_mask_L_LD = dilate(dilate(mask_L_LD))
			dilate_mask_L_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				dilate_mask_L_LD_3D[:, :, z] = dilate_mask_L_LD
			end

			pixel_size = DICOMUtils.get_pixel_size(dcm[1].meta)
			masks_3D = Array{Bool}(undef, size(dcm_array))
			for z in axes(dcm_array, 3)
				masks_3D[:, :, z] = masks
			end

			## Background
			local center_insert
			if _SIZE == "small"
				center_insert = [175, 320]
				# center_insert = [320, 175]
			elseif _SIZE == "medium"
				center_insert = [230, 370]
				# center_insert = [370, 230]
			elseif _SIZE == "large"
				center_insert = [285, 420]
				# center_insert = [420, 285]
			end
			background_mask = zeros(size(dcm_array)...)
			background_mask[
				(center_insert[1]-5):(center_insert[1]+5),
				(center_insert[2]-5):(center_insert[2]+5),
				2,
			] .= 1

			dilated_mask_L_bkg = dilate_mask_large_bkg(Bool.(background_mask))
			ring_mask_L_bkg = ring_mask_large(dilated_mask_L_bkg)

			dilated_mask_M_bkg = dilate_mask_medium_bkg(Bool.(background_mask))
			ring_mask_M_bkg = ring_mask_medium(dilated_mask_M_bkg)

			dilated_mask_S_bkg = dilate_mask_small_bkg(Bool.(background_mask))
			ring_mask_S_bkg = ring_mask_small(dilated_mask_S_bkg)

			#---------------- Material Decomposition ----------------#
			means1 = [
				mean(dcm_array[dilate_mask_L_HD_3D]),
				mean(dcm_array[dilate_mask_M_HD_3D]), mean(dcm_array[dilate_mask_S_HD_3D]), mean(dcm_array[dilate_mask_L_MD_3D]), mean(dcm_array[dilate_mask_M_MD_3D]), mean(dcm_array[dilate_mask_S_MD_3D]), mean(dcm_array[dilate_mask_L_LD_3D]), mean(dcm_array[dilate_mask_M_LD_3D]), mean(dcm_array[dilate_mask_S_LD_3D])
			]
			means1_bkg = [
				mean(dcm_array[dilated_mask_L_bkg]),
				mean(dcm_array[dilated_mask_M_bkg]),
				mean(dcm_array[dilated_mask_S_bkg])
			]
			
			## energy2
			pth2 = datadir("dcms_measurement_combined", _size, density, string(ENERGIES[2]))
			dcm2 = dcmdir_parse(pth2)
			dcm_array2 = load_dcm_array(dcm2)
			
			means2 = [
				mean(dcm_array2[dilate_mask_L_HD_3D]), mean(dcm_array2[dilate_mask_M_HD_3D]), mean(dcm_array2[dilate_mask_S_HD_3D]), mean(dcm_array2[dilate_mask_L_MD_3D]), mean(dcm_array2[dilate_mask_M_MD_3D]), mean(dcm_array2[dilate_mask_S_MD_3D]), mean(dcm_array2[dilate_mask_L_LD_3D]), mean(dcm_array2[dilate_mask_M_LD_3D]), mean(dcm_array2[dilate_mask_S_LD_3D])
			]
			means2_bkg = [
				mean(dcm_array2[dilated_mask_L_bkg]),
				mean(dcm_array2[dilated_mask_M_bkg]),
				mean(dcm_array2[dilated_mask_S_bkg])
			]

			## Calculate Predicted Densities
			calculated_intensities = hcat(means1, means2)
			predicted_densities = zeros(9)
			for i in 1:9
				predicted_densities[i] = predict_concentration(means1[i], means2[i], Array(small_param))
			end

			calculated_intensities_bkg = hcat(means1_bkg, means2_bkg)
			predicted_densities_bkg = zeros(3)
			for i in 1:3
				predicted_densities_bkg[i] = predict_concentration(means1_bkg[i], means2_bkg[i], Array(small_param))
			end

			## Choose Calcium Density
			if (density == "Density1") && (_size == "Large" || _size == "Medium" || _size == "Small")
				calcium_density = ca_dens[1]
			elseif (density == "Density2") && (_size == "Large" || _size == "Medium" || _size == "Small")
				calcium_density = ca_dens[2]
			elseif (density == "Density3") && (_size == "Large" || _size == "Medium" || _size == "Small")
				calcium_density = ca_dens[3]
			elseif (density == "Density1") && (_size == "Large1" || _size == "Medium1" || _size == "Small1")
				calcium_density = ca_dens[4]
			elseif (density == "Density2") && (_size == "Large1" || _size == "Medium1" || _size == "Small1")
				calcium_density = ca_dens[5]
			elseif (density == "Density3") && (_size == "Large1" || _size == "Medium1" || _size == "Small1")
				calcium_density = ca_dens[6]
			end

			## Calculate Mass
			voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3] * 1e-3 # cm^3
		
			vol_small, vol_medium, vol_large = count(dilate_mask_S_HD_3D) * voxel_size, count(dilate_mask_M_HD_3D) * voxel_size, count(dilate_mask_L_HD_3D) * voxel_size #cm^3
		
			vol_slice1 = [vol_large, vol_medium, vol_small]
			vols = vcat(vol_slice1, vol_slice1, vol_slice1)# cm^3
			predicted_masses = predicted_densities .* vols

			vol_bkg = count(dilated_mask_L_bkg) * voxel_size
			predicted_masses_bkg = predicted_densities_bkg * vol_bkg

			
			vol_small_gt, vol_medium_gt, vol_large_gt = π * (1/2)^2 * 3, π * (3/2)^2 * 3, π * (5/2)^2 * 3 # mm^3
		
			vol2 = [vol_large_gt, vol_medium_gt, vol_small_gt] * 1e-3 
			vols2 = vcat(vol2, vol2, vol2) # cm^3
			gt_masses = calcium_density .* vols2 .* 3

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				insert_sizes = ["Large", "Medium", "Small"],
				ground_truth_mass_hd = gt_masses[1:3],
				predicted_mass_hd = predicted_masses[1:3],
				ground_truth_mass_md = gt_masses[4:6],
				predicted_mass_md = predicted_masses[4:6],
				ground_truth_mass_ld = gt_masses[7:9],
				predicted_mass_ld = predicted_masses[7:9],
				bkg_mass = predicted_masses_bkg
			)

			push!(dfs_m, df_results)			

			#---------------- Calibration ----------------#
			local centers
			if _size == "Small" || _size == "Small1"
				centers = [187, 318]
			elseif _size == "Medium" || _size == "Medium1"
				centers = [238, 365]
			elseif _size == "Large" || _size == "Large1"
				centers = [285, 415]
			end
			calibration_rod = zeros(25, 25, size(dcm_array, 3))
			for z in axes(dcm_array, 3)
				rows, cols, depth = size(dcm_array)
				half_row, half_col = centers[1], centers[2]
				offset = 12
				row_range = half_row-offset:half_row+offset
				col_range = half_col-offset:half_col+offset	
				calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
			
			hu_calcium = mean(calibration_rod[:, :, 5])
			ρ_calcium = 0.200

			# erode_mask_cal = erode(erode(mask_L_HD))
			# erode_mask_cal_3D = Array{Bool}(undef, size(dcm_array))
			# for z in axes(dcm_array, 3)
			# 	erode_mask_cal_3D[:, :, z] = erode_mask_cal
			# end
			
			# hu_calcium = mean(dcm_array[erode_mask_cal_3D])
			# local ρ_calcium
			# if _size == "Small"
			# 	ρ_calcium = ca_dens[1][1] * 1e-3
			# elseif _size == "Medium"
			# 	ρ_calcium = ca_dens[2][1] * 1e-3
			# elseif _size == "Large"
			# 	ρ_calcium = ca_dens[3][1] * 1e-3
			# elseif _size == "Small1"
			# 	ρ_calcium = ca_dens[4][1] * 1e-3
			# elseif _size == "Medium1"
			# 	ρ_calcium = ca_dens[5][1] * 1e-3
			# elseif _size == "Large1"
			# 	ρ_calcium = ca_dens[6][1] * 1e-3
			# end

			# hu_calcium = 300
			# ρ_calcium = 0.200
			
			#---------------- Agatston ----------------#
			mass_calibration = ρ_calcium / hu_calcium
			@info mass_calibration
			header = dcm[1].meta
			pixel_size = DICOMUtils.get_pixel_size(header)
			alg = Agatston()
			kV=80

			overlayed_bkg_mask_L = create_mask(dcm_array, dilated_mask_L_bkg)
			overlayed_bkg_mask_M = create_mask(dcm_array, dilated_mask_M_bkg)
			overlayed_bkg_mask_S = create_mask(dcm_array, dilated_mask_S_bkg)

			agat_bkg, mass_bkg_large = score(
				overlayed_bkg_mask_L,
				pixel_size,
				mass_calibration,
				alg;
				kV=kV
			)
			agat_bkg, mass_bkg_medium = score(
				overlayed_bkg_mask_M,
				pixel_size,
				mass_calibration,
				alg;
				kV=kV
			)
			agat_bkg, mass_bkg_small = score(
				overlayed_bkg_mask_S,
				pixel_size,
				mass_calibration,
				alg;
				kV=kV
			)

			mass_bkg = [mass_bkg_large, mass_bkg_medium, mass_bkg_small]
		
			## High Density
			overlayed_mask_l_hd = create_mask(dcm_array, dilate_mask_L_HD_3D);
			agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, mass_calibration, alg);
			
			# Medium Density
			overlayed_mask_l_md = create_mask(dcm_array, dilate_mask_L_MD_3D);
			agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, mass_calibration, alg);
						

			# Low Density
			overlayed_mask_l_ld = create_mask(dcm_array, dilate_mask_L_LD_3D);
			agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, mass_calibration, alg)
		
			# Score Medium Inserts
			# High Density
			overlayed_mask_m_hd = create_mask(dcm_array, dilate_mask_M_HD_3D);
			agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, mass_calibration, alg)

			# Medium Density
			overlayed_mask_m_md = create_mask(dcm_array, dilate_mask_M_MD_3D);
			agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, mass_calibration, alg)

			#Low Density
			overlayed_mask_m_ld = create_mask(dcm_array, dilate_mask_M_LD_3D);
			agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, mass_calibration, alg)

			# Score Small Inserts
			# High Density
			overlayed_mask_s_hd = create_mask(dcm_array, dilate_mask_S_HD_3D);
			agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, mass_calibration, alg)

			# Medium Density
			overlayed_mask_s_md = create_mask(dcm_array, dilate_mask_S_MD_3D);
			agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, mass_calibration, alg)

			#Low Density
			overlayed_mask_s_ld = create_mask(dcm_array, dilate_mask_S_LD_3D);
			agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, mass_calibration, alg)

			# Results
			
			calcium_densities = [733, 733, 733, 411, 411, 411, 151, 151, 151]; 
			gt_masses = calcium_densities .* vols2 .* 3

			## Agatston
			calculated_agat_large = [agat_l_ld, agat_l_md, agat_l_hd]
			calculated_agat_medium = [agat_m_ld, agat_m_md, agat_m_hd]
			calculated_agat_small = [agat_s_ld, agat_s_md, agat_s_hd]

			## Mass
			predicted_mass_hd = [mass_l_hd, mass_m_hd, mass_s_hd]
			predicted_mass_md = [mass_l_md, mass_m_hd, mass_s_hd]
			predicted_mass_ld = [mass_l_ld, mass_m_ld, mass_s_ld]


			df = DataFrame(;
				phantom_size = _size,
				density = density,
				insert_sizes = ["Large", "Medium", "Small"],
				calculated_agat_large=calculated_agat_large,
				calculated_agat_medium=calculated_agat_medium,
				calculated_agat_small=calculated_agat_small,
				ground_truth_mass_hd = gt_masses[1:3],
				predicted_mass_hd = predicted_mass_hd,
				ground_truth_mass_md = gt_masses[4:6],
				predicted_mass_md=predicted_mass_md,
				ground_truth_mass_ld = gt_masses[7:9],
				predicted_mass_ld = predicted_mass_ld,
				bkg_mass = mass_bkg
			)
			
			push!(dfs_a, df)

			#---------------- Volume Fraction ----------------#

			# Background
			hu_heart_tissue_large_bkg = mean(dcm_array[ring_mask_L_bkg])
			mass_large_bkg = score(dcm_array[dilated_mask_L_bkg], hu_calcium, hu_heart_tissue_large_bkg, voxel_size, ρ_calcium, VolumeFraction())

			hu_heart_tissue_medium_bkg = mean(dcm_array[ring_mask_M_bkg])
			mass_medium_bkg = score(dcm_array[dilated_mask_M_bkg], hu_calcium, hu_heart_tissue_medium_bkg, voxel_size, ρ_calcium, VolumeFraction())

			hu_heart_tissue_small_bkg = mean(dcm_array[ring_mask_S_bkg])
			mass_small_bkg = score(dcm_array[dilated_mask_S_bkg], hu_calcium, hu_heart_tissue_small_bkg, voxel_size, ρ_calcium, VolumeFraction())

			mass_bkg = [mass_large_bkg, mass_medium_bkg, mass_small_bkg]

			# Score Large Inserts
			## High Density
			hu_heart_tissue_large_hd = mean(dcm_array[dilated_mask_L_bkg])
			mass_large_hd = score(dcm_array[dilated_mask_L_bkg], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

			## Medium Density
			hu_heart_tissue_large_md = mean(dcm_array[ring_mask_L_bkg])
			mass_large_md = score(dcm_array[dilated_mask_L_bkg], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

			## Low Density
			hu_heart_tissue_large_ld = mean(dcm_array[ring_mask_L_bkg])
			mass_large_ld = score(dcm_array[dilated_mask_L_bkg], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

			# Score Medium Inserts
			## High Density
			hu_heart_tissue_medium_hd = mean(dcm_array[ring_mask_M_bkg])
			mass_medium_hd = score(dcm_array[dilated_mask_M_bkg], hu_calcium, hu_heart_tissue_medium_hd, voxel_size, ρ_calcium, VolumeFraction())

			## Medium Density
			hu_heart_tissue_medium_md = mean(dcm_array[ring_mask_M_bkg])
			mass_medium_md = score(dcm_array[dilated_mask_M_bkg], hu_calcium, hu_heart_tissue_medium_md, voxel_size, ρ_calcium, VolumeFraction())

			## Low Density
			hu_heart_tissue_medium_ld = mean(dcm_array[ring_mask_M_bkg])
			mass_medium_ld = score(dcm_array[dilated_mask_M_bkg], hu_calcium, hu_heart_tissue_medium_ld, voxel_size, ρ_calcium, VolumeFraction())

			# Score Small Inserts
			## High Density
			hu_heart_tissue_small_hd = mean(dcm_array[ring_mask_S_bkg])
			mass_small_hd = score(dcm_array[dilated_mask_S_bkg], hu_calcium, hu_heart_tissue_large_hd, voxel_size, ρ_calcium, VolumeFraction())

			## Medium Density
			hu_heart_tissue_small_md = mean(dcm_array[ring_mask_S_bkg])
			mass_small_md = score(dcm_array[dilated_mask_S_bkg], hu_calcium, hu_heart_tissue_large_md, voxel_size, ρ_calcium, VolumeFraction())

			## Low Density
			hu_heart_tissue_small_ld = mean(dcm_array[ring_mask_S_bkg])
			mass_small_ld = score(dcm_array[dilated_mask_S_bkg], hu_calcium, hu_heart_tissue_large_ld, voxel_size, ρ_calcium, VolumeFraction())

			# Results
			predicted_mass_hd = [mass_large_hd, mass_medium_hd, mass_small_hd]
			predicted_mass_md = [mass_large_md, mass_medium_hd, mass_small_hd]
			predicted_mass_ld = [mass_large_ld, mass_medium_ld, mass_small_ld]

			df = DataFrame(;
				phantom_size = _size,
				density = density,
				insert_sizes = ["Large", "Medium", "Small"],
				ground_truth_mass_hd = gt_masses[1:3],
				predicted_mass_hd = predicted_mass_hd,
				ground_truth_mass_md = gt_masses[4:6],
				predicted_mass_md=predicted_mass_md,
				ground_truth_mass_ld = gt_masses[7:9],
				predicted_mass_ld = predicted_mass_ld,
				bkg_mass = mass_bkg
			)
			push!(dfs_v, df)
		end
	end
end

# ╔═╡ c56131ff-ef0e-413b-886a-ef2be0e2489f
begin
    new_df_m = vcat(dfs_m[1:length(dfs_m)]...)
    CSV.write(datadir("results", "material_decomposition.csv"), new_df_m)

	new_df_a = vcat(dfs_a[1:length(dfs_a)]...)
    CSV.write(datadir("results", "agatston.csv"), new_df_a)

	new_df_v = vcat(dfs_v[1:length(dfs_v)]...)
    CSV.write(datadir("results", "volume_fraction.csv"), new_df_v)
end

# ╔═╡ 0e9fabb1-05c3-4fcb-9f35-304c8782a7b6
new_df_m

# ╔═╡ 38ca9627-9ac7-4fe3-a017-592db6145281
new_df_a

# ╔═╡ 17e2adf7-df30-4fd7-a46b-07dd8a0b562b
new_df_v

# ╔═╡ Cell order:
# ╠═6c99eb9f-99ee-40d9-b43b-e220d0c991ab
# ╠═6b14198c-41bd-4f8f-8e2c-f9d31082f17c
# ╠═efdd5526-1c01-438f-a4d2-47adb669289c
# ╠═e716418b-ddee-4bc2-987a-d4e4b314d93a
# ╠═40521a9a-6e20-43e9-abb0-f0cc174bc06b
# ╠═57f0f760-2925-4aaf-a1c4-7ca882cb06d1
# ╠═f9232844-0fdf-4f86-93d9-e57560962453
# ╠═c56131ff-ef0e-413b-886a-ef2be0e2489f
# ╠═0e9fabb1-05c3-4fcb-9f35-304c8782a7b6
# ╠═38ca9627-9ac7-4fe3-a017-592db6145281
# ╠═17e2adf7-df30-4fd7-a46b-07dd8a0b562b
