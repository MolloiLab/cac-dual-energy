### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 6e6f2a53-f883-4b1c-917f-8e50219e7611
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ 0b2b31a5-29df-4344-9ce8-a1313f9e01d0
begin
    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!,rmsd
end

# ╔═╡ 8fc48986-5362-4303-aedc-cce3cb46863b
include(srcdir("helper_functions.jl")); include(srcdir("masks.jl"));

# ╔═╡ 3acf4253-394a-4ff5-be3b-32a49e0f947e
TableOfContents()

# ╔═╡ e4a717b2-4b60-49df-b948-76ddec11f35f
begin
	sizes = ["Small", "Medium", "Large", "Small1", "Medium1", "Large1"]
	densities = ["Density1", "Density2", "Density3"]
	energies = [80, 135]
end

# ╔═╡ c53b95fa-407c-4ba4-b031-5d0543261987
calcium_densities = [
	[733, 733, 733, 411, 411, 411, 151, 151, 151],
	[669, 669, 669, 370, 370, 370, 90, 90, 90],
	[552, 552, 552, 222, 222, 222, 52, 52, 52],
	[797, 797, 797, 101, 101, 101, 37, 37, 37], 
	[403, 403, 403, 48, 48, 48, 32, 32, 32],
	[199, 199, 199, 41, 41, 41, 27, 27, 27]
]

# ╔═╡ 831886a2-f222-4cc3-a5c5-87c676fa4d01
md"""
# Load Calibration Parameters
"""

# ╔═╡ ee139e88-4619-4636-8db8-b70917733643
begin
	param_base_pth = datadir("calibration_params/")
	small_pth = string(param_base_pth,"Small.csv")
	med_pth = string(param_base_pth,"Medium.csv")
	large_pth = string(param_base_pth,"Large.csv")

	small_param = DataFrame(CSV.File(small_pth))
	med_param = DataFrame(CSV.File(med_pth))
	large_param = DataFrame(CSV.File(large_pth))
end;

# ╔═╡ 3649478a-b065-4df6-9210-8a5e9fbb5010
md"""
# Loop
"""

# ╔═╡ 9ecdd5e4-8e78-47af-a26e-2fc30a929ddd
begin
	dfs_m = []
	dfs_a = []
	for _size in sizes 
		for density in densities
			@info _size, density
			
			if (_size == sizes[1] || _size == sizes[4])
				_SIZE = "small"
			elseif (_size == sizes[2] || _size == sizes[5])
				_SIZE = "medium"
			elseif (_size == sizes[3] || _size == sizes[6])
				_SIZE = "large"
			end
				
			root_new = string(datadir("julia_arrays/", _SIZE),"/")
			mask_L_HD = Array(CSV.read(string(root_new, "mask_L_HD.csv"), DataFrame; header=false))
			mask_M_HD = Array(CSV.read(string(root_new, "mask_M_HD.csv"), DataFrame; header=false))
			mask_S_HD = Array(CSV.read(string(root_new, "mask_S_HD.csv"), DataFrame; header=false))
			mask_L_MD = Array(CSV.read(string(root_new, "mask_L_MD.csv"), DataFrame; header=false))
			mask_M_MD = Array(CSV.read(string(root_new, "mask_M_MD.csv"), DataFrame; header=false))
			mask_S_MD = Array(CSV.read(string(root_new, "mask_S_MD.csv"), DataFrame; header=false))
			mask_L_LD = Array(CSV.read(string(root_new, "mask_L_LD.csv"), DataFrame; header=false))
			mask_M_LD = Array(CSV.read(string(root_new, "mask_M_LD.csv"), DataFrame; header=false))
			mask_S_LD = Array(CSV.read(string(root_new, "mask_S_LD.csv"), DataFrame; header=false))
			masks = mask_L_HD+mask_L_MD+mask_L_LD+mask_M_HD+mask_M_MD+mask_M_LD+mask_S_HD+mask_S_MD+mask_S_LD;
			
			pth = datadir("dcms_measurement_new/", _size, density, string(energies[1]))

			
			dcm = dcmdir_parse(pth)
			dcm_array = load_dcm_array(dcm)
			
			#high density masks
			dilate_mask_S_HD = dilate_mask_small(mask_S_HD)
			dilate_mask_S_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_S_HD_3D[:, :, z] = dilate_mask_S_HD
			end
			
			dilate_mask_M_HD = dilate_mask_medium(mask_M_HD)
			dilate_mask_M_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_M_HD_3D[:, :, z] = dilate_mask_M_HD
			end
			
			dilate_mask_L_HD = dilate_mask_large(mask_L_HD)
			dilate_mask_L_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_L_HD_3D[:, :, z] = dilate_mask_L_HD
			end
			
			#medium density masks
			dilate_mask_S_MD = dilate_mask_small(mask_S_MD)
			dilate_mask_S_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_S_MD_3D[:, :, z] = dilate_mask_S_MD
			end
			
			dilate_mask_M_MD = dilate_mask_medium(mask_M_MD)
			dilate_mask_M_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_M_MD_3D[:, :, z] = dilate_mask_M_MD
			end
			
			dilate_mask_L_MD = dilate_mask_large(mask_L_MD)
			dilate_mask_L_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_L_MD_3D[:, :, z] = dilate_mask_L_MD
			end
			
			#low density masks
			dilate_mask_S_LD = dilate_mask_small(mask_S_LD)
			dilate_mask_S_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_S_LD_3D[:, :, z] = dilate_mask_S_LD
			end
			
			dilate_mask_M_LD = dilate_mask_medium(mask_M_LD)
			dilate_mask_M_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_M_LD_3D[:, :, z] = dilate_mask_M_LD
			end
			
			dilate_mask_L_LD = dilate_mask_large(mask_L_LD)
			dilate_mask_L_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_L_LD_3D[:, :, z] = dilate_mask_L_LD
			end

			## Low energy
			pixel_size = DICOMUtils.get_pixel_size(dcm[1].meta)
			masks_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				masks_3D[:, :, z] = masks
			end
			means1 = [mean(dcm_array[dilate_mask_L_HD_3D]), mean(dcm_array[dilate_mask_M_HD_3D]), mean(dcm_array[dilate_mask_S_HD_3D]), mean(dcm_array[dilate_mask_L_MD_3D]), mean(dcm_array[dilate_mask_M_MD_3D]), mean(dcm_array[dilate_mask_S_MD_3D]), mean(dcm_array[dilate_mask_L_LD_3D]), mean(dcm_array[dilate_mask_M_LD_3D]), mean(dcm_array[dilate_mask_S_LD_3D])]

			## High energy
			pth2 = datadir("dcms_measurement_new/", _size, density, string(energies[2]))
			dcm2 = dcmdir_parse(pth2)
			dcm_array2 = load_dcm_array(dcm2)

			means2 = [mean(dcm_array2[dilate_mask_L_HD_3D]), mean(dcm_array2[dilate_mask_M_HD_3D]), mean(dcm_array2[dilate_mask_S_HD_3D]), mean(dcm_array2[dilate_mask_L_MD_3D]), mean(dcm_array2[dilate_mask_M_MD_3D]), mean(dcm_array2[dilate_mask_S_MD_3D]), mean(dcm_array2[dilate_mask_L_LD_3D]), mean(dcm_array2[dilate_mask_M_LD_3D]), mean(dcm_array2[dilate_mask_S_LD_3D])]

			## Calculate Predicted Densities
			calculated_intensities = hcat(means1, means2)
			predicted_densities = zeros(9)
			for i in 1:9
				predicted_densities[i] = predict_concentration(means1[i], means2[i], Array(small_param))
			end

			## Choose Calcium Density
			if (density == "Density1") && (_size == "Large" || _size == "Medium" || _size == "Small")
				calcium_density = calcium_densities[1]
			elseif (density == "Density2") && (_size == "Large" || _size == "Medium" || _size == "Small")
				calcium_density = calcium_densities[2]
			elseif (density == "Density3") && (_size == "Large" || _size == "Medium" || _size == "Small")
				calcium_density = calcium_densities[3]
			elseif (density == "Density1") && (_size == "Large1" || _size == "Medium1" || _size == "Small1")
				calcium_density = calcium_densities[4]
			elseif (density == "Density2") && (_size == "Large1" || _size == "Medium1" || _size == "Small1")
				calcium_density = calcium_densities[5]
			elseif (density == "Density3") && (_size == "Large1" || _size == "Medium1" || _size == "Small1")
				calcium_density = calcium_densities[6]
			end

			## Calculate Mass
			voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3] * 1e-3 # cm^3
		
			vol_small, vol_medium, vol_large = count(dilate_mask_S_HD_3D) * voxel_size, count(dilate_mask_M_HD_3D) * voxel_size, count(dilate_mask_L_HD_3D) * voxel_size #cm^3
		
			vol_slice1 = [vol_large, vol_medium, vol_small]
			vols = vcat(vol_slice1, vol_slice1, vol_slice1)# cm^3
			predicted_masses = predicted_densities .* vols

			
			vol_small_gt, vol_medium_gt, vol_large_gt = π * (1/2)^2 * 3, π * (3/2)^2 * 3, π * (5/2)^2 * 3 # mm^3
		
			vol2 = [vol_large_gt, vol_medium_gt, vol_small_gt] * 1e-3 
			vols2 = vcat(vol2, vol2, vol2) # cm^3
			gt_masses = calcium_density .* vols2 .* 3

			df = DataFrame(
				phantom_size = _size,
				density = density,
				insert_sizes = ["Large", "Medium", "Small"],
				ground_truth_mass_hd = gt_masses[1:3],
				predicted_mass_hd = predicted_masses[1:3],
				ground_truth_mass_md = gt_masses[4:6],
				predicted_mass_md = predicted_masses[4:6],
				ground_truth_mass_ld = gt_masses[7:9],
				predicted_mass_ld = predicted_masses[7:9],
			)

			push!(dfs_m, df)

			##Agatson 
			header = dcm[1].meta
			alg = Agatston()
			pixel_size = DICOMUtils.get_pixel_size(header)	
			avg_mass_cals = 0.0007975563520531468

			#high density

			overlayed_mask_s_hd = create_mask(dcm_array, dilate_mask_S_HD_3D);
			agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, avg_mass_cals, alg)

			overlayed_mask_m_hd = create_mask(dcm_array, dilate_mask_M_HD_3D);
			agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, avg_mass_cals, alg)

			overlayed_mask_l_hd = create_mask(dcm_array, dilate_mask_L_HD_3D);			
			agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, avg_mass_cals, alg);

			#medium density

			overlayed_mask_s_md = create_mask(dcm_array, dilate_mask_S_MD_3D);
			agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, avg_mass_cals, alg)
			
			overlayed_mask_m_md = create_mask(dcm_array, dilate_mask_M_MD_3D);
			agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, avg_mass_cals, alg)
			
			overlayed_mask_l_md = create_mask(dcm_array, dilate_mask_L_MD_3D);
			agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, avg_mass_cals, alg);

			#low density

			overlayed_mask_s_ld = create_mask(dcm_array, dilate_mask_S_LD_3D);
			agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, avg_mass_cals, alg)
		
			overlayed_mask_m_ld = create_mask(dcm_array, dilate_mask_M_LD_3D);
			agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, avg_mass_cals, alg)
			
			overlayed_mask_l_ld = create_mask(dcm_array, dilate_mask_L_LD_3D);
			agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, avg_mass_cals, alg)

			#Results

			calcium_densities1 = [733, 733, 733, 411, 411, 411, 151, 151, 151];
				
			vol_small_gt, vol_medium_gt, vol_large_gt = π * (1/2)^2 * 3, π * (3/2)^2 * 3, π * (5/2)^2 * 3 # mm^3
			vol2 = [vol_large_gt, vol_medium_gt, vol_small_gt] * 1e-3 
			vols2 = vcat(vol2, vol2, vol2) # cm^3
			gt_masses = calcium_densities1 .* vols2 .* 3
			inserts = ["Low Density", "Medium Density", "High Density"]

			##calculate agatston
			calculated_agat_large = [agat_l_ld, agat_l_md, agat_l_hd]
			calculated_agat_medium = [agat_m_ld, agat_m_md, agat_m_hd]
			calculated_agat_small = [agat_s_ld, agat_s_md, agat_s_hd]

			##Mass
			volume_gt = [7.065, 63.585, 176.625]
			calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]
			calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]
			predicted_mass_hd = [mass_l_hd, mass_m_hd, mass_s_hd]
			predicted_mass_md = [mass_l_md, mass_m_hd, mass_s_hd]
			predicted_mass_ld = [mass_l_ld, mass_m_ld, mass_s_ld]
			calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]

			scan = energies[1];
			df = DataFrame(;
				scan = scan,
				inserts=inserts,
				calculated_agat_large=calculated_agat_large,
				calculated_agat_medium=calculated_agat_medium,
				calculated_agat_small=calculated_agat_small,
				#ground_truth_mass_large=ground_truth_mass_large,
				ground_truth_mass_hd = gt_masses[1:3],
				predicted_mass_hd = predicted_mass_hd,
				#ground_truth_mass_medium=ground_truth_mass_medium,
				ground_truth_mass_md = gt_masses[4:6],
				predicted_mass_md=predicted_mass_md,
				#ground_truth_mass_small=ground_truth_mass_small,
				ground_truth_mass_ld = gt_masses[7:9],
				predicted_mass_ld = predicted_mass_ld,
				avg_mass_cals = avg_mass_cals,
			)
			
			
			push!(dfs_a, df)
		end
	end
end

# ╔═╡ abd691ff-e2e9-45b6-b249-681ee1425818
dfs_m

# ╔═╡ c31a910e-2f2e-4aab-8a18-1a39308bcd96
dfs_a

# ╔═╡ 821865e5-0456-4a70-8339-ba39974a0b88
begin
    new_df = vcat(dfs_m[1:length(dfs_m)]...)
    output_path = string(datadir("results"),"/","masses_matdecomp.csv")
	
	save(output_path,new_df)

	new_df = vcat(dfs_a[1:length(dfs_a)]...)
    output_path = string(datadir("results"),"/","masses_agat.csv")
	save(output_path,new_df)
end

# ╔═╡ d6021c28-9e79-422e-9be5-e133127b4017
output_path

# ╔═╡ Cell order:
# ╠═6e6f2a53-f883-4b1c-917f-8e50219e7611
# ╠═0b2b31a5-29df-4344-9ce8-a1313f9e01d0
# ╠═3acf4253-394a-4ff5-be3b-32a49e0f947e
# ╠═8fc48986-5362-4303-aedc-cce3cb46863b
# ╠═e4a717b2-4b60-49df-b948-76ddec11f35f
# ╠═c53b95fa-407c-4ba4-b031-5d0543261987
# ╟─831886a2-f222-4cc3-a5c5-87c676fa4d01
# ╠═ee139e88-4619-4636-8db8-b70917733643
# ╟─3649478a-b065-4df6-9210-8a5e9fbb5010
# ╠═9ecdd5e4-8e78-47af-a26e-2fc30a929ddd
# ╠═abd691ff-e2e9-45b6-b249-681ee1425818
# ╠═c31a910e-2f2e-4aab-8a18-1a39308bcd96
# ╠═d6021c28-9e79-422e-9be5-e133127b4017
# ╠═821865e5-0456-4a70-8339-ba39974a0b88
