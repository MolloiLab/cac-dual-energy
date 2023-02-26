### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 6c99eb9f-99ee-40d9-b43b-e220d0c991ab
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ efdd5526-1c01-438f-a4d2-47adb669289c
begin
    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!,rmsd
end

# ╔═╡ 6b14198c-41bd-4f8f-8e2c-f9d31082f17c
begin
	
	include(srcdir("helper_functions.jl")); 
	include(srcdir("masks.jl"));
end

# ╔═╡ e716418b-ddee-4bc2-987a-d4e4b314d93a
begin	
	ENERGIES = [80,135]
	SIZES = ["Small", "Medium", "Large", "Small1","Medium1","Large1"]
	DENSITIES = ["Density1","Density2","Density3"]
end

# ╔═╡ 40521a9a-6e20-43e9-abb0-f0cc174bc06b
#Load Calibration Parameters (Mass Decomp)
begin
	param_base_pth = datadir("calibration_params/")
	small_pth = string(param_base_pth,"Small.csv")
	med_pth = string(param_base_pth,"Medium.csv")
	large_pth = string(param_base_pth,"Large.csv")

	small_param = DataFrame(CSV.File(small_pth))
	med_param = DataFrame(CSV.File(med_pth))
	large_param = DataFrame(CSV.File(large_pth))
end;

# ╔═╡ 6afbc04e-e635-4599-a447-d947be4e0cb6
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 67d669ee-b09d-4276-a767-b8f6e90ccdb7
function predict_concentration(x, y, p)
	A = p[1] + (p[2] * x) + (p[3] * y) + (p[4] * x^2) + (p[5] * x * y) + (p[6] * y^2)
	B = 1 + (p[7] * x) + (p[8] * y)
	F = A / B
end

# ╔═╡ f9232844-0fdf-4f86-93d9-e57560962453
begin
	dfs = []
	for _size in SIZES 
		for density in DENSITIES
			@info _size, density

			#Load masks 
			if (_size == SIZES[1] || _size == SIZES[4])
				_SIZE = "small"
			elseif (_size == SIZES[2] || _size == SIZES[5])
				_SIZE = "medium"
			elseif (_size == SIZES[3] || _size == SIZES[6])
				_SIZE = "large"
			end
				
			root_new = string(datadir("julia_arrays",_SIZE),"/")
	
			mask_L_HD = Array(CSV.read(string(root_new, "mask_L_HD.csv"), 			DataFrame; header=false))
			mask_M_HD = Array(CSV.read(string(root_new, "mask_M_HD.csv"), 
			DataFrame; header=false))
			mask_S_HD = Array(CSV.read(string(root_new, "mask_S_HD.csv"), 			DataFrame; header=false))
			mask_L_MD = Array(CSV.read(string(root_new, "mask_L_MD.csv"), DataFrame; header=false))
			mask_M_MD = Array(CSV.read(string(root_new, "mask_M_MD.csv"), DataFrame; header=false))
			mask_S_MD = Array(CSV.read(string(root_new, "mask_S_MD.csv"), DataFrame; header=false))
			mask_L_LD = Array(CSV.read(string(root_new, "mask_L_LD.csv"), 			DataFrame; header=false))
			mask_M_LD = Array(CSV.read(string(root_new, "mask_M_LD.csv"), 			DataFrame; header=false))
			mask_S_LD = Array(CSV.read(string(root_new, "mask_S_LD.csv"), 			DataFrame; header=false))

			masks = mask_L_HD+mask_L_MD+mask_L_LD+mask_M_HD+mask_M_MD+mask_M_LD+mask_S_HD+mask_S_MD+mask_S_LD;

			##energy1
			pth = datadir("dcms_measurement_new/", _size, density, string(ENERGIES[1]))		
			dcm = dcmdir_parse(pth)
			dcm_array = load_dcm_array(dcm)
			
##Material Decomp
			dilate_mask_S_HD = dilate(dilate(mask_S_HD))
			dilate_mask_S_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_S_HD_3D[:, :, z] = dilate_mask_S_HD
			end
		
			dilate_mask_M_HD = dilate(dilate(mask_M_HD))
			dilate_mask_M_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_M_HD_3D[:, :, z] = dilate_mask_M_HD
			end
		
			dilate_mask_L_HD = dilate(dilate(mask_L_HD))
			dilate_mask_L_HD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_L_HD_3D[:, :, z] = dilate_mask_L_HD
			end
			
			dilate_mask_S_MD = dilate(dilate(mask_S_MD))
			dilate_mask_S_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_S_MD_3D[:, :, z] = dilate_mask_S_MD
			end
		
			dilate_mask_M_MD = dilate(dilate(mask_M_MD))
			dilate_mask_M_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_M_MD_3D[:, :, z] = dilate_mask_M_MD
			end
			
			dilate_mask_L_MD = dilate(dilate(mask_L_MD))
			dilate_mask_L_MD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_L_MD_3D[:, :, z] = dilate_mask_L_MD
			end
			
			dilate_mask_S_LD = dilate(dilate(mask_S_LD))
			dilate_mask_S_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_S_LD_3D[:, :, z] = dilate_mask_S_LD
			end
		
			dilate_mask_M_LD = dilate(dilate(mask_M_LD))
			dilate_mask_M_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_M_LD_3D[:, :, z] = dilate_mask_M_LD
			end
			
			dilate_mask_L_LD = dilate(dilate(mask_L_LD))
			dilate_mask_L_LD_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				dilate_mask_L_LD_3D[:, :, z] = dilate_mask_L_LD
			end

			pixel_size = DICOMUtils.get_pixel_size(dcm[1].meta)
			masks_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				masks_3D[:, :, z] = masks
			end
			means1 = [mean(dcm_array[dilate_mask_L_HD_3D]), mean(dcm_array[dilate_mask_M_HD_3D]), mean(dcm_array[dilate_mask_S_HD_3D]), mean(dcm_array[dilate_mask_L_MD_3D]), mean(dcm_array[dilate_mask_M_MD_3D]), mean(dcm_array[dilate_mask_S_MD_3D]), mean(dcm_array[dilate_mask_L_LD_3D]), mean(dcm_array[dilate_mask_M_LD_3D]), mean(dcm_array[dilate_mask_S_LD_3D])]
			
			##energy2
			pth2 = datadir("dcms_measurement_new/", _size, density, string(ENERGIES[2]))
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
			)

			push!(dfs, df_results)





			#Load Dicoms
			root_path = string(datadir("dcms_measurement_new", _SIZE,string(density),string(energy)));

			dcm = dcmdir_parse(root_path)
			dcm_array = load_dcm_array(dcm)
	
			dcm_path_list = dcm_list_builder(root_path);
# 				#@info dcm_path_list
			pth = dcm_path_list[1];
			scan = basename(pth);
			header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

		#Agatson Scoring
		
			#High Density
			arr = dcm_array[:, :, :];
			
			mask_L_HD_3D = Array{Bool}(undef, size(arr));
			for z in 1:size(arr, 3)
				##should I fix this >
				mask_L_HD_3D[:, :, z] = dilate(dilate(mask_L_HD))
			end;
			mean(arr[erode(erode(erode(erode(erode(mask_L_HD_3D)))))])
			##Dilated Mask
			dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D);
			
			pixel_size = DICOMUtils.get_pixel_size(header)				
			overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD);
			alg = Agatston()
			avg_mass_cals = 0.0007975563520531468
			agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, avg_mass_cals, alg);
			
			#Medium Density
			mask_L_MD_3D = Array{Bool}(undef, size(arr));
			for z in 1:size(arr, 3)
				mask_L_MD_3D[:, :, z] = mask_L_MD
			end;
			##Dilated Mask
			dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D);
			# overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask");
			overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD);
			agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, avg_mass_cals, alg);
						

			#Low Density
			mask_L_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_L_LD_3D[:, :, z] = mask_L_LD
			end
			##Dilated Mask
			dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D);
			overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD);
			agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, avg_mass_cals, alg)
		
		#Score Medium Inserts
			#High Density
			mask_M_HD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_HD_3D[:, :, z] = mask_M_HD
			end
			
			##Dilated mask
			dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D);
			overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD);
			agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, avg_mass_cals, alg)

			#Medium Density
			mask_M_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_MD_3D[:, :, z] = mask_M_MD
			end
			##Dilated mask
			dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D)
			overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD);
			agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, avg_mass_cals, alg)

			#Low Density

			mask_M_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_M_LD_3D[:, :, z] = mask_M_LD
			end

			dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D);
			overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD);
			agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, avg_mass_cals, alg)

		#Score Small Inserts
			#High Density
			mask_S_HD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_HD_3D[:, :, z] = mask_S_HD
			end
			
			##Dilated mask
			dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D);
			overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD);
			agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, avg_mass_cals, alg)

			#Medium Density
			mask_S_MD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_MD_3D[:, :, z] = mask_S_MD
			end
			dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D);
			overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD);
			agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, avg_mass_cals, alg)

			#Low Density
			mask_S_LD_3D = Array{Bool}(undef, size(arr))
			for z in 1:size(arr, 3)
				mask_S_LD_3D[:, :, z] = mask_S_LD
			end

			##Dilated mask
			dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D);
			overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD);
			agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, avg_mass_cals, alg)

		#Results
			
			calcium_densities = [733, 733, 733, 411, 411, 411, 151, 151, 151]; 
			gt_masses = calcium_densities .* vols2 .* 3
			inserts = ["Low Density", "Medium Density", "High Density"]

			##Agatston
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


			df = DataFrame(;
				scan=scan,
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
			
			
			push!(dfs,df)
		end
	end
end

# ╔═╡ Cell order:
# ╠═6c99eb9f-99ee-40d9-b43b-e220d0c991ab
# ╠═6b14198c-41bd-4f8f-8e2c-f9d31082f17c
# ╠═efdd5526-1c01-438f-a4d2-47adb669289c
# ╠═e716418b-ddee-4bc2-987a-d4e4b314d93a
# ╠═40521a9a-6e20-43e9-abb0-f0cc174bc06b
# ╠═6afbc04e-e635-4599-a447-d947be4e0cb6
# ╠═67d669ee-b09d-4276-a767-b8f6e90ccdb7
# ╠═f9232844-0fdf-4f86-93d9-e57560962453
