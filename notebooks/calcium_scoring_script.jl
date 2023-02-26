### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 45990afe-ae3e-11ed-1234-05d8e35be9d0
# ╠═╡ show_logs = false
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ cbe9b7fe-09fd-408f-80af-6b6ca8b6e03e
begin
    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!,rmsd
end

# ╔═╡ d27e74de-2d7e-4411-830e-6b996faf318c
begin
	
	include(srcdir("helper_functions.jl")); 
	include(srcdir("masks.jl"));
end 

# ╔═╡ 9dafd569-ebc8-4272-b802-8aa82e8f35bf
begin	
	ENERGIES = [80,135]
	SIZES = ["Small", "Medium", "Large", "Small1","Medium1","Large1"]
	DENSITIES = ["Density1","Density2","Density3"]
end 

# ╔═╡ 9c3c2916-33cb-4959-a096-a5ccee87d7b9
begin
	dfs = []
	for _size in SIZES 
		for density in DENSITIES
			for energy in ENERGIES
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
	
				#Load Dicoms	
				root_path = string(datadir("dcms_measurement_new", _SIZE,string(density),string(energy)));

				# @info root_path
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
				
	 			vol_small_gt, vol_medium_gt, vol_large_gt = π * (1/2)^2 * 3, π * (3/2)^2 * 3, π * (5/2)^2 * 3 # mm^3
	 			vol2 = [vol_large_gt, vol_medium_gt, vol_small_gt] * 1e-3 
	 			vols2 = vcat(vol2, vol2, vol2) # cm^3
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
end

# ╔═╡ 6f9620f8-adf3-47ba-be03-1e730bf97a83
dfs

# ╔═╡ Cell order:
# ╠═45990afe-ae3e-11ed-1234-05d8e35be9d0
# ╠═d27e74de-2d7e-4411-830e-6b996faf318c
# ╠═cbe9b7fe-09fd-408f-80af-6b6ca8b6e03e
# ╠═9dafd569-ebc8-4272-b802-8aa82e8f35bf
# ╠═9c3c2916-33cb-4959-a096-a5ccee87d7b9
# ╠═6f9620f8-adf3-47ba-be03-1e730bf97a83
