### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 05f80d21-203e-4934-8e4d-dfc843919918
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ a8646d6e-5c7b-11ed-2556-73f2b13f709d
# ╠═╡ show_logs = false
begin
    using PlutoUI, CairoMakie, Statistics, ImageMorphology, CSV, DataFrames, DICOM, DICOMUtils
	using StatsBase: quantile!
end

# ╔═╡ 87d68ed1-6d7b-4552-8680-4a8901964109
TableOfContents()

# ╔═╡ 4676198a-2063-4d98-9b11-af5acfa32245
md"""
# Helper Functions
"""

# ╔═╡ b32917fc-df2c-42eb-8e7a-cc2c564f2e22
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ 0b83ca47-1734-4172-b3c8-0290e07bbece
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=25, show_value=true)
end;

# ╔═╡ 860c6c90-b07f-472b-b364-9a28bc84434c
function overlay_mask_plot(array, mask, var, title::AbstractString)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    indices_lbl = findall(x -> x == zs[var], label_array[:, 3])

    fig = Figure()
    ax = Makie.Axis(fig[1, 1])
    ax.title = title
    heatmap!((array[:, :, zs[var]]); colormap=:grays)
    scatter!(
        label_array[:, 1][indices_lbl],
        label_array[:, 2][indices_lbl];
        markersize=1,
        color=:red,
    )
    return fig
end;

# ╔═╡ 502abd7e-c123-463b-8581-9b886c87fa29
function predict_concentration(x, y, p)
	A = p[1] + (p[2] * x) + (p[3] * y) + (p[4] * x^2) + (p[5] * x * y) + (p[6] * y^2)
	B = 1 + (p[7] * x) + (p[8] * y)
	F = A / B
end

# ╔═╡ e0b52502-2075-47c0-a465-af93db8110de
sizes_folders = ["Small", "Medium", "Large", "Small1", "Medium1", "Large1"]

# ╔═╡ 2750bb98-c026-440c-9d89-ebf49dac1676
sizes = sizes_folders

# ╔═╡ ba3b914d-758c-430d-8845-34f7b6c57f50
densities = ["Density1", "Density2", "Density3"]

#densities=["Density1"]

# ╔═╡ fad04012-8c80-45b5-9b44-7f682f4a28ff
energies = [80, 135]

# ╔═╡ 97d59370-e42c-4b5d-9551-49281db1cc84
calcium_densities = [
	[733, 733, 733, 411, 411, 411, 151, 151, 151],
	[669, 669, 669, 370, 370, 370, 90, 90, 90],
	[552, 552, 552, 222, 222, 222, 52, 52, 52],
	[797, 797, 797, 101, 101, 101, 37, 37, 37], 
	[403, 403, 403, 48, 48, 48, 32, 32, 32],
	[199, 199, 199, 41, 41, 41, 27, 27, 27]
]

# ╔═╡ c9f111b6-1e66-459f-a1ff-66b44f3ae785
md"""
# Load Calibration Parameters
"""

# ╔═╡ 565b608a-f6a8-4579-9efd-37522dc9cde5
begin
	param_base_pth = datadir("calibration_params/")
	small_pth = string(param_base_pth,"Small.csv")
	med_pth = string(param_base_pth,"Medium.csv")
	large_pth = string(param_base_pth,"Large.csv")

	small_param = DataFrame(CSV.File(small_pth))
	med_param = DataFrame(CSV.File(med_pth))
	large_param = DataFrame(CSV.File(large_pth))
end;

# ╔═╡ e093bf43-6e8a-4f89-8117-a578712ad822
md"""
# Loop
"""

# ╔═╡ 24d92c10-9463-40a2-913c-94ed3728a169
begin
	dfs = []
	for _size in sizes_folders 
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

			## Low Density
			pixel_size = DICOMUtils.get_pixel_size(dcm[1].meta)
			masks_3D = Array{Bool}(undef, size(dcm_array))
			for z in 1:size(dcm_array, 3)
				masks_3D[:, :, z] = masks
			end
			means1 = [mean(dcm_array[dilate_mask_L_HD_3D]), mean(dcm_array[dilate_mask_M_HD_3D]), mean(dcm_array[dilate_mask_S_HD_3D]), mean(dcm_array[dilate_mask_L_MD_3D]), mean(dcm_array[dilate_mask_M_MD_3D]), mean(dcm_array[dilate_mask_S_MD_3D]), mean(dcm_array[dilate_mask_L_LD_3D]), mean(dcm_array[dilate_mask_M_LD_3D]), mean(dcm_array[dilate_mask_S_LD_3D])]

			## High Density
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
		end
	end
end

# ╔═╡ f416a14c-84bc-464a-8b1b-a722699d91c3
dfs

# ╔═╡ 1172b423-4d84-442b-93af-4c73087f0a29
begin
    new_df = vcat(dfs[1:length(dfs)]...)
    output_path = datadir("results", "masses.csv")
	save(output_path, new_df)
end

# ╔═╡ Cell order:
# ╠═05f80d21-203e-4934-8e4d-dfc843919918
# ╠═a8646d6e-5c7b-11ed-2556-73f2b13f709d
# ╠═87d68ed1-6d7b-4552-8680-4a8901964109
# ╟─4676198a-2063-4d98-9b11-af5acfa32245
# ╠═b32917fc-df2c-42eb-8e7a-cc2c564f2e22
# ╟─0b83ca47-1734-4172-b3c8-0290e07bbece
# ╟─860c6c90-b07f-472b-b364-9a28bc84434c
# ╠═502abd7e-c123-463b-8581-9b886c87fa29
# ╠═e0b52502-2075-47c0-a465-af93db8110de
# ╠═2750bb98-c026-440c-9d89-ebf49dac1676
# ╠═ba3b914d-758c-430d-8845-34f7b6c57f50
# ╠═fad04012-8c80-45b5-9b44-7f682f4a28ff
# ╠═97d59370-e42c-4b5d-9551-49281db1cc84
# ╟─c9f111b6-1e66-459f-a1ff-66b44f3ae785
# ╠═565b608a-f6a8-4579-9efd-37522dc9cde5
# ╠═e093bf43-6e8a-4f89-8117-a578712ad822
# ╠═24d92c10-9463-40a2-913c-94ed3728a169
# ╠═f416a14c-84bc-464a-8b1b-a722699d91c3
# ╠═1172b423-4d84-442b-93af-4c73087f0a29
