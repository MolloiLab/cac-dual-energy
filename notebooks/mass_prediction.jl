### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ db978680-24b3-11ed-338e-d18072c03678
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate("."); Pkg.instantiate()

    using PlutoUI, CairoMakie, Statistics, ImageMorphology, CSV, DataFrames, DICOM, DICOMUtils, CalciumScoring
	using StatsBase: quantile!
end

# ╔═╡ 7f5a24e6-e9b0-4a72-a9fa-98ab01a20125
TableOfContents()

# ╔═╡ 562b6c00-5470-4437-9e50-a75bf1ddd030
function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

# ╔═╡ cf5b409e-2f02-4dec-afe4-3cf26dab0cd8
function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=25, show_value=true)
end;

# ╔═╡ 0358fdb5-f07e-4925-8fbf-8f7928565034
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

# ╔═╡ 85849873-319b-49de-81a8-7d9949b2a093
sizes_folders = ["Small", "Medium", "Large", "Small1", "Medium1", "Large1"]

# ╔═╡ 76540649-bf00-41e5-9dcc-4620069cd1d7
sizes = sizes_folders

# ╔═╡ 8ae2a44d-8cbb-43a3-b586-4a9eb037713d
densities = ["Density1", "Density2", "Density3"]

# ╔═╡ 80a181a0-390b-493b-b15d-03fce1ac5dbd
energies_num = [80, 135]

# ╔═╡ 818de283-2626-4afd-9465-ee065dfaff8d
md"""
# Load Segmentation Masks
"""

# ╔═╡ 4ec398d1-cc7e-4124-a351-fc1058d11cf6
begin
	size0 = sizes_folders[1]
	density0 = densities[1]
	energy0 = string(energies_num[1])
	energy1 = string(energies_num[2])
end

# ╔═╡ 9ca73379-3456-404c-8893-b20be9ea6acf
size0

# ╔═╡ 013acd6e-0465-46c9-83b1-97613ebdee53
begin
	
	SIZE = size0;
	energies = ["80", "135"]
	if (SIZE == sizes[1] || SIZE == sizes[4])
		_SIZE = "small"
	elseif (SIZE == sizes[2] || SIZE == sizes[5])
		_SIZE = "medium"
	elseif (SIZE == sizes[3] || SIZE == sizes[6])
		_SIZE = "large"
	end
		
	root_new = string(
		dirname(pwd()),
        "/julia_arrays/",
        _SIZE,
        "/",
    )
    mask_L_HD = Array(CSV.read(string(root_new, "mask_L_HD.csv"), DataFrame; header=false))
	mask_M_HD = Array(CSV.read(string(root_new, "mask_M_HD.csv"), DataFrame; header=false))
    mask_S_HD = Array(CSV.read(string(root_new, "mask_S_HD.csv"), DataFrame; header=false))
	mask_L_MD = Array(CSV.read(string(root_new, "mask_L_MD.csv"), DataFrame; header=false))
	mask_M_MD = Array(CSV.read(string(root_new, "mask_M_MD.csv"), DataFrame; header=false))
	mask_S_MD = Array(CSV.read(string(root_new, "mask_S_MD.csv"), DataFrame; header=false))
    mask_L_LD = Array(CSV.read(string(root_new, "mask_L_LD.csv"), DataFrame; header=false))
    mask_M_LD = Array(CSV.read(string(root_new, "mask_M_LD.csv"), DataFrame; header=false))
    mask_S_LD = Array(CSV.read(string(root_new, "mask_S_LD.csv"), DataFrame; header=false))
end;

# ╔═╡ 4526ec32-2491-41be-91c4-03f545f9733c
masks = mask_L_HD+mask_L_MD+mask_L_LD+mask_M_HD+mask_M_MD+mask_M_LD+mask_S_HD+mask_S_MD+mask_S_LD;

# ╔═╡ 4807a059-3711-4839-96c6-a3f1ce8b3c7c
md"""
# Load Calibration Parameters
"""

# ╔═╡ c268d8f1-706f-43d8-9d3d-dc4072666451
begin
	param_base_pth = string(dirname(pwd()), "/calibration_params/")
	small_pth = string(param_base_pth,"Small.csv")
	med_pth = string(param_base_pth,"Medium.csv")
	large_pth = string(param_base_pth,"Large.csv")

	small_param = DataFrame(CSV.File(small_pth))
	med_param = DataFrame(CSV.File(med_pth))
	large_param = DataFrame(CSV.File(large_pth))
end;

# ╔═╡ a7be826f-97c1-4322-a1ad-644576d9a07a
large_param

# ╔═╡ 1d878c11-fca2-49e4-bd12-e645cc8bcb70
md"""
# Calculate Intensities
"""

# ╔═╡ fb440caa-f753-4b86-b6f2-fed47e31e94d
begin
	pth = joinpath(dirname(pwd()), "dcms_measurement_new", size0, density0, energy0)
	dcm = dcmdir_parse(pth)
	dcm_array = load_dcm_array(dcm)
end;

# ╔═╡ 0f13a6fe-36a6-4139-83ce-48f456b4ce81
begin
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
end;

# ╔═╡ 0a4fe231-8806-4a1b-a3e1-db79e1bb8dc5
md"""
## Low Energy
"""

# ╔═╡ ff366ede-a201-48da-9f59-638c9907ed72
begin
	pixel_size = DICOMUtils.get_pixel_size(dcm[1].meta)
end

# ╔═╡ af62f7d2-8b32-40de-a789-99b112be3852
begin
    masks_3D = Array{Bool}(undef, size(dcm_array))
    for z in 1:size(dcm_array, 3)
        masks_3D[:, :, z] = masks
    end
end;

# ╔═╡ 25f36ebb-0b26-4f27-b42a-3ee08afb9496
@bind v1 overlay_mask_bind(masks_3D)

# ╔═╡ bef52d62-75f0-450e-83df-9f635dad2f5b
overlay_mask_plot(dcm_array, dilate(dilate(dilate_mask_L_HD_3D)), v1, "masks overlayed")

# ╔═╡ 1a254747-0f14-4731-b014-2d59ff6a759b
begin
	means1 = [mean(dcm_array[dilate_mask_L_HD_3D]), mean(dcm_array[dilate_mask_M_HD_3D]), mean(dcm_array[dilate_mask_S_HD_3D]),
	mean(dcm_array[dilate_mask_L_MD_3D]), mean(dcm_array[dilate_mask_M_MD_3D]), mean(dcm_array[dilate_mask_S_MD_3D]), mean(dcm_array[dilate_mask_L_LD_3D]), mean(dcm_array[dilate_mask_M_LD_3D]), mean(dcm_array[dilate_mask_S_LD_3D])]
end

# ╔═╡ 10273dd0-a77f-4c43-be5a-a5696aa5688b
md"""
## High Energy
"""

# ╔═╡ 6f95f46a-bed8-4845-b7b1-d9340a7eea82
begin
	pth2 = joinpath(dirname(pwd()), "dcms_measurement_new", size0, density0, energy1)
	dcm2 = dcmdir_parse(pth2)
	dcm_array2 = load_dcm_array(dcm2)

end;

# ╔═╡ bb55626b-93a1-47d4-9cbe-f80d064ea5f9
@bind v2 overlay_mask_bind(masks_3D)

# ╔═╡ 05e8bfc9-22a6-4299-848b-02e02486a2de
overlay_mask_plot(dcm_array2, dilate(dilate(dilate_mask_L_HD_3D)), v2, "masks overlayed")

# ╔═╡ 76490868-ac50-416b-a1f6-0c6a3206b4b5
begin
	means2 = [mean(dcm_array2[dilate_mask_L_HD_3D]), mean(dcm_array2[dilate_mask_M_HD_3D]), mean(dcm_array2[dilate_mask_S_HD_3D]),
	mean(dcm_array2[dilate_mask_L_MD_3D]), mean(dcm_array2[dilate_mask_M_MD_3D]), mean(dcm_array2[dilate_mask_S_MD_3D]), mean(dcm_array2[dilate_mask_L_LD_3D]), mean(dcm_array2[dilate_mask_M_LD_3D]), mean(dcm_array2[dilate_mask_S_LD_3D])]
end

# ╔═╡ fab56fe8-4dcc-4e22-99d8-e1a0f8214537
md"""
# Calculate Predicted Densities
"""

# ╔═╡ 0a3974a9-097b-4568-a7f5-5a1ff2755092
function predict_concentration(x, y, p)
	A = p[1] + (p[2] * x) + (p[3] * y) + (p[4] * x^2) + (p[5] * x * y) + (p[6] * y^2)
	B = 1 + (p[7] * x) + (p[8] * y)
	F = A / B
end

# ╔═╡ deef95e0-a973-4c0f-8d6c-c112306567ba
x = mean(dcm_array[dilate_mask_L_HD_3D])

# ╔═╡ 68f34789-b1fe-458e-a584-46a9e4466b4c
y = mean(dcm_array2[dilate_mask_L_HD_3D])

# ╔═╡ adf21604-8031-4740-a82f-e57c04556bb0
density = predict_concentration(x, y, Array(small_param)) # mg/mL

# ╔═╡ 065e98af-e61d-444b-bb61-f4d551b5d570
calculated_intensities = hcat(means1, means2)

# ╔═╡ 62991bf5-ccbe-48e3-a1e0-314b3308d6b5


# ╔═╡ ea38774f-f3db-459c-8b89-738dc8821595
begin
	predicted_densities = zeros(9)
	
	for i in 1:9
		predicted_densities[i] = score(means1[i], means2[i], Array(small_param), CalciumScoring.MaterialDecomposition())
	end
end

# ╔═╡ 1ac7e2e1-2ff1-4db1-a2cc-99aef25694e3
predicted_densities

# ╔═╡ 643432e9-efce-4cb2-a5d1-58b6cf0db565
md"""
# Calculate Mass
"""

# ╔═╡ 7ed5bc6c-960b-4d3c-8d41-ff48d251a238
md"""
## Predicted
"""

# ╔═╡ d6347374-594d-4075-8fc8-74f4a93e1e41
voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3] * 1e-3 # cm^3

# ╔═╡ 47a7bfbd-35e7-4d30-ac37-a7d57a81e1c4
vol_small, vol_medium, vol_large = count(dilate_mask_S_HD_3D) * voxel_size, count(dilate_mask_M_HD_3D) * voxel_size, count(dilate_mask_L_HD_3D) * voxel_size #cm^3

# ╔═╡ 87040875-8c97-4e60-aa4f-057a0732c9ca
begin
	vol_slice1 = [vol_large, vol_medium, vol_small]
	vols = vcat(vol_slice1, vol_slice1, vol_slice1)# cm^3
end

# ╔═╡ 3b0c241a-29bb-47e4-8993-fabbc6da0c1a
predicted_masses = predicted_densities .* vols

# ╔═╡ d2ca7e74-ccfb-4f64-bb76-4ad5ba78a907
md"""
## Ground Truth
"""

# ╔═╡ 3dc4d4f9-15fe-4f59-a1c4-ed241fe68f77
begin
	calcium_densities = [733, 733, 733, 411, 411, 411, 151, 151, 151]
	# calcium_densities = vcat(calcium_densities_slice1, calcium_densities_slice1, calcium_densities_slice1)
end

# ╔═╡ b8e54663-48f7-4f24-b7ec-b0c6b036ea85
vol_small_gt, vol_medium_gt, vol_large_gt = π * (1/2)^2 * 3, π * (3/2)^2 * 3, π * (5/2)^2 * 3 # mm^3

# ╔═╡ 0637bdd0-6afa-4ca9-83c9-615692fa375c
begin
	vol2 = [vol_large_gt, vol_medium_gt, vol_small_gt] * 1e-3 
	vols2 = vcat(vol2, vol2, vol2) # cm^3
end

# ╔═╡ 34288a18-57a9-4c0c-acfa-c47d5d1e0a22
gt_masses = calcium_densities .* vols2 .* 3

# ╔═╡ feaa710b-265d-4892-9f2e-3ae3a8246db0
md"""
# Prepare Results Dataframe
"""

# ╔═╡ 8a7416df-4149-46e1-8104-c35569657fce
df_results = DataFrame(
	insert_sizes = ["Large", "Medium", "Small"],
	ground_truth_mass_hd = gt_masses[1:3],
	predicted_mass_hd = predicted_masses[1:3],
	ground_truth_mass_md = gt_masses[4:6],
	predicted_mass_md = predicted_masses[4:6],
	ground_truth_mass_ld = gt_masses[7:9],
	predicted_mass_ld = predicted_masses[7:9],
)

# ╔═╡ Cell order:
# ╠═db978680-24b3-11ed-338e-d18072c03678
# ╠═7f5a24e6-e9b0-4a72-a9fa-98ab01a20125
# ╟─562b6c00-5470-4437-9e50-a75bf1ddd030
# ╟─cf5b409e-2f02-4dec-afe4-3cf26dab0cd8
# ╟─0358fdb5-f07e-4925-8fbf-8f7928565034
# ╠═85849873-319b-49de-81a8-7d9949b2a093
# ╠═76540649-bf00-41e5-9dcc-4620069cd1d7
# ╠═8ae2a44d-8cbb-43a3-b586-4a9eb037713d
# ╠═80a181a0-390b-493b-b15d-03fce1ac5dbd
# ╟─818de283-2626-4afd-9465-ee065dfaff8d
# ╠═4ec398d1-cc7e-4124-a351-fc1058d11cf6
# ╠═9ca73379-3456-404c-8893-b20be9ea6acf
# ╠═013acd6e-0465-46c9-83b1-97613ebdee53
# ╠═4526ec32-2491-41be-91c4-03f545f9733c
# ╟─4807a059-3711-4839-96c6-a3f1ce8b3c7c
# ╠═c268d8f1-706f-43d8-9d3d-dc4072666451
# ╠═a7be826f-97c1-4322-a1ad-644576d9a07a
# ╟─1d878c11-fca2-49e4-bd12-e645cc8bcb70
# ╠═fb440caa-f753-4b86-b6f2-fed47e31e94d
# ╠═0f13a6fe-36a6-4139-83ce-48f456b4ce81
# ╟─0a4fe231-8806-4a1b-a3e1-db79e1bb8dc5
# ╠═ff366ede-a201-48da-9f59-638c9907ed72
# ╠═af62f7d2-8b32-40de-a789-99b112be3852
# ╟─25f36ebb-0b26-4f27-b42a-3ee08afb9496
# ╠═bef52d62-75f0-450e-83df-9f635dad2f5b
# ╠═1a254747-0f14-4731-b014-2d59ff6a759b
# ╟─10273dd0-a77f-4c43-be5a-a5696aa5688b
# ╠═6f95f46a-bed8-4845-b7b1-d9340a7eea82
# ╟─bb55626b-93a1-47d4-9cbe-f80d064ea5f9
# ╠═05e8bfc9-22a6-4299-848b-02e02486a2de
# ╠═76490868-ac50-416b-a1f6-0c6a3206b4b5
# ╟─fab56fe8-4dcc-4e22-99d8-e1a0f8214537
# ╠═0a3974a9-097b-4568-a7f5-5a1ff2755092
# ╠═deef95e0-a973-4c0f-8d6c-c112306567ba
# ╠═68f34789-b1fe-458e-a584-46a9e4466b4c
# ╠═adf21604-8031-4740-a82f-e57c04556bb0
# ╠═065e98af-e61d-444b-bb61-f4d551b5d570
# ╠═62991bf5-ccbe-48e3-a1e0-314b3308d6b5
# ╠═ea38774f-f3db-459c-8b89-738dc8821595
# ╠═1ac7e2e1-2ff1-4db1-a2cc-99aef25694e3
# ╟─643432e9-efce-4cb2-a5d1-58b6cf0db565
# ╟─7ed5bc6c-960b-4d3c-8d41-ff48d251a238
# ╠═d6347374-594d-4075-8fc8-74f4a93e1e41
# ╠═47a7bfbd-35e7-4d30-ac37-a7d57a81e1c4
# ╠═87040875-8c97-4e60-aa4f-057a0732c9ca
# ╠═3b0c241a-29bb-47e4-8993-fabbc6da0c1a
# ╟─d2ca7e74-ccfb-4f64-bb76-4ad5ba78a907
# ╠═3dc4d4f9-15fe-4f59-a1c4-ed241fe68f77
# ╠═b8e54663-48f7-4f24-b7ec-b0c6b036ea85
# ╠═0637bdd0-6afa-4ca9-83c9-615692fa375c
# ╠═34288a18-57a9-4c0c-acfa-c47d5d1e0a22
# ╟─feaa710b-265d-4892-9f2e-3ae3a8246db0
# ╠═8a7416df-4149-46e1-8104-c35569657fce
