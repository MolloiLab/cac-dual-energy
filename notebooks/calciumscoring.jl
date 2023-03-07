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

# ╔═╡ d7ea5b44-b664-40bf-a551-13c1b3746747
# ╠═╡ show_logs = false
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ d4a49edf-eb3c-48b5-b4bb-b7c650362c09
# ╠═╡ show_logs = false
begin
    using PlutoUI, Statistics, CSV, DataFrames, GLM, CairoMakie, HypothesisTests, Colors, MLJBase, DICOM, DICOMUtils, CalciumScoring, ImageMorphology, ImageFiltering, Noise
    using StatsBase: quantile!,rmsd
end

# ╔═╡ 5fc80751-2964-409c-bbef-493842890831
include(srcdir("helper_functions.jl")); include(srcdir("masks.jl"));

# ╔═╡ d296911e-0d03-44c8-82c4-d898a1353f75
TableOfContents()

# ╔═╡ 395964f9-99af-4447-bead-e7a12b7b3e4b
md"""
## Load DICOMS
"""

# ╔═╡ 97ae8b59-bd68-4054-a470-12a76b359f77
begin
    FILE_NUMBER = 1
    ENERGIES = [80,135]
	SIZES = ["Small", "Medium", "Large", "Small1"]
    SIZE = SIZES[1]
	DENSITIES = ["Density1","Density2","Density3"]
    DENSITY = DENSITIES[1]
	ENERGY = ENERGIES[2]
    TYPE = "agatston"
    root_path = datadir("dcms_measurement_new", SIZE, string(DENSITY), string(ENERGY))
end

# ╔═╡ 2f8da4af-d096-4c2f-9986-dc069bb8279b
md"""
**Everything below should be automatic, just scroll through to visually inspect that things make sense**
"""

# ╔═╡ 50dd2530-443c-47af-af1a-7a05ffb2f113
root_path

# ╔═╡ 8aebe1b8-69e8-45b1-b647-052af4cca064
dcm_path_list = dcm_list_builder(root_path)

# ╔═╡ 57305a82-80ce-4007-bf55-599cbabd5a9f
pth = dcm_path_list[1]

# ╔═╡ 27fced43-34ad-455c-9fd7-4f08307d5445
scan = basename(pth)

# ╔═╡ 9af8ed02-2002-44ef-bab2-5d0193fb6437
header, dcm_array, slice_thick_ori1 = dcm_reader(pth);

# ╔═╡ 2782beb1-82b6-4f9c-97ec-4cf271b5aa80
md"""
# Load Segmentation Masks
"""

# ╔═╡ 2c408201-89b7-4dfc-9040-db9c0e26fe1b
begin
	density0 = DENSITY
	energy0 = string(ENERGIES[1])
	energy1 = string(ENERGIES[2])
end

# ╔═╡ c61575db-db65-4436-b52d-45936b54b414
begin
	
	energies = ["80", "135"]
	if (SIZE == SIZES[1] || SIZE == SIZES[4])
		_SIZE = "small"
	elseif (SIZE == SIZES[2] || SIZE == SIZES[5])
		_SIZE = "medium"
	elseif (SIZE == SIZES[3] || SIZE == SIZES[6])
		_SIZE = "large"
	end

	root_new = datadir("julia_arrays",_SIZE)


    mask_L_HD = Array(CSV.read(string(root_new,"/", "mask_L_HD.csv"), DataFrame; header=false))
	mask_M_HD = Array(CSV.read(string(root_new,"/", "mask_M_HD.csv"), DataFrame; header=false))
    mask_S_HD = Array(CSV.read(string(root_new,"/", "mask_S_HD.csv"), DataFrame; header=false))
	mask_L_MD = Array(CSV.read(string(root_new,"/", "mask_L_MD.csv"), DataFrame; header=false))
	mask_M_MD = Array(CSV.read(string(root_new,"/", "mask_M_MD.csv"), DataFrame; header=false))
	mask_S_MD = Array(CSV.read(string(root_new,"/", "mask_S_MD.csv"), DataFrame; header=false))
    mask_L_LD = Array(CSV.read(string(root_new,"/", "mask_L_LD.csv"), DataFrame; header=false))
    mask_M_LD = Array(CSV.read(string(root_new,"/", "mask_M_LD.csv"), DataFrame; header=false))
    mask_S_LD = Array(CSV.read(string(root_new,"/", "mask_S_LD.csv"), DataFrame; header=false))
end;

# ╔═╡ e5145242-078e-4062-94ba-df2d8fdd1081
md"""
# Agatston Scoring
"""

# ╔═╡ 1c9d4ca9-6606-4836-904d-61a1d048c772
md"""
## High Density
"""

# ╔═╡ 180a94de-1ee8-4996-a2a6-db05e7c367f5
arr = dcm_array[:, :, :];

# ╔═╡ 280e2210-aee5-4ca0-b49e-540cdf75f7f2
begin
    mask_L_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_HD_3D[:, :, z] = dilate(dilate(mask_L_HD))
    end
end;

# ╔═╡ 6f88c46c-b77f-4dcb-aa65-42b2183b3ab8
mean_hu_200 = mean(arr[erode(erode(erode(erode(erode(mask_L_HD_3D)))))])

# ╔═╡ fba2f6b5-c6c4-459e-8213-ba8b46a1d74e
md"""
#### Dilated mask
"""

# ╔═╡ 79b850f4-3004-4ac6-86b8-e84c90e96f4c
dilated_mask_L_HD = dilate_mask_large(mask_L_HD_3D);

# ╔═╡ 516f3f75-5d4f-4358-b76c-370785d4becb
@bind g2 overlay_mask_bind(dilated_mask_L_HD)

# ╔═╡ 7d5035b0-b277-4ba8-a040-0d9ca8aebccf
overlay_mask_plot(arr, dilated_mask_L_HD, g2, "dilated mask")

# ╔═╡ 6e2ca817-d552-42ec-b56b-c878257c8ad7
pixel_size = DICOMUtils.get_pixel_size(header)

# ╔═╡ 7572c134-d0a5-4c22-a251-c65f92aece17
overlayed_mask_l_hd = create_mask(arr, dilated_mask_L_HD);

# ╔═╡ 35de6e5a-d6c2-4f7b-accc-763f22b8742b
alg = Agatston()

# ╔═╡ 331df9c2-88a4-479b-ab69-88c907c0fd1c
avg_mass_cals = mean([
0.000931147497044372
0.000931147497044372
0.000931147497044372
0.0009211153369385057
0.0009211153369385057
0.0009211153369385057
0.0009096665764441353
0.0009096665764441353
0.0009096665764441353
0.0009031593242901845
0.0009031593242901845
0.0009031593242901845
0.0009012665464632671
0.0009012665464632671
0.0009012665464632671
0.0008848386753622677
0.0008848386753622677
0.0008848386753622677
0.0008846639896996343
0.0008846639896996343
0.0008846639896996343
0.0008830535691512911
0.0008830535691512911
0.0008830535691512911
0.0008828619418716314
0.0008828619418716314
0.0008828619418716314
0.000871091095403193
0.000871091095403193
0.000871091095403193
0.0008232363272292931
0.0008232363272292931
0.0008232363272292931
0.0008222808145711048
0.0008222808145711048
0.0008222808145711048
0.0008205738803591024
0.0008205738803591024
0.0008205738803591024
0.0008189027838316996
0.0008189027838316996
0.0008189027838316996
0.0008151719733384375
0.0008151719733384375
0.0008151719733384375
0.0008076950603540785
0.0008076950603540785
0.0008076950603540785
0.0008000344342467945
0.0008000344342467945
0.0008000344342467945
0.000795033977195256
0.000795033977195256
0.000795033977195256
0.0007848850854924014
0.0007848850854924014
0.0007848850854924014
0.0007833934362090658
0.0007833934362090658
0.0007833934362090658
0.0007821124363623475
0.0007821124363623475
0.0007821124363623475
0.0007812337722716096
0.0007812337722716096
0.0007812337722716096
0.0007780860432477867
0.0007780860432477867
0.0007780860432477867
0.0007778239848230719
0.0007778239848230719
0.0007778239848230719
0.0007755188620679913
0.0007755188620679913
0.0007755188620679913
0.0007743691127089254
0.0007743691127089254
0.0007743691127089254
0.0007542379553525793
0.0007542379553525793
0.0007542379553525793
0.0007509459790839494
0.0007509459790839494
0.0007509459790839494
0.0007503656580249782
0.0007503656580249782
0.0007503656580249782
0.0007476133079294162
0.0007476133079294162
0.0007476133079294162
0.0007292010274957114
0.0007292010274957114
0.0007292010274957114
0.0007228701895993029
0.0007228701895993029
0.0007228701895993029
0.0007226216721221612
0.0007226216721221612
0.0007226216721221612
0.0007222519504066924
0.0007222519504066924
0.0007222519504066924
0.0007190124559613215
0.0007190124559613215
0.0007190124559613215
0.0007177458439735389
0.0007177458439735389
0.0007177458439735389
0.0007141625582687574
0.0007141625582687574
0.0007141625582687574
0.000713788198460867
0.000713788198460867
0.000713788198460867
0.0007129340558334588
0.0007129340558334588
0.0007129340558334588
0.0007112866926356895
0.0007112866926356895
0.0007112866926356895
])

# ╔═╡ 44c48ad8-c4d8-4bdf-9dcc-f3a935a78d01
pixel_size

# ╔═╡ 390da412-58c0-4b47-bff9-1ea6ab6f5b67
agat_l_hd, mass_l_hd = score(overlayed_mask_l_hd, pixel_size, 0.0006, alg; kV=80)

# ╔═╡ de4f5ab5-7a01-4979-b8bb-14fed542a3c2
volume_score_hd = agat_l_hd * sum(dilated_mask_L_HD) * (pixel_size[1] * pixel_size[2] * pixel_size[3])

# ╔═╡ ca54de09-53da-4da8-918d-8fa53c037cb6
rel_mass_score_hd = volume_score_hd * mean(arr[dilated_mask_L_HD])

# ╔═╡ 0c214ec1-db86-4406-886e-d5b55f6c6a5a
c_agat = 200/mean_hu_200

# ╔═╡ 2c828133-352e-49ea-99f9-a0ef59365f5b
md"""
## Medium Density
"""

# ╔═╡ 3ef7fea8-b7e0-431f-adc8-7ad545381333
begin
    mask_L_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_MD_3D[:, :, z] = mask_L_MD
    end
end;

# ╔═╡ da7adbc2-748b-47e8-ac26-598e4ce7721d
md"""
#### Dilated mask
"""

# ╔═╡ 11b2b292-7970-4752-9756-93bbd2949512
dilated_mask_L_MD = dilate_mask_large(mask_L_MD_3D);

# ╔═╡ 3c81c83c-c099-4130-a9b2-49d1db4b2f7a
@bind h2 overlay_mask_bind(dilated_mask_L_MD)

# ╔═╡ 190de176-9314-49a1-9ac0-c13a708ef771
overlay_mask_plot(arr, dilated_mask_L_MD, h2, "dilated mask")

# ╔═╡ b0120b6c-2968-4d55-8ab1-299618b6e625
overlayed_mask_l_md = create_mask(arr, dilated_mask_L_MD);

# ╔═╡ 08243377-b7ff-4f88-9718-944e0ea9759d
agat_l_md, mass_l_md = score(overlayed_mask_l_md, pixel_size, avg_mass_cals, alg)

# ╔═╡ 1085e838-3e80-4cdb-b2b8-ad3b2469c8b2
md"""
## Low Density
"""

# ╔═╡ fa1667b2-fd44-4457-bb39-9e7b22a11300
begin
    mask_L_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_L_LD_3D[:, :, z] = mask_L_LD
    end
end;

# ╔═╡ 9ba0a07a-693c-4668-a3d2-b96da571e1cc
md"""
#### Dilated mask
"""

# ╔═╡ 8ad5e427-d08f-43e2-9fd2-f61b5c49c13d
dilated_mask_L_LD = dilate_mask_large(mask_L_LD_3D);

# ╔═╡ f4047965-83c0-46d4-ac79-05aeee2778c5
@bind i2 overlay_mask_bind(dilated_mask_L_LD)

# ╔═╡ f95a9b98-4e7a-44a6-a0c3-944dfa164816
overlay_mask_plot(arr, dilated_mask_L_LD, i2, "dilated mask")

# ╔═╡ a2149a60-cfcd-481f-9619-251641eafdef
overlayed_mask_l_ld = create_mask(arr, dilated_mask_L_LD);

# ╔═╡ 409f65e6-93cb-4be2-a230-d330bb65d295
agat_l_ld, mass_l_ld = score(overlayed_mask_l_ld, pixel_size, avg_mass_cals, alg)

# ╔═╡ ce3daabc-5a4f-46d9-9a20-9bca62568a49
md"""
# Score Medium Inserts
"""

# ╔═╡ 9217ff13-2131-4056-a812-5d62937f8e87
md"""
## High Density
"""

# ╔═╡ f52bbb55-0923-47b4-a595-59b89b15bd9c
begin
    mask_M_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_HD_3D[:, :, z] = mask_M_HD
    end
end;

# ╔═╡ 12e24be5-e7ef-4c00-93c8-fcc606eb2405
md"""
#### Dilated mask
"""

# ╔═╡ eb07f901-9972-4037-8a23-c06e78026d87
dilated_mask_M_HD = dilate_mask_medium(mask_M_HD_3D);

# ╔═╡ 6f2d6607-0bcb-48a3-8933-62d6d4ebcf91
@bind j2 overlay_mask_bind(dilated_mask_M_HD)

# ╔═╡ b0ed646c-257c-48ae-ad27-f0c66c1c5b6e
overlay_mask_plot(arr, dilated_mask_M_HD, j2, "dilated mask")

# ╔═╡ 319103cb-3b5f-4c93-bc5d-7c2b0c783dd9
overlayed_mask_m_hd = create_mask(arr, dilated_mask_M_HD);

# ╔═╡ 50b9289f-40f7-4b96-bf5f-9a20bf5b1e23
agat_m_hd, mass_m_hd = score(overlayed_mask_m_hd, pixel_size, avg_mass_cals, alg)

# ╔═╡ f1e1b051-6cd1-491c-a89b-02a192ddba45
md"""
## Medium Density
"""

# ╔═╡ 36d605a1-08e4-4157-ba20-e353921f4f98
begin
    mask_M_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_MD_3D[:, :, z] = mask_M_MD
    end
end;

# ╔═╡ 29eb892e-80eb-4734-932d-0ffc273e50e3
md"""
#### Dilated mask
"""

# ╔═╡ d82354a0-fe61-40de-9821-e6fcc9218a53
dilated_mask_M_MD = dilate_mask_medium(mask_M_MD_3D);

# ╔═╡ f623c22f-f5bd-411d-ade2-ffe266cba132
@bind k2 overlay_mask_bind(dilated_mask_M_MD)

# ╔═╡ 160bff34-b06e-4b85-bb70-9bddceb55239
overlay_mask_plot(arr, dilated_mask_M_MD, k2, "dilated mask")

# ╔═╡ 1d833b76-45bf-43da-ac14-f6fa53c96afe
overlayed_mask_m_md = create_mask(arr, dilated_mask_M_MD);

# ╔═╡ 4452c397-188b-467c-adee-ccaaafecdadd
agat_m_md, mass_m_md = score(overlayed_mask_m_md, pixel_size, avg_mass_cals, alg)

# ╔═╡ 6921c9bd-d831-4fff-9b98-291a7e211cca
md"""
## Low Density
"""

# ╔═╡ 8627a767-8528-41d0-9f0b-6a2b4d499237
begin
    mask_M_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_M_LD_3D[:, :, z] = mask_M_LD
    end
end;

# ╔═╡ d09357d8-f182-400f-8361-f1a2a9478c8f
md"""
#### Dilated mask
"""

# ╔═╡ aada4a7e-92d5-4895-95b5-bec2f34ae75e
dilated_mask_M_LD = dilate_mask_medium(mask_M_LD_3D);

# ╔═╡ 4a609549-f0d7-422d-a9ba-b2e6dc8b1c80
@bind l2 overlay_mask_bind(dilated_mask_M_LD)

# ╔═╡ 63c761ee-64a5-4458-8a1b-2a7c9176ca9e
overlay_mask_plot(arr, dilated_mask_M_LD, l2, "dilated mask")

# ╔═╡ e7a3a2d7-287e-401d-b926-020ef6d0a5fe
overlayed_mask_m_ld = create_mask(arr, dilated_mask_M_LD);

# ╔═╡ 47034956-848f-4712-a473-909d393d93c0
agat_m_ld, mass_m_ld = score(overlayed_mask_m_ld, pixel_size, avg_mass_cals, alg)

# ╔═╡ ae8c7f5a-c8f2-4fb7-a9a1-cdd13458bc18
md"""
# Score Small Inserts
"""

# ╔═╡ cd3d9057-bf98-486b-9d7e-04ac8055e478
md"""
## High Density
"""

# ╔═╡ a98736bb-c061-4459-abaf-2751ec72e69f
begin
    mask_S_HD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_HD_3D[:, :, z] = mask_S_HD
    end
end;

# ╔═╡ 6d175194-723b-4f31-9e55-94ca33bfaa89
md"""
#### Dilated mask
"""

# ╔═╡ 47576705-df35-4de5-9e87-edcc0cb9b287
dilated_mask_S_HD = dilate_mask_small(mask_S_HD_3D);

# ╔═╡ e8585e58-f73a-4230-b31f-c6781dfda574
@bind m2 overlay_mask_bind(dilated_mask_S_HD)

# ╔═╡ bd48a0b3-48e9-468f-84b2-77030a06cb0b
overlay_mask_plot(arr, dilated_mask_S_HD, m2, "dilated mask")

# ╔═╡ 52240110-ec3f-49a7-861c-71eb580a37af
overlayed_mask_s_hd = create_mask(arr, dilated_mask_S_HD);

# ╔═╡ 0db318d5-5217-4958-9246-5ee2aaadeb4a
agat_s_hd, mass_s_hd = score(overlayed_mask_s_hd, pixel_size, avg_mass_cals, alg)

# ╔═╡ 98fddc3c-edc0-48f4-b7d3-c2bb8a755fc2
md"""
## Medium Density
"""

# ╔═╡ 08d441fa-751e-4709-a227-bd9bd373762b
begin
    mask_S_MD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_MD_3D[:, :, z] = mask_S_MD
    end
end;

# ╔═╡ a4eb7491-89f2-4367-864c-02d0086ed279
md"""
#### Dilated mask
"""

# ╔═╡ 533cf0ee-1de4-469c-94ae-03008938cf1d
dilated_mask_S_MD = dilate_mask_small(mask_S_MD_3D);

# ╔═╡ 1d8cf8b2-08ba-4c2b-999e-1012a0dfd0c6
@bind n2 overlay_mask_bind(dilated_mask_S_MD)

# ╔═╡ 0e101865-1b15-48ab-bf9d-a14b842967eb
overlay_mask_plot(arr, dilated_mask_S_MD, n2, "dilated mask")

# ╔═╡ 07ea3d6b-f1bf-457c-9174-701fc7a3ac92
overlayed_mask_s_md = create_mask(arr, dilated_mask_S_MD);

# ╔═╡ 74bcae1f-839e-4f7d-a49a-53a2fb93277c
agat_s_md, mass_s_md = score(overlayed_mask_s_md, pixel_size, avg_mass_cals, alg)

# ╔═╡ d6b95753-0631-4bb3-a3ca-4caa5da85219
md"""
## Low Density
"""

# ╔═╡ 807e9253-b69d-4937-bec0-ec2781bbd61e
begin
    mask_S_LD_3D = Array{Bool}(undef, size(arr))
    for z in 1:size(arr, 3)
        mask_S_LD_3D[:, :, z] = mask_S_LD
    end
end;

# ╔═╡ b72fed4d-1680-4030-822d-5d78b7be4829
md"""
#### Dilated mask
"""

# ╔═╡ 051e1ed0-8a4f-4979-9881-8b79d62166fb
dilated_mask_S_LD = dilate_mask_small(mask_S_LD_3D);

# ╔═╡ caf5e276-36a7-4977-9b24-c2355dd4b9ad
@bind o2 overlay_mask_bind(dilated_mask_S_LD)

# ╔═╡ 06366efe-0442-4fad-9e52-c664ee21ba90
overlay_mask_plot(arr, dilated_mask_S_LD, o2, "dilated mask")

# ╔═╡ e322f5d9-d68a-46a4-b012-23cde439ec8e
overlayed_mask_s_ld = create_mask(arr, dilated_mask_S_LD);

# ╔═╡ 4cc3a60a-1837-4cf2-a75a-9a1d0ca044c3
agat_s_ld, mass_s_ld = score(overlayed_mask_s_ld, pixel_size, avg_mass_cals, alg)

# ╔═╡ d0848371-996f-4ee4-b768-887db088acd9
md"""
# Results
"""

# ╔═╡ b68fdb27-eaff-4917-a9d6-0a0e2a54548f
begin
	calcium_densities = [733, 733, 733, 411, 411, 411, 151, 151, 151]
	# calcium_densities = vcat(calcium_densities_slice1, calcium_densities_slice1, calcium_densities_slice1)
end

# ╔═╡ 439e8907-133e-4b35-aea2-73aab0c4b832
vol_small_gt, vol_medium_gt, vol_large_gt = π * (1/2)^2 * 3, π * (3/2)^2 * 3, π * (5/2)^2 * 3 # mm^3

# ╔═╡ fd4940c7-fdce-4ad5-80d3-f8163000121f
begin
	vol2 = [vol_large_gt, vol_medium_gt, vol_small_gt] * 1e-3 
	vols2 = vcat(vol2, vol2, vol2) # cm^3
end

# ╔═╡ d0f0ad96-fc0b-453a-9925-7e79e7d467f4
gt_masses = calcium_densities .* vols2 .* 3

# ╔═╡ a1c562c2-b8ca-474f-880c-fffcc29c7028
gt_masses

# ╔═╡ b644ed06-6780-47ab-af90-4110d381e728
inserts = ["Low Density", "Medium Density", "High Density"]

# ╔═╡ abe54baa-7821-466c-87bb-573833d61ae9
md"""
## Agatston
"""

# ╔═╡ d5f896e1-d971-442f-a385-672aa88e475e
calculated_agat_large = [agat_l_ld, agat_l_md, agat_l_hd]

# ╔═╡ 81ddb672-878e-4f14-92f3-c2733c670bdf
calculated_agat_medium = [agat_m_ld, agat_m_md, agat_m_hd]

# ╔═╡ 523c1905-e137-4276-9383-f33a57d63c95
calculated_agat_small = [agat_s_ld, agat_s_md, agat_s_hd]

# ╔═╡ 9bc71cbd-5404-4e84-80ee-61a904ac9e8c
md"""
## Mass
"""

# ╔═╡ 6ae264a1-44bd-4885-a13e-7977c3d081ad
volume_gt = [7.065, 63.585, 176.625]

# ╔═╡ 5c9738a9-3ec0-4e4b-87a6-176b2637d3f1
calculated_mass_large = [mass_l_ld, mass_l_md, mass_l_hd]

# ╔═╡ 7ae7a890-530f-48cb-b394-35f79d82217e
calculated_mass_medium = [mass_m_ld, mass_m_md, mass_m_hd]

# ╔═╡ 56f6bd62-c3c9-468a-ba23-a0e582d5dae9
predicted_mass_hd = [mass_l_hd, mass_m_hd, mass_s_hd]

# ╔═╡ 90058778-9cb7-4e08-b7ec-c8ae90ddfc36
predicted_mass_md = [mass_l_md, mass_m_hd, mass_s_hd]

# ╔═╡ d0e07edf-7466-4b3f-8dc9-6b44bdaab3b3
predicted_mass_ld = [mass_l_ld, mass_m_ld, mass_s_ld]

# ╔═╡ 7f3184b2-9c4f-482e-a974-57ec0c62066f
calculated_mass_small = [mass_s_ld, mass_s_md, mass_s_hd]

# ╔═╡ 6cb8493f-193a-48b9-86e5-e15347d52c1e
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

# ╔═╡ Cell order:
# ╠═d7ea5b44-b664-40bf-a551-13c1b3746747
# ╠═d4a49edf-eb3c-48b5-b4bb-b7c650362c09
# ╠═5fc80751-2964-409c-bbef-493842890831
# ╠═d296911e-0d03-44c8-82c4-d898a1353f75
# ╟─395964f9-99af-4447-bead-e7a12b7b3e4b
# ╠═97ae8b59-bd68-4054-a470-12a76b359f77
# ╟─2f8da4af-d096-4c2f-9986-dc069bb8279b
# ╠═50dd2530-443c-47af-af1a-7a05ffb2f113
# ╠═8aebe1b8-69e8-45b1-b647-052af4cca064
# ╠═57305a82-80ce-4007-bf55-599cbabd5a9f
# ╠═27fced43-34ad-455c-9fd7-4f08307d5445
# ╠═9af8ed02-2002-44ef-bab2-5d0193fb6437
# ╟─2782beb1-82b6-4f9c-97ec-4cf271b5aa80
# ╠═2c408201-89b7-4dfc-9040-db9c0e26fe1b
# ╠═c61575db-db65-4436-b52d-45936b54b414
# ╟─e5145242-078e-4062-94ba-df2d8fdd1081
# ╟─1c9d4ca9-6606-4836-904d-61a1d048c772
# ╠═180a94de-1ee8-4996-a2a6-db05e7c367f5
# ╠═280e2210-aee5-4ca0-b49e-540cdf75f7f2
# ╠═6f88c46c-b77f-4dcb-aa65-42b2183b3ab8
# ╟─fba2f6b5-c6c4-459e-8213-ba8b46a1d74e
# ╠═79b850f4-3004-4ac6-86b8-e84c90e96f4c
# ╟─516f3f75-5d4f-4358-b76c-370785d4becb
# ╠═7d5035b0-b277-4ba8-a040-0d9ca8aebccf
# ╠═6e2ca817-d552-42ec-b56b-c878257c8ad7
# ╠═7572c134-d0a5-4c22-a251-c65f92aece17
# ╠═35de6e5a-d6c2-4f7b-accc-763f22b8742b
# ╟─331df9c2-88a4-479b-ab69-88c907c0fd1c
# ╠═44c48ad8-c4d8-4bdf-9dcc-f3a935a78d01
# ╠═390da412-58c0-4b47-bff9-1ea6ab6f5b67
# ╠═de4f5ab5-7a01-4979-b8bb-14fed542a3c2
# ╠═ca54de09-53da-4da8-918d-8fa53c037cb6
# ╠═0c214ec1-db86-4406-886e-d5b55f6c6a5a
# ╠═a1c562c2-b8ca-474f-880c-fffcc29c7028
# ╟─2c828133-352e-49ea-99f9-a0ef59365f5b
# ╠═3ef7fea8-b7e0-431f-adc8-7ad545381333
# ╟─da7adbc2-748b-47e8-ac26-598e4ce7721d
# ╠═11b2b292-7970-4752-9756-93bbd2949512
# ╟─3c81c83c-c099-4130-a9b2-49d1db4b2f7a
# ╠═190de176-9314-49a1-9ac0-c13a708ef771
# ╠═b0120b6c-2968-4d55-8ab1-299618b6e625
# ╠═08243377-b7ff-4f88-9718-944e0ea9759d
# ╟─1085e838-3e80-4cdb-b2b8-ad3b2469c8b2
# ╠═fa1667b2-fd44-4457-bb39-9e7b22a11300
# ╟─9ba0a07a-693c-4668-a3d2-b96da571e1cc
# ╠═8ad5e427-d08f-43e2-9fd2-f61b5c49c13d
# ╟─f4047965-83c0-46d4-ac79-05aeee2778c5
# ╠═f95a9b98-4e7a-44a6-a0c3-944dfa164816
# ╠═a2149a60-cfcd-481f-9619-251641eafdef
# ╠═409f65e6-93cb-4be2-a230-d330bb65d295
# ╟─ce3daabc-5a4f-46d9-9a20-9bca62568a49
# ╟─9217ff13-2131-4056-a812-5d62937f8e87
# ╠═f52bbb55-0923-47b4-a595-59b89b15bd9c
# ╟─12e24be5-e7ef-4c00-93c8-fcc606eb2405
# ╠═eb07f901-9972-4037-8a23-c06e78026d87
# ╟─6f2d6607-0bcb-48a3-8933-62d6d4ebcf91
# ╠═b0ed646c-257c-48ae-ad27-f0c66c1c5b6e
# ╠═319103cb-3b5f-4c93-bc5d-7c2b0c783dd9
# ╠═50b9289f-40f7-4b96-bf5f-9a20bf5b1e23
# ╟─f1e1b051-6cd1-491c-a89b-02a192ddba45
# ╠═36d605a1-08e4-4157-ba20-e353921f4f98
# ╟─29eb892e-80eb-4734-932d-0ffc273e50e3
# ╠═d82354a0-fe61-40de-9821-e6fcc9218a53
# ╟─f623c22f-f5bd-411d-ade2-ffe266cba132
# ╠═160bff34-b06e-4b85-bb70-9bddceb55239
# ╠═1d833b76-45bf-43da-ac14-f6fa53c96afe
# ╠═4452c397-188b-467c-adee-ccaaafecdadd
# ╟─6921c9bd-d831-4fff-9b98-291a7e211cca
# ╠═8627a767-8528-41d0-9f0b-6a2b4d499237
# ╟─d09357d8-f182-400f-8361-f1a2a9478c8f
# ╠═aada4a7e-92d5-4895-95b5-bec2f34ae75e
# ╟─4a609549-f0d7-422d-a9ba-b2e6dc8b1c80
# ╠═63c761ee-64a5-4458-8a1b-2a7c9176ca9e
# ╠═e7a3a2d7-287e-401d-b926-020ef6d0a5fe
# ╠═47034956-848f-4712-a473-909d393d93c0
# ╟─ae8c7f5a-c8f2-4fb7-a9a1-cdd13458bc18
# ╟─cd3d9057-bf98-486b-9d7e-04ac8055e478
# ╠═a98736bb-c061-4459-abaf-2751ec72e69f
# ╟─6d175194-723b-4f31-9e55-94ca33bfaa89
# ╠═47576705-df35-4de5-9e87-edcc0cb9b287
# ╟─e8585e58-f73a-4230-b31f-c6781dfda574
# ╠═bd48a0b3-48e9-468f-84b2-77030a06cb0b
# ╠═52240110-ec3f-49a7-861c-71eb580a37af
# ╠═0db318d5-5217-4958-9246-5ee2aaadeb4a
# ╟─98fddc3c-edc0-48f4-b7d3-c2bb8a755fc2
# ╠═08d441fa-751e-4709-a227-bd9bd373762b
# ╟─a4eb7491-89f2-4367-864c-02d0086ed279
# ╠═533cf0ee-1de4-469c-94ae-03008938cf1d
# ╟─1d8cf8b2-08ba-4c2b-999e-1012a0dfd0c6
# ╠═0e101865-1b15-48ab-bf9d-a14b842967eb
# ╠═07ea3d6b-f1bf-457c-9174-701fc7a3ac92
# ╠═74bcae1f-839e-4f7d-a49a-53a2fb93277c
# ╟─d6b95753-0631-4bb3-a3ca-4caa5da85219
# ╠═807e9253-b69d-4937-bec0-ec2781bbd61e
# ╟─b72fed4d-1680-4030-822d-5d78b7be4829
# ╠═051e1ed0-8a4f-4979-9881-8b79d62166fb
# ╠═caf5e276-36a7-4977-9b24-c2355dd4b9ad
# ╠═06366efe-0442-4fad-9e52-c664ee21ba90
# ╠═e322f5d9-d68a-46a4-b012-23cde439ec8e
# ╠═4cc3a60a-1837-4cf2-a75a-9a1d0ca044c3
# ╟─d0848371-996f-4ee4-b768-887db088acd9
# ╠═b68fdb27-eaff-4917-a9d6-0a0e2a54548f
# ╠═439e8907-133e-4b35-aea2-73aab0c4b832
# ╠═fd4940c7-fdce-4ad5-80d3-f8163000121f
# ╠═d0f0ad96-fc0b-453a-9925-7e79e7d467f4
# ╠═b644ed06-6780-47ab-af90-4110d381e728
# ╟─abe54baa-7821-466c-87bb-573833d61ae9
# ╠═d5f896e1-d971-442f-a385-672aa88e475e
# ╠═81ddb672-878e-4f14-92f3-c2733c670bdf
# ╠═523c1905-e137-4276-9383-f33a57d63c95
# ╟─9bc71cbd-5404-4e84-80ee-61a904ac9e8c
# ╠═6ae264a1-44bd-4885-a13e-7977c3d081ad
# ╠═5c9738a9-3ec0-4e4b-87a6-176b2637d3f1
# ╠═7ae7a890-530f-48cb-b394-35f79d82217e
# ╠═56f6bd62-c3c9-468a-ba23-a0e582d5dae9
# ╠═90058778-9cb7-4e08-b7ec-c8ae90ddfc36
# ╠═d0e07edf-7466-4b3f-8dc9-6b44bdaab3b3
# ╠═7f3184b2-9c4f-482e-a974-57ec0c62066f
# ╠═6cb8493f-193a-48b9-86e5-e15347d52c1e
