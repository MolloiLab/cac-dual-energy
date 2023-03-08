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

# ╔═╡ b7334ee2-9e52-49da-9ae2-9ea31df8b171
# ╠═╡ show_logs = false
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ 01fe7259-ef09-45a0-bf1b-7b1276aaad3d
# ╠═╡ show_logs = false
begin
    using PlutoUI, CairoMakie, Statistics, ImageMorphology, CSV, DataFrames, DICOM, DICOMUtils, LsqFit, CalciumScoring, CSVFiles
	using StatsBase: quantile!, rmsd
end

# ╔═╡ 86e2dfb0-7fbb-11ec-26c0-5d55c499ebc6
TableOfContents()

# ╔═╡ 29c32fb0-f961-4794-9f9f-0a1f91898111
md"""
#  Segmentation
"""

# ╔═╡ 4cfef276-751a-4746-9b5b-7bf38bb24b79
begin
	SIZES = ["Small", "Medium", "Large"]
	#change number to select different size and save different csv
	SIZE = SIZES[3]
end

# ╔═╡ b9d34b46-5844-4830-ae59-eef44831b138
root = datadir("dcms_calibration", "$(SIZE)")

# ╔═╡ 02e01a44-9c7b-4881-9a2c-21125da82e9a
md"""
## Low Energy
"""

# ╔═╡ 624f9db8-2fe9-4a11-b90f-402d397b3709
path = datadir(root, "80")

# ╔═╡ 4666ed27-97bf-4096-a939-b8ab204697ab
begin
	dcm = dcmdir_parse(path)
	dcm_array = load_dcm_array(dcm)
end;

# ╔═╡ 77291dcf-b147-4f13-b932-8a83bec99aa7
if SIZE == "Small"
	center_insert1, center_insert2 = 187, 318
elseif SIZE == "Medium"
	center_insert1, center_insert2 = 238, 365
elseif SIZE == "Large"
	center_insert1, center_insert2 = 285, 415
end

# ╔═╡ 407d0840-b80e-4d1e-b664-3f3798d66932
let
	f = Figure()
	ax = Makie.Axis(f[1, 1])
	heatmap!(ax, dcm_array[:, :, 5]; colormap=:grays)
	scatter!(ax, 
        center_insert1:center_insert1,
        center_insert2:center_insert2;
        markersize=10,
        color=:red,
    )
	f
end

# ╔═╡ 12b261eb-d092-4565-a761-de758796c60d
begin
	calibration_rod = zeros(25, 25, size(dcm_array, 3))
	
	for z in axes(dcm_array, 3)
		rows, cols, depth = size(dcm_array)
		half_row, half_col = center_insert1, center_insert2
		offset = 12
		row_range = half_row-offset:half_row+offset
		col_range = half_col-offset:half_col+offset	
		calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
	end
end

# ╔═╡ 913e6c16-6e06-484f-9d35-4051546c01ec
@bind a PlutoUI.Slider(1:7, show_value=true)

# ╔═╡ a899ae02-ab0c-4741-a34a-db395b982cdb
heatmap(calibration_rod[:, :, a], colormap=:grays)

# ╔═╡ 30162a78-23da-42bf-aaea-005f93e84be5
mean(calibration_rod[:, :, a])

# ╔═╡ 9fb33d34-e6f5-4956-b8b1-77ce838a335c
md"""
#### Background
"""

# ╔═╡ a2fb9d8b-2168-4a0d-bf37-64d8c7ed149c
if SIZE == "Small"
	bkg_center_insert1 = 137, 318
elseif SIZE == "Medium"
	bkg_center_insert1 = 188, 365
elseif SIZE == "Large"
	bkg_center_insert1 = 235, 415
end

# ╔═╡ 589e889e-e3b7-4971-9cb9-37bbf8e9ee1b
begin
	bkg_calibration_rod = zeros(25, 25, size(dcm_array, 3))
	
	for z in axes(dcm_array, 3)
		rows, cols, depth = size(dcm_array)
		half_row, half_col = bkg_center_insert1[1], bkg_center_insert1[2]
		offset = 12
		row_range = half_row-offset:half_row+offset
		col_range = half_col-offset:half_col+offset	
		bkg_calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
	end
end

# ╔═╡ 40430578-a0e1-400d-bbd3-e7cc281cd0be
let
	f = Figure()
	ax = Makie.Axis(f[1, 1])
	heatmap!(ax, dcm_array[:, :, 5]; colormap=:grays)
	scatter!(ax, 
        bkg_center_insert1[1]:bkg_center_insert1[1],
        bkg_center_insert1[2]:bkg_center_insert1[2];
        markersize=10,
        color=:red,
    )
	f
end

# ╔═╡ 715a1532-0ab2-4a99-b886-60ed488760a7
begin
	means1 = [mean(bkg_calibration_rod)]
	for z in axes(calibration_rod, 3)
		append!(means1, mean(calibration_rod[:, :, z]))
	end
end

# ╔═╡ a5d5d874-6ea8-42b2-a691-de6881446a26
means1

# ╔═╡ 55038797-761f-453b-866f-e13f832b011c
md"""
## High Energy
"""

# ╔═╡ a746d138-8031-4fb1-aef3-5d0e1d9f733c
path2 = datadir(root, "135")

# ╔═╡ dd8431bc-5036-4ca7-9cb3-a0b5093ced18
begin
	dcm2 = dcmdir_parse(path2)
	dcm_array2 = load_dcm_array(dcm2)
end;

# ╔═╡ c525308b-c29c-4e46-b3a0-896374de8d40
let
	f = Figure()
	ax = Makie.Axis(f[1, 1])
	heatmap!(ax, dcm_array2[:, :, 1], colormap=:grays)
	scatter!(ax, 
        center_insert1:center_insert1,
        center_insert2:center_insert2;
        markersize=10,
        color=:red,
    )
	f
end

# ╔═╡ d8ebc18a-bf32-413d-8479-2722261a07a2
begin
	calibration_rod2 = zeros(25, 25, size(dcm_array, 3))
	
	for z in axes(dcm_array2, 3)
		rows, cols, depth = size(dcm_array2)
		half_row, half_col = center_insert1, center_insert2
		offset = 12
		row_range = half_row-offset:half_row+offset
		col_range = half_col-offset:half_col+offset	
		calibration_rod2[:, :, z] .= dcm_array2[row_range, col_range, z];
	end
end

# ╔═╡ cb766fc2-7f17-409d-b06b-49f53b5fc765
md"""
#### Background
"""

# ╔═╡ 8031567b-8177-47b9-855d-a687f35feeca
begin
	bkg_calibration_rod2 = zeros(25, 25, size(dcm_array2, 3))
	
	for z in axes(dcm_array2, 3)
		rows, cols, depth = size(dcm_array2)
		half_row, half_col = bkg_center_insert1[1], bkg_center_insert1[2]
		offset = 12
		row_range = half_row-offset:half_row+offset
		col_range = half_col-offset:half_col+offset	
		bkg_calibration_rod2[:, :, z] .= dcm_array2[row_range, col_range, z];
	end
end

# ╔═╡ 1882bb40-b9ed-4d02-9a61-fd9eae9fd5e8
let
	f = Figure()
	ax = Makie.Axis(f[1, 1])
	heatmap!(ax, dcm_array2[:, :, 5]; colormap=:grays)
	scatter!(ax, 
        bkg_center_insert1[1]:bkg_center_insert1[1],
        bkg_center_insert1[2]:bkg_center_insert1[2];
        markersize=10,
        color=:red,
    )
	f
end

# ╔═╡ 52855bed-17f3-4f7b-9714-937b3f93b434
begin
	means2 = [mean(bkg_calibration_rod2)]
	for z in axes(calibration_rod2, 3)
		append!(means2, mean(calibration_rod2[:, :, z]))
	end
end

# ╔═╡ d420ce86-f83b-4685-87f9-3904a4b98750
means2

# ╔═╡ 4667f340-4096-4b92-bac9-1bae5fecd953
md"""
# Fit Parameters
"""

# ╔═╡ 14fa4199-50d0-4b67-9bfd-1640a3b82889
calculated_intensities = hcat(means1, means2) # low energy, high energy

# ╔═╡ 6782d259-3d00-4293-88a7-3a7ac5aecc87
calcium_densities = [0, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 750, 800]

# ╔═╡ 269bb44e-0dce-42dd-ba9f-9e2ce8051cb7
ps = fit_calibration(calculated_intensities, calcium_densities)

# ╔═╡ ecd323a1-ceba-4924-ae8d-f8fe71fe2a7e
md"""
# Check Results
"""

# ╔═╡ 50657646-5c99-4a9d-be0c-7d39db72aad7
begin
	predicted_densities = []
	
	for i in 1:length(calcium_densities)
		append!(
			predicted_densities, 
			score(calculated_intensities[i, 1], calculated_intensities[i, 2], ps, MaterialDecomposition()
			)
		)
	end
end

# ╔═╡ b50622d4-b15f-46ee-b476-3cfb74ddf765
df = DataFrame(
	calcium_densities = calcium_densities,
	predicted_densities = predicted_densities,
	mean_intensities_low = means1,
	mean_intensities_high = means2,
)

# ╔═╡ 61c98667-82e6-4583-a328-f2e2c48474c6
rmsd(Float64.(calcium_densities), Float64.(predicted_densities))

# ╔═╡ 3b432411-31fa-48ec-b7fc-19807b193447
md"""
# Save Parameters
"""

# ╔═╡ 65a01fd9-15ee-4b03-89af-6dfe91693571
df_params = DataFrame(p = ps)

# ╔═╡ 566aad05-2ef6-4540-80a4-cd77e849126d
file_name = string("$(SIZE)",".csv")

# ╔═╡ 54b765ae-a3d8-4242-99a4-c3c60c42d6e5
save(datadir("calibration_params",file_name),df_params)

# ╔═╡ Cell order:
# ╠═b7334ee2-9e52-49da-9ae2-9ea31df8b171
# ╠═01fe7259-ef09-45a0-bf1b-7b1276aaad3d
# ╠═86e2dfb0-7fbb-11ec-26c0-5d55c499ebc6
# ╟─29c32fb0-f961-4794-9f9f-0a1f91898111
# ╠═4cfef276-751a-4746-9b5b-7bf38bb24b79
# ╠═b9d34b46-5844-4830-ae59-eef44831b138
# ╟─02e01a44-9c7b-4881-9a2c-21125da82e9a
# ╠═624f9db8-2fe9-4a11-b90f-402d397b3709
# ╠═4666ed27-97bf-4096-a939-b8ab204697ab
# ╠═77291dcf-b147-4f13-b932-8a83bec99aa7
# ╟─407d0840-b80e-4d1e-b664-3f3798d66932
# ╠═12b261eb-d092-4565-a761-de758796c60d
# ╟─913e6c16-6e06-484f-9d35-4051546c01ec
# ╠═a899ae02-ab0c-4741-a34a-db395b982cdb
# ╠═30162a78-23da-42bf-aaea-005f93e84be5
# ╟─9fb33d34-e6f5-4956-b8b1-77ce838a335c
# ╠═a2fb9d8b-2168-4a0d-bf37-64d8c7ed149c
# ╠═589e889e-e3b7-4971-9cb9-37bbf8e9ee1b
# ╠═40430578-a0e1-400d-bbd3-e7cc281cd0be
# ╠═715a1532-0ab2-4a99-b886-60ed488760a7
# ╠═a5d5d874-6ea8-42b2-a691-de6881446a26
# ╟─55038797-761f-453b-866f-e13f832b011c
# ╠═a746d138-8031-4fb1-aef3-5d0e1d9f733c
# ╠═dd8431bc-5036-4ca7-9cb3-a0b5093ced18
# ╟─c525308b-c29c-4e46-b3a0-896374de8d40
# ╠═d8ebc18a-bf32-413d-8479-2722261a07a2
# ╟─cb766fc2-7f17-409d-b06b-49f53b5fc765
# ╠═8031567b-8177-47b9-855d-a687f35feeca
# ╟─1882bb40-b9ed-4d02-9a61-fd9eae9fd5e8
# ╠═52855bed-17f3-4f7b-9714-937b3f93b434
# ╠═d420ce86-f83b-4685-87f9-3904a4b98750
# ╟─4667f340-4096-4b92-bac9-1bae5fecd953
# ╠═14fa4199-50d0-4b67-9bfd-1640a3b82889
# ╠═6782d259-3d00-4293-88a7-3a7ac5aecc87
# ╠═269bb44e-0dce-42dd-ba9f-9e2ce8051cb7
# ╟─ecd323a1-ceba-4924-ae8d-f8fe71fe2a7e
# ╠═50657646-5c99-4a9d-be0c-7d39db72aad7
# ╠═b50622d4-b15f-46ee-b476-3cfb74ddf765
# ╠═61c98667-82e6-4583-a328-f2e2c48474c6
# ╟─3b432411-31fa-48ec-b7fc-19807b193447
# ╠═65a01fd9-15ee-4b03-89af-6dfe91693571
# ╠═566aad05-2ef6-4540-80a4-cd77e849126d
# ╠═54b765ae-a3d8-4242-99a4-c3c60c42d6e5
