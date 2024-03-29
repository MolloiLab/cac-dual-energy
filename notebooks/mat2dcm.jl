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

# ╔═╡ bc7570e3-19c1-401e-95c5-9be6ac01a648
using DrWatson

# ╔═╡ 1649be2d-a180-4bb2-b3a1-3e1dce110e75
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ 77de3e09-4eb5-47aa-bf5e-8cce0aa1c62b
using PlutoUI, CairoMakie, MAT, DICOM

# ╔═╡ 2ae5a446-e3fe-431d-a9f7-2b03fc3fb11a
include(srcdir("dicom_utils.jl"));

# ╔═╡ b71ba9ef-eb74-464a-a8c4-5694deb5d08e
TableOfContents()

# ╔═╡ 89dba7f9-dfdf-4702-8173-fa94d23745cb
md"""
# Calibration
"""

# ╔═╡ 8076fcea-def3-4bd5-b7fe-040a81280543
densities_cal = [3, 6, 9, 10, 16, 25, 31, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800] # calcium densities

# ╔═╡ 576e9a4d-6a5d-49b9-9af3-518097a038db
energies_cal = [80, 135]

# ╔═╡ 5b849370-e656-4cc3-96a8-27de36d565f1
sizes_cal = [30]

# ╔═╡ 3a7e00d5-3800-4d91-9c35-ca4e1b8aeb0a
for density in densities_cal
	for energy in energies_cal
		for _size in sizes_cal
			# Convert first slice into dcm
			path1 = joinpath(datadir("mats", "cal"), string(density, "rod", energy, "kV", _size, "_1.mat"))
			var1 = matread(path1)
			img1 = var1[string("I")]
			img1 = Int16.(round.(img1))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img1
			dcm[tag"Slice Thickness"] = 0.5
			dcm[tag"Instance Number"] = 1
			dcm[tag"Rows"] = size(img1, 1)
			dcm[tag"Columns"] = size(img1, 2)

			output_dir = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energy))
			if !isdir(output_dir)
				mkpath(output_dir)
			end

			output_path1 = joinpath(output_dir, "1.dcm") 
			dcm_write(output_path1, dcm)

			# Convert second slice into dcm
			path2 = joinpath(datadir("mats", "cal"), string(density, "rod", energy, "kV", _size, "_2.mat"))
			var2 = matread(path2)
			img2 = var2[string("I")]
			img2 = Int16.(round.(img2))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img2
			dcm[tag"Slice Thickness"] = 0.5
			dcm[tag"Instance Number"] = 2
			dcm[tag"Rows"] = size(img2, 1)
			dcm[tag"Columns"] = size(img2, 2)

			output_path2 = joinpath(output_dir, "2.dcm") 
			dcm_write(output_path2, dcm)

			# Convert third slice into dcm
			path3 = joinpath(datadir("mats", "cal"), string(density, "rod", energy, "kV", _size, "_3.mat"))
			var3 = matread(path3)
			img3 = var3[string("I")]
			img3 = Int16.(round.(img3))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img3
			dcm[tag"Slice Thickness"] = 0.5
			dcm[tag"Instance Number"] = 2
			dcm[tag"Rows"] = size(img3, 1)
			dcm[tag"Columns"] = size(img3, 2)

			global output_path3_cal = joinpath(output_dir, "3.dcm") 
			dcm_write(output_path3_cal, dcm)
		end
	end
end

# ╔═╡ 3f4184af-0296-4771-bd77-f9341aa86788
md"""
## Check DICOM Calibration
"""

# ╔═╡ a2fdd6ec-f55e-49ab-9c7c-81ef411b3fe3
dcmdir_combined = dcmdir_parse(dirname(output_path3_cal));

# ╔═╡ 0cc8d2cb-a2bc-4782-8c49-492a24f7f46f
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ 1d144df8-05e4-4be3-975f-c43d6491ef7b
@bind z PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ 5c43c674-ce19-4eb7-b7f8-c8b6b214f06a
x, y = 187, 260

# ╔═╡ 3aa05603-f0af-4325-b7fd-5c7a9a7888ff
let
	f = Figure()
	ax = CairoMakie.Axis(
		f[1, 1]
	)
	heatmap!(transpose(vol_combined[:, :, z]); colormap=:grays)
	scatter!(y, x; markersize = 10, color = :red)
	f
end

# ╔═╡ af5489bd-febc-4e23-a113-291fba99e9ca
md"""
# Validation
"""

# ╔═╡ e9c7a7af-f24c-47b2-be7f-198c5992a47d
densities_val = [
	"15_18_22"
	"26_29_36"
	"52_59_73"
	"110_210_310"
	"410_610_780"
] # percentage water

# ╔═╡ 5ecbc232-7733-4f1c-a8d6-c0821badaf06
energies_val = [80, 120, 135]

# ╔═╡ e10c7eea-dfcf-4269-976e-b2765a5713ec
sizes_val = ["small", "medium", "large"]

# ╔═╡ 40d56172-80af-4d7a-91f7-63bb3acf6a45
for (i, density) in enumerate(densities_val)
	for energy in energies_val
		for _size in sizes_val
			# Convert first slice into dcm
			path1 = joinpath(datadir("mats", "val"), string(density, "energy", energy, _size, "_1.mat"))
			var1 = matread(path1)
			img1 = var1[string("I")]
			img1 = Int16.(round.(img1))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img1
			dcm[tag"Instance Number"] = 1
			dcm[tag"Rows"] = size(img1, 1)
			dcm[tag"Columns"] = size(img1, 2)

			output_dir = joinpath(datadir("dcms", "val"), string(density), string(_size), string(energy))
			if !isdir(output_dir)
				mkpath(output_dir)
			end

			output_path1 = joinpath(output_dir, "1.dcm") 
			dcm_write(output_path1, dcm)

			# Convert second slice into dcm
			path2 = joinpath(datadir("mats", "val"), string(density, "energy", energy, _size, "_2.mat"))
			var2 = matread(path2)
			img2 = var2[string("I")]
			img2 = Int16.(round.(img2))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img2
			dcm[tag"Instance Number"] = 1
			dcm[tag"Rows"] = size(img2, 1)
			dcm[tag"Columns"] = size(img2, 2)

			output_path2 = joinpath(output_dir, "2.dcm") 
			dcm_write(output_path2, dcm)

			# Convert second slice into dcm
			path3 = joinpath(datadir("mats", "val"), string(density, "energy", energy, _size, "_3.mat"))
			var3 = matread(path3)
			img3 = var3[string("I")]
			img3 = Int16.(round.(img3))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img3
			dcm[tag"Instance Number"] = 1
			dcm[tag"Rows"] = size(img3, 1)
			dcm[tag"Columns"] = size(img3, 2)

			global output_path3_val = joinpath(output_dir, "3.dcm") 
			dcm_write(output_path3_val, dcm)
		end
	end
end

# ╔═╡ 198cfcd1-e23a-4e58-a074-a9381f66d34b
md"""
## Check DICOM Validation
"""

# ╔═╡ a6b0fd0f-b576-45f3-b05b-bc075f64b66d
dcm_path = "/Users/daleblack/Library/CloudStorage/GoogleDrive-djblack@uci.edu/My Drive/dev/Molloilab/cac-dual-energy/data/dcms/val/52_59_73/large/120"

# ╔═╡ 4f045d3d-d06b-424c-b0bf-332983c57c7c
dcmdir_combined_val = dcmdir_parse(dcm_path);

# ╔═╡ 65de0193-2181-4394-bf01-e53c25ff896b
dcmdir_combined_val

# ╔═╡ 38323073-5134-47d4-a6a5-eb11729d7e51
vol_combined_val = load_dcm_array(dcmdir_combined_val);

# ╔═╡ 592521e1-ce60-404f-a8e3-0dd2b958030c
@bind z1 PlutoUI.Slider(1:size(vol_combined_val, 3); default=1, show_value=true)

# ╔═╡ edc309c4-4110-477f-b375-07db0170e8a8
heatmap(transpose(vol_combined_val[:, :, z1]); colormap=:grays)

# ╔═╡ Cell order:
# ╠═bc7570e3-19c1-401e-95c5-9be6ac01a648
# ╠═1649be2d-a180-4bb2-b3a1-3e1dce110e75
# ╠═77de3e09-4eb5-47aa-bf5e-8cce0aa1c62b
# ╠═2ae5a446-e3fe-431d-a9f7-2b03fc3fb11a
# ╠═b71ba9ef-eb74-464a-a8c4-5694deb5d08e
# ╟─89dba7f9-dfdf-4702-8173-fa94d23745cb
# ╠═8076fcea-def3-4bd5-b7fe-040a81280543
# ╠═576e9a4d-6a5d-49b9-9af3-518097a038db
# ╠═5b849370-e656-4cc3-96a8-27de36d565f1
# ╠═3a7e00d5-3800-4d91-9c35-ca4e1b8aeb0a
# ╟─3f4184af-0296-4771-bd77-f9341aa86788
# ╠═a2fdd6ec-f55e-49ab-9c7c-81ef411b3fe3
# ╠═0cc8d2cb-a2bc-4782-8c49-492a24f7f46f
# ╟─1d144df8-05e4-4be3-975f-c43d6491ef7b
# ╠═5c43c674-ce19-4eb7-b7f8-c8b6b214f06a
# ╠═3aa05603-f0af-4325-b7fd-5c7a9a7888ff
# ╟─af5489bd-febc-4e23-a113-291fba99e9ca
# ╠═e9c7a7af-f24c-47b2-be7f-198c5992a47d
# ╠═5ecbc232-7733-4f1c-a8d6-c0821badaf06
# ╠═e10c7eea-dfcf-4269-976e-b2765a5713ec
# ╠═40d56172-80af-4d7a-91f7-63bb3acf6a45
# ╟─198cfcd1-e23a-4e58-a074-a9381f66d34b
# ╠═a6b0fd0f-b576-45f3-b05b-bc075f64b66d
# ╠═4f045d3d-d06b-424c-b0bf-332983c57c7c
# ╠═65de0193-2181-4394-bf01-e53c25ff896b
# ╠═38323073-5134-47d4-a6a5-eb11729d7e51
# ╟─592521e1-ce60-404f-a8e3-0dd2b958030c
# ╟─edc309c4-4110-477f-b375-07db0170e8a8
