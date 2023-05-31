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

# ╔═╡ 0e529354-03a3-4f5b-8cd2-88c6f8f948f5
begin
	using DrWatson
	@quickactivate "cac-dual-energy"
end

# ╔═╡ d3ed7af3-72a4-46cc-92d1-854804a50d49
begin
    using PlutoUI, CairoMakie, MAT, DICOM, DICOMUtils
end

# ╔═╡ 94002848-149d-4253-84a0-08a5ff5babf3
densities = ["15_18_22","26_29_36","52_59_73"];

# ╔═╡ 76c7a586-6a6b-4a62-88e0-c7ae83068f2d
sizes = ["small","medium","large"];

# ╔═╡ 3373927e-6c47-4a17-8acd-59f5548f823d
energies = [80,135];

# ╔═╡ edb7621b-b37a-4e57-ab37-d7fe354eda13
numbers = 1:3;

# ╔═╡ be86df26-feb2-11ed-360c-53147fd02eb7
for size1 in sizes
	for density in densities
		for energy in energies
			for num in numbers

## ------------- Inserts ------------- ##
				## Path to the original dataset of .mat files
				path = datadir(
					"mat_validation_low_density",
					string(density, "energy", string(energy), size1,"_",num, ".mat")
				)
				vars1 = matread(path)
				array1 = vars1[string("I")]
				array1 = Int16.(round.(array1))
				
				## Path to known DICOM file				
				dcm_path = datadir("sample.dcm")
				
				dcm = dcm_parse(dcm_path)
				dcm[tag"Pixel Data"] = array1
				dcm[tag"Instance Number"] = num
				dcm[tag"Rows"] = size(array1, 1)
				dcm[tag"Columns"] = size(array1, 2)

				## Path to output the newly creted DICOM files
				output_path = datadir(
					"dcms_validation_low_density",
					size1,
					density,
					string(energy)
				)
				
				if !isdir(output_path)
					mkpath(output_path)
				end
				dcm_write(string(output_path, "/", num, ".dcm"), dcm)
			end
		end
    end
end

# ╔═╡ 281aab70-940d-4559-9aac-c0f4e97df500


# ╔═╡ bbe00749-f49f-42a3-bb9d-802e01df1d7e
begin
	SIZES = ["small", "medium", "large"]
	SIZE = SIZES[3]
end

# ╔═╡ 1c796e4a-dec3-447a-adf7-43b623d07767
output_path = datadir("dcms_validation_low_density", SIZE, densities[3], "80")

# ╔═╡ e02b2bdf-8a43-4f5d-98c7-372e1df46fff
dcmdir_combined = dcmdir_parse(output_path);

# ╔═╡ 31b35e0e-a7f9-420d-8203-85e7415a8562
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ 94c988e1-c75f-408a-9ee7-7da02a1f4210
@bind c PlutoUI.Slider(1:size(vol_combined, 3); default=4, show_value=true)

# ╔═╡ 977a6942-0032-4804-a772-c85b00bb9271
if SIZE == "small"
	centers = [187, 318]
elseif SIZE == "medium"
	centers = [238, 365]
elseif SIZE == "large"
	centers = [285, 415]
end

# ╔═╡ 425193d4-fe4e-4574-be7c-9d451e35566b
centers

# ╔═╡ b9c7f9bc-efa5-40ff-8b1b-142696c41915
let
	f = Figure()
	ax = Axis(f[1, 1])
	heatmap!(vol_combined[:, :, c]; colormap=:grays)
	scatter!(centers[1]:centers[1], centers[2]:centers[2], ;markersize=5)
	f
end

# ╔═╡ Cell order:
# ╠═0e529354-03a3-4f5b-8cd2-88c6f8f948f5
# ╠═d3ed7af3-72a4-46cc-92d1-854804a50d49
# ╠═94002848-149d-4253-84a0-08a5ff5babf3
# ╠═76c7a586-6a6b-4a62-88e0-c7ae83068f2d
# ╠═3373927e-6c47-4a17-8acd-59f5548f823d
# ╠═edb7621b-b37a-4e57-ab37-d7fe354eda13
# ╠═be86df26-feb2-11ed-360c-53147fd02eb7
# ╠═281aab70-940d-4559-9aac-c0f4e97df500
# ╠═bbe00749-f49f-42a3-bb9d-802e01df1d7e
# ╠═1c796e4a-dec3-447a-adf7-43b623d07767
# ╠═e02b2bdf-8a43-4f5d-98c7-372e1df46fff
# ╠═31b35e0e-a7f9-420d-8203-85e7415a8562
# ╠═94c988e1-c75f-408a-9ee7-7da02a1f4210
# ╠═977a6942-0032-4804-a772-c85b00bb9271
# ╠═425193d4-fe4e-4574-be7c-9d451e35566b
# ╠═b9c7f9bc-efa5-40ff-8b1b-142696c41915
