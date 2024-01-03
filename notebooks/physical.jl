### A Pluto.jl notebook ###
# v0.19.36

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

# ╔═╡ 27f94bff-b27c-4ff2-aecb-ce1ff005320e
# ╠═╡ show_logs = false
using Pkg; Pkg.instantiate()

# ╔═╡ 1e230606-179f-4875-aa86-ddb63f62b5ac
using DrWatson

# ╔═╡ aa37c780-d6c7-4b42-ab26-0ce4b8b2d82d
# ╠═╡ show_logs = false
@quickactivate "cac-dual-energy"

# ╔═╡ 6648546a-35b1-40df-8493-79284a46f609
using PlutoUI: bind, TableOfContents, Slider

# ╔═╡ cd5ede11-8d2b-43a9-bb6f-5768dd3e7f33
using CairoMakie: Figure, Axis, heatmap!, scatter!

# ╔═╡ 33307f19-7302-4c79-9c69-573996c678d6
using ImageMorphology: label_components, component_centroids

# ╔═╡ 64cb24ac-2f8f-4941-b978-6860226e0dc9
using LinearAlgebra: norm

# ╔═╡ 9d453f1c-dc43-4176-a7fa-4ea11c7ae3ee
using Statistics: mean

# ╔═╡ 00f3a2e5-9a32-4358-8611-b06c5b0b566b
using StatsBase: countmap

# ╔═╡ 76908202-063f-4a4c-b4aa-c1349da6adb7
using Clustering: kmeans, assignments

# ╔═╡ d9d39c0c-a1fc-4294-b5bf-ba2bfdcc3068
using DICOM

# ╔═╡ de0183b1-a794-4260-b3dc-a84ac4d889f5
include(srcdir("active_contours.jl"));

# ╔═╡ 64de6299-c69c-409d-bcb9-dd5ce09ae0b5
TableOfContents()

# ╔═╡ 5494e7c6-4f57-4e20-a2ec-964580b69740
md"""
# Load DICOM(s)
"""

# ╔═╡ 31840a2d-99de-4dbd-913d-d4e9ceb9459a
root = "/Volumes/USB DISK/b"

# ╔═╡ 0b3c923a-2440-446f-8169-1c3813012d5d
output_path = joinpath(root, "Cardiac 0.5", "14")

# ╔═╡ f68dd63d-fec7-48cb-b49d-536fa1eb0a14
dcms = dcmdir_parse(output_path)

# ╔═╡ 118433ed-6b47-4906-b8fe-1b4602aa50be
header = dcms[1].meta;

# ╔═╡ 59acf8dd-0157-44c2-a7e8-b4c45abe00b2
dcm_arr = header[tag"PixelData"];

# ╔═╡ d9030aa1-1ce1-4ed4-8863-7c681784ca2e
pixel_size = [0.5, 0.5, 0.5]

# ╔═╡ 4511769f-07af-4582-8f3d-f43ee8de48b4
md"""
## Visualize
"""

# ╔═╡ a4ef8928-7c3f-447f-b2d0-540d7e2d0fab
@bind a Slider(axes(dcm_arr, 3), default=160, show_value=true)

# ╔═╡ 7acaae9a-d31a-4e7f-9084-f97b71e9275a
let
	f = Figure()

	ax = Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, a], colormap=:grays)

	f
end

# ╔═╡ 2ff1c49a-7df0-499b-8f41-9c519dd31fb1
md"""
# Mask Heart Function
"""

# ╔═╡ 7ec1add5-7c6e-4bf6-8125-6efa5fe7c48e
function erode_mask(img, num_erosions)
	new_img = img
	i = 0
	while i < num_erosions
		new_img = erode(new_img)
		i += 1
	end
	return new_img
end

# ╔═╡ 77cdf367-1fe5-4458-a584-a21e26c3de44
function create_circle_mask(img, centroids, radius)
    # initialize mask with all zeros
    mask = zeros(size(img))

    # define the center of the circle
    x0, y0 = centroids[1], centroids[2]

    # set all pixels inside the circle to 1 and all pixels outside the circle to 0
    for x in axes(img, 1), y in axes(img, 2)
        if ((x - x0)^2 + (y - y0)^2) <= radius^2
            mask[x, y] = 1
        else
            mask[x, y] = 0
        end
    end
    return Bool.(mask)
end

# ╔═╡ 300e3a7b-dc1a-43b3-a50a-c955cebdb09b
function centroids_from_mask(mask)
	cc_labels = label_components(mask)
	largest_connected_component, _ = sort(collect(pairs(countmap(cc_labels[cc_labels .!= 0]))), by=x->x[2], rev=true)
	largest_connected_indices = findall(cc_labels .== largest_connected_component[1])

	new_mask = zeros(size(mask))
	for i in largest_connected_indices
		new_mask[i] = 1
	end
	centroids = Int.(round.(component_centroids(label_components(new_mask))[end]))
end

# ╔═╡ 5de49946-026c-48d0-a7d0-70c293032524
begin
	half_x, half_y = size(dcm_arr, 1) ÷ 2, size(dcm_arr, 2) ÷ 2
	init_circle = create_circle_mask(dcm_arr[:, :, 3], (half_x, half_y), 140)

	init_mask = BitArray(undef, size(dcm_arr))
	for z in axes(dcm_arr, 3)
		init_mask[:, :, z] = init_circle
	end

	init_mask = init_mask .* initial_level_set(size(init_mask))
end;

# ╔═╡ 3b0c49e1-b162-471f-ae2e-355c39356812
heart_cv = chan_vese(dcm_arr; init_level_set = init_mask);

# ╔═╡ 67d8627b-e486-4d99-805a-91c871e45a73
centroids = centroids_from_mask(heart_cv)

# ╔═╡ 3518dc88-328b-48a9-a4dd-a5d1ed919382
heart_mask = create_circle_mask(dcm_arr[:, :, 3], centroids, 100);

# ╔═╡ 7cea2363-611b-47b6-b0d0-a68a84127d81
idxs = getindex.(findall(isone, heart_mask), [1 2]);

# ╔═╡ b4b62a45-e7c3-4f8b-9962-ee7584469d83
md"""
## Visualize
"""

# ╔═╡ a3e7c703-f71d-4f23-bc16-6c399c210e0d
@bind z2 Slider(axes(dcm_arr, 3), default=130, show_value=true)

# ╔═╡ 50e5afe8-e276-4e56-9034-f566869b5021
let
	f = Figure()

	ax = Axis(f[1, 1])
	heatmap!(dcm_arr[:, :, z2], colormap=:grays)
	heatmap!(heart_mask, colormap = (:jet, 0.3))

	f
end

# ╔═╡ 1eb04f5a-593e-466c-bc32-ea3cffe260d4
md"""
# Inserts
"""

# ╔═╡ 9777cf9f-a34f-4db4-9a85-c04d7f0ec34a
dcm_heart = dcm_arr .* heart_mask;

# ╔═╡ c2450d27-f296-4bbb-b958-34a3778b3c2f
@bind d Slider(axes(dcm_heart, 3), default=130, show_value=true)

# ╔═╡ 27cb6246-d6f7-4c21-9844-31d259a5568c
let
	f = Figure()

	Axis(f[1, 1])
	heatmap!(dcm_heart[:, :, d], colormap=:grays)
	
	f
end

# ╔═╡ 1f91f759-a6b3-4095-9eb4-93944a28741f
md"""
## Endpoints
"""

# ╔═╡ eac3a8f2-7f1b-4e7c-b85e-a90151309602
function find_heart_endpoints(dcm_heart, air_threshold=-2000)
	# Find the indices of the elements in the array that are between 30 and 60
	selected_inds = findall(dcm_heart .<= air_threshold)
	
	# Create boolean array from cartesian indices
	local bool_arr4 = zeros(size(dcm_heart))
	for i in selected_inds
		bool_arr4[i] = 1
	end

	local cc_labels 
	while true
		# Use connected component labeling to identify and label all connected components
		cc_labels = label_components(bool_arr4)
		if length(unique(cc_labels)) == 3
			break
		end
		
		bool_arr4 = erode(bool_arr4)
		if (length(unique(cc_labels)) <= 2)
			@warn "Erosion failed; only background is left"
			break
		end
	end
	
	begining_slice = findlast(isone, cc_labels)[3] + 1
	end_slice = findfirst(x->x==2, cc_labels)[3] - 1
	return begining_slice, end_slice
end

# ╔═╡ 2d6c800c-10e6-47fb-a5da-600d4c294283
begining_slice, end_slice = find_heart_endpoints(dcm_heart)

# ╔═╡ 6f9baa30-8314-4f6c-a21f-1cb24ae924da
md"""
## Plane Fitting
"""

# ╔═╡ 194a806d-74c9-45a7-9078-f73400df421b
function kmeans_cleaning(points)
	points_t = points'
	kmeans_result = kmeans(points_t, 2)
	labels = assignments(kmeans_result)
	mean_distances = [mean([norm(points_t[:, i] - kmeans_result.centers[:, j]) for i in findall(x -> x == j, labels)]) for j in 1:2]
	main_cluster = argmin(mean_distances)
	cleaned_points = points_t[:, labels .== main_cluster]
	cleaned_points = cleaned_points'
end

# ╔═╡ b8ba8056-a03c-4e96-ac51-e206ead16257
function find_heart_plane(dcm_heart, endpoints; air_threshold = -300, std_threshold = 0.25)
	# Find the indices of the elements in the array that are less than
	# air_threshold and filter to exclude beginning and end slices
	remove = [collect(1:endpoints[1])..., collect(endpoints[2]:320)...]
	selected_inds = findall(dcm_heart .<= air_threshold)
	for r in remove
		selected_inds = [i for i in selected_inds if i.I[3] != r]
	end

	selected_inds = getindex.(selected_inds, [1 2 3])

	# # Clean points
	# centroid = mean(selected_inds, dims=1)
	
	# # Calculate the Euclidean distance of each point to the centroid
	# distances_to_centroid = [norm(selected_inds[i, :] .- centroid) for i in axes(selected_inds, 1)]
	
	# # Calculate mean and standard deviation of the distances
	# mean_distance = mean(distances_to_centroid)
	# std_distance = std(distances_to_centroid)
	
	# # Define the threshold as mean plus 0.5 standard deviations
	# threshold = mean_distance + (std_threshold * std_distance)
	
	# # Filter the points based on the threshold
	# filtered_points = selected_inds[distances_to_centroid .<= threshold, :]
	filtered_points = kmeans_cleaning(selected_inds)

	# Create boolean array from cartesian indices
	bool_arr = zeros(size(dcm_heart))
	for i in eachrow(filtered_points)
		bool_arr[i...] = 1
	end

	# Use connected component labeling to identify and label all connected components
	cc_labels = label_components(bool_arr)

	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	return filtered_points, cc_labels
end

# ╔═╡ 6217bc62-64a3-404b-8b2a-5ebab33e6e4d
idx_plane, cc_labels = find_heart_plane(dcm_heart, (begining_slice, end_slice));

# ╔═╡ b8c0bc97-bf63-4a31-b500-a39ca8e77b3d
function extract_plane_points(cc_labels)
	# Use the countmap function to count the number of occurrences of each value in the array, excluding 0
	counts = countmap(cc_labels[cc_labels .!= 0])
	
	# # Find the value with the most occurrences
	most_common_value, _ = sort(collect(pairs(counts)), by=x->x[2], rev=true)

	# # Find the indices of the most common value in the original array
	most_common_indices = findall(cc_labels .== most_common_value[1])
	most_common_indices = getindex.(most_common_indices, [1 2 3])

	pts = most_common_indices[1+100, :], most_common_indices[Int(round(end/2)), :], most_common_indices[end-100, :];
end;

# ╔═╡ c495ba88-fc99-4a78-ab85-c324cbca495b
centroids

# ╔═╡ 02f19951-5077-4744-9d18-4b99d6ec33e3
pts = extract_plane_points(cc_labels)

# ╔═╡ Cell order:
# ╠═1e230606-179f-4875-aa86-ddb63f62b5ac
# ╠═aa37c780-d6c7-4b42-ab26-0ce4b8b2d82d
# ╠═27f94bff-b27c-4ff2-aecb-ce1ff005320e
# ╠═6648546a-35b1-40df-8493-79284a46f609
# ╠═cd5ede11-8d2b-43a9-bb6f-5768dd3e7f33
# ╠═33307f19-7302-4c79-9c69-573996c678d6
# ╠═64cb24ac-2f8f-4941-b978-6860226e0dc9
# ╠═9d453f1c-dc43-4176-a7fa-4ea11c7ae3ee
# ╠═00f3a2e5-9a32-4358-8611-b06c5b0b566b
# ╠═76908202-063f-4a4c-b4aa-c1349da6adb7
# ╠═d9d39c0c-a1fc-4294-b5bf-ba2bfdcc3068
# ╠═de0183b1-a794-4260-b3dc-a84ac4d889f5
# ╠═64de6299-c69c-409d-bcb9-dd5ce09ae0b5
# ╟─5494e7c6-4f57-4e20-a2ec-964580b69740
# ╠═31840a2d-99de-4dbd-913d-d4e9ceb9459a
# ╠═0b3c923a-2440-446f-8169-1c3813012d5d
# ╠═f68dd63d-fec7-48cb-b49d-536fa1eb0a14
# ╠═59acf8dd-0157-44c2-a7e8-b4c45abe00b2
# ╠═118433ed-6b47-4906-b8fe-1b4602aa50be
# ╠═d9030aa1-1ce1-4ed4-8863-7c681784ca2e
# ╟─4511769f-07af-4582-8f3d-f43ee8de48b4
# ╟─a4ef8928-7c3f-447f-b2d0-540d7e2d0fab
# ╟─7acaae9a-d31a-4e7f-9084-f97b71e9275a
# ╟─2ff1c49a-7df0-499b-8f41-9c519dd31fb1
# ╠═7ec1add5-7c6e-4bf6-8125-6efa5fe7c48e
# ╠═77cdf367-1fe5-4458-a584-a21e26c3de44
# ╠═300e3a7b-dc1a-43b3-a50a-c955cebdb09b
# ╠═5de49946-026c-48d0-a7d0-70c293032524
# ╠═3b0c49e1-b162-471f-ae2e-355c39356812
# ╠═67d8627b-e486-4d99-805a-91c871e45a73
# ╠═3518dc88-328b-48a9-a4dd-a5d1ed919382
# ╠═7cea2363-611b-47b6-b0d0-a68a84127d81
# ╟─b4b62a45-e7c3-4f8b-9962-ee7584469d83
# ╟─a3e7c703-f71d-4f23-bc16-6c399c210e0d
# ╟─50e5afe8-e276-4e56-9034-f566869b5021
# ╟─1eb04f5a-593e-466c-bc32-ea3cffe260d4
# ╠═9777cf9f-a34f-4db4-9a85-c04d7f0ec34a
# ╟─c2450d27-f296-4bbb-b958-34a3778b3c2f
# ╟─27cb6246-d6f7-4c21-9844-31d259a5568c
# ╟─1f91f759-a6b3-4095-9eb4-93944a28741f
# ╠═eac3a8f2-7f1b-4e7c-b85e-a90151309602
# ╠═2d6c800c-10e6-47fb-a5da-600d4c294283
# ╟─6f9baa30-8314-4f6c-a21f-1cb24ae924da
# ╠═194a806d-74c9-45a7-9078-f73400df421b
# ╠═b8ba8056-a03c-4e96-ac51-e206ead16257
# ╠═6217bc62-64a3-404b-8b2a-5ebab33e6e4d
# ╠═b8c0bc97-bf63-4a31-b500-a39ca8e77b3d
# ╠═c495ba88-fc99-4a78-ab85-c324cbca495b
# ╠═02f19951-5077-4744-9d18-4b99d6ec33e3
