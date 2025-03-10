### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ f8ecd764-fd49-11ef-30ec-db27b238f7db


# ╔═╡ 05810a95-c8bc-4d52-a38e-17e951e51b3d
# ╠═╡ disabled = true
#=╠═╡
#This code is a possible way to have multiple files be selected by user by having an array of files [i]
begin
	mkpath("data")
	num_queries = 3
    filename = map(i->"DataSets/DataTest" * string(i) * ".csv", 1:num_queries)
    if split(pwd(),'/')[end] == "test"
      filename .= "../" .* filename
    end
    need_to_query = map(fn->filesize(fn)>0,filename)
end;
  ╠═╡ =#

# ╔═╡ ba507c4f-53a4-4af5-b743-d3437b2be56d


# ╔═╡ d2a288bd-e631-4418-9ec3-be0f309858d1


# ╔═╡ 86aaad14-762a-4726-9da5-028e14b1c1a5


# ╔═╡ Cell order:
# ╠═f8ecd764-fd49-11ef-30ec-db27b238f7db
# ╠═05810a95-c8bc-4d52-a38e-17e951e51b3d
# ╠═ba507c4f-53a4-4af5-b743-d3437b2be56d
# ╠═d2a288bd-e631-4418-9ec3-be0f309858d1
# ╠═86aaad14-762a-4726-9da5-028e14b1c1a5
