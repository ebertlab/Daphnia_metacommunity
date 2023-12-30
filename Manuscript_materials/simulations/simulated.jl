using CSV
using DataFrames
using Distributions
using LinearAlgebra

const jsL = CSV.read("./data/jsL", header=1, delim=' ', copycols = true)
const jsM = CSV.read("./data/jsM", header=1, delim=' ', copycols = true)
const jsP = CSV.read("./data/jsP", header=1, delim=' ', copycols = true)

const gammas = CSV.read("./data/gamma", header = 1, delim = ' ', copycols = true)

# To remove interactions
# gammas[:,:] = 0

const covariates = CSV.read("../data/environmental_pca", delim = ' ', header = 1, copycols = true)
# covariates[:,:] = 0

# Full Distance matrix
const distance_matrix = convert(Matrix, CSV.read("../data/distance_matrix", delim = ' ', header = false, copycols = true))
# Grouped distance matrix
# const distance_matrix = convert(Matrix, CSV.read("../data/distance_matrix_groups", delim = ' ', header = false))

const inits = convert(Matrix, CSV.read("./data/inits", delim = ' ', header = true, copycols = true))

const nsites = 546
const nyears = 35
const nrep = 500

const rID = trunc(Int, rand()*1000000)

println("simulation #$rID")

results = zeros(Float32, (nrep, 3,2,1,nyears))
results_bysites = zeros(Float32, (nrep, 3,2, nsites))

metacommunity = zeros(Int8, (3, 2, nsites, nyears))

for r in 1:nrep

	# Initialisation

	randomSample = rand(1:20000)

	rates = similar(jsL, 0)

	basal_rates = zeros(Float64, (3,2,2,nsites))

	push!(rates, jsL[randomSample,:])
	push!(rates, jsM[randomSample,:])
	push!(rates, jsP[randomSample,:])

	distmatrix = zeros(Float64, (3, nsites, nsites))
	scaling_cst = zeros(3)

	for s in 1:3
		# Random init 
		metacommunity[s,1,:,1] = rand(Binomial(1,.2), nsites)

		# Init as 83
		# metacommunity[s,1,:,1] = deepcopy(inits[:,s])


		distmatrix[s,:,:] = exp.( - rates[s,:alpha] * distance_matrix) - I
		scaling_cst[s] = sum(distmatrix[s,:,:]) / nsites

		basal_rates[s,1,1,:] = rates[s,:p1_e] .+ rates[s,:beta1_PC1_e] .* covariates[:,:PC1] + rates[s,:beta1_PC2_e] .* covariates[:,:PC2]
		basal_rates[s,1,2,:] = rates[s,:p1_c] .+ rates[s,:beta1_PC1_c] .* covariates[:,:PC1] + rates[s,:beta1_PC2_c] .* covariates[:,:PC2]

		basal_rates[s,2,1,:] = rates[s,:p2_e] .+ rates[s,:beta2_PC1_e] .* covariates[:,:PC1] + rates[s,:beta2_PC2_e] .* covariates[:,:PC2]
		basal_rates[s,2,2,:] = rates[s,:p2_c] .+ rates[s,:beta2_PC1_c] .* covariates[:,:PC1] + rates[s,:beta2_PC2_c] .* covariates[:,:PC2]
	end

	probabilities = zeros(Float64, (2, nsites))
	local_densities = zeros(nsites)

	s = 1
	probabilities[1,:] = exp.(basal_rates[s, 1, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[1,1,2]")] .* metacommunity[2,1,:,1] .+ gammas[randomSample,Symbol("gamma_e[1,1,3]")] .* metacommunity[3,1,:,1])
	local_densities = (metacommunity[s, 1, :, 1]' * distmatrix[s,:,:]) ./ scaling_cst[s]
	probabilities[2,:] = exp.(basal_rates[s, 1, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[1,1,2]")] .* metacommunity[2,1,:,1] .+ gammas[randomSample,Symbol("gamma_c[1,1,3]")] .* metacommunity[3,1,:,1]) .* local_densities'

	metacommunity[s, 2, :, 1] = map(x -> rand(Binomial(1,x)), metacommunity[s, 1, :, 1] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
															 (1 .- metacommunity[s, 1, :, 1]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))


	s = 2
	probabilities[1,:] = exp.(basal_rates[s, 1, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[1,2,1]")] .* metacommunity[1,1,:,1] .+ gammas[randomSample,Symbol("gamma_e[1,2,3]")] .* metacommunity[3,1,:,1])
	local_densities = (metacommunity[s, 1, :, 1]' * distmatrix[s,:,:]) ./ scaling_cst[s]
	probabilities[2,:] = exp.(basal_rates[s, 1, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[1,2,1]")] .* metacommunity[1,1,:,1] .+ gammas[randomSample,Symbol("gamma_c[1,2,3]")] .* metacommunity[3,1,:,1]) .* local_densities'

	metacommunity[s, 2, :, 1] = map(x -> rand(Binomial(1,x)), metacommunity[s, 1, :, 1] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
															 (1 .- metacommunity[s, 1, :, 1]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))

	s = 3
	probabilities[1,:] = exp.(basal_rates[s, 1, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[1,3,1]")] .* metacommunity[1,1,:,1] .+ gammas[randomSample,Symbol("gamma_e[1,3,2]")] .* metacommunity[2,1,:,1])
	local_densities = (metacommunity[s, 1, :, 1]' * distmatrix[s,:,:]) ./ scaling_cst[s]
	probabilities[2,:] = exp.(basal_rates[s, 1, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[1,3,1]")] .* metacommunity[1,1,:,1] .+ gammas[randomSample,Symbol("gamma_c[1,3,2]")] .* metacommunity[2,1,:,1]) .* local_densities'

	metacommunity[s, 2, :, 1] = map(x -> rand(Binomial(1,x)), metacommunity[s, 1, :, 1] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
															 (1 .- metacommunity[s, 1, :, 1]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))

	for y in 2:nyears

		# MAJ Sample 1

		s = 1
		probabilities[1,:] = exp.(basal_rates[s, 2, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[2,1,2]")] .* metacommunity[2,2,:,y-1] .+ gammas[randomSample,Symbol("gamma_e[2,1,3]")] .* metacommunity[3,2,:,y-1])
		local_densities = (metacommunity[s, 2, :, y-1]' * distmatrix[s,:,:]) ./ scaling_cst[s]
		probabilities[2,:] = exp.(basal_rates[s, 2, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[2,1,2]")] .* metacommunity[2,2,:,y-1] .+ gammas[randomSample,Symbol("gamma_c[2,1,3]")] .* metacommunity[3,2,:,y-1]) .* local_densities'

		metacommunity[s, 1, :, y] = map(x -> rand(Binomial(1,x)), metacommunity[s, 2, :, y-1] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
																 (1 .- metacommunity[s, 2, :, y-1]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))

		s = 2
		probabilities[1,:] = exp.(basal_rates[s, 2, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[2,2,1]")] .* metacommunity[1,2,:,y-1] .+ gammas[randomSample,Symbol("gamma_e[2,2,3]")] .* metacommunity[3,2,:,y-1])
		local_densities = (metacommunity[s, 2, :, y-1]' * distmatrix[s,:,:]) ./ scaling_cst[s]
		probabilities[2,:] = exp.(basal_rates[s, 2, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[2,2,1]")] .* metacommunity[1,2,:,y-1] .+ gammas[randomSample,Symbol("gamma_c[2,2,3]")] .* metacommunity[3,2,:,y-1]) .* local_densities'

		metacommunity[s, 1, :, y] = map(x -> rand(Binomial(1,x)), metacommunity[s, 2, :, y-1] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
																 (1 .- metacommunity[s, 2, :, y-1]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))

		s = 3
		probabilities[1,:] = exp.(basal_rates[s, 2, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[2,3,1]")] .* metacommunity[1,2,:,y-1] .+ gammas[randomSample,Symbol("gamma_e[2,3,2]")] .* metacommunity[2,2,:,y-1])
		local_densities = (metacommunity[s, 2, :, y-1]' * distmatrix[s,:,:]) ./ scaling_cst[s]
		probabilities[2,:] = exp.(basal_rates[s, 2, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[2,3,1]")] .* metacommunity[1,2,:,y-1] .+ gammas[randomSample,Symbol("gamma_c[2,3,2]")] .* metacommunity[2,2,:,y-1]) .* local_densities'

		metacommunity[s, 1, :, y] = map(x -> rand(Binomial(1,x)), metacommunity[s, 2, :, y-1] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
																 (1 .- metacommunity[s, 2, :, y-1]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))


		# MAJ Sample 2

		s = 1
		probabilities[1,:] = exp.(basal_rates[s, 1, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[1,1,2]")] .* metacommunity[2,1,:,y] .+ gammas[randomSample,Symbol("gamma_e[1,1,3]")] .* metacommunity[3,1,:,y])
		local_densities = (metacommunity[s, 1, :, y]' * distmatrix[s,:,:]) ./ scaling_cst[s]
		probabilities[2,:] = exp.(basal_rates[s, 1, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[1,1,2]")] .* metacommunity[2,1,:,y] .+ gammas[randomSample,Symbol("gamma_c[1,1,3]")] .* metacommunity[3,1,:,y]) .* local_densities'

		metacommunity[s, 2, :, y] = map(x -> rand(Binomial(1,x)), metacommunity[s, 1, :, y] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
																 (1 .- metacommunity[s, 1, :, y]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))

		s = 2
		probabilities[1,:] = exp.(basal_rates[s, 1, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[1,2,1]")] .* metacommunity[1,1,:,y] .+ gammas[randomSample,Symbol("gamma_e[1,2,3]")] .* metacommunity[3,1,:,y])
		local_densities = (metacommunity[s, 1, :, y]' * distmatrix[s,:,:]) ./ scaling_cst[s]
		probabilities[2,:] = exp.(basal_rates[s, 1, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[1,2,1]")] .* metacommunity[1,1,:,y] .+ gammas[randomSample,Symbol("gamma_c[1,2,3]")] .* metacommunity[3,1,:,y]) .* local_densities'

		metacommunity[s, 2, :, y] = map(x -> rand(Binomial(1,x)), metacommunity[s, 1, :, y] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
																 (1 .- metacommunity[s, 1, :, y]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))

		s = 3
		probabilities[1,:] = exp.(basal_rates[s, 1, 1, : ] .+  gammas[randomSample,Symbol("gamma_e[1,3,1]")] .* metacommunity[1,1,:,y] .+ gammas[randomSample,Symbol("gamma_e[1,3,2]")] .* metacommunity[2,1,:,y])
		local_densities = (metacommunity[s, 1, :, y]' * distmatrix[s,:,:]) ./ scaling_cst[s]
		probabilities[2,:] = exp.(basal_rates[s, 1, 2, : ] .+  gammas[randomSample,Symbol("gamma_c[1,3,1]")] .* metacommunity[1,1,:,y] .+ gammas[randomSample,Symbol("gamma_c[1,3,2]")] .* metacommunity[2,1,:,y]) .* local_densities'

		metacommunity[s, 2, :, y] = map(x -> rand(Binomial(1,x)), metacommunity[s, 1, :, y] .* ( 1 .- (1 .- exp.(-(probabilities[1,:].+probabilities[2,:]))).*(probabilities[1,:]./(probabilities[1,:] .+ probabilities[2,:]))) +
																 (1 .- metacommunity[s, 1, :, y]) .* ( (1 .-  exp.(-(probabilities[1,:].+probabilities[2,:]))) .* (probabilities[2,:]./(probabilities[1,:] .+ probabilities[2,:]))))

	end

	results[r,:,:,:,:] = mapslices(mean, metacommunity, dims = [3])

	# results_bysites[r,:,:,:] = mapslices(mean, metacommunity[:,:,:,floor(Int, nyears/2):nyears], dims = [4])

	# For the 35y series ; keep all
	results_bysites[r,:,:,:] = mapslices(mean, metacommunity[:,:,:,:], dims = [4])

	println("On en est au $r replicat")

end

CSV.write("./results_new/L_1_$rID",  convert(DataFrames.DataFrame, results[:,1,1,1,:]))
CSV.write("./results_new/L_2_$rID",  convert(DataFrames.DataFrame, results[:,1,2,1,:]))
CSV.write("./results_new/M_1_$rID",  convert(DataFrames.DataFrame, results[:,2,1,1,:]))
CSV.write("./results_new/M_2_$rID",  convert(DataFrames.DataFrame, results[:,2,2,1,:]))
CSV.write("./results_new/P_1_$rID",  convert(DataFrames.DataFrame, results[:,3,1,1,:]))
CSV.write("./results_new/P_2_$rID",  convert(DataFrames.DataFrame, results[:,3,2,1,:]))

CSV.write("./results_new/bysite_1_L_$rID", convert(DataFrames.DataFrame, results_bysites[:,1,1,:]))
CSV.write("./results_new/bysite_2_L_$rID", convert(DataFrames.DataFrame, results_bysites[:,1,2,:]))
CSV.write("./results_new/bysite_1_M_$rID", convert(DataFrames.DataFrame, results_bysites[:,2,1,:]))
CSV.write("./results_new/bysite_2_M_$rID", convert(DataFrames.DataFrame, results_bysites[:,2,2,:]))
CSV.write("./results_new/bysite_1_P_$rID", convert(DataFrames.DataFrame, results_bysites[:,3,1,:]))
CSV.write("./results_new/bysite_2_P_$rID", convert(DataFrames.DataFrame, results_bysites[:,3,2,:]))