using PyPlot, HDF5, MRIReco, LinearAlgebra

filename = "rawdata_brain_radial_96proj_12ch.h5"
data = permutedims(h5read(filename, "rawdata"),[3,2,1,4])
traj = permutedims(h5read(filename, "trajectory"),[3,2,1])
N = 300

#############################################
# load data and form Acquisition data object
##############################################
tr = Trajectory(reshape(traj[1:2,:,:],2,:) ./300, 96, 512, circular=false)
dat = Array{Array{Complex{Float64},2},3}(undef,1,1,1)
dat[1,1,1] = 1.e8.*reshape(data,:,12)
acqData = AcquisitionData(tr, dat, encodingSize=[N,N,1])

################################
# generate coil sensitivity maps
################################
@info "Espirit"
acqDataCart = regrid2d(acqData, (N,N); cgnr_iter=3)
sensitivity = espirit(acqDataCart,(6,6),30,eigThresh_1=0.02, eigThresh_2=0.98)

##########################
# reference reconstruction
##########################
@info "reference reco"
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:λ] = 1.e-2
params[:iterations] = 100
params[:solver] = "cgnr"
params[:senseMaps] = reshape(sensitivity, N, N, 1, 12)

img_ref = reconstruction(acqData, params).data

##############################
# undersampled reconstructions
##############################
@info "undersampled reco"
rf = [1,2,3,4]
Δ = Vector{Vector{Float64}}(undef,4)
img_cg = Vector{Array{ComplexF64,5}}(undef,4)
img_cg1 = Vector{Array{ComplexF64,5}}(undef,4)
img_sc = Vector{Array{ComplexF64,5}}(undef,4)
params[:iterations] = 20
params[:λ] = 1.e-2
for r in rf
  @info "r=$r"
  # undersample profiles
  acqDataSub = convertUndersampledData(sample_kspace(acqData, Float64.(r), "regular"))
  # SENSE reconstruction while monitoring error
  params[:reco] = "multiCoil"
  params[:solverInfo] = SolverInfo(store_solutions=true)
  img_cg[r] = reconstruction(acqDataSub, params).data
  img_cg1[r] = reshape(params[:solverInfo].x_iter[2],N,N,1,1,1)
  Δ[r] = [norm(vec(img_ref).-params[:solverInfo].x_iter[i])/norm(vec(img_ref)) for i=1:params[:iterations]+1]
  # single-channel reconstruction
  params[:reco] = "direct"
  img_sc[r] = reconstruction(acqDataSub, params).data
end

##############
# plot figures
##############
figure("brain reconstruction error")
plot(log.(Δ[1]),label="r=1")
plot(log.(Δ[2]),label="r=2")
plot(log.(Δ[3]),label="r=3")
plot(log.(Δ[4]),label="r=4")
xlabel("iteration numbers")
ylabel("log(Δ)")
xticks(collect(0:4:20))
legend()

figure("brain reconstructions",figsize=(6,6))
subplot(4,3,1); title("Single coil")
subplot(4,3,2); title("Initial")
subplot(4,3,3); title("Final")
for r in rf
  # plots
  subplot(4,3,(r-1)*3+1)
  imshow(abs.(img_sc[r][end:-1:1,end:-1:1,1,1,1]),cmap="gray")
  ylabel("r=$r")
  xticks([]);yticks([])
  subplot(4,3,(r-1)*3+2)
  imshow(abs.(img_cg1[r][end:-1:1,end:-1:1,1,1,1]),cmap="gray")
  xticks([]);yticks([])
  subplot(4,3,(r-1)*3+3)
  imshow(abs.(img_cg[r][end:-1:1,end:-1:1,1,1,1]),cmap="gray")
  xticks([]);yticks([])
end
subplots_adjust(wspace=0.05,hspace=0.05,left=0.0,bottom=0.0,right=1.0,top=0.95)
