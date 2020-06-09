using PyPlot, HDF5, MRIReco, LinearAlgebra

filename = "rawdata_heart_radial_55proj_34ch.h5"
data = permutedims(h5read(filename, "rawdata"),[3,2,1,4])
traj = permutedims(h5read(filename, "trajectory"),[3,2,1])
N = 240
Nc = 34

#############################################
# load data and form Acquisition data object
##############################################
tr = Trajectory(reshape(traj[1:2,:,:],2,:) ./ N, 55, 320, circular=false)
dat = Array{Array{Complex{Float64},2},3}(undef,1,1,1)
dat[1,1,1] = 1.e8.*reshape(data,:,Nc)
acqData = AcquisitionData(tr, dat, encodingSize=[N,N,1])

################################
# generate coil sensitivity maps
################################
@info "Espirit"
acqDataCart = regrid2d(acqData, (N,N); cgnr_iter=3)
sensitivity = espirit(acqDataCart,(6,6), 30, eigThresh_1=0.02, eigThresh_2=0.96)

##############################
# undersampled reconstructions
##############################
@info "undersampled reco"
rf = [1, 5/3, 5/2, 5]
img_cg = Vector{Array{ComplexF64,5}}(undef,4)

params = Dict{Symbol, Any}()
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:solver] = "cgnr"
params[:senseMaps] = reshape(sensitivity, N, N, 1, Nc)
params[:iterations] = 20
params[:relTol] = 1.e-3
params[:λ] = 1.e-2

for (i,r) in enumerate(rf)
  @info "r=$r"
  # undersample profiles
  acqDataSub = convertUndersampledData(sample_kspace(acqData, Float64.(r), "regular",step=1))
  # SENSE reconstruction while monitoring error
  params[:reco] = "multiCoil"
  img_cg[i] = reconstruction(acqDataSub, params).data
end

##############
# plot figures
##############

figure("cardiac reconstructions",figsize=(6,1.9))
clf()
subplot(1,4,1); title("11 shots")
subplot(1,4,2); title("22 shots")
subplot(1,4,3); title("33 shots")
subplot(1,4,4); title("55 shots")
for (i,r) in enumerate(rf)

  subplot(1,4,5-i)
  imshow(abs.(img_cg[i][end:-1:1,end:-1:1,1,1,1]),cmap="gray")
  xticks([]);yticks([])
end
subplots_adjust(wspace=0.04,hspace=0.04,left=0.01,bottom=0.02,right=0.99,top=0.85)
savefig("Fig6.png",dpi=300)
