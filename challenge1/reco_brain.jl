using PyPlot, HDF5, MRIReco

filename = "rawdata_brain_radial_96proj_12ch.h5"
data = permutedims(h5read(filename, "rawdata"),[3,2,1,4])
traj = permutedims(h5read(filename, "trajectory"),[3,2,1])
N = 300

# Convert the loaded data into an Acquisition data object
# Trajectory needs to be scaled bewtween [-0.5,0.5]. The remaining
# operations are just reshaping / changing dimension order, etc.
tr = Trajectory(reshape(traj[1:2,:,:],2,:) ./300, 96, 512, circular=false)
dat = Array{Array{Complex{Float64},2},3}(undef,1,1,1)
dat[1,1,1] = 1.e8.*reshape(data,:,12)
acqData = AcquisitionData(tr, dat, encodingSize=[N,N,1])

# single channel reconstruction of fully sampled data
params = Dict{Symbol, Any}()
params[:reco] = "standard"
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:solver] = "cgnr"
params[:λ] = 0.0
params[:iterations] = 1
Ireco = reconstruction(acqData, params)

figure("single-channel reco")
imshow(abs.(Ireco[:,:,1,1,1]),cmap="gray")

# generate coil sensitivity maps
acqDataCart = regrid2d(acqData, (N,N); cgnr_iter=3)
sensitivity = espirit(acqDataCart,(6,6),30,eigThresh_1=0.02, eigThresh_2=0.98)

# generate undersampled data
redFac = 4.0
# we directly generate an undersampled trajectory
# this saves unnecessary computations during reconstruction
acqDataSub = convertUndersampledData(sample_kspace(acqData, redFac, "regular"))

# CG-SENSE reconstruction
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:λ] = 0.0
params[:iterations] = 10
params[:solver] = "cgnr"
params[:senseMaps] = reshape(sensitivity, N, N, 1, 12)

IrecoSENSE = reconstruction(acqDataSub, params)


figure(2)
imshow(abs.(IrecoSENSE[:,:,1,1,1]),cmap="gray")
