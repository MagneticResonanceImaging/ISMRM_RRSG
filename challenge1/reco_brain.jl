using PyPlot, HDF5, MRIReco

filename = "rawdata_brain_radial_96proj_12ch.h5"
data = h5read(filename, "rawdata")
traj = h5read(filename, "trajectory")
N = 300

# Convert the loaded data into an Acquisition data object
# Trajectory needs to be scaled bewtween [-0.5,0.5]. The remaining
# operations are just reshaping / changing dimension order, etc.
tr = Trajectory(reshape(traj[:,:,1:2],:,2)' ./300, 96, 512, circular=false)
dat = Array{Array{Complex{Float64},2},3}(undef,1,1,1)
dat[1,1,1] = 1.e8.*reshape(data,12,:)'
acqData = AcquisitionData(tr, dat, encodingSize=[N,N,1])

# single channel reconstruction of fully sampled data
params = Dict{Symbol, Any}()
params[:reco] = "standard"
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:solver] = "cgnr"
params[:λ] = 0.0
params[:iterations] = 3
Ireco = reconstruction(acqData, params)

figure(1)
imshow(abs.(Ireco[:,:,1,1,7]))

# generate coil sensitivity maps
# sensitivity = estimateCoilSensitivities(Ireco, 1.e-1)
numCoils = size(Ireco,5)
tr = CartesianTrajectory(N,N)
ft = N*FFTOp(ComplexF64,(N,N))
kdata = zeros(ComplexF64,N*N,numCoils)
for c=1:numCoils
  kdata[:,c] .= ft*vec(Ireco[:,:,1,1,c])
end
acqDataCart = AcquisitionData(tr,[kdata for i=1:1,j=1:1,k=1:1], encodingSize=[N,N,1])
sensitivity = espirit(acqDataCart,(6,6),50,eigThresh_1=0.02, eigThresh_2=0.98)

redFac = 4.0 #4.0
# we directly generate an undersampled trajectory
# this saves unnecessary computations during reconstruction
acqDataSub = convertUndersampledData(sample_kspace(acqData, redFac, "regular"))


# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:λ] = 1.e-10
params[:iterations] = 10
params[:solver] = "cgnr"
params[:senseMaps] = reshape(sensitivity, N, N, 1, 12)

IrecoSENSE = reconstruction(acqDataSub, params)


figure(2)
imshow(abs.(IrecoSENSE[:,:,1,1,1]))
