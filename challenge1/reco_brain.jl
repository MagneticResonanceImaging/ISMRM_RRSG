using PyPlot, HDF5, MRIReco

filename = "rawdata_brain_radial_96proj_12ch.h5"
data = h5read(filename, "rawdata")
traj = h5read(filename, "trajectory")
N = 256

# Convert the loaded data into an Acquisition data object
# Trajectory needs to be scaled bewtween [-0.5,0.5]. The remaining
# operations are just reshaping / changing dimension order, etc.
tr = Trajectory(reshape(traj[:,:,1:2],:,2)' ./300, 96, 512, circular=false)
dat = Array{Array{Complex{Float64},2},3}(undef,1,1,1)
dat[1,1,1] = reshape(data,12,:)'
acq = AcquisitionData(tr, dat, encodingSize=[N,N,1])

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (N,N) 
Ireco = reconstruction(acq, params)

sensitivity = estimateCoilSensitivities(Ireco)

figure(1)
imshow(abs.(Ireco[:,:,1,1,7]))


# reco parameters
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (N,N)
params[:regularization] = "L2"
params[:λ] = 0.0#1.e-3
params[:iterations] = 1
params[:solver] = "cgnr"
params[:senseMaps] = reshape(sensitivity.data, N, N, 1, 12)

IrecoSENSE = reconstruction(acq, params)


figure(2)
imshow(abs.(IrecoSENSE[:,:,1,1,1]))


