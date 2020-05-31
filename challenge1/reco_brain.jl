using PyPlot, HDF5, MRIReco

filename = "rawdata_brain_radial_96proj_12ch.h5"
data = h5read(filename, "rawdata")
traj = h5read(filename, "trajectory")

# Convert the loaded data into an Acquisition data object
# Trajectory needs to be scaled bewtween [-0.5,0.5]. The remaining
# operations are just reshaping / changing dimension order, etc.
tr = Trajectory(reshape(traj[:,:,1:2],:,2)' ./300, 96, 512, circular=true)
dat = Array{Array{Complex{Float64},2},3}(undef,1,1,1)
dat[1,1,1] = reshape(data,12,:)'
acq = AcquisitionData(tr, dat, encodingSize=[256,256,1])

params = Dict{Symbol, Any}()
params[:reco] = "direct"
params[:reconSize] = (256,256) 
Ireco = reconstruction(acq, params)

figure(1)
imshow(abs.(Ireco[:,:,1,1,7]))
