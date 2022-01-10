import torch
import scipy.io as io
import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D

tmp = io.loadmat('/home/fs0/qijia/code/SimTraj/trajectorysimulation/pycode/opt.mat')
radial_tops = tmp['radial_tops']
base_endpoint = tmp['base_endpoint'].transpose()
ileaves = np.shape(radial_tops)[0]
print(base_endpoint.dtype, base_endpoint.shape)
print(radial_tops.dtype, radial_tops.shape)
print(ileaves)
assert base_endpoint.shape == (3,1)



# fig = plt.figure()
# ax = Axes3D(fig)
# ax.scatter(radial_tops[:,0],radial_tops[:,1],radial_tops[:,2])
# ax.set_xlim3d(-1,1)
# ax.set_ylim3d(-1,1)
# ax.set_zlim3d(-1,1)
# plt.show()

xdir = radial_tops[:,0]
ydir = radial_tops[:,1]
zdir = radial_tops[:,2]
Rshot = torch.zeros(ileaves,3,3,dtype=float)
utu = torch.zeros(ileaves,3,3)
usk = torch.zeros(ileaves,3,3)
base_endpoint = torch.tensor(base_endpoint, dtype=float).squeeze()
for pos in range(ileaves):
    phi = math.acos( zdir[pos] / math.sqrt(xdir[pos]**2 + ydir[pos]**2 + zdir[pos]**2) )
    theta = math.pi/2 + math.atan2( ydir[pos], xdir[pos])
    
    ux = xdir[pos]
    uy = ydir[pos]
    uz = zdir[pos]
    
    Rphi = torch.tensor([[1, 0, 0],
                            [0, math.cos(phi), -math.sin(phi)],
                            [0, math.sin(phi), math.cos(phi)]])
    Rtheta = torch.tensor([[math.cos(theta), -math.sin(theta), 0],
                            [math.sin(theta),  math.cos(theta), 0],
                            [0, 0, 1]])
    
    utu[pos,:,:] = torch.tensor([[ux**2, ux*uy, ux*uz],
                                    [ux*uy, uy**2, uy*uz],
                                    [ux*uz, uy*uz, uz**2]])
    usk[pos,:,:] = torch.tensor([[0, -uz, -uy],
                                    [uz, 0, -ux],
                                    [-uy, ux, 0]])
    Rshot[pos,:,:] = torch.matmul(Rtheta, Rphi)
endpoints = torch.matmul(Rshot, base_endpoint)

# fig = plt.figure()
# ax = Axes3D(fig)
# ax.scatter(endpoints[:,0],endpoints[:,1],endpoints[:,2])
# ax.set_xlim3d(-1,1)
# ax.set_ylim3d(-1,1)
# ax.set_zlim3d(-1,1)
# plt.show()

torch.autograd.set_detect_anomaly(True)
class optimize():
    def __init__(self, ileaves, radial_tops, base_endpoint, max_iter = 1000, lr=0.1):

        self.omega_all = torch.rand(ileaves,requires_grad=True)
        self.optimizer = torch.optim.Adam([self.omega_all],lr=lr)
        self.lr_scheduler = torch.optim.lr_scheduler.StepLR(self.optimizer, step_size=100, gamma=0.5)
        self.max_iter = max_iter
        self.base_endpoint = base_endpoint.float()
        xdir = radial_tops[:,0]
        ydir = radial_tops[:,1]
        zdir = radial_tops[:,2]
        self.ileaves = ileaves
        self.Rshot = torch.zeros(ileaves,3,3)
        self.utu = torch.zeros(ileaves,3,3)
        self.usk = torch.zeros(ileaves,3,3)
        self.R = torch.zeros(ileaves,3,3)
        self.endpoints = torch.zeros(ileaves, 3)
        for pos in range(self.ileaves):
            phi = math.acos( zdir[pos] / math.sqrt(xdir[pos]**2 + ydir[pos]**2 + zdir[pos]**2) )
            theta = math.pi/2 + math.atan2( ydir[pos], xdir[pos])
            
            ux = xdir[pos]
            uy = ydir[pos]
            uz = zdir[pos]
            
            Rphi = torch.tensor([[1, 0, 0],
                                 [0, math.cos(phi), -math.sin(phi)],
                                 [0, math.sin(phi), math.cos(phi)]])
            Rtheta = torch.tensor([[math.cos(theta), -math.sin(theta), 0],
                                   [math.sin(theta),  math.cos(theta), 0],
                                   [0, 0, 1]])
            
            self.utu[pos,:,:] = torch.tensor([[ux*ux, ux*uy, ux*uz],
                                            [ux*uy, uy*uy, uy*uz],
                                            [ux*uz, uy*uz, uz*uz]])
            self.usk[pos,:,:] = torch.tensor([[0, -uz, uy],
                                            [uz, 0, -ux],
                                            [-uy, ux, 0]])
            self.Rshot[pos,:,:] = torch.matmul(Rtheta, Rphi)

    def get_rotation_matrix(self):
        cosomega = torch.cos(2 * math.pi * self.omega_all.unsqueeze(1).unsqueeze(2))
        sinomega = torch.sin(2 * math.pi * self.omega_all.unsqueeze(1).unsqueeze(2))
        Raxis = self.utu + cosomega * (torch.eye(3).unsqueeze(0) - self.utu) + sinomega * self.usk
        self.R = torch.matmul(Raxis, self.Rshot)
        self.endpoints = torch.matmul(self.R, self.base_endpoint.squeeze())

    def calc_cost(self):
        tmp = self.endpoints.unsqueeze(1) - self.endpoints.unsqueeze(0)
        dist = torch.sqrt(torch.sum(tmp*tmp, 2)+0.001) + torch.eye(self.ileaves)
        self.cost = torch.sum( 1 / dist ) / (self.ileaves**2)


    def optim(self):
        for iter in range(self.max_iter):
            self.get_rotation_matrix()
            self.calc_cost()
            
            self.optimizer.zero_grad()
            self.cost.backward()
            with torch.no_grad():
                print('iter : ', iter, ' cost = ', self.cost.item())
            self.optimizer.step()
            self.lr_scheduler.step()

opt_obj = optimize(ileaves, radial_tops, base_endpoint, lr=1, max_iter=610)
opt_obj.optim()
io.savemat('/home/fs0/qijia/code/SimTraj/trajectorysimulation/pycode/optR.mat',{'optR':opt_obj.R.detach().permute(1,2,0).numpy(), 'endpoints':opt_obj.endpoints.detach().numpy(),'omega_all':opt_obj.endpoints.detach().squeeze().numpy()})