import torch
from collections import OrderedDict
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from pyDOE import lhs
import time
import lifelines

# set seed
seed = 20240705
np.random.seed(seed)
torch.manual_seed(seed)  

# CUDA support
if torch.cuda.is_available():
    device = torch.device('cuda')
    torch.cuda.manual_seed(seed)  
    torch.cuda.manual_seed_all(seed)
else:
    device = torch.device('cpu')

# deep neural network
class DNN(torch.nn.Module):
    def __init__(self, layers):
        super(DNN, self).__init__()

        self.depth = len(layers) - 1
        self.activation = torch.nn.Tanh() # activation tanh
        layer_list = []

        # add layers with activation
        for i in range(self.depth - 1):
            linear_layer = torch.nn.Linear(layers[i], layers[i + 1])
            # initialize weights with Xavier normal initialization
            torch.nn.init.xavier_normal_(linear_layer.weight, gain=np.sqrt(2))
            # initialize biases to zero
            torch.nn.init.constant_(linear_layer.bias, 0)
            # Append the layer and activation to the list
            layer_list.append(('layer_%d' % i, linear_layer))
            layer_list.append(('activation_%d' % i, self.activation))

        # add the final linear layer without activation
        final_linear_layer = torch.nn.Linear(layers[-2], layers[-1])
        torch.nn.init.xavier_normal_(final_linear_layer.weight, gain=np.sqrt(2))
        torch.nn.init.constant_(final_linear_layer.bias, 0)
        layer_list.append(('layer_%d' % (self.depth - 1), final_linear_layer))

        # Convert list to OrderedDict and create the sequential model
        layerDict = OrderedDict(layer_list)
        self.layers = torch.nn.Sequential(layerDict)

    def forward(self, x):
        return self.layers(x)

# the physics-informed deep neural network
class PhysicsInformedNN():
    def __init__(self, X_u_train, uM1_train_list, uTC_train_list, X_f_train, Exact_cyto_list, Exact_dose_list,
                 cyto_0_list, train_data_num, FP_layers, CC_layers, t_cyto, dt_cyto, dc, lb, ub):
        # parameters from the MSABM and PDEs model (Zheng, Y., et al., 2018.)
        self.r_TC = 2.41e-5  
        self.d_TC = 2.886e-7 
        self.dM1_TC = 7.1055e-6 
        self.p_M21 = 0.15552/8 
        self.a_C = 0.675  
        self.a_I = 0.45  
        self.K_1 = 0.5
        self.K_2 = 0.5
        self.K_M12 = 0.02
        self.Kd_M12 = 14.3
        self.A_0 = 0
        self.S_A = 1.2905e-7*3
        self.K_I = 5.7

        self.max_CSF1 = 4e-12
        self.max_EGF = 6e-13
        self.max_IGF1 = 1.2e-13
        self.S_CSF1 = 2.5e-14 / self.max_CSF1
        self.S_EGF = 5e-15 / self.max_EGF
        self.S_IGF1 = 5e-15 / self.max_IGF1
        self.d_CSF1 = 5e-4
        self.d_EGF = 2.5e-4
        self.d_IGF1 = 2.5e-4

        self.q_c = 0.8
        self.eta_c = 2e-3 
        self.mu_c = 2e-3 

        self.q_i = 0.8
        self.eta_i = 2e-3 
        self.mu_i = 2e-3 

        self.C_u = torch.tensor(X_u_train[:, 0:1], dtype=torch.float32, requires_grad=True).float().to(device)
        self.t_u = torch.tensor(X_u_train[:, 1:2], dtype=torch.float32, requires_grad=True).float().to(device)
        self.uM1 = torch.tensor(uM1_train_list, dtype=torch.float32, requires_grad=True).float().to(device)
        self.uTC = torch.tensor(uTC_train_list, dtype=torch.float32, requires_grad=True).float().to(device)
        self.C_f = torch.tensor(X_f_train[:, 0:1], dtype=torch.float32, requires_grad=True).float().to(device)
        self.t_f = torch.tensor(X_f_train[:, 1:2], dtype=torch.float32, requires_grad=True).float().to(device)
        self.cyto = torch.tensor(Exact_cyto_list, dtype=torch.float32, requires_grad=True).float().to(device)
        self.dose = torch.tensor(Exact_dose_list, dtype=torch.float32, requires_grad=True).float().to(device)
        self.cyto_0 = torch.tensor(cyto_0_list, dtype=torch.float32, requires_grad=True).float().to(device)

        self.train_data_num = train_data_num
        self.t_cyto = torch.tensor(t_cyto, dtype=torch.float32, requires_grad=True).float().to(device)

        self.dt = torch.tensor(dt_cyto, dtype=torch.float32, requires_grad=True).float().to(device)
        self.dc = torch.tensor(dc, dtype=torch.float32, requires_grad=True).float().to(device)
        self.c_max = 1/self.dc

        self.lb = torch.tensor(lb, dtype=torch.float32, requires_grad=True).float().to(device)
        self.ub = torch.tensor(ub, dtype=torch.float32, requires_grad=True).float().to(device)
   
        # fokker-planck parameters
        self.log_p_cyto = torch.empty(16, dtype = torch.float32).float().requires_grad_(True).to(device)
        self.log_p_TC = torch.empty(10, dtype = torch.float32).float().requires_grad_(True).to(device)
        self.log_p_M1 = torch.empty(8, dtype = torch.float32).float().requires_grad_(True).to(device)
        self.log_p_sigma = torch.empty(2, dtype = torch.float32).float().requires_grad_(True).to(device)

        # load trained parameters
        #additional_params = torch.load('surrogate_model/additional_params.pth')  
        #initial_log_p_cyto = additional_params['log_p_cyto']
        #initial_log_p_TC = additional_params['log_p_CC']
        #initial_log_p_M1 = additional_params['log_p_M1']
        #initial_log_p_sigma = additional_params['log_p_sigma']

        # initialize trainable parameters
        initial_log_p_cyto = torch.zeros_like(self.log_p_cyto) # parameters in odes of cytokines and drugs
        initial_log_p_TC = torch.zeros_like(self.log_p_TC) # parameters in sde of tumor cell
        initial_log_p_M1 = torch.zeros_like(self.log_p_M1) # parameters in sde of M1 macrophages
        initial_log_p_sigma = torch.zeros_like(self.log_p_sigma) # parameters of noise diffusion rate
        self.log_p_cyto  = torch.nn.Parameter(initial_log_p_cyto).to(device)
        self.log_p_TC = torch.nn.Parameter(initial_log_p_TC).to(device)
        self.log_p_M1 = torch.nn.Parameter(initial_log_p_M1).to(device)
        self.log_p_sigma = torch.nn.Parameter(initial_log_p_sigma).to(device)
        
        # DNN for Fokker-Planck equation of tumor cell 
        self.dnn_TC = DNN(FP_layers).to(device)
        # load trained parameters
        #para_TC = torch.load('surrogate_model/dnn_TC_model.pth')
        #self.dnn_TC.load_state_dict(para_TC, strict=False)   

        # DNN for M1 macrophages if necessary
        #self.dnn_M1 = DNN(FP_layers).to(device)
        #para_M1 = torch.load('surrogate_model/dnn_M1_model1.pth')
        #self.dnn_M1.load_state_dict(para_M1, strict=False)

        # DNN for SDE of tumor cell density 
        self.dnn_fTC = DNN(CC_layers).to(device)
        # load trained parameters
        #para_fTC = torch.load('surrogate_model/dnn_fTC_model.pth')      
        #self.dnn_fTC.load_state_dict(para_fTC, strict=False)

        # DNN for SDE of M1 macrophage density 
        self.dnn_fM1 = DNN(CC_layers).to(device)
        # load trained parameters
        #para_fM1 = torch.load('surrogate_model/dnn_fM1_model.pth')      
        #self.dnn_fM1.load_state_dict(para_fM1, strict=False)
        
        # set trainable parameters in optimizer
        parameters_1 = list(self.dnn_fTC.parameters()) + list(self.dnn_fM1.parameters()) + \
                            [self.log_p_TC, self.log_p_M1]
        parameters_2 = list(self.dnn_TC.parameters()) + [self.log_p_sigma]
                             #list(self.dnn_M1.parameters()) # if necessary
        
        parameters_1 = list(set(parameters_1))
        parameters_2 = list(set(parameters_2))

        # set LBFGS optimizer for traing parameters 
        # in SDE of TC, M1, cytokines and drugs
        self.optimizer_LBFGS_cyto = torch.optim.LBFGS(
            parameters_1,
            lr=1.0,
            max_iter=500, 
            max_eval=500,  
            history_size=50,
            tolerance_grad=1e-5,
            tolerance_change=1.0 * np.finfo(float).eps,
            line_search_fn="strong_wolfe" 
        )
        self.iter_cyto = 0

        # set Adam and LBFGS optimizer for training parameters 
        # in DNN of Fokker-Planck equation 
        # and SDE of noise diffusion rate
        self.optimizer_Adam = torch.optim.Adam(parameters_2, lr = 1e-3)
        self.optimizer_LBFGS = torch.optim.LBFGS(
            parameters_2,
            lr=1.0,
            max_iter=50000, 
            max_eval=50000,  
            history_size=50,
            tolerance_grad=1e-5,
            tolerance_change=1.0 * np.finfo(float).eps,
            line_search_fn="strong_wolfe" 
        )
        self.iter = 0
        
    # the function used for mapping time points to time sequence
    # due to the fact that prior data output by MSABM is a discrete sequence
    def mapping(self, f, t):
        t_ceil = torch.ceil((t)/self.dt + 0.5).long() - 1
        t_ceil = torch.clamp(t_ceil, 0, self.t_cyto.shape[0]-1)
        result = f[t_ceil, :].squeeze(1)
        return result
    
    # ODE solution of CSF1 
    def f_CSF1(self, c_CSF1_0, c_CSF1, dt, C_T):
        c = torch.zeros_like(c_CSF1)
        tt = dt*86400/2
        a = torch.exp(self.log_p_cyto[7])*self.S_CSF1*C_T       
        b = a + torch.exp(self.log_p_cyto[9])*self.d_CSF1
        c[0, :] = (c_CSF1_0 - a[0, :]/(b[0, :]))*torch.exp(-(b[0, :])*tt) + a[0, :]/(b[0, :])
        for k in range(c.shape[0]-1):
            c[k+1, :] = (c[k, :] - a[k+1, :]/(b[k+1, :]))*torch.exp(-(b[k+1, :])*tt) + a[k+1, :]/(b[k+1, :])  
        return c

    # ODE solution of CSF1     
    def f_EGF(self, c_EGF_0, c_EGF, dt, C_M1):
        c = torch.zeros_like(c_EGF)
        tt = dt*86400/2
        a = torch.exp(self.log_p_cyto[4])*self.S_EGF*(1 - C_M1)       
        b = a + torch.exp(self.log_p_cyto[6])*self.d_EGF   
        c[0, :] = (c_EGF_0 - a[0, :]/(b[0, :]))*torch.exp(-(b[0, :])*tt) + a[0, :]/(b[0, :])
        for k in range(c.shape[0]-1):
            c[k+1, :] = (c[k, :] - a[k+1, :]/(b[k+1, :]))*torch.exp(-(b[k+1, :])*tt) + a[k+1, :]/(b[k+1, :])     
        return c
    
    # ODE solution of IGF1
    def f_IGF1(self, c_IGF1_0, c_CSF1RI_0, c_IGF1, c_CSF1RI, dt, C_M1, t, dose_c):
        c = torch.zeros_like(c_IGF1)
        tt = dt*86400/2
        CSF1RI_pred = self.f_CSF1RI(c_CSF1RI_0, c_CSF1RI, dt, dose_c)
        a = torch.exp(self.log_p_cyto[1])*self.S_IGF1*(1 - C_M1)*self.A(CSF1RI_pred, t)
        b = a + torch.exp(self.log_p_cyto[3])*self.d_IGF1
        c[0, :] = (c_IGF1_0 - a[0, :]/(b[0, :]))*torch.exp(-(b[0, :])*tt) + a[0, :]/(b[0, :])
        for k in range(c.shape[0]-1):
            c[k+1, :] = (c[k, :] - a[k+1, :]/(b[k+1, :]))*torch.exp(-(b[k+1, :])*tt) + a[k+1, :]/(b[k+1, :])  
        return c
    
    # calculate the integral part in ODE of IGF1 
    def A(self, CSF1RI_pred, t):
        i_CSF1RI_pred = torch.cumsum(CSF1RI_pred, axis=0)
        i_CSF1RI = self.mapping(i_CSF1RI_pred, t)
        result = self.A_0 + self.S_A*i_CSF1RI*torch.exp(self.log_p_cyto[0])
        return result
    
    # ODE solution of CSF1R_I
    def f_CSF1RI(self, c_CSF1RI_0, c_CSF1RI, dt, dose_c):
        c = torch.zeros_like(c_CSF1RI)
        tt = dt
        a = torch.exp(self.log_p_cyto[10])*self.q_c
        b = torch.exp(self.log_p_cyto[11])*self.eta_c + torch.exp(self.log_p_cyto[12])*self.mu_c #d_M1 + d_M2 + d_M0 = 1
        c[0, :] = (c_CSF1RI_0 - a*dose_c[0, :]/(a+b))*torch.exp(-(a+b)*tt) + a*dose_c[0, :]/(a+b)   
        for k in range(c.shape[0]-1):
            c[k+1, :] = (c[k, :] - a*dose_c[k+1, :]/(a+b))*torch.exp(-(a+b)*tt) + a*dose_c[k+1, :]/(a+b)        
        return c

    # ODE solution of IGF1R_I
    def f_IGF1RI(self, c_IGF1RI_0, c_IGF1RI, dt, C_T, dose_I):
        c = torch.zeros_like(c_IGF1RI)
        tt = dt
        a = torch.exp(self.log_p_cyto[13])*self.q_i
        b = torch.exp(self.log_p_cyto[14])*self.eta_i + torch.exp(self.log_p_cyto[15])*self.mu_i * C_T
        c[0, :] = (c_IGF1RI_0 - a*dose_I[0, :]/(a+b[0, :]))*torch.exp(-(a+b[0, :])*tt) + a*dose_I[0, :]/(a+b[0, :])   
        for k in range(c.shape[0]-1):
            c[k+1, :] = (c[k, :] - a*dose_I[k+1, :]/(a+b[k+1, :]))*torch.exp(-(a+b[k+1, :])*tt) + a*dose_I[k+1, :]/(a+b[k+1, :])        
        return c

    # calculate Hill function H_1
    def H_1(self, c_EGF):
        result = c_EGF / (torch.exp(self.log_p_TC[4]) * self.K_1 + c_EGF)
        return result
    
    # calculate Hill function H_2
    def H_2(self, c_IGF1, c_IGF1RI):
        result = c_IGF1 / (torch.exp(self.log_p_TC[5]) * self.K_2 + 
                                c_IGF1 + torch.exp(self.log_p_TC[9]) *  c_IGF1RI)
        return result
    
    # calculate Hill function H_C
    def H_C(self, c_CSF1, c_CSF1RI):
        result = c_CSF1 / (torch.exp(self.log_p_M1[5]) * self.K_M12 + 
                                c_CSF1 + torch.exp(self.log_p_M1[6]) * self.Kd_M12 * c_CSF1RI)
        return result
    
    # calculate Hill function H_I
    def H_I(self, c_CSF1RI, t):
        a = self.A(c_CSF1RI, t)
        result = a/(torch.exp(self.log_p_M1[7])*self.K_I + a)
        return result
    
    # DNN for SDE of M1 macrophage density
    # Define input and output of the DNN
    def net_u_fM1(self, t, dose_c, dose_cum, dose_I):
        t = torch.clamp(t, min=self.dt, max = self.ub[1])
        t = 2.0*(t - self.dt)/(self.ub[1] - self.dt) - 1.0  
        dose_c = 2.0*(dose_c - 0) - 1.0
        dose_cum = 2.0*(dose_cum - 0)/(self.t_cyto.shape[0] - 0) - 1.0
        dose_I = 2.0*(dose_I - 0) - 1.0
        M1TD = torch.cat([t, dose_c, dose_cum, dose_I], dim=1) 
        c = self.dnn_fM1(M1TD)
        c = torch.clamp(c, 0.0, 1.0) # density (normalized)
        return c
    
    # use autograd to calculatie PINN
    def net_f_fM1(self, t, c_CSF1, c_CSF1RI, h_I, dose_c, dose_cum, dose_I):
        c = self.net_u_fM1(t, dose_c, dose_cum, dose_I)
        c_t = torch.autograd.grad(
            c, t, 
            grad_outputs=torch.ones_like(c),
            retain_graph=True,
            create_graph=True
        )[0]
        a = torch.exp(self.log_p_M1[0])*self.p_M21*(1 - c) - \
             torch.exp(self.log_p_M1[1])*(torch.exp(self.log_p_M1[3])*self.a_C*self.H_C(c_CSF1, c_CSF1RI) + \
                                         torch.exp(self.log_p_M1[4])*self.a_I*h_I) * c
        b = c * (-torch.exp(self.log_p_M1[1])*torch.exp(self.log_p_sigma[1])* \
             (torch.exp(self.log_p_M1[3])*self.a_C*self.H_C(c_CSF1, c_CSF1RI)+torch.exp(self.log_p_M1[4])* \
              h_I))
        f = c_t - a - b # PINN of dnn_fM1
        return f

    # the function to calculate TC proliferation rate
    def f_rTC(self, c_EGF, c_IGF1, c_IGF1RI):
        result = self.r_TC * (1 + torch.exp(self.log_p_TC[2]) * self.H_1(c_EGF) + torch.exp(self.log_p_TC[3]) * self.H_2(c_IGF1, c_IGF1RI))
        return result
    
    # DNN for SDE of tumor cell density
    # Define input and output of the DNN
    def net_u_fTC(self, t, dose_c, dose_cum, dose_I): 
        t = torch.clamp(t, min=self.dt, max = self.ub[1])
        t = 2.0*(t - self.dt)/(self.ub[1] - self.dt) - 1.0  
        dose_c = 2.0*(dose_c - 0) - 1.0
        dose_cum = 2.0*(dose_cum - 0)/(self.t_cyto.shape[0] - 0) - 1.0
        dose_I = 2.0*(dose_I - 0) - 1.0
        TCTD = torch.cat([t, dose_c, dose_cum, dose_I], dim=1)      
        c = self.dnn_fTC(TCTD)
        c = torch.clamp(c, 0.0, 1.0) # density (normalized)
        return c
    
    # use autograd to calculatie PINN
    def net_f_fTC(self, t, c_EGF, c_IGF1, dose_c, dose_cum, dose_I, c_IGF1RI):
        c = self.net_u_fTC(t, dose_c, dose_cum, dose_I)
        c_M1 = self.net_u_fM1(t - torch.exp(self.log_p_TC[8]), dose_c, dose_cum, dose_I)
        c_t = torch.autograd.grad(
            c, t, 
            grad_outputs=torch.ones_like(c),
            retain_graph=True,
            create_graph=True
        )[0]
        a = torch.exp(self.log_p_TC[0])*self.f_rTC(c_EGF, c_IGF1, c_IGF1RI) * c * (1 - c) - \
             torch.exp(self.log_p_TC[1])*self.d_TC*(1+torch.exp(self.log_p_TC[6])*self.dM1_TC* \
                                       c_M1)*c
        b = c* (-torch.exp(self.log_p_TC[1])*self.d_TC*torch.exp(self.log_p_TC[6])* \
             self.dM1_TC*c_M1 * torch.exp(self.log_p_sigma[0])) 
        f = c_t - a - b # PINN of dnn_fTC
        return f
    
    # DNN for Fokker-Planck equation of M1 macrophages 
    # Define input and output of the DNN   
    # (if necessary)     
    def net_u_M1(self, C, t, dose_c, dose_cum, dose_I):
        t = 2.0*(t- self.lb[1])/(self.ub[1] - self.lb[1]) - 1.0 
        C = 2.0*(C- self.lb[0])/(self.ub[0] - self.lb[0]) - 1.0 
        dose_c = 2.0*(dose_c - 0) - 1.0
        dose_cum = 2.0*(dose_cum - 0)/(self.t_cyto.shape[0] - 0) - 1.0
        dose_I = 2.0*(dose_I - 0) - 1.0
        HM1 = torch.cat([t, C, dose_c, dose_cum, dose_I], dim=1)
        p = self.dnn_M1(HM1) # probability
        return p

    # use autograd to calculatie PINN
    def net_f_M1(self, C, t, c_CSF1, c_CSF1RI, h_I, dose_c, dose_cum, dose_I):
        p = self.net_u_M1(C, t, dose_c, dose_cum, dose_I)
        p_t = torch.autograd.grad(
            p, t, 
            grad_outputs=torch.ones_like(p),
            retain_graph=True,
            create_graph=True
        )[0]
        a = torch.exp(self.log_p_M1[0])*self.p_M21*(self.c_max - C) - \
             torch.exp(self.log_p_M1[1])*(torch.exp(self.log_p_M1[3])*self.a_C*self.H_C(c_CSF1, c_CSF1RI) + \
                                         torch.exp(self.log_p_M1[4])*self.a_I*h_I) * C
        a_C = -torch.exp(self.log_p_M1[0])*self.p_M21 - \
                torch.exp(self.log_p_M1[1])*(torch.exp(self.log_p_M1[3])*self.a_C*self.H_C(c_CSF1, c_CSF1RI) + \
                                         torch.exp(self.log_p_M1[4])*self.a_I*h_I)
        b2 = C**2 * (-torch.exp(self.log_p_M1[1])*torch.exp(self.log_p_sigma[1])* \
             (torch.exp(self.log_p_M1[3])*self.a_C*self.H_C(c_CSF1, c_CSF1RI)+torch.exp(self.log_p_M1[4])* \
              self.a_I*h_I))**2 
        b2_C = 2 * C * (-torch.exp(self.log_p_M1[1])*torch.exp(self.log_p_sigma[1])* \
             (torch.exp(self.log_p_M1[3])*self.a_C*self.H_C(c_CSF1, c_CSF1RI)+torch.exp(self.log_p_M1[4])* \
              self.a_I*h_I))**2 
        b2_CC = 2  * (-torch.exp(self.log_p_M1[1])*torch.exp(self.log_p_sigma[1])* \
             (torch.exp(self.log_p_M1[3])*self.a_C*self.H_C(c_CSF1, c_CSF1RI)+torch.exp(self.log_p_M1[4])* \
              self.a_I*h_I))**2 
        p_C = torch.autograd.grad(
            p, 
            C, 
            grad_outputs=torch.ones_like(p),
            retain_graph=True,
            create_graph=True
        )[0]
        p_CC = torch.autograd.grad(
            p_C, C, 
            grad_outputs=torch.ones_like(p_C),
            retain_graph=True,
            create_graph=True
        )[0] 
        f = p_t + a_C * p + a * p_C - 0.5 * (b2_CC * p + 2 * b2_C * p_C + b2 * p_CC) # PINN of dnn_M1
        return f
    
    # DNN for Fokker-Planck equation of tumor cell
    # Define input and output of the DNN         
    def net_u_TC(self, C, t, dose_c, dose_cum, dose_I): 
        t = 2.0*(t- self.lb[1])/(self.ub[1] - self.lb[1]) - 1.0 
        C = 2.0*(C- self.lb[0])/(self.ub[0] - self.lb[0]) - 1.0
        dose_c = 2.0*(dose_c - 0) - 1.0
        dose_cum = 2.0*(dose_cum - 0)/(self.t_cyto.shape[0] - 0) - 1.0
        dose_I = 2.0*(dose_I - 0) - 1.0
        HTC = torch.cat([t, C, dose_c, dose_cum, dose_I], dim=1)
        p = self.dnn_TC(HTC)
        return p
    
    # use autograd to calculatie PINN
    def net_f_TC(self, C, t, c_EGF, c_IGF1, d_M1, dose_c, dose_cum, dose_I, c_IGF1RI):
        p = self.net_u_TC(C, t, dose_c, dose_cum, dose_I)
        p_t = torch.autograd.grad(
            p, t, 
            grad_outputs=torch.ones_like(p),
            retain_graph=True,
            create_graph=True
        )[0]
        a = torch.exp(self.log_p_TC[0])*self.f_rTC(c_EGF, c_IGF1, c_IGF1RI) * C * (self.c_max - C) - \
             torch.exp(self.log_p_TC[1])*self.d_TC*(1+torch.exp(self.log_p_TC[6])*self.dM1_TC* \
                                       d_M1)*C
        a_C = torch.exp(self.log_p_TC[0])*self.f_rTC(c_EGF, c_IGF1, c_IGF1RI) * (self.c_max - 2 * C) - \
                torch.exp(self.log_p_TC[1])*self.d_TC*(1+torch.exp(self.log_p_TC[6])*self.dM1_TC* \
                                                                              d_M1*C)
        b2 = C**2 * (-torch.exp(self.log_p_TC[1])*self.d_TC*torch.exp(self.log_p_TC[6])* \
             self.dM1_TC*d_M1 * torch.exp(self.log_p_sigma[0]))**2 
        b2_C = 2 * C * (-torch.exp(self.log_p_TC[1])*self.d_TC*torch.exp(self.log_p_TC[6])* \
             self.dM1_TC*d_M1 * torch.exp(self.log_p_sigma[0]))**2 
        b2_CC = 2 * (-torch.exp(self.log_p_TC[1])*self.d_TC*torch.exp(self.log_p_TC[6])* \
             self.dM1_TC*d_M1 * torch.exp(self.log_p_sigma[0]))**2 
        p_C = torch.autograd.grad(
            p, C, 
            grad_outputs=torch.ones_like(p),
            retain_graph=True,
            create_graph=True
        )[0]
        p_CC = torch.autograd.grad(
            p_C, C, 
            grad_outputs=torch.ones_like(p_C),
            retain_graph=True,
            create_graph=True
        )[0]     
        f = p_t + a_C * p + a * p_C - 0.5 * (b2_CC * p + 2 * b2_C * p_C + b2 * p_CC) # PINN of dnn_TC
        return f

    # save model
    def save(self):
        torch.save(self.dnn_M1.state_dict(), 'surrogate_model/dnn_M1_model.pth')
        print("dnn_M1 saved")
        torch.save(self.dnn_TC.state_dict(), 'surrogate_model/dnn_TC_model.pth')
        print("dnn_TC saved")
        torch.save(self.dnn_fM1.state_dict(), 'surrogate_model/dnn_fM1_model.pth')
        print("dnn_fM1 saved")
        torch.save(self.dnn_fTC.state_dict(), 'surrogate_model/dnn_fTC_model.pth')
        print("dnn_fTC saved")
        torch.save({'log_p_TC': self.log_p_TC, 'log_p_M1': self.log_p_M1, 'log_p_cyto': self.log_p_cyto,
                    'log_p_sigma': self.log_p_sigma}, 'surrogate_model/additional_params.pth')
        print("log_p saved")

    # LBFGS loss function for dnn_TC (and dnn_M1 if necessary)
    def loss_func(self):
        loss_LBFGS = torch.zeros(self.train_data_num) 
        self.optimizer_LBFGS.zero_grad()
        
        for i in range(self.train_data_num):
            #uM1_exact = self.uM1[i] # if necessary
            uTC_exact = self.uTC[i]

            cyto_i = self.cyto[i]
            d_M1 = cyto_i[:, 0:1]
            #c_CSF1 = cyto_i[:, 1:2]
            c_IGF1 = cyto_i[:, 2:3]
            c_EGF = cyto_i[:, 3:4]
            c_CSF1RI = cyto_i[:, 4:5]
            #d_TC = cyto_i[:, 5:6]
            c_IGF1RI = cyto_i[:, 6:7]

            #fCSF1 = self.mapping(c_CSF1, self.t_f)
            #fCSF1RI = self.mapping(c_CSF1RI, self.t_f)
            fIGF1RI = self.mapping(c_IGF1RI, self.t_f)
            fEGF = self.mapping(c_EGF, self.t_f)
            fIGF1 = self.mapping(c_IGF1, self.t_f)
            #h_I = self.H_I(c_CSF1RI, self.t_f)
            f_dM1 = self.mapping(d_M1, self.t_f)

            #dose_i = self.dose[i]
            #dose_i_CSF1RI = dose_i[:, 0:1]
            dose_u_CSF1RI = self.mapping(c_CSF1RI, self.t_u)
            dose_f_CSF1RI = self.mapping(c_CSF1RI, self.t_f)
            #dose_i_IGF1RI = dose_i[:, 1:2]
            dose_u_IGF1RI = self.mapping(c_IGF1RI, self.t_u)
            dose_f_IGF1RI = self.mapping(c_IGF1RI, self.t_f)

            dose_cum = torch.cumsum(c_CSF1RI, dim = 0)
            dose_u_cum = self.mapping(dose_cum, self.t_u)
            dose_f_cum = self.mapping(dose_cum, self.t_f)

            #uM1_pred = self.net_u_M1(self.C_u, self.t_u, dose_u_CSF1RI, dose_u_cum, dose_u_IGF1RI) # if necessary
            #fM1_pred = self.net_f_M1(self.C_f, self.t_f, fCSF1, fCSF1RI, h_I, dose_f_CSF1RI, dose_f_cum, dose_f_IGF1RI) # if necessary
            uTC_pred = self.net_u_TC(self.C_u, self.t_u, dose_u_CSF1RI, dose_u_cum, dose_u_IGF1RI)
            fTC_pred = self.net_f_TC(self.C_f, self.t_f, fEGF, fIGF1, f_dM1, dose_f_CSF1RI, dose_f_cum, dose_f_IGF1RI, fIGF1RI)

            #loss_uM1 = torch.mean((uM1_exact - uM1_pred) ** 2) # if necessary
            #loss_fM1 = torch.mean(fM1_pred ** 2) # if necessary
            loss_uTC = torch.mean((uTC_exact - uTC_pred) ** 2)
            loss_fTC = torch.mean(fTC_pred ** 2)

            loss_LBFGS[i] = loss_uTC + loss_fTC # + loss_uM1 + loss_fM1 # if necessary
            if self.iter % 1 == 0:
                print(
                    'LBFGS Iter %d, LBFGS_Loss[%d]: %.5e' 
                    % (self.iter, i, loss_LBFGS[i].item())
                )
                #print(
                #    'Loss_uTC: %.5e, Loss_fTC: %.5e' 
                #    % (loss_uTC.item(), loss_fTC.item())
                #)
                #print(
                #    'Loss_uM1: %.5e, Loss_fM1: %.5e, Loss_uTC: %.5e, Loss_fTC: %.5e' 
                #    % (loss_uM1.item(), loss_fM1.item(), loss_uTC.item(), loss_fTC.item())
                #)

        # select top 10 largest loss
        top_losses, _ = torch.topk(loss_LBFGS, 10)
        top_loss_sum = torch.sum(top_losses)

        print(
            'LBFGS Iter %d, top_loss_sum: %.5e' 
            % (self.iter, top_loss_sum.item())
        )
        print("-----------------------------------------------------")                     
        top_loss_sum.backward(retain_graph=True)
        self.iter += 1
        if self.iter % 200 == 0:
            self.save()
        return top_loss_sum

    # use Adam and LBFGS optimizer to train dnn_TC (and dnn_M1 if necessary)
    def train(self):
        #self.dnn_M1.train() # if necessary
        self.dnn_TC.train()

        print("before adam self.p_cyto: ", torch.exp(self.log_p_cyto))
        print("before adam self.p_TC: ", torch.exp(self.log_p_TC))
        print("before adam self.p_M1: ", torch.exp(self.log_p_M1))
        print("before adam self.p_sigma: ", torch.exp(self.log_p_sigma))

        # use Adam optimizer first
        for k in range(5000):
            loss_Adam = torch.zeros(self.train_data_num)
            for i in range(self.train_data_num):
                #uM1_exact = self.uM1[i]
                uTC_exact = self.uTC[i]

                cyto_i = self.cyto[i]
                d_M1 = cyto_i[:, 0:1]
                #c_CSF1 = cyto_i[:, 1:2]
                c_IGF1 = cyto_i[:, 2:3]
                c_EGF = cyto_i[:, 3:4]
                c_CSF1RI = cyto_i[:, 4:5]
                #d_TC = cyto_i[:, 5:6]
                c_IGF1RI = cyto_i[:, 6:7]

                #fCSF1 = self.mapping(c_CSF1, self.t_f)
                #fCSF1RI = self.mapping(c_CSF1RI, self.t_f)
                fIGF1RI = self.mapping(c_IGF1RI, self.t_f)
                fEGF = self.mapping(c_EGF, self.t_f)
                fIGF1 = self.mapping(c_IGF1, self.t_f)
                #h_I = self.H_I(c_CSF1RI, self.t_f)
                f_dM1 = self.mapping(d_M1, self.t_f)

                #dose_i = self.dose[i]
                #dose_i_CSF1RI = dose_i[:, 0:1]
                dose_u_CSF1RI = self.mapping(c_CSF1RI, self.t_u)
                dose_f_CSF1RI = self.mapping(c_CSF1RI, self.t_f)
                #dose_i_IGF1RI = dose_i[:, 1:2]
                dose_u_IGF1RI = self.mapping(c_IGF1RI, self.t_u)
                dose_f_IGF1RI = self.mapping(c_IGF1RI, self.t_f)

                dose_cum = torch.cumsum(c_CSF1RI, dim = 0)
                dose_u_cum = self.mapping(dose_cum, self.t_u)
                dose_f_cum = self.mapping(dose_cum, self.t_f)

                #uM1_pred = self.net_u_M1(self.C_u, self.t_u, dose_u_CSF1RI, dose_u_cum, dose_u_IGF1RI) # if necessary
                #fM1_pred = self.net_f_M1(self.C_f, self.t_f, fCSF1, fCSF1RI, h_I, dose_f_CSF1RI, dose_f_cum, dose_f_IGF1RI) # if necessary
                uTC_pred = self.net_u_TC(self.C_u, self.t_u, dose_u_CSF1RI, dose_u_cum, dose_u_IGF1RI)
                fTC_pred = self.net_f_TC(self.C_f, self.t_f, fEGF, fIGF1, f_dM1, dose_f_CSF1RI, dose_f_cum, dose_f_IGF1RI, fIGF1RI)

                #loss_uM1 = torch.mean((uM1_exact - uM1_pred) ** 2) # if necessary
                #loss_fM1 = torch.mean(fM1_pred ** 2) # if necessary
                loss_uTC = torch.mean((uTC_exact - uTC_pred) ** 2)
                loss_fTC = torch.mean(fTC_pred ** 2)

                loss_Adam[i] = loss_uTC + loss_fTC # + loss_uM1 + loss_fM1 # if necessary

                if k % 100 == 0:
                    print(
                        'Adam Epoches %d, Adam_Loss[%d]: %.5e' 
                        % (k, i, loss_Adam[i].item())
                    )
                    #print(
                    #    'Loss_uTC: %.5e, Loss_fTC: %.5e' 
                    #    % (loss_uTC.item(), loss_fTC.item())
                    #)
                    #print(
                    #    'Loss_uM1: %.5e, Loss_fM1: %.5e, Loss_uTC: %.5e, Loss_fTC: %.5e' 
                    #    % (loss_uM1.item(), loss_fM1.item(), loss_uTC.item(), loss_fTC.item())
                    #)

            # select top 10 largest loss
            top_losses, _ = torch.topk(loss_Adam, 10)
            top_loss_sum = torch.sum(top_losses)

            print('Adam Epoches %d, top_loss_sum: %.5e' % (k, top_loss_sum.item()))
            print("-----------------------------------------------------")
            if top_loss_sum.item() < 0.1:
                print("Adam break")
                break
            if k % 200 == 0:
                self.save()
            self.optimizer_Adam.zero_grad()
            top_loss_sum.backward(retain_graph=True)
            self.optimizer_Adam.step()

        self.save()    

        print("after Adam self.p_cyto: ", torch.exp(self.log_p_cyto))
        print("after Adam self.p_CC: ", torch.exp(self.log_p_TC))
        print("after Adam self.p_M1: ", torch.exp(self.log_p_M1))
        print("after Adam self.p_sigma: ", torch.exp(self.log_p_sigma))
       
        # use LBFGS optimizer last
        self.optimizer_LBFGS.step(self.loss_func)

    # LBFGS loss function for parameters 
    # in SDE of TC, M1, cytokines and drugs
    def loss_func_cyto(self):
    
        loss_LBFGS = torch.zeros(self.train_data_num) 
        self.optimizer_LBFGS_cyto.zero_grad()
        
        for i in range(self.train_data_num):

            cyto_i = self.cyto[i]
            d_M1 = cyto_i[:, 0:1]
            c_CSF1 = cyto_i[:, 1:2]
            c_IGF1 = cyto_i[:, 2:3]
            c_EGF = cyto_i[:, 3:4]
            c_CSF1RI = cyto_i[:, 4:5]
            d_TC = cyto_i[:, 5:6]
            c_IGF1RI = cyto_i[:, 6:7]

            cyto_0_i = self.cyto_0[i]
            dose_i = self.dose[i]
            dose_i_CSF1RI = dose_i[:, 0:1]
            dose_i_IGF1RI = dose_i[:, 1:2]
            dose_cum = torch.cumsum(dose_i_CSF1RI, dim = 0)

            d_M1_0 = cyto_0_i[0]
            c_CSF1_0 = cyto_0_i[1]
            c_IGF1_0 = cyto_0_i[2]
            c_EGF_0 = cyto_0_i[3]
            c_CSF1RI_0 = cyto_0_i[4]
            d_TC_0 = cyto_0_i[5]
            c_IGF1RI_0 = cyto_0_i[6]

            f_EGF_pred = self.f_EGF(c_EGF_0, c_EGF, self.dt, d_M1)
            f_IGF1_pred = self.f_IGF1(c_IGF1_0, c_CSF1RI_0, c_IGF1, c_CSF1RI, self.dt,
                                        d_M1, self.t_cyto, dose_i_CSF1RI)
            f_CSF1_pred = self.f_CSF1(c_CSF1_0, c_CSF1, self.dt, d_TC)
            f_CSF1RI_pred = self.f_CSF1RI(c_CSF1RI_0, c_CSF1RI, self.dt, dose_i_CSF1RI)
            f_IGF1RI_pred = self.f_IGF1RI(c_IGF1RI_0, c_IGF1RI, self.dt, d_TC, dose_i_IGF1RI)

            fCSF1 = self.mapping(c_CSF1, self.t_cyto)
            fCSF1RI = self.mapping(c_CSF1RI, self.t_cyto)
            fIGF1RI = self.mapping(c_IGF1RI, self.t_cyto)
            fEGF = self.mapping(c_EGF, self.t_cyto)
            fIGF1 = self.mapping(c_IGF1, self.t_cyto)
            h_I = self.H_I(c_CSF1RI, self.t_cyto)
            f_dM1 = self.mapping(d_M1, self.t_cyto)

            uM1_pred = self.net_u_fM1(self.t_cyto, dose_i_CSF1RI, dose_cum, dose_i_IGF1RI)
            fM1_pred = self.net_f_fM1(self.t_cyto, fCSF1, fCSF1RI, h_I, dose_i_CSF1RI, dose_cum, dose_i_IGF1RI)
            uTC_pred = self.net_u_fTC(self.t_cyto, dose_i_CSF1RI, dose_cum, dose_i_IGF1RI)
            fTC_pred = self.net_f_fTC(self.t_cyto, fEGF, fIGF1, dose_i_CSF1RI, dose_cum, dose_i_IGF1RI, fIGF1RI)

            loss_uM1 = torch.mean((d_M1 - uM1_pred) ** 2)
            loss_fM1 = torch.mean(fM1_pred ** 2)
            loss_uTC = torch.mean((d_TC - uTC_pred) ** 2)
            loss_fTC= torch.mean(fTC_pred ** 2)

            loss_fIGF1 = torch.mean((c_IGF1 - f_IGF1_pred) ** 2)
            loss_fEGF = torch.mean((c_EGF - f_EGF_pred) ** 2)
            loss_fCSF1 = torch.mean((c_CSF1 - f_CSF1_pred) ** 2)
            loss_fCSF1RI = torch.mean((c_CSF1RI - f_CSF1RI_pred) ** 2) 
            loss_fIGF1RI = torch.mean((c_IGF1RI - f_IGF1RI_pred) ** 2)

            # train loss_fIGF1 + loss_fCSF1 + loss_fEGF + loss_fCSF1RI + loss_fIGF1RI first
            # and then train loss_uM1 + loss_fM1 + loss_uTC + loss_fTC
            loss_LBFGS[i] = loss_fIGF1 + loss_fCSF1 + loss_fEGF + loss_fCSF1RI + loss_fIGF1RI 
                            # + loss_uM1 + loss_fM1 + loss_uTC + loss_fTC
               
            if self.iter_cyto % 100 == 0:
                print(
                    'LBFGS Iter_cyto %d, LBFGS_Loss[%d]: %.5e' 
                    % (self.iter_cyto, i, loss_LBFGS[i].item())
                )
                #print(
                #    'Loss_uTC: %.5e, Loss_fTC: %.5e' 
                #    % (loss_uTC.item(), loss_fTC.item())
                #)
                #print(
                #    'Loss_uM1: %.5e, Loss_fM1: %.5e, Loss_uTC: %.5e, Loss_fTC: %.5e' 
                #    % (loss_uM1.item(), loss_fM1.item(), loss_uTC.item(), loss_fTC.item())
                #)
                #print(
                #    'Loss_fCSF1: %.5e, Loss_fEGF: %.5e, Loss_fIGF1: %.5e, Loss_fCSF1RI: %.5e, Loss_fIGF1RI: %.5e' 
                #    % (loss_fCSF1.item(), loss_fEGF.item(), loss_fIGF1.item(), loss_fCSF1RI.item(), loss_fIGF1RI.item())
                #)
 
        loss_LBFGS_sum = torch.sum(loss_LBFGS) 
                                
        if(loss_LBFGS_sum > 1e6 or loss_LBFGS_sum == np.nan):
            exit()
        
        loss_LBFGS_sum.backward(retain_graph=True)
        if self.iter_cyto % 200 == 0:
            self.save()
        self.iter_cyto += 1

        return loss_LBFGS_sum

    # use LBFGS optimizer to train parameters 
    # in SDE of TC, M1, cytokines and drugs   
    def train_cyto(self):
        self.dnn_fM1.train()
        self.dnn_fTC.train()
        self.optimizer_LBFGS_cyto.step(self.loss_func_cyto)

    # predict cytokines and drugs
    def predict_cyto(self, train_data_num):
        cyto_i = self.cyto[train_data_num]
        cyto_0_i = self.cyto_0[train_data_num]
        dose_i = self.dose[train_data_num]
        dose_i_CSF1RI = dose_i[:, 0:1]
        dose_i_IGF1RI = dose_i[:, 1:2]

        d_M1 = cyto_i[:, 0:1]
        c_CSF1 = cyto_i[:, 1:2]
        c_IGF1 = cyto_i[:, 2:3]
        c_EGF = cyto_i[:, 3:4]
        c_CSF1RI = cyto_i[:, 4:5]
        d_TC = cyto_i[:, 5:6]
        c_IGF1RI = cyto_i[:, 6:7]

        d_M1_0 = cyto_0_i[0]
        c_CSF1_0 = cyto_0_i[1]
        c_IGF1_0 = cyto_0_i[2]
        c_EGF_0 = cyto_0_i[3]
        c_CSF1RI_0 = cyto_0_i[4]
        d_TC_0 = cyto_0_i[5]
        c_IGF1RI_0 = cyto_0_i[6]

        f_EGF_pred = self.f_EGF(c_EGF_0, c_EGF, self.dt, d_M1)
        f_IGF1_pred = self.f_IGF1(c_IGF1_0, c_CSF1RI_0, c_IGF1, c_CSF1RI, self.dt,
                                    d_M1, self.t_cyto, dose_i_CSF1RI)
        f_CSF1_pred = self.f_CSF1(c_CSF1_0, c_CSF1, self.dt, d_TC)
        f_CSF1RI_pred = self.f_CSF1RI(c_CSF1RI_0, c_CSF1RI, self.dt, dose_i_CSF1RI)
        f_IGF1RI_pred = self.f_IGF1RI(c_IGF1RI_0, c_IGF1RI, self.dt, d_TC, dose_i_IGF1RI)

        f_CSF1_pred = f_CSF1_pred.detach().cpu().numpy()
        f_EGF_pred = f_EGF_pred.detach().cpu().numpy()
        f_IGF1_pred = f_IGF1_pred.detach().cpu().numpy()
        f_CSF1RI_pred = f_CSF1RI_pred.detach().cpu().numpy()
        f_IGF1RI_pred = f_IGF1RI_pred.detach().cpu().numpy()

        return f_CSF1_pred, f_EGF_pred, f_IGF1_pred, f_CSF1RI_pred, f_IGF1RI_pred

    # predict TC and M1 density
    def predict_cell(self, train_data_i):
        self.dnn_fM1.eval()
        self.dnn_fTC.eval()

        dose_i = self.dose[train_data_i]
        dose_i_CSF1RI = dose_i[:, 0:1]
        dose_i_IGF1RI = dose_i[:, 1:2]
        dose_cum = torch.cumsum(dose_i_CSF1RI, dim = 0)
          
        fTC_pred = self.net_u_fTC(self.t_cyto, dose_i_CSF1RI, dose_cum, dose_i_IGF1RI)
        fM1_pred = self.net_u_fM1(self.t_cyto, dose_i_CSF1RI, dose_cum, dose_i_IGF1RI)

        fTC_pred = fTC_pred.detach().cpu().numpy()
        fM1_pred = fM1_pred.detach().cpu().numpy()

        return fTC_pred, fM1_pred

    # predict Fokker-Planck equation of TC (and M1 if necessary)
    def predict_FP(self, train_data_i, C_v, t_v):
        #self.dnn_M1.eval() # if necessary
        self.dnn_TC.eval()  

        C_v = torch.tensor(C_v, dtype=torch.float32, requires_grad=True).float().to(device)
        t_v = torch.tensor(t_v, dtype=torch.float32, requires_grad=True).float().to(device)

        cyto_i = self.cyto[train_data_i]
        c_CSF1RI = cyto_i[:, 4:5]
        c_IGF1RI = cyto_i[:, 6:7]
        dose_v_CSF1RI = self.mapping(c_CSF1RI, t_v)
        dose_v_IGF1RI = self.mapping(c_IGF1RI, t_v)
        dose_cum = torch.cumsum(c_CSF1RI, dim = 0)
        dose_v_cum  = self.mapping(dose_cum, t_v)

        uCC_pred = self.net_u_TC(C_v, t_v, dose_v_CSF1RI, dose_v_cum, dose_v_IGF1RI)
        uM1_pred = self.net_u_M1(C_v, t_v, dose_v_CSF1RI, dose_v_cum, dose_v_IGF1RI)

        uCC_pred = uCC_pred.detach().cpu().numpy()
        uM1_pred = uM1_pred.detach().cpu().numpy()

        return uCC_pred, uM1_pred

# convert survival rate sequence to time & event
def convert_to_time_event(survival_rate):
    time = []
    event = []
    survival_num = 100
    for i in range(len(survival_rate) - 1):
        cut_num = survival_num - int(np.floor((survival_rate[i]) * 100))
        if cut_num < 0:
            cut_num = 0
        survival_num = survival_num - cut_num
        for j in range(cut_num):
            time.append(i + 1)
            event.append(1)
    for i in range(survival_num):
            time.append(len(survival_rate))
            event.append(0)
    return time, event

# add the missing time points
def supplement_survival_function(t_max, survival_df):
    column_label = survival_df.columns[0]
    seq = np.arange(0, 201, 5) 

    if t_max not in survival_df.index:
        survival_df = pd.concat([survival_df, pd.DataFrame({column_label: [0]}, index=[t_max])])

    missing_points = [t for t in seq if t not in survival_df.index]
    times = survival_df.index.values
    survival_values = survival_df[column_label].values
    interpolated_values = []

    for mp in missing_points:
        idx = np.searchsorted(times, mp, side='right') - 1
        interpolated_values.append(survival_values[idx])

    new_rows = pd.DataFrame({column_label: interpolated_values}, index=missing_points)
    survival_df = pd.concat([survival_df, new_rows]).sort_index()
    return survival_df

if __name__ == '__main__':
    start_time = time.time() 

    N_f = 200 # number of ramdomly selected points
    FP_layers = [5, 40, 120, 250, 500, 1000, 1000, 600, 300, 150, 1]
    CC_layers = [4, 40, 40, 40, 40, 40, 40, 40, 40, 1]

    t_max = 200
    dt = 5 
    t = np.arange(dt, t_max + dt, dt) # for training
    dt_cyto = 1/48
    t_cyto = np.arange(dt_cyto, t_max + dt_cyto, dt_cyto) # for training
    
      
    c_min, c_max = 0, 1
    dc = 0.05 
    c = np.arange(c_min + dc, c_max + dc, dc) # for training
    

    denominator = 10 # 1-10, help make predicting and plotting more precise
    dtt = dt/denominator
    dcc = dc/denominator
    tt = np.arange(dt, t_max + dtt, dtt) # for predicting and plotting
    cc = np.arange(c_min + dc, c_max + dcc, dcc) # for predicting and plotting

    t = np.array(t)[:, None]
    tt = np.array(tt)[:, None]
    t_cyto = np.array(t_cyto)[:, None]
    c = np.array(c)[:, None]
    cc = np.array(cc)[:, None]

    C, T = np.meshgrid(c,t)
    C_star = np.hstack((C.flatten()[:,None], T.flatten()[:,None]))
    X_u_train = C_star # prior data points 

    lb = C_star.min(0)
    ub = C_star.max(0)  
    X_f_train = lb + (ub-lb)*lhs(2, N_f)
    X_f_train = np.vstack((X_f_train, X_u_train))

    Exact_TC_list = []
    Exact_M1_list = []
    Exact_cyto_list = []
    Exact_dose_list = []
    uM1_train_list = []
    uTC_train_list = []
    cyto_0_list = []
    
    train_data_num = 20
    for i in range(train_data_num):
        Exact_TC_filepath = f'PINN_load_data/data_{i}/CC_pdf_20_40.csv'
        Exact_TC = pd.read_csv(Exact_TC_filepath, header=None) #20*40
        Exact_TC = Exact_TC.T * 1/dc #40*20
        Exact_TC_list.append(Exact_TC.values)

        Exact_M1_filepath = f'PINN_load_data/data_{i}/M1_pdf_20_40.csv'
        Exact_M1 = pd.read_csv(Exact_M1_filepath, header=None) #20*40
        Exact_M1 = Exact_M1.T * 1/dc #40*20
        Exact_M1_list.append(Exact_M1.values)

        #[d_M1; CSF1; IGF1; EGF; in vivo CSF1R_I; d_TC; in vivo IGF1R_I]
        Exact_cyto_filepath = f'PINN_load_data/data_{i}/r_cyto_9600.csv'
        Exact_cyto = pd.read_csv(Exact_cyto_filepath, header=None) #7*9600
        Exact_cyto = Exact_cyto.T #9600*7
        Exact_cyto_list.append(Exact_cyto.values)

        #[dose CSF1R_I; dose IGF1R_I]
        Exact_dose_filepath = f'PINN_load_data/data_{i}/r_dose_9600.csv'
        Exact_dose = pd.read_csv(Exact_dose_filepath, header=None) #2*9600
        Exact_dose = Exact_dose.T #9600*2
        Exact_dose_list.append(Exact_dose.values)

        # inital values [d_M1, CSF1, IGF1, EGF, in vivo CSF1R_I, d_TC, in vivo IGF1R_I]
        cyto_0_list.append([35/450, 0.3128, 0.0, 0.5, 0.0, 0.58, 0.0])
        
        u_starM1 = Exact_M1.values.flatten()[:,None]
        uM1_train_list.append(u_starM1)              
        u_starCC = Exact_TC.values.flatten()[:,None]
        uTC_train_list.append(u_starCC)
    
    # PINN model for training and predicting
    model = PhysicsInformedNN(X_u_train, uM1_train_list, uTC_train_list, X_f_train, Exact_cyto_list, Exact_dose_list, 
                              cyto_0_list, train_data_num, FP_layers, CC_layers, t_cyto, dt_cyto, dc, lb, ub)

    # train parameters 
    model.train_cyto()
    model.save()
    
    # plot
    figure_save_path = "PINN_figures"
    for i in range(train_data_num):
        f_CSF1, f_EGF, f_IGF1, f_CSF1RI, f_IGF1RI= model.predict_cyto(i)
        f_CC, f_M1 = model.predict_cell(i)

        # plot in vivo cytokines and drugs in Fig 7C
        cyto_i = Exact_cyto_list[i].T

        fig = plt.figure(figsize=(14, 12))
        ax = fig.add_subplot(111)
        gs1 = gridspec.GridSpec(2, 3)
        gs1.update(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.3, hspace=0.3)

        ax = plt.subplot(gs1[0, 0])
        ax.plot(t_cyto, cyto_i[3,:], 'b-', linewidth = 2, label = 'Exact_EGF')       
        ax.plot(t_cyto, f_EGF, 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$t$')
        ax.set_ylabel('$c$')    
        ax.set_title('$EGF$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,t_max])
        ax.set_ylim([0.0,1])
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        ax.legend()
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)
            
        ax = plt.subplot(gs1[0, 1])
        ax.plot(t_cyto, cyto_i[2,:], 'b-', linewidth = 2, label = 'Exact_IGF1')       
        ax.plot(t_cyto, f_IGF1, 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$t$')
        ax.set_ylabel('$c$')    
        ax.set_title('$IGF1$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,t_max])
        ax.set_ylim([0.0,1])
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15) 

        ax = plt.subplot(gs1[1, 0])
        ax.plot(t_cyto, cyto_i[1,:], 'b-', linewidth = 2, label = 'Exact_CSF1')       
        ax.plot(t_cyto, f_CSF1, 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$t$')
        ax.set_ylabel('$c$')    
        ax.set_title('$CSF1$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,t_max])
        ax.set_ylim([0.0,1])
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        ax = plt.subplot(gs1[1, 1])
        ax.plot(t_cyto, cyto_i[4,:], 'b-', linewidth = 2, label = 'Exact_CSF1RI')     
        ax.plot(t_cyto, f_CSF1RI, 'r--', linewidth = 2, label = 'Prediction')  
        ax.set_xlabel('$t$')
        ax.set_ylabel('$c$')    
        ax.set_title('$CSF1RI$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,t_max])
        ax.set_ylim([0.0,1])
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        ax = plt.subplot(gs1[0, 2])
        ax.plot(t_cyto, cyto_i[6,:], 'b-', linewidth = 2, label = 'Exact_IGF1RI')       
        ax.plot(t_cyto, f_IGF1RI, 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$t$')
        ax.set_ylabel('$c$')    
        ax.set_title('$IGF1RI$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,t_max])
        ax.set_ylim([0.0,1])
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        ax = plt.subplot(gs1[1, 2])
        ax.plot(t_cyto, cyto_i[5,:], 'b-', linewidth = 2, label = 'Exact_TC')       
        ax.plot(t_cyto, f_CC, 'r--', linewidth = 2, label = 'Predictive_TC')
        ax.plot(t_cyto, cyto_i[0,:], 'b-', linewidth = 2, label = 'Exact_M1')       
        ax.plot(t_cyto, f_M1, 'g--', linewidth = 2, label = 'Predictive_M1')
        ax.set_xlabel('$t$')
        ax.set_ylabel('$c$')    
        ax.set_title('$TC & M1$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,t_max])
        ax.set_ylim([0.0,1])
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)


        save_path = f'{figure_save_path}/rCyto_{i}.pdf'
        plt.savefig(save_path, format='pdf')
        print(f'cyto_{i} done')

        # plot drug dose
        dose_i = Exact_dose_list[i].T

        fig = plt.figure(figsize=(14, 12))
        ax = fig.add_subplot(111)
        gs1 = gridspec.GridSpec(2, 1)

        ax = plt.subplot(gs1[0, 0])
        gs1.update(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.3, hspace=0.3)
        ax.plot(t_cyto, dose_i[0,:], 'b-', linewidth = 4, label = 'CSF1RI_dose')       
        ax.set_xlabel('$t$')
        ax.set_ylabel('$c$')    
        ax.set_title('$CSF1RI dose$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,t_max])
        ax.set_ylim([0.0,1])
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        ax = plt.subplot(gs1[1, 0])
        gs1.update(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.3, hspace=0.3)
        ax.plot(t_cyto, dose_i[1,:], 'g-', linewidth = 4, label = 'IGF1RI_dose')       
        ax.set_xlabel('$t$')
        ax.set_ylabel('$c$')    
        ax.set_title('$IGF1RI dose$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,t_max])
        ax.set_ylim([0.0,1])
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                        ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        save_path = f'{figure_save_path}/rdose_{i}.pdf'
        plt.savefig(save_path, format='pdf')
        print(f'dose_{i} done')            

    # points for prediction 
    Cv, Tv = np.meshgrid(cc, tt)
    Cv_star = np.hstack((Cv.flatten()[:,None], Tv.flatten()[:,None]))

    # train dnn for TC (and M1 if necessary)
    model.train()       
    model.save()

    # lists for statistics
    p_value_list = [] # log-rank test p-value prediction vs excat
    p_basal_list = [] # log-rank test p-value prediction vs basal
    mst_list = [] # median survival time
    tqst_list = [] # three quarter survival time
    nqst_list = [] # 90 percent survival time 
    KLD_list = [] # KL divergence
    KS_list = [] # KS statistic
    KSp_list = [] # KS statistic p-value
    rmse_list = [] # root mean square error


    for i in range(train_data_num):  
        # u_predCC, u_predM1 = model.predict_FP(i, Cv_star[:, 0:1], Cv_star[:, 1:2]) # more denser prediction points
        u_predCC, u_predM1 = model.predict_FP(i, X_u_train[:, 0:1], X_u_train[:, 1:2]) # use training points as prediction points

        ''' # plot dnn_M1 results (if necessary)

        # plot Exact_M1
        Exact_M1 = Exact_M1_list[i]
        fig = plt.figure(figsize=(9, 5))
        ax = fig.add_subplot(111)
        print("Exact_M1.shape: ", Exact_M1.shape)
        h = ax.imshow(Exact_M1.T, interpolation='nearest', cmap='rainbow', 
                    extent=[t.min(), t.max(), c.min(), c.max()], 
                    origin='lower', aspect='auto')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.10)
        cbar = fig.colorbar(h, cax=cax)
        cbar.ax.tick_params(labelsize=15)            
        ax.plot(
            X_u_train[:,1], 
            X_u_train[:,0], 
            'kx', label = 'Data (%d points)' % (X_u_train.shape[0]), 
            markersize = 1,  # marker size doubled
            clip_on = False,
            alpha=1.0
        )    
        ax.set_xlabel('$t$', size=20)
        ax.set_ylabel('$C$', size=20)
        ax.legend(
            loc='upper center', 
            bbox_to_anchor=(0.9, -0.05), 
            ncol=5, 
            frameon=False, 
            prop={'size': 15}
        )
        ax.set_title('$p(t,C)$', fontsize = 20) # font size doubled
        ax.tick_params(labelsize=15)
        save_path = f'{figure_save_path}/Exact_M1_{i}.pdf'
        plt.savefig(save_path, format='pdf')
        
        # plot Pred_M1  
        U_pred = scipy.interpolate.griddata(Cv_star, u_predM1.flatten(), (C, T), method='cubic')
        U_pred[U_pred < 0] = 0
        print("U_pred.shape: ", U_pred.shape)
        fig = plt.figure(figsize=(9, 5))
        ax = fig.add_subplot(111)
        
        h = ax.imshow(U_pred.T, interpolation='nearest', cmap='rainbow', 
                    extent=[t.min(), t.max(), c.min(), c.max()], 
                    origin='lower', aspect='auto')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.10)
        cbar = fig.colorbar(h, cax=cax)
        cbar.ax.tick_params(labelsize=15)  
        ax.set_xlabel('$t$', size=20)
        ax.set_ylabel('$C$', size=20)

        ax.set_title('$p(t,C)$', fontsize = 20) # font size doubled
        ax.tick_params(labelsize=15)
        save_path = f'{figure_save_path}/Pred_M1_{i}.pdf'
        plt.savefig(save_path, format='pdf')

    
        # plot Pred_M1_silce
        fig = plt.figure(figsize=(14, 12))
        ax = fig.add_subplot(111)
        
        gs1 = gridspec.GridSpec(2, 3)
        gs1.update(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.3, hspace=0.3)
        
        ax = plt.subplot(gs1[0, 0])
        ax.plot(c, Exact_M1[30,:], 'b-', linewidth = 2, label = 'Exact_M1')       
        ax.plot(c, U_pred[30,:], 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$p(t,C)$')    
        ax.set_title('$t = 30$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([0.0,c_max])
        ax.set_ylim([0.0,c_max/dc])  
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)
        
        ax = plt.subplot(gs1[0, 1])
        ax.plot(c, Exact_M1[60,:], 'b-', linewidth = 2, label = 'Exact_M1')       
        ax.plot(c, U_pred[60,:], 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$u(t,C)$')
        ax.axis('square')
        ax.set_xlim([0.0,c_max])
        ax.set_ylim([0.0,c_max/dc])  
        ax.set_title('$t = 50$', fontsize = 15)
        ax.legend(
            loc='upper center', 
            bbox_to_anchor=(0.5, -0.15), 
            ncol=5, 
            frameon=False, 
            prop={'size': 15}
        )
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)
        
        ax = plt.subplot(gs1[0, 2])
        ax.plot(c, Exact_M1[90,:], 'b-', linewidth = 2, label = 'Exact_M1')       
        ax.plot(c, U_pred[90,:], 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$u(t,C)$')
        ax.axis('square')
        ax.set_xlim([0.0,c_max])
        ax.set_ylim([0.0,c_max/dc])   
        ax.set_title('$t = 75$', fontsize = 15)
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        ax = plt.subplot(gs1[1, 0])
        ax.plot(c, Exact_M1[120,:], 'b-', linewidth = 2, label = 'Exact_M1')       
        ax.plot(c, U_pred[120,:], 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$u(t,C)$')
        ax.axis('square')
        ax.set_xlim([0.0,c_max])
        ax.set_ylim([0.0,c_max/dc])   
        ax.set_title('$t = 120$', fontsize = 15)
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        ax = plt.subplot(gs1[1, 1])
        ax.plot(c, Exact_M1[150,:], 'b-', linewidth = 2, label = 'Exact_M1')       
        ax.plot(c, U_pred[150,:], 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$u(t,C)$')
        ax.axis('square')
        ax.set_xlim([0.0,c_max])
        ax.set_ylim([0.0,c_max/dc])   
        ax.set_title('$t = 150$', fontsize = 15)
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        ax = plt.subplot(gs1[1, 2])
        ax.plot(c, Exact_M1[180,:], 'b-', linewidth = 2, label = 'Exact_M1')       
        ax.plot(c, U_pred[180,:], 'r--', linewidth = 2, label = 'Prediction')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$u(t,C)$')
        ax.axis('square')
        ax.set_xlim([0.0,c_max])
        ax.set_ylim([0.0,c_max/dc])   
        ax.set_title('$t = 180$', fontsize = 15)
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        save_path = f'{figure_save_path}/Pred_M1_slice_{i}.pdf'
        plt.savefig(save_path, format='pdf')        
        '''
        
        # plot Exact_TC in Fig 7A
        Exact_TC = Exact_TC_list[i]

        fig = plt.figure(figsize=(9, 5))
        ax = fig.add_subplot(111)
        h = ax.imshow(Exact_TC.T, interpolation='nearest', cmap='rainbow', 
                    extent=[t.min(), t.max(), c.min(), c.max()], 
                    origin='lower', aspect='auto')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.10)
        cbar = fig.colorbar(h, cax=cax)
        cbar.ax.tick_params(labelsize=15)         
        ax.plot(
            X_u_train[:,1], 
            X_u_train[:,0], 
            'kx', label = 'Data (%d points)' % (X_u_train.shape[0]), 
            markersize = 1,  # marker size doubled
            clip_on = False,
            alpha=1.0
        )      
        ax.set_xlabel('$t$', size=20)
        ax.set_ylabel('$C$', size=20)
        ax.legend(
            loc='upper center', 
            bbox_to_anchor=(0.9, -0.05), 
            ncol=5, 
            frameon=False, 
            prop={'size': 15}
        )
        ax.set_title('$p(t,C)$', fontsize = 20) # font size doubled
        ax.tick_params(labelsize=15)
        save_path = f'{figure_save_path}/Exact_TC_{i}.pdf'
        plt.savefig(save_path, format='pdf')
        
        # plot Pred_CC in Fig 7A
        u_predCC_v, u_predM1_v = model.predict_FP(i, Cv_star[:, 0:1], Cv_star[:, 1:2])
        # more denser points using cubic interpolation 
        U_pred_v = scipy.interpolate.griddata(Cv_star, u_predCC_v.flatten(), (Cv, Tv), method='cubic')
        U_pred_v[U_pred_v < 0] = 0
        print("U_pred.shape: ", U_pred_v.shape)
        fig = plt.figure(figsize=(9, 5))
        ax = fig.add_subplot(111)    
        h = ax.imshow(U_pred_v.T, interpolation='nearest', cmap='rainbow', 
                    extent=[t.min(), t.max(), c.min(), c.max()], 
                    origin='lower', aspect='auto')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.10)
        cbar = fig.colorbar(h, cax=cax)
        cbar.ax.tick_params(labelsize=15) 
        ax.set_xlabel('$t$', size=20)
        ax.set_ylabel('$C$', size=20)
        ax.set_title('$p(t,C)$', fontsize = 20) # font size doubled
        ax.tick_params(labelsize=15)
        ax.plot(
            X_u_train[:,1], 
            X_u_train[:,0], 
            'kx', label = 'exact data (%d points)' % (X_u_train.shape[0]), 
            markersize = 2,  # marker size doubled
            clip_on = False,
            alpha=1.0
        )    

        save_path = f'{figure_save_path}/Pred_CC_{i}.pdf'
        plt.savefig(save_path, format='pdf')
        
        # plot Pred_CC_slice and statistics in Fig 7B
        U_pred = scipy.interpolate.griddata(C_star, u_predCC.flatten(), (C, T), method='cubic') 
        U_pred[U_pred < 0] = 0      
        fig = plt.figure(figsize=(14, 12))
        ax = fig.add_subplot(111)
        
        gs1 = gridspec.GridSpec(2, 2)
        gs1.update(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.3, hspace=0.3)
        
        ax = plt.subplot(gs1[0, 0])
        ax.plot(c, Exact_TC[9,:], 'b-', linewidth = 2, label = 'Exact_TC')       
        ax.plot(c, U_pred[9,:], 'r-', linewidth = 2, label = 'prediction')
        # MSE
        mse = np.mean((U_pred[9,:] - Exact_TC[9,:]) ** 2)
        # RMSE
        rmse = np.sqrt(mse)
        #relative_mse = mse / np.mean(Exact_TC[9,:])
        # KL divergence
        kl_divergence = np.sum(scipy.special.rel_entr(U_pred[9,:]+1e-7, Exact_TC[9,:]+1e-7))
        # R^2
        ss_res = np.sum((Exact_TC[9,:] - U_pred[9,:]) ** 2)
        ss_tot = np.sum((Exact_TC[9,:] - np.mean(Exact_TC[9,:])) ** 2)
        r_square = 1 - (ss_res / ss_tot)

        ax.text(0.05, 0.95, f'KL divergence: {kl_divergence:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.05, 0.9, f'R square: {r_square:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        rmse = np.sqrt(mse)
        ax.text(0.05, 0.85, f'RMSE: {rmse:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$p(t,C)$')    
        ax.set_title('$t = 50$', fontsize = 15)
        ax.axis('square')
        ax.set_xlim([dc,c_max])
        ax.set_ylim([0.0,c_max/dc])  
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)
        
        ax = plt.subplot(gs1[0, 1])
        ax.plot(c, Exact_TC[19,:], 'b-', linewidth = 2, label = 'Exact_TC')       
        ax.plot(c, U_pred[19,:], 'r-', linewidth = 2, label = 'prediction')
        # MSE
        mse = np.mean((U_pred[19,:] - Exact_TC[19,:]) ** 2)
        # RMSE
        rmse = np.sqrt(mse)
        #relative_mse = mse / np.mean(Exact_TC[19,:])
        # KL divergence
        kl_divergence = np.sum(scipy.special.rel_entr(U_pred[19,:]+1e-7, Exact_TC[19,:]+1e-7))
        # R^2
        ss_res = np.sum((Exact_TC[19,:] - U_pred[19,:]) ** 2)
        ss_tot = np.sum((Exact_TC[19,:] - np.mean(Exact_TC[19,:])) ** 2)
        r_square = 1 - (ss_res / ss_tot)

        ax.text(0.05, 0.95, f'KL divergence: {kl_divergence:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.05, 0.9, f'R square: {r_square:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.05, 0.85, f'RMSE: {rmse:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$u(t,C)$')
        ax.axis('square')
        ax.set_xlim([dc,c_max])
        ax.set_ylim([0.0,c_max/dc])  
        ax.set_title('$t = 100$', fontsize = 15)
        ax.legend(
            loc='upper center', 
            bbox_to_anchor=(0.5, -0.15), 
            ncol=5, 
            frameon=False, 
            prop={'size': 15}
        )
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)
        
        ax = plt.subplot(gs1[1, 0])
        ax.plot(c, Exact_TC[29,:], 'b-', linewidth = 2, label = 'Exact_TC')       
        ax.plot(c, U_pred[29,:], 'r-', linewidth = 2, label = 'prediction')
        # MSE
        mse = np.mean((U_pred[29,:] - Exact_TC[29,:]) ** 2)
        # RMSE
        rmse = np.sqrt(mse)
        #relative_mse = mse / np.mean(Exact_TC[29,:])
        # KL divergence
        kl_divergence = np.sum(scipy.special.rel_entr(U_pred[29,:]+1e-7, Exact_TC[29,:]+1e-7))
        # R^2
        ss_res = np.sum((Exact_TC[29,:] - U_pred[29,:]) ** 2)
        ss_tot = np.sum((Exact_TC[29,:] - np.mean(Exact_TC[29,:])) ** 2)
        r_square = 1 - (ss_res / ss_tot)

        ax.text(0.05, 0.95, f'KL divergence: {kl_divergence:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.05, 0.9, f'R square: {r_square:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.05, 0.85, f'RMSE: {rmse:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$u(t,C)$')
        ax.axis('square')
        ax.set_xlim([dc,c_max])
        ax.set_ylim([0.0,c_max/dc])  
        ax.set_title('$t = 150$', fontsize = 15)
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)

        ax = plt.subplot(gs1[1, 1])
        ax.plot(c, Exact_TC[39,:], 'b-', linewidth = 2, label = 'Exact_TC')       
        ax.plot(c, U_pred[39,:], 'r-', linewidth = 2, label = 'prediction')
        # MSE
        mse = np.mean((U_pred[39,:] - Exact_TC[39,:]) ** 2)
        # RMSE
        rmse = np.sqrt(mse)
        #relative_mse = mse / np.mean(Exact_TC[39,:])
        # KL divergence
        kl_divergence = np.sum(scipy.special.rel_entr(U_pred[39,:]+1e-7, Exact_TC[39,:]+1e-7))
        # R^2
        ss_res = np.sum((Exact_TC[39,:] - U_pred[39,:]) ** 2)
        ss_tot = np.sum((Exact_TC[39,:] - np.mean(Exact_TC[39,:])) ** 2)
        r_square = 1 - (ss_res / ss_tot)

        ax.text(0.05, 0.95, f'KL divergence: {kl_divergence:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.05, 0.9, f'R square: {r_square:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.text(0.05, 0.85, f'RMSE: {rmse:.4g}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
        ax.set_xlabel('$C$')
        ax.set_ylabel('$u(t,C)$')
        ax.axis('square')
        ax.set_xlim([dc,c_max])
        ax.set_ylim([0.0,c_max/dc])  
        ax.set_title('$t = 200$', fontsize = 15)
        aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        ax.set_aspect(aspect_ratio)    
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(15)


        save_path = f'{figure_save_path}/Pred_CC_slice_{i}.pdf'
        plt.savefig(save_path, format='pdf')
        print(f'predCC_{i} done')      
        
        # calculate statistics and plot KM curves
        U_pred = scipy.interpolate.griddata(C_star, u_predCC.flatten(), (C, T), method='cubic')
        U_pred[U_pred < 0] = 0  
        Exact_TC = Exact_TC_list[i]
        Basal_CC = Exact_TC_list[0]
        
        Exact_survival_rate = 1 - np.sum(Exact_TC.T[-3:, :], axis = 0)*dc
        Pred_survival_rate = 1 - np.sum(U_pred.T[-3:, :], axis = 0)/np.sum(U_pred.T, axis = 0)
        Basal_survival_rate = 1 - np.sum(Basal_CC.T[-3:, :], axis = 0)*dc

        Exact_survival_rate[0] = 1
        Basal_survival_rate[0] = 1
        Pred_survival_rate[0] = 1

        for j in range(1, len(Exact_survival_rate)):
            if Exact_survival_rate[j] > Exact_survival_rate[j-1]:
                Exact_survival_rate[j] = Exact_survival_rate[j-1]
            if Basal_survival_rate[j] > Basal_survival_rate[j-1]:
                Basal_survival_rate[j] = Basal_survival_rate[j-1]
        for j in range(1, len(Pred_survival_rate)):
            if Pred_survival_rate[j] > Pred_survival_rate[j-1]:
                Pred_survival_rate[j] = Pred_survival_rate[j-1]

        exact_time, exact_event = convert_to_time_event(Exact_survival_rate)
        pred_time, pred_event = convert_to_time_event(Pred_survival_rate)
        basal_time, basal_event = convert_to_time_event(Basal_survival_rate)

        exact_time = [t * dt for t in exact_time]
        pred_time = [t * dt for t in pred_time]
        basal_time = [t * dt for t in basal_time]

        ks_statistic, pvalue = scipy.stats.ks_2samp(pred_time, exact_time)
        
        KS_list.append("{:.4g}".format(ks_statistic))
        KSp_list.append("{:.4g}".format(pvalue))

        # frameworks
        data_exact = pd.DataFrame({'time': exact_time, 'event': exact_event, 'group': ['Exact'] * len(exact_time)})
        data_pred = pd.DataFrame({'time': pred_time, 'event': pred_event, 'group': ['Pred'] * len(pred_time)})
        data_basal = pd.DataFrame({'time': basal_time, 'event': basal_event, 'group': ['Basal'] * len(basal_time)})

        df = pd.concat([data_exact, data_pred, data_basal])

        # groups
        group_exact = df[df['group'] == 'Exact']
        group_pred = df[df['group'] == 'Pred']
        group_basal = df[df['group'] == 'Basal']

        # KM fitting
        kmf_exact = lifelines.KaplanMeierFitter()
        kmf_exact.fit(group_exact['time'], event_observed=group_exact['event'], label='Exact')

        kmf_pred = lifelines.KaplanMeierFitter()
        kmf_pred.fit(group_pred['time'], event_observed=group_pred['event'], label='Pred')

        kmf_basal = lifelines.KaplanMeierFitter()
        kmf_basal.fit(group_basal['time'], event_observed=group_basal['event'], label='Basal')

        # log-rank test
        results_pred_exact = lifelines.statistics.logrank_test(group_pred['time'], group_exact['time'], event_observed_A=group_pred['event'], event_observed_B=group_exact['event'])
        results_exact_basal = lifelines.statistics.logrank_test(group_exact['time'], group_basal['time'], event_observed_A=group_exact['event'], event_observed_B=group_basal['event'])
        results_pred_basal = lifelines.statistics.logrank_test(group_pred['time'], group_basal['time'], event_observed_A=group_pred['event'], event_observed_B=group_basal['event'])

        # print
        #print(f'Exact vs Pred - p-value: {results_exact_pred.p_value}, test statistic: {results_exact_pred.test_statistic}')
        #print(f'Exact vs Basal - p-value: {results_exact_basal.p_value}, test statistic: {results_exact_basal.test_statistic}')
        #print(f'Pred vs Basal - p-value: {results_pred_basal.p_value}, test statistic: {results_pred_basal.test_statistic}')
        
        # survival functions
        survival_exact = kmf_exact.survival_function_
        survival_pred = kmf_pred.survival_function_

        # supplement missing time points
        survival_exact_interpolated = supplement_survival_function(t_max, survival_exact)
        survival_pred_interpolated = supplement_survival_function(t_max, survival_pred)

        # used to calculate p.d.f. later
        time_points_exact = survival_exact_interpolated.index.values
        survival_values_exact = survival_exact_interpolated.iloc[:, 0].values

        time_points_pred = survival_pred_interpolated.index.values
        survival_values_pred = survival_pred_interpolated.iloc[:, 0].values

        # MSE
        mse = np.mean((survival_values_pred - survival_values_exact) ** 2)
        # RMSE
        rmse = np.sqrt(mse)
        #relative_mse = mse / np.mean(survival_values_exact)
        rmse_list.append("{:.4g}".format(rmse))

        # 75% and 90% survival time
        survival_75_time = survival_pred_interpolated[survival_pred_interpolated['Pred'] <= 0.75].index.min()
        survival_90_time = survival_pred_interpolated[survival_pred_interpolated['Pred'] <= 0.9].index.min()

        p_value_list.append("{:.4g}".format(results_pred_exact.p_value))
        p_basal_list.append("{:.4g}".format(results_pred_basal.p_value))
        mst_list.append("{:.4g}".format(kmf_pred.median_survival_time_))
        tqst_list.append("{:.4g}".format(survival_75_time))
        nqst_list.append("{:.4g}".format(survival_90_time))

        # plot KM curves and note statistics in Fig 7D
        fig = plt.figure(figsize=(14, 12))
        kmf_exact.plot_survival_function(ci_show=False)
        kmf_pred.plot_survival_function(ci_show=False)
        kmf_basal.plot_survival_function(ci_show=False)

        # add notes
        plt.annotate(f'pred vs exact log-rank-test p-value: {results_pred_exact.p_value:.4g}', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=12, color='blue', horizontalalignment='right')
        plt.annotate(f'exact vs basal log-rank-test p-value: {results_exact_basal.p_value:.4g}', xy=(0.95, 0.1), xycoords='axes fraction', fontsize=12, color='blue', horizontalalignment='right')
        plt.annotate(f'pred vs basal log-rank-test p-value: {results_pred_basal.p_value:.4g}', xy=(0.95, 0.15), xycoords='axes fraction', fontsize=12, color='orange', horizontalalignment='right')
        plt.annotate(f'pred median survival time: {kmf_pred.median_survival_time_:.4g}', xy=(0.95, 0.2), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
        plt.annotate(f'pred vs exact RMSE: {rmse:.4g}', xy=(0.95, 0.25), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
        plt.annotate(f'pred vs excat KS: {ks_statistic:.4g}', xy=(0.95, 0.3), xycoords='axes fraction', fontsize=12, color='red', horizontalalignment='right')
        plt.annotate(f'pred vs exact KS p-value: {pvalue:.4g}', xy=(0.95, 0.35), xycoords='axes fraction', fontsize=12, color='red', horizontalalignment='right')

        plt.xlabel('time (days)')
        plt.ylabel('survival probability')
        plt.legend()
        plt.xlim(0, t_max)
        plt.ylim(0, 1)

        # save figure
        save_path = f'{figure_save_path}/KM_{i}.pdf'
        plt.savefig(save_path, format='pdf')
        print(f'KM_{i} done')

        # calculate gradients
        pdf_exact = -np.gradient(survival_values_exact, time_points_exact)
        pdf_pred = -np.gradient(survival_values_pred, time_points_pred)

        # avoid denominator being zero
        epsilon = 1e-7
        pdf_exact += epsilon
        pdf_pred += epsilon

        # integral
        excat_integral_result = scipy.integrate.simps(pdf_exact, time_points_exact)
        pdf_exact /= excat_integral_result
        pred_integral_result = scipy.integrate.simps(pdf_pred, time_points_pred)
        pdf_pred /= pred_integral_result

        # KL divergences
        kl_divergence = np.sum(scipy.special.rel_entr(pdf_pred, pdf_exact))
        KLD_list.append("{:.4g}".format(kl_divergence))

        # plot survival distribution in Fig 7E
        plt.figure(figsize=(10, 6))
        plt.plot(time_points_exact, pdf_exact, label='exact p.d.f', drawstyle='steps-post')
        plt.plot(time_points_pred, pdf_pred, label='pred p.d.f', drawstyle='steps-post')

        plt.title('Probability density functions of pred & exact')
        plt.xlabel('time (days)')
        plt.ylabel('density')
        plt.legend()
        # add notes
        plt.text(0.05, 0.7, f'KL Divergence: {kl_divergence:.4g}', transform=plt.gca().transAxes)

        # save figures
        save_path = f'{figure_save_path}/KM_pdf_{i}.pdf'
        plt.savefig(save_path, format='pdf')
        print(f'KM_pdf_{i} done')

    # print statistics
    print("log-rank-test p value", p_value_list) 
    print("p basal value", p_basal_list) 
    print("median survival time", mst_list)  
    print("three quarter survival time", tqst_list)
    print("90 percent survival time", nqst_list)
    print("KL divergence", KLD_list)
    print("KS statistic", KS_list)
    print("KS statistic p-value", KSp_list)
    print("RMSE", rmse_list)

# REFERENCE

# M. Raissi, P.Perdikaris, G.E.Ksrniadakis, Physics-informed neural networks: 
# A deep learning framework for solving forward and inverse problems involving 
# nonlinear partial differential equations. Journal of Computational Physics 
# 378, 686-707 (2019).

