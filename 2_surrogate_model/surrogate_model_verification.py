import torch
from collections import OrderedDict
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import lifelines

# set seed
seed = 20240712
np.random.seed(seed)
torch.manual_seed(seed)  

# CUDA support
if torch.cuda.is_available():
    device = torch.device('cuda')
    torch.cuda.manual_seed(seed)  
    torch.cuda.manual_seed_all(seed)
else:
    device = torch.device('cpu')

save_folder = f'test_and_verification'

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

# PINN based surrogate model used for reinfocement learning
class RL_dNN():
    def __init__(self, FP_layers, lb, ub, dc, dt_cyto): 
        # ode parameters
        self.q_c = 0.8
        self.eta_c = 2e-3 
        self.mu_c = 2e-3 
        self.q_i = 0.8
        self.eta_i = 2e-3 
        self.mu_i = 2e-3         
        # fokker-planck parameters
        #load parameters
        additional_params = torch.load('surrogate_model/additional_params.pth')  
        initial_log_p_cyto = additional_params['log_p_cyto']
        self.log_p_cyto  = torch.nn.Parameter(initial_log_p_cyto).to(device)        
        # deep neural networks
        self.dnn_TC = DNN(FP_layers).to(device)
        para_TC = torch.load('surrogate_model/dnn_TC_model1.pth')
        self.dnn_TC.load_state_dict(para_TC, strict=False)        
        
        # other parameters
        self.lb = torch.tensor(lb, dtype=torch.float32, requires_grad=True).float().to(device)
        self.ub = torch.tensor(ub, dtype=torch.float32, requires_grad=True).float().to(device)
        self.c_max = torch.tensor(1/dc, dtype=torch.float32, requires_grad=True).float().to(device)
        self.dt_cyto = dt_cyto
        self.CSF1RI_cum_max = 200 / self.dt_cyto
  
    # DNN for Fokker-Planck equation of tumor cell
    def net_u_CC(self, C, t, dose_c, dose_cum, dose_I): 
        t = 2.0*(t- self.lb[1])/(self.ub[1] - self.lb[1]) - 1.0 
        C = 2.0*(C- self.lb[0])/(self.ub[0] - self.lb[0]) - 1.0
        dose_c = 2.0*(dose_c - 0) - 1.0
        dose_cum = 2.0*(dose_cum - 0)/self.CSF1RI_cum_max - 1.0
        dose_I = 2.0*(dose_I - 0) - 1.0
        HTC = torch.cat([t, C, dose_c, dose_cum, dose_I], dim=1)
        p = self.dnn_TC(HTC)
        return p

    # ODE solution of CSF1R_I
    def f_CSF1RI(self, c_CSF1RI_0, dose_c, interval):
        itera_cyto = int(interval/self.dt_cyto)
        c = torch.empty(itera_cyto)
        a = torch.exp(self.log_p_cyto[10])*self.q_c
        b = torch.exp(self.log_p_cyto[11])*self.eta_c + torch.exp(self.log_p_cyto[12])*self.mu_c #d_M1 + d_M2 + d_M0 = 1 
        c[0] = (c_CSF1RI_0 - a*dose_c/(a+b))*torch.exp(-(a+b)*self.dt_cyto) + a*dose_c/(a+b) 
        for k in range(itera_cyto-1):
            c[k+1] = (c[k] - a*dose_c/(a+b))*torch.exp(-(a+b)*self.dt_cyto) + a*dose_c/(a+b)        
        return c

    # ODE solution of IGF1R_I
    def f_IGF1RI(self, c_IGF1RI_0, dose_I, C_T, interval):
        itera_cyto = int(interval/self.dt_cyto)
        c = torch.empty(itera_cyto)
        a = torch.exp(self.log_p_cyto[13])*self.q_i
        b = torch.exp(self.log_p_cyto[14])*self.eta_i + torch.exp(self.log_p_cyto[15])*self.mu_i * C_T
        c[0] = (c_IGF1RI_0 - a*dose_I/(a+b))*torch.exp(-(a+b)*self.dt_cyto) + a*dose_I/(a+b)
        for k in range(itera_cyto-1):
            c[k+1] = (c[k] - a*dose_I/(a+b))*torch.exp(-(a+b)*self.dt_cyto) + a*dose_I/(a+b)        
        return c
    
    # predict Fokker-Planck equation of TC for RL
    def predict_DRL_FP(self, c, t, dose_c, CSF1RI_cum, dose_I, c_CSF1RI_0, c_IGF1RI_0, c_T_0, interval):
        self.dnn_TC.eval()  

        c = torch.tensor(c, dtype=torch.float32, requires_grad=True).float().to(device)
        t = torch.tensor(t, dtype=torch.float32, requires_grad=True).float().to(device)
        t_k = t * torch.ones_like(c)

        c_CSF1RI_0 = torch.tensor(c_CSF1RI_0, dtype=torch.float32, requires_grad=True).float().to(device)
        c_IGF1RI_0 = torch.tensor(c_IGF1RI_0, dtype=torch.float32, requires_grad=True).float().to(device)
        c_T_0 = torch.tensor(c_T_0, dtype=torch.float32, requires_grad=True).float().to(device)

        dose_c = torch.tensor(dose_c, dtype=torch.float32, requires_grad=True).float().to(device)
        CSF1RI_cum = torch.tensor(CSF1RI_cum, dtype=torch.float32, requires_grad=True).float().to(device)
        dose_I = torch.tensor(dose_I, dtype=torch.float32, requires_grad=True).float().to(device)
        
        c_CSF1RI = self.f_CSF1RI(c_CSF1RI_0, dose_c, interval)
        c_IGF1RI = self.f_IGF1RI(c_IGF1RI_0, dose_I, c_T_0, interval)
        CSF1RI_cum = CSF1RI_cum + torch.sum(c_CSF1RI)

        c_CSF1RI_k = c_CSF1RI[-1] * torch.ones_like(c)
        CSF1RI_cum_k = CSF1RI_cum * torch.ones_like(c)
        c_IGF1RI_k = c_IGF1RI[-1] * torch.ones_like(c)

        uTC_pred = self.net_u_CC(c, t_k, c_CSF1RI_k, CSF1RI_cum_k, c_IGF1RI_k)

        uTC_pred = uTC_pred.detach().cpu().numpy()
        CSF1RI_cum = CSF1RI_cum.detach().cpu().numpy()
        c_CSF1RI_k = c_CSF1RI_k[-1, 0].detach().cpu().numpy()
        c_IGF1RI_k = c_IGF1RI_k[-1, 0].detach().cpu().numpy()

        uTC_pred[uTC_pred < 0] = 0

        return uTC_pred, CSF1RI_cum, c_CSF1RI_k, c_IGF1RI_k

# initialize parameters of policy net and value net
def normalized_columns_initializer(weights, std):
    with torch.no_grad():
        out = torch.randn_like(weights)
        out *= std / torch.sqrt(out.pow(2).sum(dim=0, keepdim=True))
        return out

# long short-term memory (LSTM) net
class LSTMNet(torch.nn.Module):
    def __init__(self, input_size, lstm_hidden_dim, num_layers):
        super(LSTMNet, self).__init__()
        self.lstm = torch.nn.LSTM(input_size, lstm_hidden_dim, num_layers, batch_first=True)

    def forward(self, x):
        batch_size, seq_len, _ = x.shape
        lstm_out, _ = self.lstm(x)
        hidden_size = lstm_out.shape[2]
        lstm_out_flattened = lstm_out.view(batch_size * seq_len, hidden_size)
        return lstm_out_flattened
 
# policy net, actor               
class PolicyNet(torch.nn.Module):
    def __init__(self, lstm_hidden_dim, action_size):
        super(PolicyNet, self).__init__()
        # fully connected layers
        self.fc_layers = torch.nn.Sequential(
            torch.nn.Linear(lstm_hidden_dim, 128),
            torch.nn.ReLU(),
            torch.nn.Linear(128, 64),
            torch.nn.ReLU(),
            torch.nn.Linear(64, 32),
            torch.nn.ReLU(),
            torch.nn.Linear(32, 16),
            torch.nn.ReLU(),
            torch.nn.Linear(16, 8),
            torch.nn.ReLU(),
            torch.nn.Linear(8, action_size)
        )
        for m in self.fc_layers:
            if isinstance(m, torch.nn.Linear):
                m.weight.data = normalized_columns_initializer(m.weight.data, std=0.01)

    def forward(self, x):
        fc_out = self.fc_layers(x)
        policy_dist = torch.nn.functional.softmax(fc_out, dim=1)
        return policy_dist

# value net, critic 
class ValueNet(torch.nn.Module):
    def __init__(self, lstm_hidden_dim):
        super(ValueNet, self).__init__()
        # fully connected layers
        self.fc_layers = torch.nn.Sequential(
            torch.nn.Linear(lstm_hidden_dim, 128),
            torch.nn.ReLU(),
            torch.nn.Linear(128, 64),
            torch.nn.ReLU(),
            torch.nn.Linear(64, 32),
            torch.nn.ReLU(),
            torch.nn.Linear(32, 16),
            torch.nn.ReLU(),
            torch.nn.Linear(16, 8),
            torch.nn.ReLU(),
            torch.nn.Linear(8, 1)
        )
        for m in self.fc_layers:
            if isinstance(m, torch.nn.Linear):
                m.weight.data = normalized_columns_initializer(m.weight.data, std=1.0)

    def forward(self, x):
        value_estimate = self.fc_layers(x)
        return value_estimate.squeeze(1)

# RL environmen       
class Env:
    def __init__(self, model, initial):  
        self.state = [] # state of RL environment
        self.action_c_space = np.array([0.0, 1.0]) # action space of CSF1R_I dose 
        self.action_I_space = np.array([0.0, 1.0]) # action space of IGF1R_I dose   
        self.surrogat_model = model 
        self.initial = initial
        self.predict_CSF1RI_list = []
        self.predict_CSF1RI_list.append(self.initial[0])
        self.predict_IGF1RI_list = []
        self.predict_IGF1RI_list.append(self.initial[1])
        self.predict_TC_list = []
        self.predict_TC_list.append(self.initial[2])
        self.predict_CSF1RI_cum = 0
        self.CSF1RI_dose_cum = 0
        self.last_state_survival = 0
        self.time = 0

    # reset RL envrionment
    def reset(self):  
        self.state = []
        self.predict_CSF1RI_list = []
        self.predict_CSF1RI_list.append(self.initial[0])
        self.predict_IGF1RI_list = []
        self.predict_IGF1RI_list.append(self.initial[1])
        self.predict_TC_list = []
        self.predict_TC_list.append(self.initial[2])
        self.predict_CSF1RI_cum = 0
        self.CSF1RI_dose_cum = 0
        self.last_state_survival = 0
        self.time = 0

    # simulate each step in RL envrionment and calculate reward
    def step(self, action_c, action_I, epi, predict_dt, cc, lamda):  
        last_state = self.state.copy()
        done = False
        self.time += predict_dt
        dose_c = self.action_c_space[action_c]
        dose_I = self.action_I_space[action_I]
        self.CSF1RI_dose_cum += dose_c * predict_dt
        c_CSF1RI_0 = self.predict_CSF1RI_list[-1]
        c_IGF1RI_0 = self.predict_IGF1RI_list[-1]
        c_T_0 = self.predict_TC_list[-1]

        # surrogate model 
        rl_FP_TC, c_CSF1RI_cum, c_CSF1RI, c_IGF1RI = self.surrogat_model.predict_DRL_FP(cc, self.time, dose_c, self.predict_CSF1RI_cum, dose_I, 
                                                   c_CSF1RI_0, c_IGF1RI_0, c_T_0, predict_dt)

        # normalization
        rl_FP_TC_norm = rl_FP_TC/(np.sum(rl_FP_TC) + 1e-7)

        self.predict_CSF1RI_cum = c_CSF1RI_cum
        self.predict_CSF1RI_list.append(c_CSF1RI)
        self.predict_IGF1RI_list.append(c_IGF1RI)
        self.predict_TC_list.append(np.sum(rl_FP_TC_norm*cc))

        prob_sum = np.sum(rl_FP_TC)
        survival_prob = np.sum(rl_FP_TC[:-3, :]) / (prob_sum + 1e-7)    
        death_prob = 1 - survival_prob
        cure_prob = np.sum(rl_FP_TC[:1, :]) / (prob_sum + 1e-7)

        # punish no treatment in high death probability case
        dose_punish = (self.action_c_space.shape[0] - 1) + 1e-7

        # set thresholds
        death_threshold = 0.2
        cure_threshold = 0.99

        # calculate the change in survival probability
        state_survival = (survival_prob - (1 - death_threshold)) / death_threshold
        if epi == 0:
            survival_prob_delta = 0
        else:
            survival_prob_delta = state_survival - self.last_state_survival
        
        # update state of RL environment
        self.state.append([state_survival, survival_prob_delta])
        self.last_state_survival = state_survival
        now_state = self.state.copy()

        # reward function
        reward = 0.1
        if death_prob >= death_threshold:
            reward = -0.1
            done = True
        elif cure_prob >= cure_threshold:
            reward = reward + 1.0
            done = True
        if not done:
            if death_prob >= death_threshold * 0.5:
                reward = reward - 0.1 * (dose_punish - action_c)
            elif action_c == 0:
                reward = reward + 0.05

            if survival_prob <= 0.9:
                reward = reward - 0.1 * (dose_punish - action_I)
            elif action_I == 0:
                reward = reward + 0.05

            if survival_prob > 0.9:
                reward = reward + lamda * (self.time - self.CSF1RI_dose_cum)    

        return last_state, reward, done, now_state

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


def main():
    FP_layers = [5, 40, 120, 250, 500, 1000, 1000, 600, 300, 150, 1]
    predict_dt = 1/48
    dt = predict_dt
    dt_cyto = 1/48
    t_max = 200
    
    t = np.arange(dt, t_max + dt, dt)
    t_cyto = np.arange(dt_cyto, t_max + dt_cyto, dt_cyto)

    dc = 0.01
    c_min, c_max = 0, 1
    c = np.arange(c_min + dc, c_max + dc, dc)

    t = np.array(t)[:, None]
    t_cyto = np.array(t_cyto)[:, None]
    c = np.array(c)[:, None]

    c_show = np.insert(c, 0, 0) 
    t_show = np.insert(t, 0, 0)
    T_show, C_show = np.meshgrid(t_show, c_show)

    initial_TC = 0.58
    initial_CSF1RI = 0
    initial_IGF1RI = 0
    C, T = np.meshgrid(c,t)
    C_star = np.hstack((C.flatten()[:,None], T.flatten()[:,None]))

    lb = [dc, predict_dt]
    ub = [1.0, 200]   

    surrogate_model = RL_dNN(FP_layers, lb, ub, dc, dt_cyto)
    env = Env(surrogate_model, [initial_CSF1RI, initial_IGF1RI, initial_TC])

    action_list = []

    # 'switch' treatment in experiment
    actions_list = []
    for i in range(7*4 * 48):
        action_list.append([1, 0])
    for i in range(172 * 48):
        action_list.append([0, 1])
    actions_list.append(action_list)

    # 'add' treatment in experiment
    # the same as 1-week interval RL treatment
    action_list = []
    for i in range(7*4 * 48):
        action_list.append([1, 0])
    for i in range(172 * 48):
        action_list.append([1, 1])
    actions_list.append(action_list)
    
    # surrogate model prediction vs 'add' treatment in experiment
    predict_actions = actions_list[1]
    
    predict_CSF1RI_cum = 0
    predict_CC_list = []
    predict_dose_c_list = []
    predict_dose_I_list = []
    predict_CSF1RI_list = []
    predict_CSF1RI_list.append(initial_CSF1RI)
    predict_IGF1RI_list = []
    predict_IGF1RI_list.append(initial_IGF1RI)
    predict_TC_list = []
    predict_TC_list.append(initial_TC)

    length_time = len(predict_actions)
    for k in range(length_time):

        predict_dose_c = env.action_c_space[predict_actions[k][0]]   
        predict_dose_I = env.action_I_space[predict_actions[k][1]]
        predict_tt = (k+1) * predict_dt
        c_CSF1RI_0 = predict_CSF1RI_list[-1]
        c_IGF1RI_0 = predict_IGF1RI_list[-1]
        c_T_0 = predict_TC_list[-1]
        rl_FP_TC, c_CSF1RI_cum, c_CSF1RI, c_IGF1RI = surrogate_model.predict_DRL_FP(c, predict_tt, predict_dose_c, predict_CSF1RI_cum, predict_dose_I,
                                                                c_CSF1RI_0, c_IGF1RI_0, c_T_0, predict_dt)
        predict_CSF1RI_cum = c_CSF1RI_cum
        rl_FP_TC_norm = rl_FP_TC/(np.sum(rl_FP_TC) + 1e-7)

        predict_CC_list.append(rl_FP_TC_norm * 1/dc)
        predict_CSF1RI_list.append(c_CSF1RI)
        predict_IGF1RI_list.append(c_IGF1RI)
        predict_TC_list.append(np.sum(rl_FP_TC_norm*c))
        predict_dose_c_list.append(predict_dose_c)
        predict_dose_I_list.append(predict_dose_I)

    predict_CClist_flat = np.concatenate(predict_CC_list).ravel()
    dose_c = np.array(predict_dose_c_list)
    dose_I = np.array(predict_dose_I_list)
    predict_CSF1RI = np.array(predict_CSF1RI_list)
    predict_IGF1RI = np.array(predict_IGF1RI_list)
    predict_TC = np.array(predict_TC_list)

    # plot and save Pred_TC_add
    U_pred = scipy.interpolate.griddata(C_star, predict_CClist_flat, (C, T), method='cubic')
    U_pred[U_pred < 0] = 0
    U_show = U_pred.T
    fig, ax = plt.subplots(figsize=(9, 5))
    h = ax.pcolormesh(T_show, C_show, U_show, shading='auto', cmap='rainbow')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.10)
    cbar = fig.colorbar(h, cax=cax)
    cbar.ax.tick_params(labelsize=15)
    ax.set_xlabel('$t$', size=20)
    ax.set_ylabel('$C$', size=20)
    ax.set_title('$p(t,C)$', fontsize=20)
    ax.tick_params(labelsize=15)
    t_ticks = np.arange(0, t_max + 1e-7, 20)
    ax.set_xticks(t_ticks)
    ax.tick_params(axis='x', labelsize=6)
    save_path = f'{save_folder}/Pred_TC_add.pdf'
    plt.savefig(save_path, format='pdf')

    # plot and save 'add' treatment
    # Fig 8A
    fig = plt.figure(figsize=(5, 12))
    ax = fig.add_subplot(111)
    gs1 = gridspec.GridSpec(2, 1)
    gs1.update(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.3, hspace=0.3)
    
    ax = plt.subplot(gs1[0, 0])
    ax.plot(t_cyto, dose_c, 'b-', linewidth = 2, label = 'CSF1RI dose')  
    ax.plot(t_show, predict_CSF1RI, 'k-', linewidth = 2, label = 'invivo CSF1RI')   
    ax.set_xlabel('$t$')
    ax.set_ylabel('$c$')    
    ax.set_title('$CSF1RI dose$', fontsize = 15)
    ax.axis('square')
    ax.set_xlim([0.0,t_max])
    ax.set_ylim([0.0,1])
    aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.set_aspect(aspect_ratio)
    ax.set_xticks(t_ticks)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    plt.legend()

    ax = plt.subplot(gs1[1, 0])
    ax.plot(t_cyto, dose_I, 'g-', linewidth = 2, label = 'IGF1RI dose')  
    ax.plot(t_show, predict_IGF1RI, 'k-', linewidth = 2, label = 'invivo IGF1RI')     
    ax.plot(t_show, predict_TC, 'r-', linewidth = 2, label = 'TC expectation')  
    ax.set_xlabel('$t$')
    ax.set_ylabel('$c$')    
    ax.set_title('$IGF1RI dose$', fontsize = 15)
    ax.axis('square')
    ax.set_xlim([0.0,t_max])
    ax.set_ylim([0.0,1])
    aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.set_aspect(aspect_ratio)
    ax.set_xticks(t_ticks)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    plt.legend()
    save_path = f'{save_folder}/dose_add.pdf'
    plt.savefig(save_path, format='pdf')
    print(f'dose_add done')  


    # plot and save KM curves
    # log-rank test on surrogate model result
    # compared with experiment and MSABM result
    # Fig 8C & D
    experiemnt_add_filepath = f'experiment_data/experiment_CSF1RI_add_IGF1RI.csv'
    experiemnt_data = pd.read_csv(experiemnt_add_filepath, header=None) # 9 mice
    add_data = experiemnt_data.values

    ABM_add_filepath = f'MSABM_data/CC_add_pdf_20_100.csv'
    ABM_data = pd.read_csv(ABM_add_filepath, header=None) # p.d.f. of 100 virtual mice * 200 days
    ABM_add_data = ABM_data.values # p.d.f. of 100 virtual mice * 196 days

    Pred_survival_rate = (1 - np.sum(U_show[-14:, :], axis = 0)/np.sum(U_show, axis = 0))
    ABM_survival_rate = (1 - np.sum(ABM_add_data[-3:, :], axis = 0)/np.sum(ABM_add_data, axis = 0))
    
    ABM_survival_rate[0] = 1
    Pred_survival_rate[0] = 1
    for j in range(1, len(ABM_survival_rate)):
        if j <= 14:
            ABM_survival_rate[j] = 1
        if ABM_survival_rate[j] > ABM_survival_rate[j-1]:
            ABM_survival_rate[j] = ABM_survival_rate[j-1]

    for j in range(1, len(Pred_survival_rate)):
        if j <= int(28 / predict_dt):
            Pred_survival_rate[j] = 1
        if Pred_survival_rate[j] > Pred_survival_rate[j-1]:
            Pred_survival_rate[j] = Pred_survival_rate[j-1]

    pred_time, pred_event = convert_to_time_event(Pred_survival_rate)
    ABM_time, ABM_event = convert_to_time_event(ABM_survival_rate)

    pred_time = [t * predict_dt for t in pred_time]
    ABM_time = [t * 2 for t in ABM_time]
    exact_time = predict_dt * np.ceil(add_data[:, 0]/predict_dt)
    ks_statistic_ABM, pvalue_ABM = scipy.stats.ks_2samp(ABM_time, exact_time)
    ks_statistic_FP, pvalue_FP = scipy.stats.ks_2samp(pred_time, exact_time)

    data_pred = pd.DataFrame({'time': pred_time, 'event': pred_event, 'group': ['Pred'] * len(pred_time)})
    data_exp = pd.DataFrame({'time': exact_time, 'event': add_data[:, 1], 'group': ['Exp'] * len(add_data[:, 0])})
    data_ABM = pd.DataFrame({'time': ABM_time, 'event': ABM_event, 'group': ['ABM'] * len(ABM_time)})
    df = pd.concat([data_pred, data_exp, data_ABM])

    group_pred = df[df['group'] == 'Pred']
    group_exp = df[df['group'] == 'Exp']
    group_ABM = df[df['group'] == 'ABM']

    kmf_pred = lifelines.KaplanMeierFitter()
    kmf_pred.fit(group_pred['time'], event_observed=group_pred['event'], label='Pred')

    kmf_exp = lifelines.KaplanMeierFitter()
    kmf_exp.fit(group_exp['time'], event_observed=group_exp['event'], label='Exp')

    kmf_ABM = lifelines.KaplanMeierFitter()
    kmf_ABM.fit(group_ABM['time'], event_observed=group_ABM['event'], label='ABM')
 
    results_pred_ABM = lifelines.statistics.logrank_test(group_pred['time'], group_ABM['time'], event_observed_A=group_pred['event'], event_observed_B=group_ABM['event'])
    results_ABM_exp = lifelines.statistics.logrank_test(group_ABM['time'], group_exp['time'], event_observed_A=group_ABM['event'], event_observed_B=group_exp['event'])
    results_pred_exp = lifelines.statistics.logrank_test(group_pred['time'], group_exp['time'], event_observed_A=group_pred['event'], event_observed_B=group_exp['event'])
    
    survival_exp = kmf_exp.survival_function_
    survival_ABM = kmf_ABM.survival_function_
    survival_pred = kmf_pred.survival_function_

    survival_exp_interpolated = supplement_survival_function(t_max, predict_dt, survival_exp)
    survival_ABM_interpolated = supplement_survival_function(t_max, predict_dt, survival_ABM)
    survival_pred_interpolated = supplement_survival_function(t_max, predict_dt, survival_pred)

    kmf_pred.plot_survival_function(ci_show=False)
    kmf_exp.plot_survival_function(ci_show=False)
    kmf_ABM.plot_survival_function(ci_show=False)

    survival_90_time_pred = survival_pred_interpolated[survival_pred_interpolated['Pred'] <= 0.9].index.min()
    survival_90_time_exp = survival_exp_interpolated[survival_exp_interpolated['Exp'] <= 0.9].index.min()
    survival_90_time_ABM = survival_ABM_interpolated[survival_ABM_interpolated['ABM'] <= 0.9].index.min()
    
    survival_values_exp = survival_exp_interpolated.iloc[:, 0].values
    survival_values_pred = survival_pred_interpolated.iloc[:, 0].values
    survival_values_ABM = survival_ABM_interpolated.iloc[:, 0].values

    # RMSE
    mse_FP = np.mean((survival_values_pred - survival_values_exp) ** 2)
    rmse_FP = np.sqrt(mse_FP)
    #relative_mse_FP = "{:.4g}".format(mse_FP / np.mean(survival_values_exp))

    mse_ABM = np.mean((survival_values_ABM - survival_values_exp) ** 2)
    rmse_ABM = np.sqrt(mse_ABM)
    #relative_mse_ABM = "{:.4g}".format(mse_ABM / np.mean(survival_values_exp))

    # plot and save
    fig = plt.figure(figsize=(14, 12))
    # add notes
    plt.annotate(f'ABM 90% survival time: {survival_90_time_ABM:.4g}', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'FP model 90% survival time: {survival_90_time_pred:.4g}', xy=(0.95, 0.1), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'exp 90% survival time: {survival_90_time_exp:.4g}', xy=(0.95, 0.15), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'ABM vs exp RMSE: {rmse_ABM}', xy=(0.95, 0.2), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')    
    plt.annotate(f'FP model vs exp RMSE: {rmse_FP}', xy=(0.95, 0.25), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'ABM vs exp KS: {ks_statistic_ABM:.4g}', xy=(0.95, 0.4), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'ABM vs exp log-rank test p-value: {results_ABM_exp.p_value:.4g}', xy=(0.95, 0.45), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'FP model vs exp KS: {ks_statistic_FP:.4g}', xy=(0.95, 0.5), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')    
    plt.annotate(f'FP model vs exp log-rank test p-value: {results_pred_exp.p_value:.4g}', xy=(0.95, 0.55), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')

    plt.xlabel('time (days)')
    plt.ylabel('survival probability')
    plt.legend()
    plt.xlim(0, t_max)
    plt.ylim(0, 1)
    plt.xticks(t_ticks)
    save_path = f'{save_folder}/KM_add.pdf'
    plt.savefig(save_path, format='pdf')
    print(f'KM done')   

    # surrogate model prediction vs 'switch' treatment in experiment
    predict_actions = actions_list[0]

    predict_CSF1RI_cum = 0
    predict_CC_list = []
    predict_dose_c_list = []
    predict_dose_I_list = []
    predict_CSF1RI_list = []
    predict_CSF1RI_list.append(initial_CSF1RI)
    predict_IGF1RI_list = []
    predict_IGF1RI_list.append(initial_IGF1RI)
    predict_TC_list = []
    predict_TC_list.append(initial_TC)
    
    length_time = len(predict_actions)
    for k in range(length_time):

        predict_dose_c = env.action_c_space[predict_actions[k][0]]   
        predict_dose_I = env.action_I_space[predict_actions[k][1]]
        predict_tt = (k+1) * predict_dt
        c_CSF1RI_0 = predict_CSF1RI_list[-1]
        c_IGF1RI_0 = predict_IGF1RI_list[-1]
        c_T_0 = predict_TC_list[-1]
        rl_FP_TC, c_CSF1RI_cum, c_CSF1RI, c_IGF1RI = surrogate_model.predict_DRL_FP(c, predict_tt, predict_dose_c, predict_CSF1RI_cum, predict_dose_I,
                                                                c_CSF1RI_0, c_IGF1RI_0, c_T_0, predict_dt)
        predict_CSF1RI_cum = c_CSF1RI_cum
        rl_FP_TC_norm = rl_FP_TC/(np.sum(rl_FP_TC) + 1e-7)
        predict_CC_list.append(rl_FP_TC_norm * 1/dc)
        predict_CSF1RI_list.append(c_CSF1RI)
        predict_IGF1RI_list.append(c_IGF1RI)
        predict_TC_list.append(np.sum(rl_FP_TC_norm*c))
        predict_dose_c_list.append(predict_dose_c)
        predict_dose_I_list.append(predict_dose_I)

    predict_CClist_flat = np.concatenate(predict_CC_list).ravel()
    dose_c = np.array(predict_dose_c_list)
    dose_I = np.array(predict_dose_I_list)
    predict_CSF1RI = np.array(predict_CSF1RI_list)
    predict_IGF1RI = np.array(predict_IGF1RI_list)
    predict_TC = np.array(predict_TC_list)

    # plot and save Pred_TC_add
    # Fig 7I
    U_pred = scipy.interpolate.griddata(C_star, predict_CClist_flat, (C, T), method='cubic')
    U_pred[U_pred < 0] = 0
    T_show, C_show = np.meshgrid(t_show, c_show)
    U_show = U_pred.T

    fig, ax = plt.subplots(figsize=(9, 5))
    h = ax.pcolormesh(T_show, C_show, U_show, shading='auto', cmap='rainbow')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.10)
    cbar = fig.colorbar(h, cax=cax)
    cbar.ax.tick_params(labelsize=15)
    ax.set_xlabel('$t$', size=20)
    ax.set_ylabel('$C$', size=20)
    ax.set_title('$p(t,C)$', fontsize=20)
    ax.tick_params(labelsize=15)
    t_ticks = np.arange(0, t_max + 1e-7, 20)
    ax.set_xticks(t_ticks)
    ax.tick_params(axis='x', labelsize=6)
    save_path = f'{save_folder}/Pred_TC_switch.pdf'
    plt.savefig(save_path, format='pdf')

    # plot and save 'switch' treatment
    # Fig 7F
    fig = plt.figure(figsize=(5, 12))
    ax = fig.add_subplot(111)
    gs1 = gridspec.GridSpec(2, 1)
    gs1.update(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.3, hspace=0.3)

    ax = plt.subplot(gs1[0, 0])
    ax.plot(t_cyto, dose_c, 'b-', linewidth = 2, label = 'CSF1RI dose')  
    ax.plot(t_show, predict_CSF1RI, 'k-', linewidth = 2, label = 'invivo CSF1RI')  
    ax.set_xlabel('$t$')
    ax.set_ylabel('$c$')    
    ax.set_title('$CSF1RI dose$', fontsize = 15)
    ax.axis('square')
    ax.set_xlim([0.0,t_max])
    ax.set_ylim([0.0,1])
    aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.set_aspect(aspect_ratio)
    ax.set_xticks(t_ticks)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    plt.legend()

    ax = plt.subplot(gs1[1, 0])
    ax.plot(t_cyto, dose_I, 'g-', linewidth = 2, label = 'IGF1RI dose')  
    ax.plot(t_show, predict_IGF1RI, 'k-', linewidth = 2, label = 'invivo IGF1RI')     
    ax.plot(t_show, predict_TC, 'r-', linewidth = 2, label = 'TC expectation')  
    ax.set_xlabel('$t$')
    ax.set_ylabel('$c$')    
    ax.set_title('$IGF1RI dose$', fontsize = 15)
    ax.axis('square')
    ax.set_xlim([0.0,t_max])
    ax.set_ylim([0.0,1])
    aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.set_aspect(aspect_ratio)
    ax.set_xticks(t_ticks)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    plt.legend()
    save_path = f'{save_folder}/dose_switch.pdf'
    plt.savefig(save_path, format='pdf')
    print(f'dose_switch done')  

    # plot and save KM curves
    # log-rank test on surrogate model result
    # compared with experiment and MSABM result
    # Fig 7H & J
    experiemnt_switch_filepath = f'experiment_data/experiment_CSF1RI_switch_IGF1RI.csv'
    experiemnt_data = pd.read_csv(experiemnt_switch_filepath, header=None) #100*200
    switch_data = experiemnt_data.values

    ABM_switch_filepath = f'MSABM_data/CC_switch_pdf_20_100.csv'
    ABM_data = pd.read_csv(ABM_switch_filepath, header=None) #100*200
    ABM_switch_data = ABM_data.values

    U_pred = scipy.interpolate.griddata(C_star, predict_CClist_flat, (C, T), method='cubic')
    U_pred[U_pred < 0] = 0
    Pred_survival_rate = (1 - np.sum(U_show[-15:, :], axis = 0)/np.sum(U_show, axis = 0))
    ABM_survival_rate = (1 - np.sum(ABM_switch_data[-3:, :], axis = 0)/np.sum(ABM_switch_data, axis = 0))
    Pred_survival_rate = Pred_survival_rate[:int(t_max/2/predict_dt)]
    ABM_survival_rate = ABM_survival_rate[:50]
    ABM_survival_rate[0] = 1
    Pred_survival_rate[0] = 1
    for j in range(1, len(ABM_survival_rate)):
        if j <= 14:
            ABM_survival_rate[j] = 1
        if ABM_survival_rate[j] > ABM_survival_rate[j-1]:
            ABM_survival_rate[j] = ABM_survival_rate[j-1]

    for j in range(1, len(Pred_survival_rate)):
        if j <= int(28/predict_dt):
            Pred_survival_rate[j] = 1
        if Pred_survival_rate[j] > Pred_survival_rate[j-1]:
            Pred_survival_rate[j] = Pred_survival_rate[j-1]

    pred_time, pred_event = convert_to_time_event(Pred_survival_rate)
    ABM_time, ABM_event = convert_to_time_event(ABM_survival_rate)

    pred_time = [t * predict_dt for t in pred_time]
    ABM_time = [t * 2 for t in ABM_time]
    exact_time = predict_dt * np.ceil(switch_data[:, 0]/predict_dt)
    ks_statistic_ABM, pvalue_ABM = scipy.stats.ks_2samp(ABM_time, exact_time)
    ks_statistic_FP, pvalue_FP = scipy.stats.ks_2samp(pred_time, exact_time)

    data_pred = pd.DataFrame({'time': pred_time, 'event': pred_event, 'group': ['Pred'] * len(pred_time)})
    data_exp = pd.DataFrame({'time': exact_time, 'event': switch_data[:, 1], 'group': ['Exp'] * len(switch_data[:, 0])})
    data_ABM = pd.DataFrame({'time': ABM_time, 'event': ABM_event, 'group': ['ABM'] * len(ABM_time)})

    df = pd.concat([data_pred, data_exp, data_ABM])

    group_pred = df[df['group'] == 'Pred']
    group_exp = df[df['group'] == 'Exp']
    group_ABM = df[df['group'] == 'ABM']

    kmf_pred = lifelines.KaplanMeierFitter()
    kmf_pred.fit(group_pred['time'], event_observed=group_pred['event'], label='Pred')

    kmf_exp = lifelines.KaplanMeierFitter()
    kmf_exp.fit(group_exp['time'], event_observed=group_exp['event'], label='Exp')

    kmf_ABM = lifelines.KaplanMeierFitter()
    kmf_ABM.fit(group_ABM['time'], event_observed=group_ABM['event'], label='ABM')
 
    # log-rank test
    results_pred_ABM = lifelines.statistics.logrank_test(group_pred['time'], group_ABM['time'], event_observed_A=group_pred['event'], event_observed_B=group_ABM['event'])
    results_ABM_exp = lifelines.statistics.logrank_test(group_ABM['time'], group_exp['time'], event_observed_A=group_ABM['event'], event_observed_B=group_exp['event'])
    results_pred_exp = lifelines.statistics.logrank_test(group_pred['time'], group_exp['time'], event_observed_A=group_pred['event'], event_observed_B=group_exp['event'])
    
    survival_exp = kmf_exp.survival_function_
    survival_ABM = kmf_ABM.survival_function_
    survival_pred = kmf_pred.survival_function_

    survival_exp_interpolated = supplement_survival_function(t_max/2, predict_dt, survival_exp)
    survival_ABM_interpolated = supplement_survival_function(t_max/2, predict_dt, survival_ABM)
    survival_pred_interpolated = supplement_survival_function(t_max/2, predict_dt, survival_pred)

    kmf_pred.plot_survival_function(ci_show=False)
    kmf_exp.plot_survival_function(ci_show=False)
    kmf_ABM.plot_survival_function(ci_show=False)

    survival_90_time_pred = survival_pred_interpolated[survival_pred_interpolated['Pred'] <= 0.9].index.min()
    survival_90_time_exp = survival_exp_interpolated[survival_exp_interpolated['Exp'] <= 0.9].index.min()
    survival_90_time_ABM = survival_ABM_interpolated[survival_ABM_interpolated['ABM'] <= 0.9].index.min()
    
    survival_values_exp = survival_exp_interpolated.iloc[:, 0].values
    survival_values_pred = survival_pred_interpolated.iloc[:, 0].values
    survival_values_ABM = survival_ABM_interpolated.iloc[:, 0].values

    # RMSE
    mse_FP = np.mean((survival_values_pred - survival_values_exp) ** 2)
    rmse_FP = np.sqrt(mse_FP)
    #relative_mse_FP = "{:.4g}".format(mse_FP / np.mean(survival_values_exp))

    mse_ABM = np.mean((survival_values_ABM - survival_values_exp) ** 2)
    rmse_ABM = np.sqrt(mse_ABM)
    #relative_mse_ABM = "{:.4g}".format(mse_ABM / np.mean(survival_values_exp))

    # plot 
    fig = plt.figure(figsize=(14, 12))
    # add notes
    plt.annotate(f'ABM 90% survival time: {survival_90_time_ABM:.4g}', xy=(0.95, 0.05), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'FP model 90% survival time: {survival_90_time_pred:.4g}', xy=(0.95, 0.1), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'exp 90% survival time: {survival_90_time_exp:.4g}', xy=(0.95, 0.15), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'ABM vs exp RMSE: {rmse_ABM}', xy=(0.95, 0.2), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')    
    plt.annotate(f'FP model vs exp RMSE: {rmse_FP}', xy=(0.95, 0.25), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')

    plt.annotate(f'ABM vs exp KS: {ks_statistic_ABM:.4g}', xy=(0.95, 0.4), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'ABM vs exp log-rank test p-value: {results_ABM_exp.p_value:.4g}', xy=(0.95, 0.45), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')
    plt.annotate(f'FP model vs exp KS: {ks_statistic_FP:.4g}', xy=(0.95, 0.5), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')    
    plt.annotate(f'FP model vs exp log-rank test p-value: {results_pred_exp.p_value:.4g}', xy=(0.95, 0.55), xycoords='axes fraction', fontsize=12, color='black', horizontalalignment='right')

    plt.xlabel('time (days)')
    plt.ylabel('survival probability')
    plt.legend()
    plt.xlim(0, t_max/2)
    plt.ylim(0, 1)
    t_ticks = np.arange(0, t_max/2 + 1e-7, 10)
    plt.xticks(t_ticks)
    save_path = f'{save_folder}/KM_switch.pdf'
    plt.savefig(save_path, format='pdf')
    print(f'KM done')   

    
    # most effective RL treatment
    action_list = []
    for i in range(7*4 * 48):
        action_list.append([1, 0])
    for i in range(7*16 * 48):
        action_list.append([1, 1])
    for i in range(7*8 * 48):
        action_list.append([0, 1])
    actions_list.append(action_list)
    t_max = 7 * 28
    t = np.arange(dt, t_max + dt, dt)
    t_cyto = np.arange(dt_cyto, t_max + dt_cyto, dt_cyto)

    dc = 0.01
    c_min, c_max = 0, 1
    c = np.arange(c_min + dc, c_max + dc, dc)

    t = np.array(t)[:, None]
    t_cyto = np.array(t_cyto)[:, None]
    c = np.array(c)[:, None]

    c_show = np.insert(c, 0, 0) 
    t_show = np.insert(t, 0, 0)
    T_show, C_show = np.meshgrid(t_show, c_show)

    initial_TC = 0.58
    initial_CSF1RI = 0
    initial_IGF1RI = 0
    C, T = np.meshgrid(c,t)
    C_star = np.hstack((C.flatten()[:,None], T.flatten()[:,None]))
    predict_actions = actions_list[2]
    
    predict_CSF1RI_cum = 0
    predict_CC_list = []
    predict_dose_c_list = []
    predict_dose_I_list = []
    predict_CSF1RI_list = []
    predict_CSF1RI_list.append(initial_CSF1RI)
    predict_IGF1RI_list = []
    predict_IGF1RI_list.append(initial_IGF1RI)
    predict_TC_list = []
    predict_TC_list.append(initial_TC)

    length_time = len(predict_actions)
    for k in range(length_time):

        predict_dose_c = env.action_c_space[predict_actions[k][0]]   
        predict_dose_I = env.action_I_space[predict_actions[k][1]]
        predict_tt = (k+1) * predict_dt
        c_CSF1RI_0 = predict_CSF1RI_list[-1]
        c_IGF1RI_0 = predict_IGF1RI_list[-1]
        c_T_0 = predict_TC_list[-1]
        rl_FP_TC, c_CSF1RI_cum, c_CSF1RI, c_IGF1RI = surrogate_model.predict_DRL_FP(c, predict_tt, predict_dose_c, predict_CSF1RI_cum, predict_dose_I,
                                                                c_CSF1RI_0, c_IGF1RI_0, c_T_0, predict_dt)
        predict_CSF1RI_cum = c_CSF1RI_cum
        rl_FP_TC_norm = rl_FP_TC/(np.sum(rl_FP_TC) + 1e-7)
        predict_CC_list.append(rl_FP_TC_norm * 1/dc)
        predict_CSF1RI_list.append(c_CSF1RI)
        predict_IGF1RI_list.append(c_IGF1RI)
        predict_TC_list.append(np.sum(rl_FP_TC_norm*c))
        predict_dose_c_list.append(predict_dose_c)
        predict_dose_I_list.append(predict_dose_I)

    predict_CClist_flat = np.concatenate(predict_CC_list).ravel()
    dose_c = np.array(predict_dose_c_list)
    dose_I = np.array(predict_dose_I_list)
    predict_CSF1RI = np.array(predict_CSF1RI_list)
    predict_IGF1RI = np.array(predict_IGF1RI_list)
    predict_TC = np.array(predict_TC_list)

    # plot and save surrogate model result with 
    # the most effective RL treatment in Fig 8F
    U_pred = scipy.interpolate.griddata(C_star, predict_CClist_flat, (C, T), method='cubic')
    U_pred[U_pred < 0] = 0
    U_show = U_pred.T

    fig, ax = plt.subplots(figsize=(9, 5))
    h = ax.pcolormesh(T_show, C_show, U_show, shading='auto', cmap='rainbow')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.10)
    cbar = fig.colorbar(h, cax=cax)
    cbar.ax.tick_params(labelsize=15)
    ax.set_xlabel('$t$', size=20)
    ax.set_ylabel('$C$', size=20)
    ax.set_title('$p(t,C)$', fontsize=20)
    ax.tick_params(labelsize=15)
    t_ticks = np.arange(0, t_max + 1e-7, 7)
    ax.set_xticks(t_ticks)
    ax.tick_params(axis='x', labelsize=6)
    save_path = f'{save_folder}/Pred_TC_most_effective.pdf'
    plt.savefig(save_path, format='pdf')

    # plot and save the most effective RL treatment
    fig = plt.figure(figsize=(5, 12))
    ax = fig.add_subplot(111)
    gs1 = gridspec.GridSpec(2, 1)
    gs1.update(top=0.9, bottom=0.1, left=0.05, right=0.95, wspace=0.3, hspace=0.3)
    
    ax = plt.subplot(gs1[0, 0])
    ax.plot(t_cyto, dose_c, 'b-', linewidth = 2, label = 'CSF1RI dose')  
    ax.plot(t_show, predict_CSF1RI, 'k-', linewidth = 2, label = 'invivo CSF1RI')  
            
    ax.set_xlabel('$t$')
    ax.set_ylabel('$c$')    
    ax.set_title('$CSF1RI dose$', fontsize = 15)
    ax.axis('square')
    ax.set_xlim([0.0,t_max])
    ax.set_ylim([0.0,1])
    aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.set_aspect(aspect_ratio)
    ax.set_xticks(t_ticks)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    plt.legend()

    ax = plt.subplot(gs1[1, 0])
    ax.plot(t_cyto, dose_I, 'g-', linewidth = 2, label = 'IGF1RI dose')  
    ax.plot(t_show, predict_IGF1RI, 'k-', linewidth = 2, label = 'invivo IGF1RI')     
    ax.plot(t_show, predict_TC, 'r-', linewidth = 2, label = 'TC expectation')  
    ax.set_xlabel('$t$')
    ax.set_ylabel('$c$')    
    ax.set_title('$IGF1RI dose$', fontsize = 15)
    ax.axis('square')
    ax.set_xlim([0.0,t_max])
    ax.set_ylim([0.0,1])
    aspect_ratio = 1 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.set_aspect(aspect_ratio)
    ax.set_xticks(t_ticks)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(15)
    plt.legend()
    save_path = f'{save_folder}/dose_most_effective.pdf'
    plt.savefig(save_path, format='pdf')
    print(f'dose_most_effective done')  


if __name__ == "__main__":
    main()
