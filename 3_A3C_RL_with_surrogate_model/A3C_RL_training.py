import torch
from collections import OrderedDict
import numpy as np
import warnings
import time
import multiprocessing as mp
import os

# set initial seed
seed = 20240712

# CUDA support
if torch.cuda.is_available():
    device = torch.device('cuda')
else:
    device = torch.device('cpu')

# set prediction dt 7, 14, 21 or 28
predict_dt = 28
# save folder
save_folder = f'A3C_RL_train_{predict_dt}'

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

        # Fokker-Planck parameters
        # load parameters
        additional_params = torch.load('surrogate_model/additional_params.pth')  
        initial_log_p_cyto = additional_params['log_p_cyto']
        self.log_p_cyto  = torch.nn.Parameter(initial_log_p_cyto).to(device)        
        # deep neural networks
        self.dnn_TC = DNN(FP_layers).to(device)
        # initialize DNN parameters
        para_TC = torch.load('surrogate_model/dnn_TC_model.pth') # trained with 20-40 data-constraint
        #para_TC = torch.load('surrogate_model/dnn_TC_model_20_100.pth') # trained with 20-100 data-constraint
        self.dnn_TC.load_state_dict(para_TC, strict=False)     
        
        # other parameters
        self.lb = torch.tensor(lb, dtype=torch.float32, requires_grad=True).float().to(device)
        self.ub = torch.tensor(ub, dtype=torch.float32, requires_grad=True).float().to(device)
        self.c_max = torch.tensor(1/dc, dtype=torch.float32, requires_grad=True).float().to(device)
        self.dt_cyto = dt_cyto
        self.CSF1RI_cum_max = 200 / self.dt_cyto  

    # DNN for Fokker-Planck equation of tumor cell
    def net_u_TC(self, C, t, dose_c, dose_cum, dose_I): 
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
    def predict_RL_FP(self, c, t, dose_c, CSF1RI_cum, dose_I, c_CSF1RI_0, c_IGF1RI_0, c_T_0, interval):
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

        uTC_pred = self.net_u_TC(c, t_k, c_CSF1RI_k, CSF1RI_cum_k, c_IGF1RI_k)

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
        self.surrogate_model = model 
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
        last_state = self.state.copy() # copy state
        done = False # termination condtion
        self.time += predict_dt # record time
        dose_c = self.action_c_space[action_c] # get action
        dose_I = self.action_I_space[action_I] # get action
        self.CSF1RI_dose_cum += dose_c * predict_dt # record CSF1R_I accummulation
        c_CSF1RI_0 = self.predict_CSF1RI_list[-1] # get current invivo CSF1R_I concentration
        c_IGF1RI_0 = self.predict_IGF1RI_list[-1] # get current invivo IGF1R_I concentration
        c_T_0 = self.predict_TC_list[-1] # get current TC density

        # surrogate model 
        rl_FP_TC, c_CSF1RI_cum, c_CSF1RI, c_IGF1RI = self.surrogate_model.predict_RL_FP(cc, self.time, dose_c, self.predict_CSF1RI_cum, dose_I, 
                                                   c_CSF1RI_0, c_IGF1RI_0, c_T_0, predict_dt)

        # normalization
        rl_FP_TC_norm = rl_FP_TC/(np.sum(rl_FP_TC) + 1e-7)
        
        self.predict_CSF1RI_cum = c_CSF1RI_cum 
        self.predict_CSF1RI_list.append(c_CSF1RI)
        self.predict_IGF1RI_list.append(c_IGF1RI)
        self.predict_TC_list.append(np.sum(rl_FP_TC_norm*cc)) # expectation TC density

        prob_sum = np.sum(rl_FP_TC)
        survival_prob = np.sum(rl_FP_TC[:-3, :]) / (prob_sum + 1e-7) # survival probabiltiy
        death_prob = 1 - survival_prob # death probabiltiy
        cure_prob = np.sum(rl_FP_TC[:1, :]) / (prob_sum + 1e-7) # cure probabiltiy

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
        reward = 0.1 # basal reward
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

# local net
class A3C_RL_local:
    def __init__(self, local_lstm, local_actor_c, local_actor_I, local_critic, beta = 0.05, gamma=0.9999):
        self.lstm = local_lstm
        self.actor_c = local_actor_c
        self.actor_I = local_actor_I
        self.critic = local_critic
        self.beta = beta # entropy control coefficient
        self.gamma = gamma # discounted
        

    def compute_gradients(self, memory, current_epoch):
        self.lstm.train()
        self.actor_c.train()
        self.actor_I.train()
        self.critic.train()

        length_of_memory = len(memory)
        last_states, actions_c, actions_I, rewards, dones, now_states = zip(*memory)

        last_states_tensor = [torch.tensor(state, dtype=torch.float32, requires_grad=True).to(device) for state in last_states]
        now_states_tensor = [torch.tensor(state, dtype=torch.float32, requires_grad=True).to(device) for state in now_states]

        input_last_states_tensor = torch.stack(last_states_tensor).permute(1, 0, 2)
        input_now_states_tensor = torch.stack(now_states_tensor).permute(1, 0, 2)
        
        actions_c = torch.tensor(actions_c, dtype=torch.float32).float().to(device)
        actions_I = torch.tensor(actions_I, dtype=torch.float32).float().to(device)

        discounted_rewards = torch.empty(length_of_memory, dtype=torch.float32).float().requires_grad_(True).to(device)

        for i in reversed(range(length_of_memory)):
            with torch.no_grad():
                if dones[i]:
                    discounted_rewards[i] = rewards[i]
                else:
                    discounted_rewards[i] = rewards[i] + self.gamma * discounted_rewards[i + 1]

        td_targets = torch.empty(length_of_memory, dtype=torch.float32).float().requires_grad_(True).to(device)
        advantages = torch.empty(length_of_memory, dtype=torch.float32).float().requires_grad_(True).to(device)
        discounted_advantages = torch.empty(length_of_memory, dtype=torch.float32).float().requires_grad_(True).to(device)

        current_values_lstm = self.lstm(input_last_states_tensor)
        current_values = self.critic(current_values_lstm)

        with torch.no_grad():
            next_values_lstm = self.lstm(input_now_states_tensor)
            next_values = self.critic(next_values_lstm)

        for i in reversed(range(length_of_memory)):
            with torch.no_grad():
                if dones[i]:
                    td_targets[i] = rewards[i]
                else:
                    td_targets[i] = rewards[i] + self.gamma * next_values[i]

                advantages[i] = td_targets[i] - current_values[i]
                if dones[i]:
                    discounted_advantages[i] = advantages[i]
                else:
                    discounted_advantages[i] = advantages[i] + self.gamma * discounted_advantages[i + 1]

        critic_loss = torch.nn.functional.mse_loss(discounted_rewards, current_values, reduction='sum')

        policy_dists_c = self.actor_c(next_values_lstm)
        policy_dists_I = self.actor_I(next_values_lstm)

        dists_c = torch.distributions.Categorical(policy_dists_c)
        dists_I = torch.distributions.Categorical(policy_dists_I)
        entropy_c = dists_c.entropy().sum()
        entropy_I = dists_I.entropy().sum()
        action_c_log_probs = dists_c.log_prob(actions_c)
        action_I_log_probs = dists_I.log_prob(actions_I)

        actor_c_loss = -torch.sum(action_c_log_probs * discounted_advantages) - self.beta * entropy_c
        actor_I_loss = -torch.sum(action_I_log_probs * discounted_advantages) - self.beta * entropy_I

        loss = critic_loss + actor_c_loss + actor_I_loss
        loss.backward()

        # print
        '''
        if current_epoch % 100 == 99:
            print(
                'global epoches %d, critic_loss: %.5e' 
                % (current_epoch, critic_loss.item())
            )
            print(
                'actor_c_loss: %.5e, actor_I_loss: %.5e, entropy_c: %.5e, entropy_I: %.5e' 
                % (actor_c_loss.item(), actor_I_loss.item(), entropy_c.item(), entropy_I.item())
            )
            print("input_last_states_tensor", input_last_states_tensor.detach().cpu().numpy())
            print("current_values_lstm", current_values_lstm.detach().cpu().numpy())
            print("current_values", current_values.detach().cpu().numpy())
            print("rewards", [f"{reward:.2f}" for reward in rewards])
            print("discounted_rewards", discounted_rewards)
            print("action_c_probs", torch.exp(action_c_log_probs))
            print("action_I_probs", torch.exp(action_I_log_probs))
        '''

# run multiple local A3C nets
def A3C_RL_train_worker(surrogate_model, env_initial_conditions, state_size, cc, lstm_hidden_dim, num_layers, predict_time, predict_dt, 
                        pre_therapy_time, lamda, global_lstm, global_actor_c, global_actor_I, global_critic, critic_optimizer,
                        actor_optimizer, roll_out, global_epoch_num, max_epochs, worker_seed):
    # set seed
    np.random.seed(worker_seed)
    torch.manual_seed(worker_seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(worker_seed)
        torch.cuda.manual_seed_all(worker_seed)

    # create local RL environment
    env = Env(surrogate_model, env_initial_conditions)

    str_lamda = str(lamda).replace('.', '_')

    # get number of actions
    action_c_size = env.action_c_space.shape[0]
    action_I_size = env.action_I_space.shape[0]

    # initialize local nets using global net parameters
    local_lstm = LSTMNet(state_size, lstm_hidden_dim, num_layers).to(device)
    local_lstm.load_state_dict(global_lstm.state_dict())

    local_actor_c = PolicyNet(lstm_hidden_dim, action_c_size).to(device)
    local_actor_c.load_state_dict(global_actor_c.state_dict())

    local_actor_I = PolicyNet(lstm_hidden_dim, action_I_size).to(device)
    local_actor_I.load_state_dict(global_actor_I.state_dict())

    local_critic = ValueNet(lstm_hidden_dim).to(device)
    local_critic.load_state_dict(global_critic.state_dict())

    # initialize local A3C net
    local_model = A3C_RL_local(local_lstm, local_actor_c, local_actor_I, local_critic)

    # set epsilon greedy (if necessary)
    epsilon = 0.0    
    epsilon_decay = 0.95

    while True:

        with global_epoch_num.get_lock():
            current_epoch = global_epoch_num.value
            # save global A3C net
            if current_epoch % 1000 == 0:
                torch.save(global_critic.state_dict(), f'{save_folder}/global_critic_{str_lamda}.pth')
                print("global_critic saved")
                torch.save(global_actor_c.state_dict(), f'{save_folder}/global_actor_c_{str_lamda}.pth')
                print("global_actor_c saved")
                torch.save(global_actor_I.state_dict(), f'{save_folder}/global_actor_I_{str_lamda}.pth')
                print("global_actor_I saved")
                torch.save(global_lstm.state_dict(), f'{save_folder}/global_LSTM_{str_lamda}.pth')
                print("global_LSTM saved")
            if current_epoch >= max_epochs:
                torch.save(global_critic.state_dict(), f'{save_folder}/global_critic_{str_lamda}.pth')
                print("global_critic saved")
                torch.save(global_actor_c.state_dict(), f'{save_folder}/global_actor_c_{str_lamda}.pth')
                print("global_actor_c saved")
                torch.save(global_actor_I.state_dict(), f'{save_folder}/global_actor_I_{str_lamda}.pth')
                print("global_actor_I saved")
                torch.save(global_lstm.state_dict(), f'{save_folder}/global_LSTM_{str_lamda}.pth')
                print("global_LSTM saved")
                break
            # count global epoches
            global_epoch_num.value += 1

        # initialization
        memory = []
        env.reset()
        actions_taken = []
        done = False
        epi = 0 

        while not done:
            
            if epi == 0: # preprocess treatment
                action_c = action_c_size - 1
                action_I = 0
                last_state, reward, done, now_state = env.step(action_c, action_I, epi, pre_therapy_time, cc, lamda)
                actions_taken.append([action_c, action_I]) 
                epi += 1
            else:
                rollout_memory = []  # save rollout memory
                steps = 0

                # get actions
                local_actor_c.eval()
                local_actor_I.eval()
                state_tensor = torch.tensor(now_state, dtype=torch.float32, requires_grad=True).float().to(device).unsqueeze(0)
                state_tensor_lstm = local_lstm(state_tensor)

                policy_dist_c = local_actor_c(state_tensor_lstm[-1:,:])
                dist_c = torch.distributions.Categorical(policy_dist_c)

                policy_dist_I = local_actor_I(state_tensor_lstm[-1:,:])
                dist_I = torch.distributions.Categorical(policy_dist_I)
                
                action_c = dist_c.sample().item()
                action_I = dist_I.sample().item()

                # set rollout (if necessary) default = 1
                while steps < roll_out and not done:
                    # run each step
                    last_state, reward, done, now_state = env.step(action_c, action_I, epi, predict_dt, cc, lamda)
                    # record actions
                    actions_taken.append([action_c, action_I])
                    if epi >= ((predict_time - pre_therapy_time)/predict_dt):
                        done = True                   
                    epi += 1 # count RL environment steps
                    steps += 1  # count rollout steps  
                    rollout_memory.append((last_state[-1:], action_c, action_I, reward, done, now_state[-1:]))
                memory.extend(rollout_memory)      

        if len(memory) > 0:
            # Reset gradients to zero
            critic_optimizer.zero_grad()
            actor_optimizer.zero_grad()

            # compute gradients on local A3C net
            local_model.compute_gradients(memory, current_epoch)

            # apply the local nets gradients to the global nets
            for local_param, global_param in zip(local_lstm.parameters(), global_lstm.parameters()):
                global_param._grad = local_param.grad

            for local_param, global_param in zip(local_actor_c.parameters(), global_actor_c.parameters()):
                global_param._grad = local_param.grad

            for local_param, global_param in zip(local_actor_I.parameters(), global_actor_I.parameters()):
                global_param._grad = local_param.grad

            for local_param, global_param in zip(local_critic.parameters(), global_critic.parameters()):
                global_param._grad = local_param.grad

            # clip
            critic_global_params = list(global_lstm.parameters()) + list(global_critic.parameters())
            actor_global_params = list(global_actor_c.parameters()) + list(global_actor_I.parameters())                     
            torch.nn.utils.clip_grad_norm_(critic_global_params, max_norm=10.0)
            torch.nn.utils.clip_grad_norm_(actor_global_params, max_norm=10.0)

            # update global optimzers
            critic_optimizer.step()
            actor_optimizer.step()

            # Synchronize local nets parameters
            local_lstm.load_state_dict(global_lstm.state_dict())
            local_actor_c.load_state_dict(global_actor_c.state_dict())
            local_actor_I.load_state_dict(global_actor_I.state_dict())
            local_critic.load_state_dict(global_critic.state_dict())

        print(f"Epoch: {current_epoch} Actions taken: {actions_taken}")
        
        # epsilon greedy
        if epsilon > 0.1:
            epsilon *= epsilon_decay

def main():
    lamda = 0.005 # control coefficient in reward funtion
    str_lamda = str(lamda).replace('.', '_')

    FP_layers = [5, 40, 120, 250, 500, 1000, 1000, 600, 300, 150, 1]

    dt = 7 # 7 days
    dt_cyto = 1/48 # 1/48 day = 0.5 hour
    t_max = 28*7 # total simulation time 49 weeks
    t = np.arange(dt, t_max + dt, dt)

    dc = 0.05 # normalized density
    c_min, c_max = 0, 1 # normalized density
    c = np.arange(c_min + dc, c_max + dc, dc)

    t = np.array(t)[:, None]
    c = np.array(c)[:, None]

    # initial conditions
    initial_TC = 0.58
    initial_CSF1RI = 0
    initial_IGF1RI = 0
    env_initial_conditions = [initial_CSF1RI, initial_IGF1RI, initial_TC]

    # lower bound and upper bound
    lb = [0.05, 5]
    ub = [1.0, 200]  

    state_size = 2 # RL environment state size
    lstm_hidden_dim = state_size * 2 # LSTM hidden layer size
    num_layers = 1 # LSTM number of layer
    
    predict_time = t_max # prediction time
    pre_therapy_time = 28 # preprocess therapy time
    roll_out = 1 # rollout default = 1
    
    surrogate_model = RL_dNN(FP_layers, lb, ub, dc, dt_cyto) # initialize surrogate model
    env = Env(surrogate_model, env_initial_conditions) # create RL environment

    # action space size
    action_c_size = env.action_c_space.shape[0]
    action_I_size = env.action_I_space.shape[0]

    # create folder for saving
    figures_folder = f'{save_folder}/{str_lamda}'
    if not os.path.exists(figures_folder):
        os.makedirs(figures_folder)

    # record start time
    start_time = time.time()

    # initialize global A3C net
    global_lstm = LSTMNet(state_size, lstm_hidden_dim, num_layers).to(device)
    para_global_lstm = torch.load(f'{save_folder}/global_LSTM_{str_lamda}.pth')      
    global_lstm.load_state_dict(para_global_lstm, strict=False)

    global_actor_c = PolicyNet(lstm_hidden_dim, action_c_size).to(device)
    para_global_actor_c = torch.load(f'{save_folder}/global_actor_c_{str_lamda}.pth')      
    global_actor_c.load_state_dict(para_global_actor_c, strict=False)

    global_actor_I = PolicyNet(lstm_hidden_dim, action_I_size).to(device)
    para_global_actor_I = torch.load(f'{save_folder}/global_actor_I_{str_lamda}.pth')      
    global_actor_I.load_state_dict(para_global_actor_I, strict=False)

    global_critic = ValueNet(lstm_hidden_dim).to(device)
    para_global_critic = torch.load(f'{save_folder}/global_critic_{str_lamda}.pth')      
    global_critic.load_state_dict(para_global_critic, strict=False)
    
    # share parameters
    global_lstm.share_memory()
    global_actor_c.share_memory()
    global_actor_I.share_memory()
    global_critic.share_memory()

    # trainable parameters
    critic_parameters = list(global_lstm.parameters()) + list(global_critic.parameters())
    actor_parameters = list(global_actor_c.parameters()) + list(global_actor_I.parameters())

    # set optimizers
    critic_optimizer = torch.optim.Adam(critic_parameters, lr=1e-4, weight_decay=1e-5)
    actor_optimizer = torch.optim.Adam(actor_parameters, lr=1e-4, weight_decay=1e-5)

    # count global epoches
    global_epoch_num = mp.Value('i', 0)
    max_epochs = 20000  # maximum global training epoches

    # start local workers
    workers = []
    #for i in range(mp.cpu_count()):
    for i in range(30):
        worker_seed = seed + i
        worker = mp.Process(target = A3C_RL_train_worker, args = (surrogate_model, env_initial_conditions, state_size, c, lstm_hidden_dim, num_layers, 
                                    predict_time, predict_dt, pre_therapy_time, lamda, global_lstm, global_actor_c, global_actor_I, global_critic, 
                                    critic_optimizer, actor_optimizer, roll_out, global_epoch_num, max_epochs, worker_seed))
        worker.start()
        workers.append(worker)
        time.sleep(0.2)

    for worker in workers:
        worker.join()

    # record end time
    end_time = time.time()

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")        


if __name__ == "__main__":
    mp.set_start_method('spawn')
    main()

'''
    REFERENCE

    K. Gallagher, M. A. R. Strobl, D. S. Park, F. C. Spoendlin, R. A. Gatenby, P. K. Maini, 
    A. R. A. Anderson, Mathematical Model-Driven Deep Learning Enables Personalized Adaptive 
    Therapy. Cancer Res 84, 1929-1941 (2024).

    V. Mnih, A. P. Badia, M. Mirza, A. Graves, T. Harley, T. P. Lillicrap, D. Silver, 
    K. Kavukcuoglu, Asynchronous Methods for Deep Reinforcement Learning, paper presented 
    at the 33rd International Conference on International Conference on Machine Learning, 
    New York, NY, 11 June 2016.
'''
