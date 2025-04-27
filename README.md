# Multiscale mathematical model-informed reinforcement learning (M4RL) 

## Article short title ##

M4RL for dynamic optimization of cancer treatment

## Abstract ##

Dynamic tumor-microenvironment interactions greatly affect growth and drug resistance, highlighting the importance and challenge of developing mathematical models to optimize treatment schedules. Here, we describe a multiscale mathematical model-informed-reinforcement learning (M4RL) framework to simulate dynamic tumor-microenvironment interactions and optimize drug combination scheduling. We first develop a multiscale agent-based model (MSABM) for a critical biological scenario where interactions between tumor-associated macrophages (TAMs) and tumor cells (TCs) underlie immunotherapy resistance in glioblastoma. Next, we learn a surrogate model based on Fokker-Planck equations for the MSABM using a physics-informed neural network approach. We then design a surrogate model-based reinforcement learning method, using an efficient parallel actor-critic algorithm, to predict optimal scheduling of combination treatments. The most effective regimen of dynamic combination of CSF1R inhibitor (targeting TAMs) and IGF1R inhibitor (targeting TCs) is identified and then verified using spatial transcriptomics data. Overall, the M4RL framework introduces a computational approach for characterizing tumor-microenvironment interactions and optimizing dynamic scheduling of drug combinations. 


## Usage instructions ##

1. Folder [**1_MSABM**](https://github.com/aewsion/M4RL/tree/main/1_MSABM) contains supporting codes for MSABM [**tumor_growth_case**](https://github.com/aewsion/M4RL/tree/main/1_MSABM/tumor_growth_case), [**no_treatment_case**](https://github.com/aewsion/M4RL/tree/main/1_MSABM/no_treatment_case) and [**continuous_CSF1R_I_treatment_case**](https://github.com/aewsion/M4RL/tree/main/1_MSABM/continuous_CSF1R_I_treatment_case).<br><br>
To simulate MSABM cases, first, download and save the required folder. Then, open a terminal in the folder's directory and compile the C++ project using the commands `cmake .` and `make`. Finally, execute the compiled binary by running the bash script `sh run_multiple_simulations.sh`. The results of the simulations will be collected in the folder `YOUR_INPUT_TREATMENT_case_output`. <br><br>
An example of the results of **tumor_growth_case** have been shown in [S1_video_tumor_growth_case](https://github.com/aewsion/M4RL/blob/main/supplement_videos/S1_video_tumor_growth_case.mp4). <br>
An example of the results of **no_treatment_case** have been shown in [S2_video_no_treatment_case](https://github.com/aewsion/M4RL/blob/main/supplement_videos/S2_video_no_treatment_case.mp4). <br>
Two examples of the results of **continuous_CSF1R_I_treatment_case** have been shown in [S3_video_responder_in_continuous_CSF1R_I_treatment_case](https://github.com/aewsion/M4RL/blob/main/supplement_videos/S3_video_responder_in_continuous_CSF1R_I_treatment_case.mp4) and [S4_video_non-responder_in_continuous_CSF1R_I_treatment_case](https://github.com/aewsion/M4RL/blob/main/supplement_videos/S4_video_non-responder_in_continuous_CSF1R_I_treatment_case.mp4). <br><br>

2. Folder [**2_surrogate_model**](https://github.com/aewsion/M4RL/tree/main/2_surrogate_model) contains supporting codes for [**surrogate_model_training**](https://github.com/aewsion/M4RL/blob/main/2_surrogate_model/surrogate_model_training.py) and [**surrogate_model_verification**](https://github.com/aewsion/M4RL/blob/main/2_surrogate_model/surrogate_model_verification.py), along with the sub folder [**PINN_data_constraints**](https://github.com/aewsion/M4RL/tree/main/2_surrogate_model/PINN_data_constraints) for surrogate model training and experimental data [(Daniela F. Quail et al., 2016)](https://www.science.org/doi/10.1126/science.aad3018) for surrogate model verification. The subfolder [**experiment_data**](https://github.com/aewsion/M4RL/tree/main/2_surrogate_model/experiment_data) contains the experimental survival results of 'switch treatment' and 'add treatment' (consistent with 1-week interval RL treatment), and the subfolder [**MSABM_data**](https://github.com/aewsion/M4RL/tree/main/2_surrogate_model/MSABM_data) contains the MSABM prediction results of 'switch treatment' and 'add treatment'.<br><br>
Use the command `python3 surrogate_model_training.py` for model training and generating related plots and `python3 surrogate_model_verification.py` for model verification.<br><br>

3. Folder [**3_A3C_RL_with_surrogate_model**](https://github.com/aewsion/M4RL/tree/main/3_A3C_RL_with_surrogate_model) contains supporting codes for [**A3C_RL_training**](https://github.com/aewsion/M4RL/blob/main/3_A3C_RL_with_surrogate_model/A3C_RL_training.py) and [**A3C_RL_evaluation**](https://github.com/aewsion/M4RL/blob/main/3_A3C_RL_with_surrogate_model/A3C_RL_evaluation.py). The subfolder [**surrogate_model**](https://github.com/aewsion/M4RL/tree/main/3_A3C_RL_with_surrogate_model/surrogate_model) , which contains the pre-trained surrogate model (trained from [**2_surrogate_model**](https://github.com/aewsion/M4RL/tree/main/2_surrogate_model)) for asynchronous advantage actor-critic reinforcement learning (A3C_RL). Results of training are shown in subfolder [**A3C_RL_train**](https://github.com/aewsion/M4RL/tree/main/3_A3C_RL_with_surrogate_model/A3C_RL_train). Results of evaluation and log-rank test on optimal RL treatment compared with combination treatment (S4 Figure) are shown in [**A3C_eval_vs_combination**](https://github.com/aewsion/M4RL/tree/main/3_A3C_RL_with_surrogate_model/A3C_eval_vs_combination). <br><br>
Use the command `python3 A3C_RL_training.py` for training and `python3 A3C_RL_evaluation.py` for evaluation and plotting.<br><br>

4. Folder [**4_testing_treatments_with_ST_based_TME**](https://github.com/aewsion/M4RL/tree/main/4_testing_treatments_with_ST_based_TME) contains supporting codes for MSABM simulation with ST data-based TME (saved in [**stGBM_data**](https://github.com/aewsion/M4RL/tree/main/4_testing_treatments_with_ST_based_TME/stGBM_data)).<br><br>
To simulate MSABM cases, we also need to compile the C++ project using the commands `cmake .` and `make`. Finally, execute the compiled binary by running the bash script `sh run_multiple_simulations.sh`. The results of the simulations will be collected in the folder `YOUR_INPUT_TREATMENT_case_output`. <br><br>
To change the treatment thearpy in MSABM, we need to open [diffusible.cpp](https://github.com/aewsion/M4RL/blob/main/4_testing_treatments_with_ST_based_TME/src/diffusibles.cpp) and navigate lines 172-200. Then, select or customize the treatment as needed.
```cpp
    // optimal RL treatment
    if(time < 140*24){
        CSF1R_I[i][j] = CSF1R_I[i][j] + 0.8 * (1.0 - CSF1R_I[i][j]);
    }
    else{
        CSF1R_I[i][j] = CSF1R_I[i][j] + 0.8 * (0 - CSF1R_I[i][j]);
    }
    if(time >= 28*24){
        IGF1R_I[i][j] = IGF1R_I[i][j] + 0.8 * (1.0 - IGF1R_I[i][j]);
    }
    // 'CSF1R_I only' treatment    
    /*
        *CSF1R_I[i][j] = CSF1R_I[i][j] + 0.8 * (0.7 - CSF1R_I[i][j]);
        */

    // 'IGF1R_I only' treatment
    /*
        *if(time >= 28*24){
        *   IGF1R_I[i][j] = IGF1R_I[i][j] + 0.8 * (1.0 - IGF1R_I[i][j]);
        *}
        */

    // 'CSF1R_I & IGF1R_I' treatment
    /*
        *CSF1R_I[i][j] = CSF1R_I[i][j] + 0.8 * (0.7 - CSF1R_I[i][j]);
        *if(time >= 28*24){
        *   IGF1R_I[i][j] = IGF1R_I[i][j] + 0.8 * (1.0 - IGF1R_I[i][j]);
        *}
        */ 
```
5. Folder [**5_stMLnet_analysis_with_ST_based_data(UKF_262)**](https://github.com/aewsion/M4RL/tree/main/5_stMLnet_analysis_with_ST_based_data(UKF_262)) contains supporting codes for stMLnet analysis with ST data (UKF_262).<br><br>
