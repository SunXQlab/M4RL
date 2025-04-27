#include <iostream>

#include "diffusibles.h"

Diffusibles::Diffusibles(double DX) {
    // representation, unit, references of parameters
    // are listed in S2 Table Cytokines part & Drugs part
    dx = DX;
    dt = 2;
    dt_drug = 2;
    time = 0;

    D_CSF1 = 3e-7; 
    D_IGF1 = 3e-7;
    D_EGF = 3e-7;

    D_CSF1R_I = 2e-10;
    D_IGF1R_I = 2e-10; 
    D_EGFR_I = 2e-10; // if necessary

    q_CSF1R_I = 0.8; 
    q_IGF1R_I = 0.8; 
    q_EGFR_I = 0.8; // if necessary

    S_CSF1 = 2.5e-14; // 2.5e-26 mol per second per cell = 2.5e-14 mol/L/s
    S_IGF1 = 5e-15; // 2.5e-26 mol per second per cell; M2 cell diameter = 20um
    S_EGF = 5e-15; // same as S_IGF1

    d_CSF1 = 5e-4; 
    d_IGF1 = 2.5e-4; 
    d_EGF = 2.5e-4; 

    u_CSF1R_I = 2e-7;
    u_IGF1R_I = 2e-7;
    u_EGFR_I = 2e-7;

    d_CSF1R_I = 2e-7; 
    d_IGF1R_I = 2e-7;
    d_EGFR_I = 2e-7;

    drug_time = 0;
    drug_feed_time = (205)*24*60*60; // no_treatment_case here
    // set drug_feed_time = 0-200 days for continuous treatment
    drug_suspend_time = (205)*24*60*60;

    S_A = 1.2905e-7*3;
    A_0 = 0.0;

    max_IGF1 = 1.2e-13;
    max_EGF = 6e-13;
    max_CSF1 = 4e-12;

    if_CSF1R_I = 0;
    if_IGF1R_I = 0;
    if_EGFR_I = 0;   

    for(int i=0; i<100; ++i){
        for(int j=0; j<100; ++j){
            CSF1[i][j] = 0;
            EGF[i][j] = 0;
            IGF1[i][j] = 0;
            I[i][j] = 0;

            CSF1R_I[i][j] = 0;
            IGF1R_I[i][j] = 0;
            EGFR_I[i][j] = 0;
            
            AKT[i][j] = 0;
            ERK[i][j] = 0;
            IGF1R[i][j] = 0;
            EGFR[i][j] = 0;

            prob[i][j] = 0;
            prob_0[i][j] = 0;
        }
    }
}

// diffusion
void Diffusibles::diffusion(CellGrids &cg, double tstep) {

    // set zero, only for recording
    // no calculation here
    for(int i=0; i<100; ++i){
        for(int j=0; j<100; ++j){
            AKT[i][j] = 0;
            ERK[i][j] = 0;
            IGF1R[i][j] = 0;
            EGFR[i][j] = 0;
            prob[i][j] = 0;
            prob_0[i][j] = 0;
        }
    }

    // if feeding CSF1R_I
    double S_a = 0;
    if(time * 60 * 60 >= drug_feed_time && if_CSF1R_I == 1){
        drug_time = drug_time + tstep * 60 * 60;
        S_a = S_A;
    }
    //std::cout<<"drug_time = "<<drug_time<<std::endl;

    // break condition
    double maxDif;
    double c;

    // steps of forward Eular method
    double t = 60 * 60 * tstep / dt; 
    double t_drug = 60 * 60 * tstep / dt_drug;

    // finite difference & forward Eular method
    for (int q = 0; q < t; q++) {
        maxDif = 0;
        for (int i = 1; i < 99; i++) {
            for (int j = 1; j < 99; j++) {
                double CSF1_0 = CSF1[i][j];
                double EGF_0 = EGF[i][j];       

                CSF1[i][j] = CSF1[i][j] + dt * ((S_CSF1*cg.tc[i][j]) * (1 - CSF1[i][j]/max_CSF1)) - dt * d_CSF1 * CSF1[i][j]
                        + (dt * D_CSF1 / (dx * dx)) * (CSF1[i + 1][j] + CSF1[i - 1][j]+ CSF1[i][j + 1] + CSF1[i][j - 1]- 4 * CSF1[i][j]);
                                                        
                EGF[i][j] = EGF[i][j] + dt * ((S_EGF*cg.m2[i][j]) * (1 - EGF[i][j]/max_EGF)) - dt * d_EGF * EGF[i][j]
                        + (dt * D_EGF / (dx * dx)) * (EGF[i + 1][j] + EGF[i - 1][j] + EGF[i][j + 1] + EGF[i][j - 1] - 4 * EGF[i][j]);
         
                if(CSF1_0 > 0){
                    c = (CSF1[i][j] - CSF1_0)/CSF1_0;
                    if(c > maxDif){maxDif=c;}
                }
                if(EGF_0 > 0){
                    c = (EGF[i][j] - EGF_0)/EGF_0;
                    if(c > maxDif){maxDif=c;}
                }
                // check positivity
                if(EGF[i][j] < 0 || CSF1[i][j] < 0){
                    std::cout << "EGF " << EGF[i][j] << std::endl;
                    std::cout << "CSF1 " << CSF1[i][j] << std::endl;
                    throw std::runtime_error("Diffusibles::diffusion error 1");
                }
            }            
        }

        // zero flux boundary condition
        for(int i = 1; i < 99; i++){
            CSF1[i][0] = CSF1[i][2]; 
            CSF1[i][99] = CSF1[i][97];
            EGF[i][0] = EGF[i][2];
            EGF[i][99] = EGF[i][97];
        }
        for(int j = 1; j < 99; j++){
            CSF1[0][j] = CSF1[2][j];
            CSF1[99][j] = CSF1[97][j];
            EGF[0][j] = EGF[2][j];
            EGF[99][j] = EGF[97][j];
        }

        // break condition
        if(maxDif < 0.00001 && q > 5){
            std::cout<<"diff break q1 = "<<q<<std::endl;
            break;}      
    }

    // finite difference & forward Eular method
    for (int q = 0; q < t_drug; q++) {
        // feed drug 
        if(time*60*60 >= drug_feed_time && time*60*60 < drug_suspend_time){
            if_CSF1R_I = 1;
            for(int i=1; i<99; ++i){
                for(int j=1; j<99; ++j){
                    if(cg.vas[i][j] == 1){
                        CSF1R_I[i][j] = CSF1R_I[i][j] + q_CSF1R_I * (0.0 - CSF1R_I[i][j]);
                    }      
                }
            }
        }

        for (int i = 1; i < 99; i++) {
            for (int j = 1; j < 99; j++) {
                double IGF1_0 = IGF1[i][j];
                double CSF1R_I_0 = CSF1R_I[i][j];
                double IGF1R_I_0 = IGF1R_I[i][j];
                double EGFR_I_0 = EGFR_I[i][j];
                for(int k=0; k<dt_drug/dt; k++){
                    IGF1[i][j] = IGF1[i][j] + dt * ((S_IGF1*cg.m2[i][j]) * (A_0 + I[i][j])
                            * (1 - IGF1[i][j]/max_IGF1)) - dt * d_IGF1 * IGF1[i][j]
                            + (dt * D_IGF1 / (dx * dx)) * (IGF1[i + 1][j] + IGF1[i - 1][j]+ IGF1[i][j + 1] + IGF1[i][j - 1]- 4 * IGF1[i][j]);        
                }

                // integral parts in IGF1 concentration caculation
                I[i][j] = I[i][j] + dt_drug * S_a * CSF1R_I[i][j];
         
                CSF1R_I[i][j] = CSF1R_I[i][j] + ((dt_drug) * D_CSF1R_I / (dx * dx)) 
                         * (CSF1R_I[i + 1][j] + CSF1R_I[i - 1][j]+ CSF1R_I[i][j + 1] + CSF1R_I[i][j - 1]- 4 * CSF1R_I[i][j])
                         - (dt_drug) * (d_CSF1R_I + (cg.m0[i][j]+cg.m1[i][j]+cg.m2[i][j]) * u_CSF1R_I) * CSF1R_I[i][j];

                IGF1R_I[i][j] = IGF1R_I[i][j] + ((dt_drug) * D_IGF1R_I / (dx * dx))
                         * (IGF1R_I[i + 1][j] + IGF1R_I[i - 1][j]+ IGF1R_I[i][j + 1] + IGF1R_I[i][j - 1]- 4 * IGF1R_I[i][j])
                         - (dt_drug) * (d_IGF1R_I + cg.tc[i][j] * u_IGF1R_I) * IGF1R_I[i][j];

                EGFR_I[i][j] = EGFR_I[i][j] + ((dt_drug) * D_EGFR_I / (dx * dx))
                         * (EGFR_I[i + 1][j] + EGFR_I[i - 1][j]+ EGFR_I[i][j + 1] + EGFR_I[i][j - 1]- 4 * EGFR_I[i][j])
                         - (dt_drug) * (d_EGFR_I + cg.tc[i][j] * u_EGFR_I) * EGFR_I[i][j];
                
                if(IGF1[i][j] < 0 || CSF1R_I[i][j] < 0 || IGF1R_I[i][j] < 0 || EGFR_I[i][j] < 0){
                    std::cout << "IGF1 " << IGF1[i][j] << std::endl;
                    std::cout << "CSF1R_I "<< CSF1R_I[i][j] << std::endl;
                    std::cout << "CSF1R_I "<< CSF1R_I[i][j] << std::endl;
                    std::cout << "IGF1R_I "<< IGF1R_I[i][j] << std::endl;
                    std::cout << "EGFR_I "<< EGFR_I[i][j] << std::endl;
                    throw std::runtime_error("Diffusibles::diffusion error 2");
                }
            }
        }

        // zero flux boundary condition
        for(int i = 1; i < 99; i++){
            IGF1[i][0] = IGF1[i][2];
            IGF1[i][99] = IGF1[i][97];
            CSF1R_I[i][0] = CSF1R_I[i][2];
            CSF1R_I[i][99] = CSF1R_I[i][97];
            IGF1R_I[i][0] = IGF1R_I[i][2];
            IGF1R_I[i][99] = IGF1R_I[i][97];
            EGFR_I[i][0] = EGFR_I[i][2];
            EGFR_I[i][99] = EGFR_I[i][97];
        }
        for(int j = 1; j < 99; j++){
            IGF1[0][j] = IGF1[2][j];
            IGF1[99][j] = IGF1[97][j];
            CSF1R_I[0][j] = CSF1R_I[2][j];
            CSF1R_I[99][j] = CSF1R_I[97][j];
            IGF1R_I[0][j] = IGF1R_I[2][j];
            IGF1R_I[99][j] = IGF1R_I[97][j];
            EGFR_I[0][j] = EGFR_I[2][j];
            EGFR_I[99][j] = EGFR_I[97][j];
        }    
    }

    // record simulation time
    time = time + tstep;          
}



            
            
    

