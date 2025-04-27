#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "environment.h"

std::random_device rd;

Environment::Environment(double stepSize, std::string folder, int set, double M0RecProb, double TCProl, double k_M12,
                        double kd_M12, double k_I, double p_tc, double p_d, int pha, int rand,
                        double a_c, double a_i, double a_erk){   
    // directory to save time courses
    saveDir ="./"+folder+"/set_"+std::to_string(set); 

    dx = 0.0015; // 15um, diameter of cell grid

    p_TC = p_tc;
    p_D = p_d;
    Pha = pha;
    tstep = stepSize; // 0.5 hour

    // key parameters
    a_C = a_c;
    a_I = a_i;
    a_ERK = a_erk;
   
    M0recProb = M0RecProb; 
    TCProlTime = TCProl;

    // Maximum environment capacity number
    TC_max = 2154;
    M_max = 450;
    DC_max = 1000;

    // Michaelis constants of TAMs
    K_M12 = k_M12;
    Kd_M12 = kd_M12;
    K_I = k_I;

    // initial value of signal pathways in tumor cell
    EGFRi = 0.906171095;
    IGF1Ri = 0.027626518;
    AKTi = 0.174679186;
    ERKi = 0.496913113;

    // an array for identifying tumor cell state
    for(int i=0; i<100; i++){
        for(int j=0; j<100; j++){
            tcs[i][j] = 0;
        }
    }

    // end time
    endTime = 0;
}

// plot spatial states of the simulation
void Environment::plot(CellGrids &cg, Diffusibles &diff, int s){

    // identify tumor cell state
    for(int i=0; i<100; i++){
        for(int j=0; j<100; j++){
            tcs[i][j] = 0;
        }
    }
    for(auto & cell : tc_list){
        int i = cell.location[0];
        int j = cell.location[1];
        if(cell.state == "alive"){
            tcs[i][j] = 1*cg.tc[i][j];
        }
        else if(cell.state == "quiescent"){
            tcs[i][j] = 7*cg.tc[i][j];
        }
        else{
            std::cout<<"Environment tumor cell state = "<<cell.state<<std::endl;
        }
    }

    // directory to save spatial information
    std::string str = "mkdir -p "+saveDir+"/spatial/"+ std::to_string(s);
    std::system(str.c_str());

    // save cell grid
    std::ofstream myfile; 
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/cg.csv");
    for(int i=0; i<100; ++i){
        myfile << tcs[i][0] + 2*cg.m0[i][0] + 3*cg.m1[i][0] + 4*cg.m2[i][0] 
            + 5*cg.vas[i][0] + 6*cg.dead[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << tcs[i][j] + 2*cg.m0[i][j] + 3*cg.m1[i][j] 
            + 4*cg.m2[i][j] + 5*cg.vas[i][j]  + 6*cg.dead[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save CSF1
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/CSF1.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.CSF1[i][0]/diff.max_CSF1;
        for(int j=1; j<100; ++j){
            myfile << "," << diff.CSF1[i][j]/diff.max_CSF1;
        }
        myfile << std::endl;
    }
    myfile.close();

    // save IGF1
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/IGF1.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.IGF1[i][0]/diff.max_IGF1;
        for(int j=1; j<100; ++j){
            myfile << "," << diff.IGF1[i][j]/diff.max_IGF1;
        }
        myfile << std::endl;
    }
    myfile.close();

    // save EGF
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/EGF.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.EGF[i][0]/diff.max_EGF;
        for(int j=1; j<100; ++j){
            myfile << "," << diff.EGF[i][j]/diff.max_EGF;
        }
        myfile << std::endl;
    }
    myfile.close();

    // save AKT
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/AKT.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.AKT[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.AKT[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save ERK
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/ERK.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.ERK[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.ERK[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save IGF1R
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/IGF1R.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.IGF1R[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.IGF1R[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save EGFR
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/EGFR.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.EGFR[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.EGFR[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save CSF1R_I
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/CSF1R_I.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.CSF1R_I[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.CSF1R_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save IGF1R_I
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/IGF1R_I.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.IGF1R_I[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.IGF1R_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save EGFR_I (if necessary)
    /*
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/EGFR_I.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.EGFR_I[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.EGFR_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    */

    // save tumor cell proliferation probability (if necessary)
    /*
    // proliferation probability restricted by maximum capacity
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/prob.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.prob[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.prob[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // proliferation probability without maximum capacity restriction
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/prob0.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.prob_0[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.prob_0[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    */

    // save vascular cells
    if(s == 0){
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/vas.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.vas[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.vas[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    }
    
    // reset zero
    for(int i=0; i<100; i++){
        for(int j=0; j<100; j++){
            diff.AKT[i][j] = 0;
            diff.ERK[i][j] = 0;
            diff.IGF1R[i][j] = 0;
            diff.EGFR[i][j] = 0;
            diff.prob[i][j] = 0;
            diff.prob_0[i][j] = 0;
        }
    }
    
}

// save time courses of the simulation
void Environment::save(Diffusibles &diff){
    std::fstream myfile;

    // save tumor cell number sequence
    myfile.open(saveDir+"/tc_num.csv", std::ios_base::app);
    myfile << tc_num[0];
    for(int i=1; i<tc_num.size(); i++){
        myfile << "," << tc_num[i];
    }
    myfile << std::endl;
    myfile.close();

    // calculate maximum number of cells
    int maxTumor = 0;
    int maxM0 = 0;
    int maxM1 = 0;
    int maxM2 = 0;
    for(int i=0; i<tc_num.size(); ++i){
        if(tc_num[i] > maxTumor){maxTumor = tc_num[i];}
        if(m0_num[i] > maxM0){maxM0 = m0_num[i];}
        if(m1_num[i] > maxM1){maxM1 = m1_num[i];}
        if(m2_num[i] > maxM2){maxM2 = m2_num[i];} 
    }

    // save end time information & maximun number of cells
    myfile.open(saveDir+"/endtime_information.csv", std::ios_base::app);
    myfile << times[times.size()-1] << ",";
    myfile << tc_num[tc_num.size()-1] << ",";
    myfile << maxTumor << ",";
    myfile << maxM1 << ",";
    myfile << maxM2 << ",";
    myfile << maxM0 << ",";
    myfile << std::endl;
    myfile.close();

    // save M0 macrophages number sequence
    myfile.open(saveDir+"/m0_num.csv", std::ios_base::app);
    myfile << m0_num[0];
    for(int i=1; i<m0_num.size(); i++){
        myfile << "," << m0_num[i];
    }
    myfile << std::endl;
    myfile.close();

    // save M1 macrophages number sequence
    myfile.open(saveDir+"/m1_num.csv", std::ios_base::app);
    myfile << m1_num[0];
    for(int i=1; i<m1_num.size(); i++){
        myfile << "," << m1_num[i];
    }
    myfile << std::endl;
    myfile.close();

    // save M2 macrophages number sequence
    myfile.open(saveDir+"/m2_num.csv", std::ios_base::app);
    myfile << m2_num[0];
    for(int i=1; i<m2_num.size(); i++){
        myfile << "," << m2_num[i];
    }
    myfile << std::endl;
    myfile.close();
    
    // save TAMs total number sequence
    myfile.open(saveDir+"/m_num.csv", std::ios_base::app);
    myfile << m0_num[0]+m1_num[0]+m2_num[0];
    for(int i=1; i<m0_num.size(); i++){
        int m_num = m0_num[i]+m1_num[i]+m2_num[i];
        myfile << "," << m_num;
    }
    myfile << std::endl;
    myfile.close();

    // save M1 macrophages phagocytosis number sequence
    myfile.open(saveDir+"/phago_num.csv", std::ios_base::app);
    myfile << phago_num[0];
    for(int i=1; i<phago_num.size(); i++){
        if(i%2 == 0){
            myfile << "," << phago_num[i];
        }
    }
    myfile << std::endl;
    myfile.close();

    // save average concentration sequence of IGF1
    myfile.open(saveDir+"/avgIGF1.csv", std::ios_base::app);
    myfile << IGF1_avg[0]/diff.max_IGF1;
    for(int i=1; i<IGF1_avg.size(); ++i){
	myfile << "," << IGF1_avg[i]/diff.max_IGF1;
    }
    myfile << std::endl;
    myfile.close();

    // save maximum concentration sequence of IGF1
    myfile.open(saveDir+"/maxIGF1.csv", std::ios_base::app);
    myfile << IGF1_max[0]/diff.max_IGF1;
    for(int i=1; i<IGF1_max.size(); ++i){
	myfile << "," << IGF1_max[i]/diff.max_IGF1;
    }
    myfile << std::endl;
    myfile.close();

    // save average concentration sequence of CSF1
    myfile.open(saveDir+"/avgCSF1.csv", std::ios_base::app);
    myfile << CSF1_avg[0]/diff.max_CSF1;
    for(int i=1; i<CSF1_avg.size(); ++i){
	myfile << "," << CSF1_avg[i]/diff.max_CSF1;
    }
    myfile << std::endl;
    myfile.close();

    // save maximum concentration sequence of CSF1
    myfile.open(saveDir+"/maxCSF1.csv", std::ios_base::app);
    myfile << CSF1_max[0]/diff.max_CSF1;
    for(int i=1; i<CSF1_max.size(); ++i){
	myfile << "," << CSF1_max[i]/diff.max_CSF1;
    }
    myfile << std::endl;
    myfile.close();

    // save average concentration sequence of EGF
    myfile.open(saveDir+"/avgEGF.csv", std::ios_base::app);
    myfile << EGF_avg[0]/diff.max_EGF;
    for(int i=1; i<EGF_avg.size(); ++i){
	myfile << "," << EGF_avg[i]/diff.max_EGF;
    }
    myfile << std::endl;
    myfile.close();

    // save maximum concentration sequence of EGF
    myfile.open(saveDir+"/maxEGF.csv", std::ios_base::app);
    myfile << EGF_max[0]/diff.max_EGF;
    for(int i=1; i<EGF_max.size(); ++i){
	myfile << "," << EGF_max[i]/diff.max_EGF;
    }
    myfile << std::endl;
    myfile.close();

    // save average concentration sequence of CSF1R_I
    myfile.open(saveDir+"/avgCSF1RI.csv", std::ios_base::app);
    myfile << CSF1R_I_avg[0];
    for(int i=1; i<CSF1R_I_avg.size(); ++i){
	myfile << "," << CSF1R_I_avg[i];
    }
    myfile << std::endl;
    myfile.close();

    // save maximum concentration sequence of CSF1R_I
    myfile.open(saveDir+"/maxCSF1RI.csv", std::ios_base::app);
    myfile << CSF1R_I_max[0];
    for(int i=1; i<CSF1R_I_max.size(); ++i){
	myfile << "," << CSF1R_I_max[i];
    }
    myfile << std::endl;
    myfile.close();

    // save average concentration sequence of AKT
    myfile.open(saveDir+"/avgAKT.csv", std::ios_base::app);
    myfile << AKT_avg[0];
    for(int i=1; i<AKT_avg.size(); ++i){
	myfile << "," << AKT_avg[i];
    }
    myfile << std::endl;
    myfile.close();

    // save maximum concentration sequence of AKT
    myfile.open(saveDir+"/maxAKT.csv", std::ios_base::app);
    myfile << AKT_max[0];
    for(int i=1; i<AKT_max.size(); ++i){
	myfile << "," << AKT_max[i];
    }
    myfile << std::endl;
    myfile.close();

    // save average concentration sequence of ERK
    myfile.open(saveDir+"/avgERK.csv", std::ios_base::app);
    myfile << ERK_avg[0];
    for(int i=1; i<ERK_avg.size(); ++i){
	myfile << "," << ERK_avg[i];
    }
    myfile << std::endl;
    myfile.close();

    // save maximum concentration sequence of ERK
    myfile.open(saveDir+"/maxERK.csv", std::ios_base::app);
    myfile << ERK_max[0];
    for(int i=1; i<ERK_max.size(); ++i){
	myfile << "," << ERK_max[i];
    }
    myfile << std::endl;
    myfile.close();

    // save average concentration sequence of ERK
    myfile.open(saveDir+"/avgIGF1R.csv", std::ios_base::app);
    myfile << IGF1R_avg[0];
    for(int i=1; i<IGF1R_avg.size(); ++i){
	myfile << "," << IGF1R_avg[i];
    }
    myfile << std::endl;
    myfile.close();

    // save maximum concentration sequence of IGF1R
    myfile.open(saveDir+"/maxIGF1R.csv", std::ios_base::app);
    myfile << IGF1R_max[0];
    for(int i=1; i<IGF1R_max.size(); ++i){
	myfile << "," << IGF1R_max[i];
    }
    myfile << std::endl;
    myfile.close();

    // save average concentration sequence of EGFR
    myfile.open(saveDir+"/avgEGFR.csv", std::ios_base::app);
    myfile << EGFR_avg[0];
    for(int i=1; i<EGFR_avg.size(); ++i){
	myfile << "," << EGFR_avg[i];
    }
    myfile << std::endl;
    myfile.close();

    // save maximum concentration sequence of EGFR
    myfile.open(saveDir+"/maxEGFR.csv", std::ios_base::app);
    myfile << EGFR_max[0];
    for(int i=1; i<EGFR_max.size(); ++i){
	myfile << "," << EGFR_max[i];
    }
    myfile << std::endl;
    myfile.close();
}

// clean dead cells from cell grid and cell lists
void Environment::clean(CellGrids &cg, Diffusibles &diff){ 

    // reset allcells array
    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){
            if(cg.tc[i][j] == 1||cg.m0[i][j]==1||cg.m1[i][j]==1||cg.m2[i][j]==1){
                cg.allcells[i][j]=1;
            }    
        }
    }
    
    std::vector<Tumor> new_Tumor;
    int tc_alive = 0; // number of alive tumor cell
    int phago = 0; // number of tumor cell phagocytosed by M1
    int tc_dead = 0; // number of naturally dead tumor cell (age > lifespan)
    std::cout<<"tclist.size="<<tc_list.size()<<" mplist.size="<<mp_list.size()
        <<" dclist.size="<<dc_list.size()<<std::endl;
    
    for(auto & cell : tc_list){
        int i = cell.location[0];
        int j = cell.location[1];

        if(cell.state == "dead"){
            cg.allcells[i][j] = 1;
            cg.tc[i][j] = 0;
            cg.m1[i][j] = 0;
            cg.m2[i][j] = 0;
            cg.m0[i][j] = 0;
            cg.dead[i][j] = 1;
            diff.AKT[i][j] = 0;
            diff.ERK[i][j] = 0;
            diff.IGF1R[i][j] = 0;
            diff.EGFR[i][j] = 0;
            dc_list.push_back(DeadCell({i,j},"tc"));
            tc_dead++;
            continue;
        }      
        else if((cell.state == "alive" || cell.state == "quiescent")&& cg.tc[i][j] == 0){
            cg.allcells[i][j] = 1;
            cg.m1[i][j] = 1;
            cg.m2[i][j] = 0;
            cg.m0[i][j] = 0;
            diff.AKT[i][j] = 0;
            diff.ERK[i][j] = 0;
            diff.IGF1R[i][j] = 0;
            diff.EGFR[i][j] = 0;
            phago++;         
            continue;
        }       
        else if(cell.state == "alive" || cell.state == "quiescent"){
            cg.allcells[i][j]=1;
            cg.tc[i][j]=1;
            cg.m0[i][j]=0;
            cg.m1[i][j]=0;
            cg.m2[i][j]=0;
            new_Tumor.push_back(cell);
            tc_alive++;
        }     
    }
    
    int m_dead=0; // number of natrually dead TAMs (age > lifespan)
    int num_m0=0, num_m1=0, num_m2=0; // number of alive TAMs
    std::vector<Macrophage> new_mp;
    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){
            if(cg.m1[i][j] == 1 || cg.m2[i][j] == 1 || cg.m0[i][j] == 1){
                cg.m1[i][j] = 0;
                cg.m0[i][j] = 0;
                cg.m2[i][j] = 0;
                cg.allcells[i][j] = 0;
            }
        }
    }
    for(auto & cell : mp_list){
        int i = cell.location[0];
        int j = cell.location[1];
        bool error = false;
        for(auto & new_cell : new_mp){
            if(new_cell.location[0]==i && new_cell.location[1]==j){
                //std::cout << "Environ::clean macro! i="<<i<<" j="<<j<<" state="<<new_cell.state<<std::endl;
                error = true;
            }
        }
        if(error == true){continue;}
        if(cell.state == "dead"){
            cg.allcells[i][j] = 1;
            cg.dead[i][j] = 1;
            cg.m0[i][j] = 0;
            cg.m1[i][j] = 0;
            cg.m2[i][j] = 0;
            cg.tc[i][j] = 0;
            dc_list.push_back(DeadCell({i,j},"M"));
            m_dead++;
        }
        else if(cell.state == "M1"){
            cg.m1[i][j] = 1;
            cg.allcells[i][j] = 1;
            cg.m0[i][j] = 0;
            cg.m2[i][j] = 0;
            cg.tc[i][j] = 0;
            new_mp.push_back(cell);
            num_m1++;            
            }
        else if(cell.state == "M2"){
            cg.m2[i][j] = 1;
            cg.allcells[i][j] = 1;
            cg.m0[i][j] = 0;
            cg.m1[i][j] = 0;
            cg.tc[i][j] = 0;
            for(auto & new_cell : new_mp){
                if(new_cell.location[0]==i && new_cell.location[1]==j){
                    std::cout << "Environ::clean macro! i="<<i<<" j="<<j<<" state="<<new_cell.state<<std::endl;
                }
            }
            new_mp.push_back(cell);
            num_m2++;  
            }
        else if(cell.state == "M0" ){
            cg.m0[i][j] = 1;
            cg.allcells[i][j] = 1;
            cg.m1[i][j] = 0;
            cg.m2[i][j] = 0;
            cg.tc[i][j] = 0;
            for(auto & new_cell : new_mp){
                if(new_cell.location[0]==i && new_cell.location[1]==j){
                    std::cout << "Environ::clean macro! i="<<i<<" j="<<j<<" state="<<new_cell.state<<std::endl;
                }
            }
            new_mp.push_back(cell);
            num_m0++;
            }
        else{
            std::cout<<cell.location[0]<<","<<cell.location[1]<<" "<<cell.state
            <<" m0="<<cg.m0[i][j]<<" m1="<<cg.m1[i][j]<<" m2="<<cg.m2[i][j]<<std::endl;
            throw std::runtime_error("clean");
        } 
    }

    int phago_dead = 0; // dead cell phagotosed by M1
    int removal = 0; // dead cell removal from cell grid
    std::vector<DeadCell> new_dc;
    for(auto & cell : dc_list){
        int i = cell.location[0];
        int j = cell.location[1];
        if(cg.dead[i][j] == 1 && cell.source != "removal"){
            new_dc.push_back(cell);
        }
        else if(cg.dead[i][j] == 0){
            phago_dead++;
        }
        else if(cg.dead[i][j] == 1 && cell.source == "removal"){
            cg.dead[i][j] = 0;
            cg.allcells[i][j] = 0;
            removal++;
        }
    }
    tc_list = new_Tumor;
    tc_list.reserve(TC_max);
    mp_list = new_mp;
    mp_list.reserve(M_max);
    dc_list = new_dc;
    dc_list.reserve(DC_max);
    if(phago_num.size() > 0){
        phago_num.push_back(phago + phago_num[phago_num.size()-1]);
    }
    else{
        phago_num.push_back(phago);
    }

    // count cell number
    int sumc = 0;
    int sum0 = 0;
    int sum1 = 0;
    int sum2 = 0;
    int sumd = 0;
    for(int i = 0; i < 100; i++){
        for(int j = 0; j < 100; j++){
            sumc += cg.tc[i][j];
            sum0 += cg.m0[i][j];
            sum1 += cg.m1[i][j];
            sum2 += cg.m2[i][j];
            sumd += cg.dead[i][j];
        }
    }
    std::cout<<"tc_dead="<<tc_dead<<" phago="<<phago<<" m_dead="<<m_dead<<" removal = "<<removal
        <<" phago_dead="<<phago_dead<<" m0="<<num_m0<<" m1="<<num_m1<<" m2="<<num_m2<<std::endl;
    std::cout<<"new tclist.size="<<tc_list.size()<<" new mplist.size="<<mp_list.size()
        <<" new dclist.size="<<dc_list.size()<<std::endl;
    std::cout<<"cg.tc = "<<sumc<<" cg.m0 = "<<sum0<<" cg.m1 = "<<sum1<<" cg.m2 = "<<sum2
        <<" cg.dead = "<<sumd<<std::endl;
    std::cout<<" "<<std::endl;
  
}

// recruitment of M0
void Environment::recruitm0(CellGrids &cg, double rec) {

    // randomly place them in available sites based on vascular cell location
    std::mt19937 g(rd());
    std::uniform_real_distribution<double> dis(0.0,1.0);
    double p;
    double k;
    int num_m=0;
    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){
            p = dis(g); 
            if(cg.vas[i][j]==1){
                std::array<int, 4> ix={0,1,0,-1};
                std::array<int, 4> jx={1,0,-1,0};
                double probs[4]={0,0,0,0};
                double sum=0;
                for(int q=0; q<4; q++){
                    probs[q] = 1 - cg.allcells[i+ix[q]][j+jx[q]];
                    sum += probs[q];
                }
                if(sum == 0){continue;}
                double norm_probs[4];
                for(int q=0;q<4;q++){norm_probs[q] = probs[q]/(sum);} 
                for(int q=1; q<4; q++){norm_probs[q] = norm_probs[q] + norm_probs[q-1];} 
                int choice = 0;
                for(double norm_prob : norm_probs){
                    if(p > norm_prob){choice++;}
                }
                int ni = i + ix[choice];
                int nj = j + jx[choice];
                int m = mp_list.size();
                int c = tc_list.size();
                k = dis(g);
                if(k<(rec*(M0recProb)*tstep*60*60*(1-m/M_max)) ){  
                    num_m++;
                    cg.allcells[ni][nj] = 1;
                    cg.m0[ni][nj] = 1;
                    mp_list.push_back(Macrophage({ni,nj},1,"M0", K_M12, Kd_M12, K_I, Pha, a_I, a_C));
                } 
            }           
        }
    }
    std::cout<<"recruit m0="<<num_m<<std::endl;
}

// TAMs simulation
void Environment::run_Macrophage(CellGrids &cg, Diffusibles &diff, double depletion) {

    // shuffle the TAM list
    std::vector<int> run;
    for(int q=0; q<mp_list.size(); q++){
        run.push_back(q);
    }
    std::mt19937 g(rd());
    std::shuffle(run.begin(), run.end(), g);

    // run each TAM agent
    for(int L1 : run){
        
        // check cell grid error
        for(int i=1; i<99; i++){
            for(int j=1; j<99; j++){
                if(cg.tc[i][j]+cg.m0[i][j]+cg.m2[i][j]+cg.m1[i][j]+cg.vas[i][j]+cg.dead[i][j]!=cg.allcells[i][j]){
                    std::cout<<"L1 = "<<L1<<std::endl;
                    std::cout<<cg.tc[i][j]+cg.m2[i][j]+cg.m1[i][j]<<" "<<"i="<<i<<" "<<"j="<<j<<" "<<"allcells="<<
                    cg.allcells[i][j]<<" m0="<<cg.m0[i][j]<<" m1="<<cg.m1[i][j]<<
                    " m2="<<cg.m2[i][j]<<" tc="<<cg.tc[i][j]<<" vas="<<cg.vas[i][j]<<" dead="<<cg.dead[i][j]<<std::endl;
                    throw std::runtime_error("Environ::run_Macro");
                }
            }
        }
        
        // reset TAMs grid
        if(mp_list[L1].state == "M1" ){
            cg.m1[mp_list[L1].location[0]][mp_list[L1].location[1]] = 1;
            cg.m0[mp_list[L1].location[0]][mp_list[L1].location[1]] = 0;
            cg.m2[mp_list[L1].location[0]][mp_list[L1].location[1]] = 0;
        }
        else if(mp_list[L1].state == "M2" ){
            cg.m2[mp_list[L1].location[0]][mp_list[L1].location[1]] = 1;
            cg.m0[mp_list[L1].location[0]][mp_list[L1].location[1]] = 0;
            cg.m1[mp_list[L1].location[0]][mp_list[L1].location[1]] = 0;
        }
        else if(mp_list[L1].state == "M0" ){
            cg.m0[mp_list[L1].location[0]][mp_list[L1].location[1]] = 1;
            cg.m2[mp_list[L1].location[0]][mp_list[L1].location[1]] = 0;
            cg.m1[mp_list[L1].location[0]][mp_list[L1].location[1]] = 0;
        }
        else{
            // check TAM grid error
            std::cout<<mp_list[L1].location[0]<<","<<mp_list[L1].location[1]<<" "<<mp_list[L1].state
            <<" m0="<<cg.m0[mp_list[L1].location[0]][mp_list[L1].location[1]]<<" m1="<<cg.m1[mp_list[L1].location[0]][mp_list[L1].location[1]]<<
            " m2="<<cg.m2[mp_list[L1].location[0]][mp_list[L1].location[1]]<<std::endl;
        }

        // TAM simulation
        mp_list[L1].simulaion(tstep, cg, diff, depletion);     
    }

}

// tumor cells simulation
void Environment::run_Tumor(CellGrids &cg, Diffusibles &diff){

    // shuffle the tumor cell list
    std::vector<int> run;
    for(int q=0; q<tc_list.size(); q++){
        run.push_back(q);
    }
    std::shuffle(run.begin(), run.end(), rd);

    // run each tumor cell agent
    for(int L2 : run){
        tc_list[L2].simulation(tstep, cg, tc_list, diff, TCProlTime);
    }
}

// dead cells simulation
void Environment::run_DeadCell(CellGrids &cg){

    // shuffle the dead cell list
    std::vector<int> run;
    for(int q=0; q<dc_list.size(); q++){
        run.push_back(q);
    }
    std::shuffle(run.begin(), run.end(), rd);
    
    // run each dead cell agent
    for(int L3 : run){
        dc_list[L3].simulation(tstep, cg);
    }
}

// initialize cells for tumor growth simulation
void Environment::initializeCells(CellGrids &cg, Diffusibles &diff) {
    
    std::uniform_real_distribution<double> randType(0.0,1.0);
    std::mt19937 g(rd());

    // place 4 initial tumor cells in the center of the TME
    int z = 0;
    for(int i=49; i<51; i++){
        for(int j=49; j<51; j++){
                if(cg.allcells[i][j] == 0){
                cg.tc[i][j] = 1;
                cg.allcells[i][j] = 1;
                diff.AKT[i][j] = AKTi;
                diff.ERK[i][j] = ERKi;
                diff.IGF1R[i][j] = IGF1Ri;
                diff.EGFR[i][j] = EGFRi;
                tc_list.push_back(Tumor({i,j}, TCProlTime, z, 1, p_TC, p_D,
                                                EGFRi, IGF1Ri, AKTi, ERKi, a_ERK, TC_max));
                z++;
                }  
        }
    }

    // place initial TAMs randomly in the center of the TME
    double r = 0.5*M_max/TC_max;    
    for(int i=35; i<65; i++){
        for(int j=35; j<65; j++){
            double rand = randType(g);
            if(cg.allcells[i][j]==0 && rand < 0.02*r){
                cg.m0[i][j] = 1;
                cg.allcells[i][j] = 1;
                if(rand < 0.0*r){
                    mp_list.push_back(Macrophage({i,j}, 1,"M1", K_M12, Kd_M12, K_I, Pha, a_I, a_C));
                    cg.m0[i][j] = 0;
                    cg.m1[i][j] = 1;
                }
                else{
                    mp_list.push_back(Macrophage({i,j}, 1,"M2", K_M12, Kd_M12, K_I, Pha, a_I, a_C));
                    cg.m0[i][j] = 0;
                    cg.m2[i][j] = 1;
                }
            }
        }
    }
    // place initial TAMs randomly throughout of the TME
    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){       
            double rand = randType(g);
            if(cg.allcells[i][j]==0 && rand < 0.01*r){
                cg.m0[i][j] = 1;
                cg.allcells[i][j] = 1;
                if(rand < 0.0*r){
                    mp_list.push_back(Macrophage({i,j}, 1,"M1", K_M12, Kd_M12, K_I, Pha, a_I, a_C));
                    cg.m0[i][j] = 0;
                    cg.m1[i][j] = 1;
                }
                else{
                    mp_list.push_back(Macrophage({i,j}, 1,"M2", K_M12, Kd_M12, K_I, Pha, a_I, a_C));
                    cg.m0[i][j] = 0;
                    cg.m2[i][j] = 1;
                }
            }
        }
    }
}

// record time courses
void Environment::updateTimeCourses(int s, CellGrids &cg, Diffusibles &diff) {

    int sumc = 0; // number of tumor cell
    int sum0 = 0; // number of M0
    int sum1 = 0, sum11 = 0; // number of M1
    int sum2 = 0; // number of M2

    double avgIGF1 = 0;
    double avgEGF = 0;
    double avgCSF1 = 0;
    double avgCSF1R_I = 0;
    double avgIGF1R_I = 0;
    double avgEGFR_I = 0;
    double avgAKT = 0;
    double avgERK = 0;
    double avgIGF1R = 0;
    double avgEGFR = 0;
    double avgProb = 0;
    double avgProb_0 = 0;

    double maxIGF1 = 0;
    double maxEGF = 0;
    double maxCSF1 = 0;
    double maxCSF1R_I = 0;
    double maxIGF1R_I = 0;
    double maxEGFR_I = 0;
    double maxAKT = 0;
    double maxERK = 0;
    double maxIGF1R = 0;
    double maxEGFR = 0;
    double maxProb = 0;
    double maxProb_0 = 0;

    for(int i=1;i<99;i++){
        for(int j=1;j<99;j++){
            sumc += cg.tc[i][j];
            sum0 += cg.m0[i][j];
            sum1 += cg.m1[i][j];
            if(cg.m1[i][j] == 1){
                sum11++;
            }
            sum2 += cg.m2[i][j];

            avgIGF1 += diff.IGF1[i][j];
            avgEGF += diff.EGF[i][j];
            avgCSF1 += diff.CSF1[i][j];
            avgCSF1R_I += diff.CSF1R_I[i][j];
            avgIGF1R_I += diff.IGF1R_I[i][j];
            avgEGFR_I += diff.EGFR_I[i][j];
            avgAKT += diff.AKT[i][j];
            avgERK += diff.ERK[i][j];
            avgIGF1R += diff.IGF1R[i][j];
            avgEGFR += diff.EGFR[i][j];

            if(diff.CSF1[i][j] > maxCSF1){maxCSF1=diff.CSF1[i][j];}
            if(diff.IGF1[i][j] > maxIGF1){maxIGF1=diff.IGF1[i][j];}
            if(diff.EGF[i][j] > maxEGF){maxEGF=diff.EGF[i][j];}
            if(diff.CSF1R_I[i][j] > maxCSF1R_I){maxCSF1R_I=diff.CSF1R_I[i][j];}
            if(diff.IGF1R_I[i][j] > maxIGF1R_I){maxIGF1R_I=diff.IGF1R_I[i][j];}
            if(diff.EGFR_I[i][j] > maxEGFR_I){maxEGFR_I=diff.EGFR_I[i][j];}
            if(diff.AKT[i][j] > maxAKT){maxAKT=diff.AKT[i][j];}
            if(diff.ERK[i][j] > maxERK){maxERK=diff.ERK[i][j];}
            if(diff.IGF1R[i][j] > maxIGF1R){maxIGF1R=diff.IGF1R[i][j];}
            if(diff.EGFR[i][j] > maxEGFR){maxEGFR=diff.EGFR[i][j];}
        }
    }

    avgIGF1 = avgIGF1/(100*100);
    avgEGF = avgEGF/(100*100);
    avgCSF1 = avgCSF1/(100*100);
    avgCSF1R_I = avgCSF1R_I/(100*100);
    avgIGF1R_I = avgIGF1R_I/(100*100);
    avgEGFR_I = avgEGFR_I/(100*100);

    if(sumc > 0){
        avgAKT = avgAKT/sumc;
        avgERK = avgERK/sumc;
        avgIGF1R = avgIGF1R/sumc;
        avgEGFR = avgEGFR/sumc;
    }

    times.push_back(tstep*(s+1)/24);
    tc_num.push_back(sumc);
    m0_num.push_back(sum0);
    m1_num.push_back(sum1);
    m2_num.push_back(sum2);

    IGF1_max.push_back(maxIGF1);
    IGF1_avg.push_back(avgIGF1);
    EGF_max.push_back(maxEGF);
    EGF_avg.push_back(avgEGF);
    CSF1_max.push_back(maxCSF1);
    CSF1_avg.push_back(avgCSF1);

    CSF1R_I_avg.push_back(avgCSF1R_I);
    CSF1R_I_max.push_back(maxCSF1R_I);
    IGF1R_I_avg.push_back(avgIGF1R_I);
    IGF1R_I_max.push_back(maxIGF1R_I);
    EGFR_I_avg.push_back(avgEGFR_I);
    EGFR_I_max.push_back(maxEGFR_I);

    AKT_max.push_back(maxAKT);
    AKT_avg.push_back(avgAKT);
    ERK_avg.push_back(avgERK);
    ERK_max.push_back(maxERK);
    IGF1R_avg.push_back(avgIGF1R);
    IGF1R_max.push_back(maxIGF1R);
    EGFR_avg.push_back(avgEGFR);
    EGFR_max.push_back(maxEGFR);

    endTime = (s+1)*(tstep/24);
}

// check if the cell gird matches the cell list
void Environment::checkError(int s, CellGrids &cg) {

    int sumc = tc_num[s];
    int sum0 = m0_num[s];
    int sum1 = m1_num[s];
    int sum2 = m2_num[s];
    int sumv = 0;
    int sumd = 0;
    int allSum = 0;
    for(int i=1; i<99; ++i){
        for(int j=1; j<99; ++j){
            if(cg.allcells[i][j] > 1){
                throw std::runtime_error("environment::checkError 1");
            }
            if(cg.vas[i][j] == 1){
                sumv++;
            }
            if(cg.dead[i][j] == 1){
                sumd++;
            }
            allSum += cg.allcells[i][j];
        }
    }

    int macrophages = 0;
    for(auto & cell : mp_list){
        if(cell.state=="M0" || cell.state=="M1" || cell.state=="M2"){
            macrophages++;
        }
        else{
            std::cout<<"environment::checkError 2"<<std::endl;
        }
    }

    int canceralive = 0;
    int cancerquiescent = 0;
    for(auto & cell : tc_list){
        if(cell.state == "alive"){
            canceralive++;
        }
        else if(cell.state == "quiescent"){
            cancerquiescent++;
        }
    }

    if(canceralive+cancerquiescent != sumc){
        std::cout <<"tc_num error:"<< "canceralive = "<<canceralive <<" cancerquiescent = "<<cancerquiescent<< " sumc = " << sumc << " tc_list.size() = " << tc_list.size() << " allsum = " << allSum << std::endl;
        throw std::runtime_error("environment::checkError 3");
    }
    if(macrophages != (sum0 + sum1 + sum2)){
        std::cout <<"macro_num error:"<< macrophages << " " << sum0 << " " << sum1 << " " << sum2 << std::endl;
        throw std::runtime_error("environment::checkError 4");
    }
    if(sumc+sum0+sum1+sum2+sumv+sumd != allSum){
        std::cout << "allSum: " << allSum << std::endl;
        std::cout << "tc: " << sumc << std::endl;
        std::cout << "m0: " << sum0 << std::endl;
        std::cout << "macros: " << macrophages << std::endl;
        std::cout << s << std::endl;
        std::cout << std::endl;

        for(int i=1; i<99; ++i){
            for(int j=1; j<99; ++j){
                if(cg.allcells[i][j]==1 && cg.tc[i][j]+cg.m0[i][j]+cg.m1[i][j] +cg.m2[i][j]!=1){
                    std::cout << i << " " << j << " " << cg.tc[i][j]+cg.m0[i][j] <<  std::endl;
                    for(auto & cell : tc_list){
                        if(cell.location[0]==i && cell.location[1]==j){
                            std::cout << "cancer!\n";
                        }
                    }
                    for(auto & cell : mp_list){
                        if(cell.location[0]==i && cell.location[1]==j){
                            std::cout << "macro!\n";
                        }
                    }
                }
            }
        }
        throw std::runtime_error("environment::checkError 5");
    }
}

// print simulation information after a time step
void Environment::printStep(Diffusibles &diff) {

    int s = tc_num.size() - 1;
    double k1 = (tc_num[s]*1.0/TC_max);
    double k2 = (tc_list.size()*1.0/TC_max);
    std::cout
        << "Sim time: " << (s+1)*tstep/24 << "\n"
        << "Num tumor cells: " << tc_num[s] << " " << tc_list.size() << "\n"
        << "Density tumor cells (normalized): " << k1 << " " << k2 << "\n"
        << "Num M0: " << m0_num[s] <<" " <<"\n"
        << "Num M1: " << m1_num[s] <<" " <<"\n"
        << "Num M2: " << m2_num[s] <<" " <<"\n"
        << "Avg/Max IGF1: " << IGF1_avg[s]/diff.max_IGF1 << " | " << IGF1_max[s]/diff.max_IGF1 << "\n"
        << "Avg/Max EGF: " << EGF_avg[s]/diff.max_EGF << " | " << EGF_max[s]/diff.max_EGF << "\n"
        << "Avg/Max CSF1: " << CSF1_avg[s]/diff.max_CSF1 << " | " << CSF1_max[s]/diff.max_CSF1 << "\n"
        << "Avg/Max CSF1R_I: " << CSF1R_I_avg[s] << " | " << CSF1R_I_max[s] << "\n"
        << "Avg/Max IGF1R_I: " << IGF1R_I_avg[s] << " | " << IGF1R_I_max[s] << "\n";             
}

// save current virtual TME to csv files
void Environment::saveall(CellGrids &cg, Diffusibles &diff, int s){

    // directory to save virtual TME
    std::string str = "mkdir -p "+saveDir+"/saveall/"+ std::to_string(s);
    std::system(str.c_str());
    std::ofstream myfile;

    // save tumor cell list information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/tc_list.csv");
    int z = tc_list.size();

    myfile << "location[0]" <<","
        <<"location[1]"<<","
        <<"index"<<","
        <<"state"<<","
        <<"age"<<","
        <<"if_migra"<<","
        <<"mature"<<","
        <<"div_time"<<","
        <<"maxDiv"<<","
        <<"nDiv"<<","
        <<"cell_cycle"<<","
        <<"life_span"<<","
        <<"prob"<<","
        <<"prob_0"<<","
        <<"EGFR"<<","
        <<"IGF1R"<<","
        <<"AKT"<<","
        <<"ERK"<<","
        <<"a_ERK"<<","
        <<"a_IKT"<<std::endl;

    for(int i=0; i<z; i++){
        myfile <<tc_list[i].location[0]<< ","
            <<tc_list[i].location[1]<< ","
            <<tc_list[i].idx<< ","
            <<tc_list[i].state<< ","
            <<tc_list[i].age<< ","
            <<tc_list[i].if_migra<< ","
            <<tc_list[i].mature<< ","
            <<tc_list[i].div_time<< ","
            <<tc_list[i].maxDiv<< ","
            <<tc_list[i].nDiv<< ","
            <<tc_list[i].cell_cycle<< ","
            <<tc_list[i].life_span<< ","
            <<tc_list[i].prob<< ","
            <<tc_list[i].prob_0<< ","
            <<tc_list[i].EGFR<< ","
            <<tc_list[i].IGF1R<< ","
            <<tc_list[i].AKT<< ","
            <<tc_list[i].ERK<< ","
            <<tc_list[i].a_ERK<< ","
            <<tc_list[i].a_AKT<< std::endl;
    }    
    myfile.close();

    // save macrophage list information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/mp_list.csv");
    int k = mp_list.size();
    myfile << "location[0]" <<","
        <<"location[1]"<<","
        <<"state"<<","
        <<"age"<<","
        <<"now_phagocytose"<<","
        <<"max_phagocytose"<<","
        <<"CSF1R"<<","
        <<"trans"<<","
        <<"trans_num"<<","
        <<"rest_time"<<","
        <<"life_span"<<","
        <<"a_C"<<","
        <<"a_I"<<std::endl;

    for(int i=0; i<k; i++){
        myfile <<mp_list[i].location[0]<< ","
            <<mp_list[i].location[1]<< ","
            <<mp_list[i].state<< ","
            <<mp_list[i].age<< ","
            <<mp_list[i].now_phagocytose<< ","
            <<mp_list[i].max_phagocytose<< ","
            <<mp_list[i].CSF1R<< ","
            <<mp_list[i].trans<< ","
            <<mp_list[i].trans_num<< ","
            <<mp_list[i].rest_time<< ","
            <<mp_list[i].life_span<< ","
            <<mp_list[i].a_C<< ","
            <<mp_list[i].a_I<< std::endl;
    }   
    myfile.close();

    // save dead cell list information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/dc_list.csv");
    int p = dc_list.size();
    myfile << "location[0]" <<","
        <<"location[1]"<<","
        <<"source"<<","
        <<"removal_time"<<std::endl;

    for(int i=0; i<p; i++){
        myfile <<dc_list[i].location[0]<< ","
            <<dc_list[i].location[1]<< ","
            <<dc_list[i].source<<"," 
            <<dc_list[i].removal_time<< std::endl;
    }   
    myfile.close();

    // save CSF1 spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/CSF1.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.CSF1[i][0]/diff.max_CSF1;
        for(int j=1; j<100; ++j){
            myfile << "," << diff.CSF1[i][j]/diff.max_CSF1;
        }
        myfile << std::endl;
    }
    myfile.close();

    // save IGF1 spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/IGF1.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.IGF1[i][0]/diff.max_IGF1;
        for(int j=1; j<100; ++j){
            myfile << "," << diff.IGF1[i][j]/diff.max_IGF1;
        }
        myfile << std::endl;
    }
    myfile.close();

    // save integral parts in IGF1 concentration caculating
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/I.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.I[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save EGF spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/EGF.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.EGF[i][0]/diff.max_EGF;
        for(int j=1; j<100; ++j){
            myfile << "," << diff.EGF[i][j]/diff.max_EGF;
        }
        myfile << std::endl;
    }
    myfile.close();

    // save CSF1R_I spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/CSF1R_I.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.CSF1R_I[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.CSF1R_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save IGF1R_I spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/IGF1R_I.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.IGF1R_I[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.IGF1R_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save EGFR_I spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/EGFR_I.csv");
    for(int i=0; i<100; ++i){
        myfile << diff.EGFR_I[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << diff.EGFR_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    // save vascular cell spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/vas.csv");
    for(int i=0; i<100; ++i){
        myfile << cg.vas[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << cg.vas[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
}

// load virtual TME from csv files
void Environment::load(){
    
    std::string line; 
    std::string address = "./initial_virtual_TME"; // tumor cell density 0.58
    std::ifstream datatc( address + "/tc_list.csv");
    std::getline(datatc, line);

    // load tumor cells
    while(std::getline(datatc, line)){
        std::istringstream sin(line);      
        std::string tmp1;
        std::vector<std::string> tmp2;
        while(std::getline(sin, tmp1, ',')){           
            tmp2.push_back(tmp1);
        }
        tc_matrix.push_back(tmp2);  
    }
    std::cout<<"tc_matrix loaded"<<std::endl;  

    // load TAMs
    std::ifstream datamp(address + "/mp_list.csv");
    std::getline(datamp, line);
    while(std::getline(datamp, line)){
        std::istringstream sin(line);      
        std::string tmp1;
        std::vector<std::string> tmp2;
        while(std::getline(sin, tmp1, ',')){       
            tmp2.push_back(tmp1);
        }
        mp_matrix.push_back(tmp2);  
    }
    std::cout<<"mp_matrix loaded"<<std::endl; 

    // load dead cells
    std::ifstream datadc(address + "/dc_list.csv");
    std::getline(datadc, line);
    while(std::getline(datadc, line)){ 
        std::istringstream sin(line);      
        std::string tmp1;
        std::vector<std::string> tmp2;
        while(std::getline(sin, tmp1, ',')){       
            tmp2.push_back(tmp1);
        }
        dc_matrix.push_back(tmp2);  
    }
    std::cout<<"dc_matrix loaded"<<std::endl; 

    // load CSF1
    std::ifstream data_C(address + "/CSF1.csv");
    while(std::getline(data_C, line)){
        std::istringstream sin(line);      
        std::string tmp5;
        std::vector<double> tmp6;
        while(std::getline(sin, tmp5, ',')){          
            tmp6.push_back(std::stod(tmp5));
        }
        CSF1_matrix.push_back(tmp6);  
    }
    std::cout<<"CSF1_matrix loaded"<<std::endl; 

    // load IGF1
    std::ifstream data_IGF1(address + "/IGF1.csv");
    while(std::getline(data_IGF1, line)){ 
        std::istringstream sin(line);      
        std::string tmp7;
        std::vector<double> tmp8;
        while(std::getline(sin, tmp7, ',')){         
            tmp8.push_back(std::stod(tmp7));
        }
        IGF1_matrix.push_back(tmp8);  
    }
    std::cout<<"IGF1_matrix loaded"<<std::endl; 

    // load EGF
    std::ifstream data_EGF(address + "/EGF.csv");
    while(std::getline(data_EGF, line)){ 
        std::istringstream sin(line);      
        std::string tmp9;
        std::vector<double> tmp10;
        while(std::getline(sin, tmp9, ',')){       
            tmp10.push_back(std::stod(tmp9));
        }
        EGF_matrix.push_back(tmp10);  
    }
    std::cout<<"EGF_matrix loaded"<<std::endl; 

    // load integral parts in IGF1 concentration caculating
    std::ifstream data_I(address + "/I.csv");
    while(std::getline(data_I, line)){ 
        std::istringstream sin(line);      
        std::string tmp11;
        std::vector<double> tmp12;
        while(std::getline(sin, tmp11, ',')){       
            tmp12.push_back(std::stod(tmp11));
        }
        I_matrix.push_back(tmp12);  
    }
    std::cout<<"I_matrix loaded"<<std::endl;  

    // load CSF1R_I
    std::ifstream data_CR_I(address + "/CSF1R_I.csv");
    while(std::getline(data_CR_I, line)){ 
        std::istringstream sin(line);      
        std::string tmp13;
        std::vector<double> tmp14;
        while(std::getline(sin, tmp13, ',')){           
            tmp14.push_back(std::stod(tmp13));
        }
        CSF1R_I_matrix.push_back(tmp14);  
    }
    std::cout<<"CSF1RI_matrix loaded"<<std::endl; 

    // load vascular cells
    std::ifstream data_vas(address + "/vas.csv");
    while(std::getline(data_vas, line)){ 
        std::istringstream sin(line);      
        std::string tmp13;
        std::vector<int> tmp14;
        while(std::getline(sin, tmp13, ',')){           
            tmp14.push_back(std::stoi(tmp13));
        }
        vas_matrix.push_back(tmp14);  
    }
    std::cout<<"vas_matrix loaded"<<std::endl; 
}

// initialize virtual TME
void Environment::initialization(CellGrids &cg, Diffusibles &diff){

    std::cout<<"initialization start"<<std::endl;

    // initialize vascular cells in cell grid
    for(int i=1; i<99; i++){
        for(int j=1; j<99; j++){   
            cg.vas[i][j] = vas_matrix[i][j];
            cg.allcells[i][j] = vas_matrix[i][j];
        }
    }

    // initialize cytokines and drugs
    for(int i=0; i<100; i++){
        for(int j=0; j<100; j++){    
            diff.CSF1[i][j] = CSF1_matrix[i][j] * diff.max_CSF1;
            diff.EGF[i][j] = EGF_matrix[i][j] * diff.max_EGF;
            diff.IGF1[i][j] = IGF1_matrix[i][j] * diff.max_IGF1;
            diff.CSF1R_I[i][j] = CSF1R_I_matrix[i][j];
            diff.I[i][j] = I_matrix[i][j];
        }
    }
    std::cout<<"vas initialized; diff initialized"<<std::endl;

    // initialize tumor cells
    for(int i=0; i<tc_matrix.size(); i++){

        // add tumor cell to the list
        tc_list.push_back(Tumor({std::stoi(tc_matrix[i][0]),std::stoi(tc_matrix[i][1])},
                                std::stod(tc_matrix[i][6]),
                                std::stoi(tc_matrix[i][2]),
                                0,p_TC, p_D, 
                                std::stod(tc_matrix[i][14]),
                                std::stod(tc_matrix[i][15]),
                                std::stod(tc_matrix[i][16]),
                                std::stod(tc_matrix[i][17]),
                                a_ERK, TC_max));

        // initialize tumor cell attributes
        tc_list[i].state = tc_matrix[i][3];
        tc_list[i].age = std::stod(tc_matrix[i][4]);
        tc_list[i].if_migra = std::stoi(tc_matrix[i][5]);
        tc_list[i].div_time = std::stod(tc_matrix[i][7]);
        tc_list[i].maxDiv = std::stoi(tc_matrix[i][8]);
        tc_list[i].nDiv = std::stoi(tc_matrix[i][9]);
        tc_list[i].cell_cycle = std::stod(tc_matrix[i][10]);
        tc_list[i].life_span = std::stod(tc_matrix[i][11]);
        tc_list[i].prob = std::stod(tc_matrix[i][12]);
        tc_list[i].prob_0 = 0;
        tc_list[i].a_AKT = std::stod(tc_matrix[i][19]);

        // check loaction available
        if(cg.allcells[std::stoi(tc_matrix[i][0])][std::stoi(tc_matrix[i][1])] != 0){
            std::cout<<"initialization error1"<<std::endl;
        }

        // add tumor cell to cell grid
        cg.tc[std::stoi(tc_matrix[i][0])][std::stoi(tc_matrix[i][1])] = 1;
        cg.allcells[std::stoi(tc_matrix[i][0])][std::stoi(tc_matrix[i][1])] = 1; 
    }

    // number of M0, M1, M2
    int sum0 = 0;
    int sum1 = 0;
    int sum2 = 0;

    // initialize TAMs 
    for(int i=0; i<mp_matrix.size(); i++){
        // add TAM to the list        
        mp_list.push_back(Macrophage({std::stoi(mp_matrix[i][0]),std::stoi(mp_matrix[i][1])},
                                1,mp_matrix[i][2], K_M12, Kd_M12, K_I, Pha, a_I, a_C));

        // initialize TAM attributes
        mp_list[i].age = std::stod(mp_matrix[i][3]);
        mp_list[i].now_phagocytose = std::stoi(mp_matrix[i][4]);
        mp_list[i].max_phagocytose = std::stoi(mp_matrix[i][5]);
        mp_list[i].CSF1R = std::stoi(mp_matrix[i][6]);
        mp_list[i].trans = std::stoi(mp_matrix[i][7]);
        mp_list[i].trans_num = std::stoi(mp_matrix[i][8]);
        mp_list[i].rest_time = std::stod(mp_matrix[i][9]);
        mp_list[i].life_span = std::stod(mp_matrix[i][10]);

        // check location available
        if(cg.allcells[std::stoi(mp_matrix[i][0])][std::stoi(mp_matrix[i][1])] != 0){
            std::cout<<"initialization error2"<<std::endl;
        }

        // add TAM to cell grid
        cg.allcells[std::stoi(mp_matrix[i][0])][std::stoi(mp_matrix[i][1])] = 1; 
        if(mp_matrix[i][2] == "M0"){
            sum0++;
            //std::cout<<"mp "<<i<<" "<<mp_matrix[i][2]<<std::endl;
            cg.m0[std::stoi(mp_matrix[i][0])][std::stoi(mp_matrix[i][1])] = 1; 
        }
        else if(mp_matrix[i][2] == "M1"){
            sum1++;
            cg.m1[std::stoi(mp_matrix[i][0])][std::stoi(mp_matrix[i][1])] = 1; 
        }
        else if(mp_matrix[i][2] == "M2"){
            sum2++;
            cg.m2[std::stoi(mp_matrix[i][0])][std::stoi(mp_matrix[i][1])] = 1; 
        }
        else{
            std::cout<<"mp_matrix error"<<std::endl;
        }

        // check one grid one cell
        if(cg.m0[std::stoi(mp_matrix[i][0])][std::stoi(mp_matrix[i][1])]
            + cg.m1[std::stoi(mp_matrix[i][0])][std::stoi(mp_matrix[i][1])]
            + cg.m2[std::stoi(mp_matrix[i][0])][std::stoi(mp_matrix[i][1])]
            > 1){
            std::cout<<"initialization error3"<<std::endl;    
        }    
    }
    std::cout<<"mp_matrix size = "<<mp_matrix.size()<<" cg.m0 = "<<sum0<<" cg.m1 = "<<sum1<<" cg.m2 = "<<sum2<<std::endl;

    // initialize dead cells
    for(int i=0; i<dc_matrix.size(); i++){
        dc_list.push_back(DeadCell({std::stoi(dc_matrix[i][0]),std::stoi(dc_matrix[i][1])}, dc_matrix[i][2]));
        dc_list[i].removal_time = std::stod(dc_matrix[i][3]);
        cg.dead[std::stoi(dc_matrix[i][0])][std::stoi(dc_matrix[i][1])] = 1;
        cg.allcells[std::stoi(dc_matrix[i][0])][std::stoi(dc_matrix[i][1])] = 1;
    }
}

// run a simulation
void Environment::simulate(double days){
  
    // create cell gird
    CellGrids cg;
    std::cout<<"cell grid created \n"<<std::endl;

    // create diffusibles for cytokines and drugs
    Diffusibles diff(dx);
    std::cout<<"Diffusibles created \n"<<std::endl;

    // initialize cells in the cell gird
    initializeCells(cg, diff);
    std::cout<<"initialize cells for tumor growth \n"<<std::endl;

    // plot initialization virtual TME spatial state
    plot(cg, diff, 100000);

    double depletion = 0; // depletion of TAMs by drug
    double rec = 1; // '1': recuit M0
                    // '0': no recuitment

    int save_time = 0; // save virtual TME

    // simulation in each time step
    for(int s=0; s<days*24/tstep; s++){
        std::cout<< "-------------------------"<<std::endl;
        std::cout<<"time step "<<s<<" start"<<std::endl;
        diff.diffusion(cg, tstep);
        std::cout<<"diffusion"<<std::endl;

        recruitm0(cg, rec);
        std::cout<<"recruit_M0"<<std::endl;

        run_Macrophage(cg, diff, depletion);
        std::cout<<"run_Macro"<<std::endl;
        clean(cg,diff);

        run_Tumor(cg, diff);
        std::cout<<"run_Tmuor"<<std::endl;
        clean(cg,diff);

        run_DeadCell(cg);
        std::cout<<"run_DeadCell"<<std::endl;
        clean(cg, diff);

        updateTimeCourses(s, cg, diff);
        checkError(s, cg);       
        printStep(diff);

        // plot spatial states
        if(fmod(s*tstep, 480*tstep) == 0){           
            plot(cg, diff, s);
        }
        
        // save virtual TME
        if(tc_list.size() > 0.58*TC_max && save_time <= 2){
            std::cout<<"saveall"<<std::endl;
            saveall(cg, diff, s);
            save_time++;
        }
        if(tc_list.size() > 0.7*TC_max && save_time <= 4){
            std::cout<<"saveall"<<std::endl;
            saveall(cg, diff, s);
            save_time++;
        }
        if(tc_list.size() > 0.85*TC_max && save_time <= 6){
            std::cout<<"saveall"<<std::endl;
            saveall(cg, diff, s);
            save_time++;
        }
        if(tc_list.size() > 0.92*TC_max && save_time <= 8){
            std::cout<<"saveall"<<std::endl;
            saveall(cg, diff, s);
            save_time++;
        }
        if(tc_list.size() > 0.95*TC_max && save_time <= 10){
            std::cout<<"saveall"<<std::endl;
            saveall(cg, diff, s);
            save_time++;
        }

        // break condition
        if(tc_num[s] == 0 || tc_num[s] >= TC_max){
            std::cout<<"Break!"<<std::endl;
            break;}
    }
    // plot terminal spatial state
    plot(cg, diff, 100001);

    // save number of cells and average/maximum concentration of cytokines
    save(diff);
}

/*
    REFERENCES

    D. F. Quail, R. L. Bowman, L. Akkari, M. L. Quick, A. J. Schuhmacher, J. T. Huse, 
    E. C. Holland, J. C. Sutton, J. A. Joyce, The tumor microenvironment underlies 
    acquired resistance to CSF-1R inhibition in gliomas. Science 352, aad3018 (2016).

    Y. Zheng, J. Bao, Q. Zhao, T. Zhou, X. Sun, A spatio-temporal model of macrophage-mediated 
    drug resistance in glioma immunotherapy. Mol Cancer Ther 17, 814-824 (2018).

    C. G. Cess, S. D. Finley, Multi-scale modeling of macrophage-T cell interactions 
    within the tumor microenvironment. PLoS Comput Biol 16, e1008519 (2020).

    A. R. A. Anderson, A hybrid mathematical model of solid tumour invasion: the 
    importance of cell adhesion. Mathematical medicine and biology : a journal of 
    the IMA 22, 163-186 (2005).

    A. R. A. Anderson, “A Hybrid Discrete-continuum Technique for Individual-based 
    Migration Models” in Polymer and Cell Dynamics, W. Alt, M. Chaplain, M. Griebel, 
    J. Lenz, Eds. (Birkhäuser, Basel, Switzerland, 2003).
 */
