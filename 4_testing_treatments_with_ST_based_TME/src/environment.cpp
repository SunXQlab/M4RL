#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "environment.h"

std::random_device rd;

Environment::Environment(double stepSize, std::string folder, std::string st_folder, int set, double M0RecProb, double TCProl,
                         double k_M12, double kd_M12, double k_I, double p_tc, double p_d, int pha, int rand, 
                         double a_c, double a_i, double a_erk, int grid_size) {
  // directory to save time courses
    saveDir ="./"+folder+"/set_"+std::to_string(set);
    load_folder = st_folder;

    dx = 0.006; // 60um, diameter of cell spot based on ST data

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
    G_size = grid_size;
    
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
    
    std::cout<<"G_size = "<<G_size<<std::endl;
    std::cout<<"default TC_max = "<<TC_max<<std::endl;
    std::cout<<"default M_max = "<<M_max<<std::endl;
    std::cout<<"default DC_max = "<<DC_max<<std::endl;
}

// plot spatial states of the simulation
void Environment::plot(CellGrids &cg, Diffusibles &diff, int s){ 
  
    // identify tumor cell state
    for(int i=0; i<G_size; i++){
        for(int j=0; j<G_size; j++){
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
    for(int i=0; i<G_size; ++i){
        myfile << tcs[i][0] + 2*cg.m0[i][0] + 3*cg.m1[i][0] + 4*cg.m2[i][0] 
            + 5*cg.vas[i][0] + 6*cg.dead[i][0];
        for(int j=1; j<100; ++j){
            myfile << "," << tcs[i][j] + 2*cg.m0[i][j] + 3*cg.m1[i][j] 
            + 4*cg.m2[i][j] + 5*cg.vas[i][j]  + 6*cg.dead[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/CSF1.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.CSF1[i][0]/diff.max_CSF1;
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.CSF1[i][j]/diff.max_CSF1;
        }
        myfile << std::endl;
    }
    myfile.close();
    
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/IGF1.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.IGF1[i][0]/diff.max_IGF1;
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.IGF1[i][j]/diff.max_IGF1;
        }
        myfile << std::endl;
    }
    myfile.close();
    
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/EGF.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.EGF[i][0]/diff.max_EGF;
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.EGF[i][j]/diff.max_EGF;
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/AKT.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.AKT[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.AKT[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/ERK.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.ERK[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.ERK[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/IGF1R.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.IGF1R[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.IGF1R[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/EGFR.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.EGFR[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.EGFR[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/CSF1R_I.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.CSF1R_I[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.CSF1R_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/IGF1R_I.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.IGF1R_I[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.IGF1R_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    /*
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/EGFR_I.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.EGFR_I[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.EGFR_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/prob.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.prob[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.prob[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();

    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/prob0.csv");
    for(int i=0; i<G_size; ++i){//  ӡ100*100  grid
        myfile << diff.prob_0[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.prob_0[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    */
    
    if(s == 0){
    myfile.open(saveDir+"/spatial/"+std::to_string(s)+"/vas.csv");
    for(int i=0; i<G_size; ++i){
        myfile << cg.vas[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << cg.vas[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    }
    // reset zero    
    for(int i=0; i<G_size; i++){
        for(int j=0; j<G_size; j++){
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
    for(int i=1; i<G_size-1; i++){
        for(int j=1; j<G_size-1; j++){
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
            dc_list.push_back(DeadCell({i,j},"TC"));
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
            tc_alive++;
            new_Tumor.push_back(cell);
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
    for(int i = 0; i < G_size; i++){
        for(int j = 0; j < G_size; j++){
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
    for(int i=1; i<G_size-1; i++){
        for(int j=1; j<G_size-1; j++){
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
        for(int i=1; i<G_size-1; i++){
            for(int j=1; j<G_size-1; j++){
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
    
    for(int i=1;i<G_size-1;i++){
        for(int j=1;j<G_size-1;j++){
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
    for(int i=1; i<G_size-1; ++i){
        for(int j=1; j<G_size-1; ++j){
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

        for(int i=1; i<G_size-1; ++i){
            for(int j=1; j<G_size-1; ++j){
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
    for(int i=0; i<G_size; ++i){
        myfile << diff.CSF1[i][0]/diff.max_CSF1;
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.CSF1[i][j]/diff.max_CSF1;
        }
        myfile << std::endl;
    }
    myfile.close();
    
    // save IGF1 spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/IGF1.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.IGF1[i][0]/diff.max_IGF1;
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.IGF1[i][j]/diff.max_IGF1;
        }
        myfile << std::endl;
    }
    myfile.close();
    
    // save integral parts in IGF1 concentration caculating
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/I.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.I[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    
    // save EGF spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/EGF.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.EGF[i][0]/diff.max_EGF;
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.EGF[i][j]/diff.max_EGF;
        }
        myfile << std::endl;
    }
    myfile.close();
    
    // save CSF1R_I spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/CSF1R_I.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.CSF1R_I[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.CSF1R_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    
    // save IGF1R_I spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/IGF1R_I.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.IGF1R_I[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.IGF1R_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    
    // save EGFR_I spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/EGFR_I.csv");
    for(int i=0; i<G_size; ++i){
        myfile << diff.EGFR_I[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << diff.EGFR_I[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
    
    // save vascular cell spatial information
    myfile.open(saveDir+"/saveall/"+std::to_string(s)+"/vas.csv");
    for(int i=0; i<G_size; ++i){
        myfile << cg.vas[i][0];
        for(int j=1; j<G_size; ++j){
            myfile << "," << cg.vas[i][j];
        }
        myfile << std::endl;
    }
    myfile.close();
}

// load virtual TME from ST data-based csv files
void Environment::load(){
    std::string line; 
    std::string address = "./stGBM_data/UKF_" + load_folder;
    std::cout<<"load start: "<<address<<std::endl;
    
    std::ifstream data_cell(address + "/After_interpolation_cellgrid_" 
                              + std::to_string(G_size) +".csv");
    while(std::getline(data_cell, line)){
        std::istringstream sin(line);      
        std::string tmp1;
        std::vector<int> tmp2;
        while(std::getline(sin, tmp1, ',')){
            tmp2.push_back(std::stoi(tmp1));
        }
        cell_matrix.push_back(tmp2);
    }
    std::cout<<"cell_matrix loaded"<<std::endl; 
  
    std::ifstream data_CSF1(address + "/Interpolation_CSF1_" 
                              + std::to_string(G_size) +".csv");
    while(std::getline(data_CSF1, line)){
        std::istringstream sin(line);      
        std::string tmp5;
        std::vector<double> tmp6;
        while(std::getline(sin, tmp5, ',')){
            tmp6.push_back(std::stod(tmp5));
        }
        CSF1_matrix.push_back(tmp6);  
    }
    std::cout<<"CSF1_matrix loaded"<<std::endl; 
    
    std::ifstream data_IGF1(address + "/Interpolation_IGF1_"
                               + std::to_string(G_size) +".csv");
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
    
    std::ifstream data_EGF(address + "/Interpolation_EGF_" 
                             + std::to_string(G_size) +".csv");
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
}

// initialize virtual TME
void Environment::initialization(CellGrids &cg, Diffusibles &diff){
    std::cout<<"initialization start"<<std::endl;
    int tc = 0;
    int m1 = 0;
    int m2 = 0;
    int vas = 0;
    int dead = 0;
    for(int i=1; i<G_size-1; i++){
        for(int j=1; j<G_size-1; j++){
            if(cell_matrix[i][j] == 1){
                tc = tc + 1;
            }
            else if(cell_matrix[i][j] == 3){
                m1 = m1 + 1;
            }
            else if(cell_matrix[i][j] == 4){
                m2 = m2 + 1;
            }
            else if(cell_matrix[i][j] == 5){
                vas = vas + 1;
            }
            else if(cell_matrix[i][j] == 6){
                dead = dead + 1;
            }
        }
    }
    std::cout << "TC number: " << tc << std::endl;
    std::cout << "M1 number: " << m1 << std::endl;
    std::cout << "M2 number: " << m2 << std::endl;
    std::cout << "VS number: " << vas << std::endl;
    std::cout << "DC number: " << dead << std::endl;
    
    TC_max = std::ceil(tc * 1.2);
    if((m1 + m2) > M_max){
        M_max = std::ceil((m1 + m2) * 1.1);
    }
    std::cout<<"TC_max reset: "<< TC_max <<" M_max reset: "<< M_max <<std::endl;
    tc_list.reserve(TC_max);
    mp_list.reserve(M_max);
    dc_list.reserve(DC_max);
    
    std::mt19937 g(rd()); 
    std::uniform_real_distribution<> dis(0.0,1.0);
    for(int i=1; i<G_size-1; i++){
        for(int j=1; j<G_size-1; j++){
            if(cell_matrix[i][j] == 1 && cg.allcells[i][j] == 0){
                cg.tc[i][j] = 1;
                cg.allcells[i][j] = 1;
                tc_list.push_back(Tumor({i,j},
                                        TCProlTime,
                                         tc,
                                         0,p_TC, p_D,
                                         EGFRi,
                                         IGF1Ri,
                                         AKTi,
                                         ERKi,
                                         a_ERK, TC_max));
                
                tc_list[tc-1].age = 24.0 * 8 * dis(g);
                tc_list[tc-1].div_time = TCProlTime * dis(g);
                tc_list[tc-1].nDiv = std::floor(8 * dis(g));
              }
              else if(cell_matrix[i][j] == 3 && cg.allcells[i][j] == 0){
                  cg.m1[i][j] = 1;
                  cg.allcells[i][j] = 1;
                  mp_list.push_back(Macrophage({i,j}, 1,"M1", K_M12, Kd_M12, K_I, Pha, a_I, a_C));
                  mp_list[m1+m2-1].now_phagocytose = std::floor(2 * dis(g));
              }
              else if(cell_matrix[i][j] == 4 && cg.allcells[i][j] == 0){
                  cg.m2[i][j] = 1;
                  cg.allcells[i][j] = 1;
                  mp_list.push_back(Macrophage({i,j}, 1,"M2", K_M12, Kd_M12, K_I, Pha, a_I, a_C));
                  mp_list[m1+m2-1].now_phagocytose = std::floor(2 * dis(g));
              }
              else if(cell_matrix[i][j] == 5 && cg.allcells[i][j] == 0){
                  cg.vas[i][j] = 1;
                  cg.allcells[i][j] = 1;
              }
              else if(cell_matrix[i][j] == 6 && cg.allcells[i][j] == 0){
                  cg.dead[i][j] = 1;
                  cg.allcells[i][j] = 1;
                  dc_list.push_back(DeadCell({i,j}, "TC"));
                  dc_list[dead-1].removal_time = 24.0 * dis(g);
              }
          }
      }
      int quiescent_num = 0;
      for(int k=0; k<tc_list.size(); k++){
          int i = tc_list[k].location[0];
          int j = tc_list[k].location[1];
          int z = 0;
          std::vector<int> ix;
          std::vector<int> jx;
          for(int il=-1; il<2; il++){
              for(int jl=-1; jl<2; jl++){
                  ix.push_back(il);
                  jx.push_back(jl);
                  z++;
              }
          }
          double probs[z];
          double sum=0;
          for(int q=0; q<z; q++){
              probs[q] = (1 - cg.allcells[i+ix[q]][j+jx[q]]);    
              sum += probs[q];
          }
          if(sum == 0){
              tc_list[k].state = "quiescent";
              tc_list[k].prob = 1.0;
              tc_list[k].div_time = TCProlTime + 1e-7;
              quiescent_num++;
          }
      }
      std::cout << "quiescent_num : " << quiescent_num << std::endl;
      std::cout<<"cell grids loaded; cell lists loaded"<<std::endl;
      for(int i=0; i<G_size; i++){
          for(int j=0; j<G_size; j++){
            if(!std::isnan(CSF1_matrix[i][j])){
                diff.CSF1[i][j] = CSF1_matrix[i][j] * diff.max_CSF1;
            }
            if(!std::isnan(EGF_matrix[i][j])){
                diff.EGF[i][j] = EGF_matrix[i][j] * diff.max_EGF;
            } 
            if(!std::isnan(IGF1_matrix[i][j])){
                diff.IGF1[i][j] = IGF1_matrix[i][j] * diff.max_IGF1;
            }       
        }
    }
    std::cout<<"diff loaded"<<std::endl;
}

// run a simulation
void Environment::simulate(double days){

    // create cell gird
    CellGrids cg(G_size);
    std::cout<<"cell grid created \n"<<std::endl;
    
    // create diffusibles for cytokines and drugs
    Diffusibles diff(dx, G_size);
    std::cout<<"Diffusibles created \n"<<std::endl;
    
    // load virtual TME from ST data-based csv files
    load();
    std::cout<<"virtual TME loaded \n"<<std::endl;
    
    // initialize virtual TME
    initialization(cg, diff);
    std::cout<<"virtual TME initialized \n"<<std::endl;
    
    // plot initialization virtual TME spatial state
    plot(cg, diff, 100000);
    
    double depletion = 0; // depletion of TAMs by drug
    double rec = 1; // '1': recuit M0
                    // '0': no recuitment
    
    //int save_time = 0; // save virtual TME (if necessary)
    
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
        
        // break condition
        if(tc_num[s] == 0 || (tc_num[s] > 0.9 * TC_max && s > 50*48 )){
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
 
    V. M. Ravi, P. Will, J. Kueckelhaus, N. Sun, K. Joseph, H. Salié, L. Vollmer, 
    U. Kuliesiute, J. von Ehr, J. K. Benotmane, N. Neidert, M. Follo, F. Scherer, 
    J. M. Goeldner, S. P. Behringer, P. Franco, M. Khiat, J. Zhang, U. G. Hofmann, 
    C. Fung, F. L. Ricklefs, K. Lamszus, M. Boerries, M. Ku, J. Beck, R. Sankowski, 
    M. Schwabenland, M. Prinz, U. Schüller, S. Killmer, B. Bengsch, A. K. Walch, 
    D. Delev, O. Schnell, D. H. Heiland, Spatially resolved multi-omics deciphers 
    bidirectional tumor-host interdependence in glioblastoma. Cancer Cell 40, 
    639-655.e613 (2022).

    J. Kueckelhaus, S. Frerich, J. Kada-Benotmane, C. Koupourtidou, J. Ninkovic, 
    M. Dichgans, J. Beck, O. Schnell, D. H. Heiland, Inferring histology-associated 
    gene expression gradients in spatial transcriptomic studies. Nature Communications 
    15, 7280 (2024).
 */
