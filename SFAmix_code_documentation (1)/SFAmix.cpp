#include <stdlib.h>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>

#include <iostream>
#include <Eigen/Dense>
#include "myHeader.cpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
//#include <unsupported/Eigen/MatrixFunctions>


using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

//const double PI  =3.141592653589793238462;
/*
Chuan Gao C++
*/

/* 
usage
./SFAmix --nf 20 --y ./data/Y_1.txt --out result --sep tab
*/

int main(int argc,char *argv[]){
    
    // declare variables 
    int nf=50,nf2=0,s_n=0,d_y=0,i_seed=0,n_itr=1000000;
    double a=0.5,b=0.5,c=0.5,d=0.5,g=0.5,h=0.5,alpha=1,beta=1;
    
    string file_y,dir_out,sep;
    stringstream ss;

    sep="tab";
    
    // read in argument
    string s_df="--nf",s_y="--y",s_out="--out",s_sep="--sep",s_a="--a",s_b="--b",s_c="--c",s_d="--d",s_seed="--seed",s_itr="--itr";
    for(int i=0;i<argc;i++){
        if(s_df.compare(argv[i])==0){nf=atoi(argv[i+1]);}
        if(s_y.compare(argv[i])==0){file_y=argv[i+1];}
        if(s_out.compare(argv[i])==0){dir_out=argv[i+1];}
        if(s_sep.compare(argv[i])==0){sep=argv[i+1];}
        //if(s_a.compare(argv[i])==0){a=atof(argv[i+1]);}
        //if(s_b.compare(argv[i])==0){b=atof(argv[i+1]);}
        //if(s_c.compare(argv[i])==0){c=atof(argv[i+1]);}
        //if(s_d.compare(argv[i])==0){d=atof(argv[i+1]);}
        //if(s_seed.compare(argv[i])==0){i_seed=atoi(argv[i+1]);}
        //if(s_itr.compare(argv[i])==0){n_itr=atoi(argv[i+1]);}
    }
    
    
    // calculate the sample size and the gene numbers 
    string line;
    string field;
    ifstream f;
    
    f.open(file_y.c_str());
    if (! f.is_open())
    {
        printf("Gene expression file open failed\n");exit(0);
    }
    getline(f,line);
    s_n++;
    istringstream iss(line);
    if(sep.compare("space")==0){
        while(getline(iss,field,' ')){d_y++;}
    }else if(sep.compare("tab")==0){
        while(getline(iss,field,'\t')){d_y++;}
    }else{
        cout << "Please specify a valid separator." << endl << endl;
    }
    while(getline(f,line)){s_n++;}

    // write command into file for later references, also write the dimension of the gene expression matrix 
    ss.str("");
    ss.clear();
    ss << dir_out << "/command.txt";
    ofstream f_com (ss.str().c_str());
    if (f_com.is_open()){
        for(int i=0;i<argc;i++){
            f_com << argv[i] << " ";
        }
        f_com << endl;
    }
    f_com << endl << "Y_TMP has dimension of " << s_n << " by " << d_y << endl << endl;
    
    // read in the Y matrix 
    MatrixXd Y_TMP=MatrixXd::Constant(s_n,d_y,0);
    
    f.clear();
    f.seekg (0, ios_base::beg);
    int i=0,j=0;
    while(getline(f,line)){
        istringstream iss(line);
        j=0;
        if(sep.compare("space")==0){
            while(getline(iss,field,' ')){
                Y_TMP(i,j)=atof(field.c_str());
                j++;
            }
            i++;
        }else{
            while(getline(iss,field,'\t')){
                Y_TMP(i,j)=atof(field.c_str());
                j++;
            }
            i++;
        }
    }
    f.close();
    
    // chunk Y to a smaller sub matrix to test if necessary
    // s_n=500;
    MatrixXd Y=MatrixXd::Constant(s_n,d_y,0);
    //Y=Y_TMP.block(0,0,s_n,d_y);
    Y=Y_TMP.transpose();
    //f_com << "Submatrix of Y has dimension of " << s_n << " by " << d_y << endl;
    f_com.close();
 
    // The orginal model is set up as Y=\Lambda X, to paper set it up as Y=X \Lambda, So the model is transposed
    int d_y_bak=d_y;
    d_y=s_n;
    s_n=d_y_bak;
    
    // Declare variables independent of factor number to prepare for the EM algorithm
    
    VectorXd psi_v = VectorXd::Constant(s_n,1);
    MatrixXd PSI=MatrixXd::Constant(s_n,s_n,0);
    MatrixXd PSI_INV=MatrixXd::Constant(s_n,s_n,0);
    MatrixXd LX=MatrixXd::Constant(s_n,d_y,0);
    MatrixXd YLX=MatrixXd::Constant(s_n,d_y,0);
    MatrixXd YLXi=VectorXd::Constant(s_n,0);
    MatrixXd YLX2=MatrixXd::Constant(s_n,d_y,0);
    VectorXd TOP=VectorXd::Constant(s_n,0);
    
    // Declare variables that are dependent of the factor number
    nf2 = nf;
    int nt=nf;
    
    
    MatrixXd EX=MatrixXd::Constant(nt,d_y,0);
    MatrixXd TEX=MatrixXd::Constant(d_y,nt,0);
    MatrixXd VX=MatrixXd::Constant(nt,nt,0);
    MatrixXd EXX=MatrixXd::Constant(nt,nt,0);
    MatrixXd LP=MatrixXd::Constant(nt,s_n,0);
    MatrixXd ID=MatrixXd::Constant(nt,nt,0);
    VectorXd id_v = VectorXd::Constant(nt,1);
    
    MatrixXd LAM=MatrixXd::Constant(s_n,nt,0);
    MatrixXd LAM_BAK=MatrixXd::Constant(s_n,nt,0);
    MatrixXd THETA=MatrixXd::Constant(s_n,nf,0);
    MatrixXd DELTA=MatrixXd::Constant(s_n,nf,0);
    VectorXd PHI = VectorXd::Constant(nf,0);
    VectorXd TAU = VectorXd::Constant(nf,0);
    double nu = 1;
    double ETA = 1;
    double GAMMA = 1;
    
    VectorXd count_lam = VectorXd::Constant(nf,0);
    VectorXd index = VectorXd::Constant(nf,0);
    
    
    MatrixXd LAM_TOP=MatrixXd::Constant(s_n,nf,0);
    MatrixXd LAM_BOT=MatrixXd::Constant(s_n,nf,0);
    
    double nmix=2;
    double zi = double(1)/nmix;
    MatrixXd Z = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logZ = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGV = MatrixXd::Constant(nmix,1,log(0.5));
    MatrixXd Zm = MatrixXd::Constant(nmix,nf,log(zi));

	long seed;
	gsl_rng *r;  // random number generator
	r=gsl_rng_alloc(gsl_rng_mt19937);
	
	seed = time (NULL) * getpid();    
	gsl_rng_set (r, seed);                  // set seed
	
	//gsl_rng_set (r, i_seed);
	
    // Initialize parameters 
    ID.diagonal()=id_v;
   
    PSI.diagonal() = psi_v;
    inv_psi(PSI,PSI_INV,s_n);

    // fill in the lambda matrix  
    for (int i=0; i<s_n; i++) {
        for(int j=0;j<nt;j++){
            LAM(i,j)=gsl_ran_gaussian(r,1);
        }
    }
    
    for(int i=0;i<s_n;i++){
        for(int j=0;j<nf;j++){
            THETA(i,j)=1;
            DELTA(i,j)=1;
        }
    }
    
    // continue initializing 
    for(int i=0;i<nf;i++){
        PHI(i)=1;
        TAU(i)=1;
    }
    VectorXd lam_count_v = VectorXd::Constant(n_itr,0);

    
    for(int itr=0;itr<(n_itr-1);itr++){
        // Expectation step 
        LP=LAM.transpose()*PSI_INV;
        VX=(LP*LAM+ID).inverse();
        EX=VX*LP*Y;
        EXX=EX*EX.transpose()+VX*d_y;
        TEX=EX.transpose();
        
        // Mazimization step 
        double alpha_sum=TAU.sum();
        ETA=double((d*nf+g-1))/(GAMMA+alpha_sum);
        GAMMA=double(g+h)/(ETA+nu);
        for(int i=0;i<nf;i++){
            TAU(i)=double(c+d)/(PHI(i)+ETA);
        }
        for(int i=0;i<nf;i++){
            double sum_c=s_n*b*Z(0,i)+c-1-0.5*s_n*Z(1,i);
            double at = 2*(TAU(i)+Z(0,i)*(DELTA.col(i).sum()));
            double lam_sum=LAM.col(i).dot(LAM.col(i));
            double bt = Z(1,i)*lam_sum;
            PHI(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
        }
        
        
        // count the number of loadings that are set to all zero
        count_lam.setZero();
        for(int i=0;i<nf;i++){
            for(int j=0;j<s_n;j++){
                if(LAM(j,i)!=0){
                    count_lam(i) +=  1;
                }
            }
        }
        // Count the number of loadings that are active, either all zero or PHI_k zero will kill 
        int count_nonzero = 0;
        for(int i=0;i<nf;i++){
            if(count_lam(i)!=0&&PHI(i)!=0){
                index(count_nonzero)=i;
                count_nonzero++;
            }
        }
        
        // remove inactive factors, loadings etc. and assign to new matrix 
        if(count_nonzero != nf){
            nf=count_nonzero;
            nt=nf;
            MatrixXd EX2=MatrixXd::Constant(nt,d_y,0);
            MatrixXd TEX2=MatrixXd::Constant(d_y,nt,0);
            MatrixXd VX2=MatrixXd::Constant(nt,nt,0);
            MatrixXd EXX2=MatrixXd::Constant(nt,nt,0);
            MatrixXd LP2=MatrixXd::Constant(nt,s_n,0);
            MatrixXd ID2=MatrixXd::Constant(nt,nt,0);
            VectorXd id_v2 = VectorXd::Constant(nt,1);
            MatrixXd LAM2=MatrixXd::Constant(s_n,nt,0);
            MatrixXd LAM_BAK2=MatrixXd::Constant(s_n,nt,0);
            MatrixXd THETA2=MatrixXd::Constant(s_n,nf,0);
            MatrixXd DELTA2=MatrixXd::Constant(s_n,nf,0);
            VectorXd PHI2 = VectorXd::Constant(nf,0);
            VectorXd TAU2 = VectorXd::Constant(nf,0);
            MatrixXd Z2 = MatrixXd::Constant(nmix,nf,zi);
            MatrixXd logZ2 = MatrixXd::Constant(nmix,nf,log(zi));
            VectorXd count_lam2 = VectorXd::Constant(nf,0);
            VectorXd index2 = VectorXd::Constant(nf,0);
            MatrixXd LAM_TOP2=MatrixXd::Constant(s_n,nt,0);
            MatrixXd LAM_BOT2=MatrixXd::Constant(s_n,nt,0);
            ID2.diagonal()=id_v2;
            
            for(int i=0;i<nf;i++){
                EX2.row(i)=EX.row(index(i));
                TEX2=EX2.transpose();
                for(int j=0;j<nf;j++){
                    VX2(i,j)=VX(index(i),index(j));
                    EXX2(i,j)=EXX(index(i),index(j));
                }
                LP2.row(i)=LP.row(index(i));
                LAM2.col(i)=LAM.col(index(i));
                LAM_BAK2.col(i)=LAM_BAK.col(index(i));
                THETA2.col(i)=THETA.col(index(i));
                DELTA2.col(i)=DELTA.col(index(i));
                PHI2(i)=PHI(index(i));
                TAU2(i)=TAU(index(i));
                Z2.col(i)=Z.col(index(i));
                logZ2.col(i)=logZ.col(index(i));
                LAM_TOP2.col(i)=LAM_TOP.col(index(i));
                LAM_BOT2.col(i)=LAM_BOT.col(index(i));
                count_lam2(i)=count_lam(index(i));
                index2(i)=index(i);
            }
            
            // Assign the new parameters back 
		    EX=EX2;
            TEX=TEX2;
            VX=VX2;
            EXX=EXX2;
            LP=LP2;
            ID=ID2;
            LAM=LAM2;
            LAM_BAK=LAM_BAK2;
            THETA=THETA2;
            DELTA=DELTA2;
            PHI=PHI2;
            TAU=TAU2;
            Z=Z2;
            logZ=logZ2;
            LAM_TOP=LAM_TOP2;
            LAM_BOT=LAM_BOT2;
            count_lam=count_lam2;
            index=index2;
            
        }
        
        
        // continue updating parameters 
        for(int i=0;i<s_n;i++){
            for(int j=0;j<nf;j++){
                DELTA(i,j)=double((a+b))/(THETA(i,j)+PHI(j));
            }
        }
        for(int i=0;i<s_n;i++){
            for(int j=0;j<nf;j++){
                double a23=(2*a-3);
                THETA(i,j)=double(a23+sqrt(a23*a23+8*LAM(i,j)*LAM(i,j)*DELTA(i,j)))/4/DELTA(i,j);
            }
        }

        
        // Need to update separately because 0/0=NA
        
        YLX=Y*TEX-LAM*EXX;
        for(int i=0;i<nf;i++){
            YLX=YLX+LAM.col(i)*(EXX.row(i));
            for(int j=0;j<s_n;j++){
                if(Z(0,i)==0){
                    LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*Z(1,i)/PHI(i));
                }
                else if(Z(1,i)==0){
                    LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*Z(0,i)/THETA(j,i));
                }
                else{
                    LAM(j,i) = YLX(j,i)/(EXX(i,i)+PSI(j,j)*(Z(1,i)/PHI(i)+Z(0,i)/THETA(j,i)));
                }
            }
            YLX=YLX-LAM.col(i)*(EXX.row(i));
        }
        
         // If theta=0, then LAM=0
        for(int i=0;i<s_n;i++){
            for(int j=0;j<nf;j++){
                if(THETA(i,j)==0){
                    LAM(i,j)=0;
                }
            }
        }
        
        // logZ 
        for(int i=0;i<nf;i++){
            logZ(0,i)=LOGV(0,0);
            logZ(1,i)=LOGV(1,0);
            for(int j=0;j<s_n;j++){
                logZ(0,i)=logZ(0,i)+log_norm(LAM(j,i),0,THETA(j,i))+log_gamma(THETA(j,i),a,DELTA(j,i))+log_gamma(DELTA(j,i),b,PHI(i));
                logZ(1,i)=logZ(1,i)+log_norm(LAM(j,i),0,PHI(i));
            }
        }
        
        // Z 
        for(int i=0;i<nf;i++){
            Z(0,i)=double(1)/(1+exp(logZ(1,i)-logZ(0,i)));
            Z(1,i)=1-Z(0,i);
        }
        
        double ps1=alpha;
        double ps2=beta;
        
        // Probability 
        for(int i=0;i<nf;i++){
            ps1=ps1+Z(0,i);
            ps2=ps2+Z(1,i);
        }
        double dgama = gsl_sf_psi(ps1+ps2);
        LOGV(0,0)=gsl_sf_psi(ps1)-dgama;
        LOGV(1,0)=gsl_sf_psi(ps2)-dgama;
        
        
        // PSI 
        LX=LAM*EX;
        for(int i=0;i<s_n;i++){
            PSI(i,i)=Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i));
        }
        PSI=PSI/d_y;
        inv_psi(PSI,PSI_INV,s_n);
        
        // LAM active 
        int lam_count=0;
        for(int i=0;i<s_n;i++){
            for(int j=0;j<nf;j++){
                if(LAM(i,j)!=0){
                    lam_count++;
                }
            }
        }
        
        lam_count_v(itr+1)=lam_count;
        int interval=200;
        
        // claim convergence if the number of values is stable for 200 iterations.
        if(itr>interval){
            if((itr%interval==0)|| (lam_count_v(itr+1)-lam_count_v(itr-interval)==0)){
                int lam_count=0;
                for(int i=0;i<s_n;i++){
                    for(int j=0;j<nf;j++){
                        if(LAM(i,j)!=0){
                            lam_count++;
                        }
                    }
                }
                cout << "Iteration " << itr << endl;

                ss.str("");
                ss.clear();
                ss << dir_out << "/itr";
                ofstream f_itr (ss.str().c_str());
                if (f_itr.is_open()){
                    f_itr << itr << endl;
                }
                f_itr.close();
    
                ss.str("");
                ss.clear();
                ss << dir_out << "/LAM";
                ofstream f_lam (ss.str().c_str());
                if (f_lam.is_open()){
		  f_lam << LAM.transpose() << endl;
                }
                f_lam.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/Z";
                ofstream f_Z (ss.str().c_str());
                if (f_Z.is_open()){
                    f_Z << Z << endl;
                }
                f_Z.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/LOGV";
                ofstream f_V (ss.str().c_str());
                if (f_V.is_open()){
                    f_V << LOGV << endl;
                }
                f_V.close();
                
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/EX";
                ofstream f_EX (ss.str().c_str());
                if (f_EX.is_open()){
		  f_EX << EX.transpose() << endl;
                }
                f_EX.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/EXX";
                ofstream f_EXX (ss.str().c_str());
                if (f_EXX.is_open()){
                    f_EXX << EXX << endl;
                }
                f_EXX.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/THETA";
                ofstream f_THETA (ss.str().c_str());
                if (f_THETA.is_open()){
		  f_THETA << THETA.transpose() << endl;
                }
                f_THETA.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/DELTA";
                ofstream f_DELTA (ss.str().c_str());
                if (f_DELTA.is_open()){
		  f_DELTA << DELTA.transpose() << endl;
                }
                f_DELTA.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/PHI";
                ofstream f_PHI (ss.str().c_str());
                if (f_PHI.is_open()){
                    f_PHI << PHI << endl;
                }
                f_PHI.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/TAU";
                ofstream f_TAU (ss.str().c_str());
                if (f_TAU.is_open()){
                    f_TAU << TAU << endl;
                }
                f_TAU.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/ETA";
                ofstream f_ETA (ss.str().c_str());
                if (f_ETA.is_open()){
                    f_ETA << ETA << endl;
                }
                f_ETA.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/GAMMA";
                ofstream f_GAMMA (ss.str().c_str());
                if (f_GAMMA.is_open()){
                    f_GAMMA << GAMMA << endl;
                }
                f_GAMMA.close();
                
                ss.str("");
                ss.clear();
                ss << dir_out << "/PSI";
                ofstream f_PSI (ss.str().c_str());
                if (f_PSI.is_open()){
                    f_PSI << PSI.diagonal() << endl;
                }
                f_PSI.close();
                
            }
        }
        
        if(itr>interval && lam_count_v(itr+1)-lam_count_v(itr-interval)==0){
            ss.str("");
            ss.clear();
            ss << dir_out << "/final";
            ofstream f_FINAL (ss.str().c_str());
            if (f_FINAL.is_open()){
                f_FINAL << "done" << endl;
            }
            f_FINAL.close();
            exit(0);
        }
        
    }
}
