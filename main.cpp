#include<nr3.h>
#include<iostream>
#include<cmath>
#include<nr3.h>
#include<quadrature.h>
#include<roots.h>
#include<gamma.h>
#include<gauss_wgts.h>
#include<vector>
#include<fstream>
#include<ctime>
#include<mins.h>
#include<Amoeba.h>
#include<Ran.h>
#include "/Users/jers0730/Graduate School Work/Computational Projects/Source Code/Scythe2/Scythe_Double_Matrix.h"
#include "/Users/jers0730/Graduate School Work/Computational Projects/Source Code/Scythe2/Scythe_Optimize.h"
#include "/Users/jers0730/Graduate School Work/Computational Projects/Source Code/Scythe2/Scythe_Simulate.h"



using namespace SCYTHE;
using namespace std;

int Nvariables; // an indicator that tell us the number of variables you plan to use to characterize the states
int NvarWars; // an indicator that tells us the number of variables  you plan to use to charaterize the probability of a tie
int Constant; // an indicator that tells us whether or not to use a constant
int WhichData; // an indicator identifying which of the five datasets you want to use
int nsim; // an integer that tels us how many simulations we want to run
Doub wari; // an indicator telling us which war we are fighting
Doub r; // a variable for the ratio between log-likelihoods
Doub p; // a variable for the uniform random number that we draw
Doub sd_jump; // a variable that establishes the variance of the jumping distribution
int Estimation; // an indicator that simply tells us whether or no we are doing this unstructured or what
int nChains; //an indicator telling us how many chains to run
int dummy; //a system pause dummy
Doub pi = 3.141593; //

/************************************************
//Here is the normal pdf which will be our prior:
************************************************/
Doub normal_pdf(Doub x, Doub mu, Doub s2){
	return 1/(sqrt(2*pi*s2))*exp(-pow(x-mu,2)/(2*s2));
}

//this will be a function that computes the likelihood for any given vector of parameters!
//In general, I will need to change this depending on the actual explanatory variables that I include an whatnot!  This is set up so that I MUST include regressors in ability.
//It still might be interesting to do an unstructured estimation, so we'll think about that.  I suspect there is some clustering because my data is unbalanced, and that is causing trouble!
Doub llik(Matrix b, Matrix DataStates, Matrix DataWars, Matrix Y){
	

	Doub Priorlik = 0;
	for(int j=0;j<rows(b);j++){
		Priorlik = Priorlik + log(normal_pdf(b(j,0),0,20));
	}
	
	//Matrix SKILL(78,2);
	Matrix p1(78,1);
	Matrix p2(78,1);
	Matrix pt(78,1);
	
	for(int j=0;j<78;j++){
		
		//first we compute coalition skill given the parameters we have chosen...
		wari = DataWars(j,0);
		Matrix I(275,1);
	
		for(int i=0;i<275;i++){
			if(DataStates(i,0)==wari){I(i,0) = 1;}else{I(i,0) = 0;}
		}
		
		Matrix Wari;
		Wari = selif(DataStates,I);
		
		Matrix IA;
		IA = Wari(_,1);
		//IA.print();
		
		Matrix IB;
		IB = 1-IA;
		//IB.print();
		
		Matrix sideA;
		Matrix sideB;
		sideA = selif(Wari,IA);
		sideB = selif(Wari,IB);
		//sideA.print();
		//sideB.print();
		
		Matrix sideAskill;
		Matrix sideBskill;
		sideAskill = b(0,0)*sideA(_,(0+2));
		sideBskill = b(0,0)*sideB(_,(0+2));
		for(int i=1;i<(cols(DataStates)-2);i++){
			sideAskill = sideAskill + b(i,0)*sideA(_,(i+2));
			sideBskill = sideBskill + b(i,0)*sideB(_,(i+2));			
		}
		
		//exp(sideAskill).print();
		
		Doub suma = 0;
		Doub sumb = 0;
		for(int i=0;i<rows(sideAskill);i++){suma += exp(sideAskill)(i,0);}
		for(int i=0;i<rows(sideBskill);i++){sumb += exp(sideBskill)(i,0);}
		
		//Next we need to compute our theta for the probability of a tie:
		Doub etai = 0;
		int count = 1;
		for(int i=(cols(DataStates)-2);i<(cols(DataStates)-2+cols(DataWars)-1);i++){
			//cout << DataWars(j,count);
			etai = etai + b(i,0)*DataWars(j,count);
			count = count + 1;
			//cout << etai;
		}
		Doub thetai;
		thetai = exp(etai);
		//cout << thetai << endl;
		p1(j,0) = suma/(suma + sumb + thetai*sqrt(suma*sumb));
		p2(j,0) = sumb/(suma + sumb + thetai*sqrt(suma*sumb));
		pt(j,0) = (thetai*sqrt(suma*sumb))/(suma + sumb + thetai*sqrt(suma*sumb));
		
		//cout << suma << "  " << sumb << "  " << thetai << endl;
		
	}
	

	
	Doub lik = Priorlik;
	for(int i=0;i<78;i++){
		
		lik = lik + Y(i,2)*log(p1(i,0)) + Y(i,1)*log(p2(i,0)) + Y(i,0)*log(pt(i,0));
		
	}
	
	return lik;
	
}




Matrix MH(int nsim, Doub sd_jump,int Chain,Matrix DataStates,Matrix DataWars,Matrix Y,Matrix Variables,Matrix VariablesWars){


		/********************************************************************************
		//Now all I should need to do is to draw the draws, and I'll have an estimator!!!
		//initialize the matrix
		********************************************************************************/
		Matrix beta1((cols(DataStates)-2+cols(DataWars)-1),nsim);
		Matrix beta_candidate((cols(DataStates)-2+cols(DataWars)-1),1);
		Matrix beta_current((cols(DataStates)-2+cols(DataWars)-1),1);	
	
	

		//set the starting estimates for the parameters of the model:
		for(int i=0;i<(cols(DataStates)-2+cols(DataWars)-1);i++){beta1(i,0) = rnorm();}

		for(int i=1;i<nsim;i++){
		
			beta_current = beta1(_,(i-1));
			cout << "Chain " << Chain << ", iteration: " << i << endl;
			beta_current.print();
			for(int j=0;j<(cols(DataStates)-2+cols(DataWars)-1);j++){
				beta_candidate = beta_current;
				beta_candidate(j,0) = rnorm(1,1,beta_candidate(j,0),sd_jump)(0,0);
				//beta_candidate.print();
				//cout << llik(beta_current,DataStates,DataWars,Y) << endl;
				r = exp(llik(beta_candidate,DataStates,DataWars,Y)-llik(beta_current,DataStates,DataWars,Y));
				//cout << r << endl;
				p = runif();
				if(r>=p){beta1(j,i) = beta_candidate(j,0);}else{beta1(j,i) = beta_current(j,0);}
			
			}
		
		}
		//beta.print();
	
		return beta1;

}

Doub llikU(Matrix b_current,Matrix btheta_current,Matrix Data, Matrix DataWars, Matrix Y){

	//Matrix SKILL(78,2);
	Matrix p1(77,1);
	Matrix p2(77,1);
	Matrix pt(77,1);
	
	for(int j=0;j<77;j++){
		
		//first we compute coalition skill given the parameters we have chosen...
		wari = DataWars(j,0);
		Matrix I(272,1);
	
		for(int i=0;i<272;i++){
			if(Data(i,0)==wari){I(i,0) = 1;}else{I(i,0) = 0;}
		}
		
		Matrix Wari;
		Wari = selif(Data,I);
		//Wari.print();
		
		
		Matrix IA;
		IA = Wari(_,2);
		//IA.print();
		
		Matrix IB;
		IB = 1-IA;
		//IB.print();
		
		
		Matrix sideA;
		Matrix sideB;
		sideA = selif(Wari,IA);
		sideB = selif(Wari,IB);
		//sideA.print();
		//sideB.print();
		
		
		Matrix sideAskill(1,1);
		Matrix sideBskill(1,1);
		//cout << sideA(0,1) << endl;
		sideAskill(0,0) = exp(b_current(sideA(0,1)-1,0));
		sideBskill(0,0) = exp(b_current(sideB(0,1)-1,0));
		//sideAskill.print();
		//sideBskill.print();
		
		for(int i=1;i<(rows(sideA));i++){
			sideAskill = sideAskill + exp(b_current(sideA(i,1)-1,0));
			
		}
		//sideAskill.print();
		
		for(int i=1;i<(rows(sideB));i++){
			sideBskill = sideBskill + exp(b_current(sideB(i,1)-1,0));

		}		
		//sideBskill.print();		
		
		//Next we need to compute our theta for the probability of a tie:
		Doub etai = 0;
		int count = 1;
		for(int i=0;i<(cols(DataWars)-1);i++){
			//cout << DataWars(j,count);
			etai = etai + btheta_current(i,0)*DataWars(j,count);
			count = count + 1;
			//cout << etai;
		}
		Doub thetai;
		thetai = exp(etai);
		//cout << thetai << endl;
		
		p1(j,0) = sideAskill(0,0)/(sideAskill(0,0) + sideBskill(0,0) + thetai*sqrt(sideAskill(0,0)*sideBskill(0,0)));
		p2(j,0) = sideBskill(0,0)/(sideAskill(0,0) + sideBskill(0,0) + thetai*sqrt(sideAskill(0,0)*sideBskill(0,0)));
		pt(j,0) = (thetai*sqrt(sideAskill(0,0)*sideBskill(0,0)))/(sideAskill(0,0) + sideBskill(0,0) + thetai*sqrt(sideAskill(0,0)*sideBskill(0,0)));
		
		//cout << suma << "  " << sumb << "  " << thetai << endl;
		
	}

	//(p1+p2+pt).print();

	Doub lik = 0;
	for(int i=0;i<77;i++){
		
		lik = lik + Y(i,2)*log(p1(i,0)) + Y(i,1)*log(p2(i,0)) + Y(i,0)*log(pt(i,0));
		
	}

	return lik;

}


/******************************************************
//EXAMPLE CODE: How to draw random normals and uniforms
******************************************************/
	//rnorm(now many, 1,mean,variance)
	//Matrix N(10,1);
	//N = rnorm(10,1,0,1);
	//Matrix N(5,1);
	//cout << rnorm(0.0,1.0);
	//Matrix Adam(10,1);
	//Adam = rnorm(10,1,0,1);
	//Adam.print();
	
	//Adam = runif(10,1);
	//Adam.print();
	//N.print();
	

	
	//Matrix q(10,1);
	//q = rnorm(10,1,0,1);
	//q.print();
/*****************************************************/


int main (int argc, char *argv[]) {
	
	//Matrix y(1,1);
	//y(0,0) = rnorm(1,1,0,1)(0,0);
	//y.print();
	
	cout << "Enter 1 for an unstructured estimation, and anything else to pick variables for a structured estimation: ";
	cin >> Estimation;
	
	if(Estimation!=1){
	
		/**********************************************************
		//This is code for a structured Estimation using covariates
		/****************************************************************************************************
		//The first thing I'm doing is to go through the data and pick which variables I want to run this on:
		****************************************************************************************************/
		cout << "The column values for the Data (observations by participant) Matrix are:" << endl;
		cout << "0) War Number" << endl;
		cout << "1) Country Code" << endl;
		cout << "2) The year the war began" << endl;
		cout << "3) The year the war ended" << endl;
		cout << "4) Deaths in battle" << endl;
		cout << "5) Whether or not the state in question was on side A of the war" << endl;
		cout << "6) Whether the state in question was an initator" << endl;
		cout << "7) Pre war population for the state in question" << endl;
		cout << "8) Pre war armaments for the state in question" << endl;
		cout << "9) Iron and steel production (in units of 100000 tons...)" << endl;
		cout << "10) Military expenditure (in units of 10 billion pounds/dollars...)" << endl;
		cout << "11) Military personnel (in millions of men...)" << endl;
		cout << "12) Energy consumption (in billions of coal ton equivalents...)" << endl;
		cout << "13) Total population (in billions of people...)" << endl;
		cout << "14) CINC score" << endl;
		cout << endl << "In an attempt to facilitate empirical analysis of these and other historical trends, Polity IV includes constructed annual measures for both institutionalized democracy (DEMOC) and autocracy (AUTOC), as many polities exhibit mixed qualities of both of these distinct authority patterns. The measures are composite indices derived from the coded values of authority characteristic component variables (variables 3.1- 3.6, see below) according to the formulas, originally designed by Gurr, provided below." << endl;
		cout << endl << "15) The Polity IV measure of democracy.  Democracy is conceived as three essential, interdependent elements. One is the presence of institutions and procedures through which citizens can express effective preferences about alternative policies and leaders. Second is the existence of institutionalized constraints on the exercise of power by the executive. Third is the guarantee of civil liberties to all citizens in their daily lives and in acts of political participation. The operational indicator of democracy is derived from codings of the competitiveness of political participation (variable 2.6), the openness and competitiveness of executive recruitment (variables 2.3 and 2.2), and constraints on the chief executive (variable 2.4)." << endl;
		cout << endl << "16) The Polity IV measure of autocracy.  We use the more neutral term Autocracy and define it operationally in terms of the presence of a distinctive set of political characteristics. In mature form, autocracies sharply restrict or suppress competitive political participation. Their chief executives are chosen in a regularized process of selection within the political elite, and once in office they exercise power with few institutional constraints. Most modern autocracies also exercise a high degree of directiveness over social and economic activity, but we regard this as a function of political ideology and choice, not a defining property of autocracy." << endl;
		cout << endl <<"17) The old POLITY variable itself. The POLITY score is computed by subtracting the AUTOC score from the DEMOC score" << endl;
		cout << endl << "The Polity IV dataset contains three indicators of the structural characteristics by which chief executives are recruited: (1) the extent of institutionalization of executive transfers, XRREG; (2) the competitiveness of executive selection, XRCOMP; and (3) the openness of executive recruitment, XROPEN." << endl;
		cout << endl << "18) XRREG - Regulation of Chief Executive Recruitment: Regulation refers to the extent to which a polity has institutionalized procedures for transferring executive power. (page 19 of polity codebook)" << endl;
		cout << endl << "19) XRCOMP - Competitiveness of Executive Recruitment: Competitiveness refers to the extent that prevailing modes of advancement give subordinates equal opportunities to become superordinates. (see page 20 of the polity codebook)" << endl;
		cout << endl << "20) XROPEN - Openness of Executive Recruitment: Recruitment of the chief executive is open to the extent that all the politically active population has an opportunity, in principle, to attain the position through a regularized process.  NOT INDEPENDENT OF XRREG.  (see page 20 of polity codebook)" << endl;
		cout << endl << "21) XCONST - Executive Constraints (Decision Rules): Operationally, this variable refers to the extent of institutionalized constraints on the decisionmaking powers of chief executives, whether individuals or collectivities." << endl;
		cout << endl << "22) PARREG - Regulation of Participation: Participation is regulated to the extent that there are binding rules on when, whether, and how political preferences are expressed. One-party states and Western democracies both regulate participation but they do so in different ways, the former by channeling participation through a single party structure, with sharp limits on diversity of opinion; the latter by allowing relatively stable and enduring groups to compete nonviolently for political influence." << endl;
		cout << endl << "23) PARCOMP - The Competitiveness of Participation: The competitiveness of participation refers to the extent to which alternative preferences for policy and leadership can be pursued in the political arena. (NOT INDEPENDENT OF PARREG, see page 25 of polity codebook)" << endl;
		cout << endl << "24) EXREC* - Executive Recruitment: Concept variable combines information presented in three component variables: XRREG, XRCOMP, and XROPEN. (see polity codebook page 27.)" << endl;
		cout << endl << "25) EXCONST* - Executive Constraints: Concept variable is identical to XCONST (variable 3.4) above. See Addendum B for detailed explanations of the codes used. (see polity codebook page 28.)" << endl;
		cout << endl << "26) POLCOMP* - Political Competition: Concept variable combines information presented in two component variables: PARREG and PARCOMP. (see polity codebook page 28.)" << endl;
		cout << endl << "27) Alliance - the number of alliances each state had during the period of war"<< endl;
		cout << endl << "28) Polity*Initiator" << endl;
		cout << "29) Polity*Target" << endl;
		cout << "30) Polity*Military Expenditure" << endl;
		cout << "31) Polity*Military Personnel" << endl;
		cout << "32) Quality (money per guy...)" << endl;		
		cout << "33) The data from the first principle component*"<< endl;
		cout << "34) The data from the second principle component*"<< endl;
		cout << "35) The data from the third principle component"<< endl;
		cout << "36) The data from the fourth principle component"<< endl;
		cout << "37) The data from the fifth principle component"<< endl;
		cout << "38) relative Iron and steel production" << endl;
		cout << "39) relative Military expenditure" << endl;
		cout << "40) relative Military personnel" << endl;
		cout << "41) relative Energy consumption" << endl;
		cout << "42) relative Total population" << endl;
		cout << "43) relative urban population" << endl;		
		cout << "44) Exconst*Initiator" << endl;
		cout << "45) Exconst*Target" << endl;
		cout << "46) End: Iron and steel production (in units of 100000 tons...)" << endl;
		cout << "47) End: Military expenditure (in units of 10 billion pounds/dollars...)" << endl;
		cout << "48) End: Military personnel (in millions of men...)" << endl;
		cout << "49) End: Energy consumption (in billions of coal ton equivalents...)" << endl;
		cout << "50) End: Total population (in billions of people...)" << endl;
															
		
		cout << endl << endl << "How many explanatory variables do you want to use to parameterize ability?" << endl;
		cin >> Nvariables;
	
		Matrix Variables(Nvariables,1);
		for(int i=0;i<Nvariables;i++){
			cout << endl << endl << "Choose the column of the matrix that you want to use as a regressor (CHOOSE COLUMNS 0 AND 5 AS THE FIRST TWO): ";
			cin >> Variables(i,0);
			cout << Nvariables-i << " left...";
		}
		cout << "You chose: ";
		for(int i=0;i<Nvariables;i++){ cout << Variables(i,0) << " ";}
	
		cout << endl << endl << "Great!  Next, the column values for the Data by war are:" << endl;
		cout << "0) War Number"<< endl;
		cout << "1) The year the war began" << endl;
		cout << "2) The year the war ended" << endl;
		cout << "3) The country code for the largest state on one side..." << endl;
		cout << "4) The country code for the largest state on the other side..." << endl;	
		cout << "5) The duration of the war in months..." << endl;	
		cout << "6) The logged duration..." << endl;	
		cout << "7) An indicator of whether the war is censored (0) or not." << endl;	
		cout << "8) Strategy: offensive attrition, defense maneuver" << endl;
		cout << "9) Strategy: offensive attrition, defensive attrition" << endl;
		cout << "10) Strategy: offensive attrition, defensive punishment" << endl;
		cout << "11) Strategy: offensive manuever, defensive manuever" << endl;
		cout << "12) Strategy: offensive manuever, defensive attrition" << endl;
		cout << "13) Strategy: offensive manuever, defensive punishment" << endl;
		cout << "14) Strategy: offensive punishment, defensive manuever (this is always zero...)" << endl;
		cout << "15) Strategy: offensive punishment, defensive attrition" << endl;
		cout << "16) Strategy: offensive punishment, defensive punishment" << endl;
		cout << "17) Terrain coding, from 0 as completely open to 1 as impassible" << endl;
		//cout << "18) Whether or not the fighting took place in the Western Hemisphere" << endl;
		//cout << "19) Whether or not the fighting took place in Europe" << endl;
		//cout << "20) Whether or not the fighting took place in Africa" << endl;
		//cout << "21) Whether or not the fighting took place in the middle east" << endl;
		//cout << "22) Whether or not there was a land war in Asia..." << endl;
		//cout << "23) Whether or not the fighing took place in Oceania" << endl;

	
	
		cout << "How many explantory variables do you want to use to parameterize ties in war?" << endl;
		cin >> NvarWars;
	
		Matrix VariablesWars(NvarWars,1);	
		for(int i =0;i<NvarWars;i++){
			cout << endl << endl << "Choose the column of the matrix that you want to use as a regressor (CHOOSE COLUMN 0 AS THE FIRST): ";		
			cin >> VariablesWars(i,0);
			cout << NvarWars-i << "left...";
		}
	
		cout << "You chose: ";
		for(int i=0;i<NvarWars;i++){ cout << VariablesWars(i,0) << " ";}
	
		cout << endl << "Do you want to use a constant?  If so, enter 1: ";
		cin >> Constant;
	
		cout << endl << "How many simulations do we want to run?";
		cin >> nsim;
	
		cout << endl << "What do you want to set the jumping distribution at? ";
		cin >> sd_jump;
		
		cout << endl << "How many chains do you want to run?";
		cin >> nChains;
	
		cout << endl << "Finally, just tell me which data set you want to use (recall that we imputed five of them):";
		cin >> WhichData;
	
		cout << endl << "DONT FORGET TO CHANGE THE FILENAME TO SOMETHING APPROPRIATE!!!  Enter 0 to continue: ";
		cin >> dummy;
	
		/***************************************************************************************************************************************************************
		First thing to do is to read in the old Data - this includes data on wars, on duration of wars and strategy, regime characteristics, and finally on capabilities
		***************************************************************************************************************************************************************/
		Matrix Data;
		if(WhichData==1){
			Data = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Data1.txt",275,51);
			//Data1.print();
		}
		if(WhichData==2){
			Data = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Data2.txt",275,51);
			//Data2.print();
		}
		if(WhichData==3){
			Data = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Data3.txt",275,51);
			//Data3.print();
		}
		if(WhichData==4){
			Data = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Data4.txt",275,51);
			//Data4.print();
		}
		if(WhichData==5){
			Data = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Data5.txt",275,51);
			//Data5.print();
		}
		Data.print();
		
		cin >> Estimation;
		
		Matrix Y;
		Y = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Y.txt",78,3);
		//Y.print();

		Matrix Duration;
		Duration = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Duration.txt",78,18);
		//Duration.print();
	
		
		
		//for(int j=0;j<275;j++){cout << Data(j,0) << endl;}
		
		
		
		
		/******************************************************************************************************************************
		//Here we select the columns of the matrix that I am going to use to explain state ability in war and probabilities of a tie!!!
		******************************************************************************************************************************/
		Matrix DataStates(275,Nvariables);
		for(int j=0;j<275;j++){
			for(int i=0;i<Nvariables;i++){
				DataStates(j,i) = Data(j,Variables(i,0));
				//cout << Variables(i,0) << " ";
				//cout << DataStates(j,i) << "     ";
			}
			//cout << endl;
		}
		//DataStates.print();
		
		//cin >> Estimation;
	
		Matrix DataWars(78,NvarWars);
		for(int j=0;j<78;j++){
			for(int i=0;i<NvarWars;i++){
				DataWars(j,i) = Duration(j,VariablesWars(i,0));
				//cout << DataWars(j,i) << "     ";
			}
			//cout << endl;
		}	
		if(Constant==1){
		
			Matrix C(78,1);
			C = C(_,0) + 1;
			DataWars = cbind(DataWars,C);
			//DataWars.print();
		}
	
		
		
		
		ofstream out1("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Chainsmodelpolitics1.R");
		
		
		do{
			Matrix beta((cols(DataStates)-2+cols(DataWars)-1),nsim);
			beta = MH(nsim,sd_jump,nChains,DataStates,DataWars,Y,Variables,VariablesWars);
			

			out1 << "#BEGIN CHAIN!!!!!!" << endl;

			for(int j=0;j<(cols(DataStates)-2+cols(DataWars)-1)-1;j++){
		
				if(j<cols(DataStates)-2){
					for(int i=0;i<nsim;i++){
						if(i==0){out1 << "#This is state variable: " << Variables((j+2),0) << endl << "SVCH" << nChains << "." << Variables((j+2),0) << "<-c(" << beta(j,i) << ",";}
						if(i==(nsim-1)){out1 << beta(j,i) << ")";}
						if(i>0&i<(nsim-1)){out1 << beta(j,i) << ",";}
					}
				}else{
					for(int i=0;i<nsim;i++){
						if(i==0){out1 << "#This is war variable: " << VariablesWars(j-(cols(DataStates)-3),0) << endl << "WVCH1" << nChains << "." << VariablesWars(j-(cols(DataStates)-3),0) << "<-c(" << beta(j,i) << ",";}
						if(i==(nsim-1)){out1 << beta(j,i) << ")";}
						if(i>0&i<(nsim-1)){out1 << beta(j,i) << ",";}
					}
				}
				out1 << endl << endl << endl;
			}
	
			for(int i=0;i<nsim;i++){
				if(i==0){out1 << "#This is the Constant: " << endl << "ConstantCH" << nChains << "<-c(" << beta(((cols(DataStates)-2+cols(DataWars)-1)-1),i) << ",";}
				if(i==(nsim-1)){out1 << beta(((cols(DataStates)-2+cols(DataWars)-1)-1),i) << ")";}
				if(i>0&i<(nsim-1)){out1 << beta(((cols(DataStates)-2+cols(DataWars)-1)-1),i) << ",";}
			}
			out1 << endl << endl << endl;
			nChains = nChains-1;
		}while(nChains>0);




	}else{
	
		/****************************************************************
		//Here we have code for the unstructured estimation of the model!
		****************************************************************/
		
		cout << endl << endl << "Great!  The column values for the Data by war are:" << endl;
		cout << "0) War Number"<< endl;
		cout << "1) The year the war began" << endl;
		cout << "2) The year the war ended" << endl;
		cout << "3) The country code for the largest state on one side..." << endl;
		cout << "4) The country code for the largest state on the other side..." << endl;	
		cout << "5) The duration of the war in months..." << endl;	
		cout << "6) The logged duration..." << endl;	
		cout << "7) An indicator of whether the war is censored (0) or not." << endl;	
		cout << "8) Strategy: offensive attrition, defense maneuver" << endl;
		cout << "9) Strategy: offensive attrition, defensive attrition" << endl;
		cout << "10) Strategy: offensive attrition, defensive punishment" << endl;
		cout << "11) Strategy: offensive manuever, defensive manuever" << endl;
		cout << "12) Strategy: offensive manuever, defensive attrition" << endl;
		cout << "13) Strategy: offensive manuever, defensive punishment" << endl;
		cout << "14) Strategy: offensive punishment, defensive manuever" << endl;
		cout << "15) Strategy: offensive punishment, defensive attrition" << endl;
		cout << "16) Strategy: offensive punishment, defensive punishment" << endl;
		cout << "17) Terrain coding, from 0 as completely open to 1 as impassible" << endl;
	
		cout << "How many explantory variables do you want to use to parameterize ties in war?" << endl;
		cin >> NvarWars;
	
		Matrix VariablesWars(NvarWars,1);	
		for(int i =0;i<NvarWars;i++){
			cout << endl << endl << "Choose the column of the matrix that you want to use as a regressor (CHOOSE COLUMN 0 AS THE FIRST): ";		
			cin >> VariablesWars(i,0);
			cout << NvarWars-i << "left...";
		}
	
		cout << "You chose: ";
		for(int i=0;i<NvarWars;i++){ cout << VariablesWars(i,0) << " ";}
	
		cout << endl << "Do you want to use a constant?  If so, enter 1: ";
		cin >> Constant;
	
		cout << endl << "How many simulations do we want to run?";
		cin >> nsim;
	
		cout << endl << "What do you want to set the jumping distribution at? ";
		cin >> sd_jump;

		
		
		Matrix Data;
		Data = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Dataunstructured.txt",272,3);	
		
		Matrix Y;
		Y = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Yunstructured.txt",77,3);		

		Matrix Duration;
		Duration = read_into_matrix("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/Durationunstructured.txt",77,18);			
	
		Matrix DataWars(78,NvarWars);
		for(int j=0;j<78;j++){
			for(int i=0;i<NvarWars;i++){
				DataWars(j,i) = Duration(j,VariablesWars(i,0));
				//cout << DataWars(j,i) << "     ";
			}
			//cout << endl;
		}	
		if(Constant==1){
		
			Matrix C(78,1);
			C = C(_,0) + 1;
			DataWars = cbind(DataWars,C);
			//DataWars.print();
		}
	
		
		//Data.print();
		/**************************************************************************************************************************************************************
		//Here I set up the matrix of parameters; the first rows(unique(Data(_,1)))-1 are the state abilities, while the last NvarWars-1 are the coefficients for theta
		**************************************************************************************************************************************************************/
		Matrix States;
		States = unique(Data(_,1));
		
		//next thing I need is container matrices for the draws...
		Matrix BETA(rows(States),nsim);
		Matrix BETAtheta(cols(DataWars)-1,nsim);
		for(int i=0;i<rows(States);i++){
			BETA(i,0) = 10*rnorm(1,1,0,1)(0,0);
		}
		
		//I am choosing, as before, to set the United Kingdom as the baseline category
		BETA(16,0) = 0;
		for(int i=0;i<cols(DataWars)-1;i++){
			BETAtheta(i,0) = 10*rnorm(1,1,0,1)(0,0);
		}
		//BETA.print();
	
		//these are the current vectors of parameters!!!
		Matrix beta_current(max(States),1);
		Matrix betatheta_current(cols(DataWars)-1,1);
		Matrix beta_candidate(max(States),1);
		Matrix betatheta_candidate(cols(DataWars)-1,1);
	
		/****************************************************************
		//All thats really left to do is to run this and draw my draws!!!
		****************************************************************/
		for(int i=1;i<nsim;i++){
			
		
			for(int j=0;j<rows(States);j++){
				//cout << States(j,0) << " ";
				//cout << BETA(j,(i-1)) << endl;
				beta_current(States(j,0)-1,0) = BETA(j,(i-1));
			}
			//beta_current.print();
			
			for(int j=0;j<cols(DataWars)-1;j++){
				betatheta_current(j,0) = BETAtheta(j,(i-1));
			}
			//betatheta_current.print();
			
			//this updates the state specific ability parameters:
			for(int j=0;j<rows(States);j++){
				
				beta_candidate = beta_current;
			
			
				//This is what we change
				//including the US seems to cause problems...
				if(j!=0&j!=1&j!=2&j!=3&j!=4&j!=6&j!=7){beta_candidate(States(j,0)-1,0) = 0;}else{beta_candidate(States(j,0)-1,0) = rnorm(1,1,beta_candidate(States(j,0)-1,0),sd_jump)(0,0);}
				
				//cout << beta_candidate(States(j,0)-1,0) << " " << beta_current(States(j,0)-1,0) << endl;
				//cout << llikU(beta_candidate,betatheta_current,Data,DataWars,Y) << endl;
				//cout << llikU(beta_current,betatheta_current,Data,DataWars,Y) << endl;
				r = exp(llikU(beta_candidate,betatheta_current,Data,DataWars,Y)-llikU(beta_current,betatheta_current,Data,DataWars,Y));
				//cout << r << endl;
				p = runif();
				if(r>=p){
					BETA(j,i) = beta_candidate(States(j,0)-1,0);
					beta_current = beta_candidate;
				}else{BETA(j,i) = beta_current(States(j,0)-1,0);}
			
			}
		
			//next we need to update the covariate estimates for theta:
			
			for(int j=0;j<cols(DataWars)-1;j++){
			
				betatheta_candidate = betatheta_current;
				betatheta_candidate(j,0) = rnorm(1,1,betatheta_candidate(j,0),sd_jump)(0,0);
				//cout << betatheta_candidate(j,0) << " " << betatheta_current(j,0) << endl;
				r = exp(llikU(beta_current,betatheta_candidate,Data,DataWars,Y)-llikU(beta_current,betatheta_current,Data,DataWars,Y));				
				//cout << r << endl;
				p = runif();
				if(r>=p){
					BETAtheta(j,i) = betatheta_candidate(j,0);
					beta_current = beta_candidate;
				}else{BETAtheta(j,i) = betatheta_current(j,0);}				
				if(j==0){cout << "Iteration: " << i << ", ";}
			}
		

		}
	
		
		/**************************************************************************
		//Last thing to do: set this up to write this stuff out for post proc in R:
		**************************************************************************/
		ofstream out1("/Users/jers0730/Graduate School Work/571 Quantitative IR/571 Paper/R Stuff/ChainsU.R");
		for(int j=0;j<rows(States);j++){
		
				for(int i=0;i<nsim;i++){
					if(i==0){out1 << "#This is state ability parameter for: " << States(j,0) << endl << "S." << States(j,0) << "<-" << "c(" << BETA(j,i) << ",";}
					if(i==(nsim-1)){out1 << BETA(j,i) << ")";}
					if(i>0&i<(nsim-1)){out1 << BETA(j,i) << ",";}
				}
			
			out1 << endl << endl << endl;
		}
		
		if(0<cols(DataWars)-2){
			for(int j=0;j<cols(DataWars)-2;j++){
			
				for(int i=0;i<nsim;i++){
					if(i==0){out1 << "#This is the war parameter: " << VariablesWars(j+1,0) << endl << "W." << VariablesWars(j+1,0) << "<-" << "c(" << BETAtheta(j,i) << ",";}
					if(i==(nsim-1)){out1 << BETAtheta(j,i) << ")";}
					if(i>0&i<(nsim-1)){out1 << BETAtheta(j,i) << ",";}
				}
				out1 << endl << endl << endl;
			}
		}
		for(int i=0;i<nsim;i++){
			if(i==0){out1 << "#This is the constant: " << endl << "Constant <-" << "c(" << BETAtheta(cols(DataWars)-2,i) << ",";}
			if(i==(nsim-1)){out1 << BETAtheta(cols(DataWars)-2,i) << ")";}
			if(i>0&i<(nsim-1)){out1 << BETAtheta(cols(DataWars)-2,i) << ",";}
		}
		//BETAtheta.print();
		//cout << BETAtheta(cols(DataWars)-1,0);
		
		
		
		
		
		
	
	
	}

	
	return 0;
}
