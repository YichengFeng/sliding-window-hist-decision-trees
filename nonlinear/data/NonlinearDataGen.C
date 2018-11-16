#include "stdio.h"
#include "stdlib.h"
#include <vector>
using namespace std;

//Histogram variable 1 (V1) has 4 bins.

int NObs;

struct V1
{
	float bin[4];
	float cond; //condition: cond = sqrt((V1.bin1-0.3)^2 + (V1.bin2-0.3)^2);
	int   label;
};

bool sortV1ByCond(const V1 &lhs, const V1 &rhs) { return abs(lhs.cond-0.3) < abs(rhs.cond-0.3); }


void NonlinearPattenHist(const char* outName)
{
	//V1.cond = sqrt((V1.bin1-0.3)^2 + (V1.bin2-0.3)^2) < 0.3.
	vector<V1> var1(NObs);

	//randomly assign values
	for(int i=0; i<NObs; i++){
		float sum = 0;
		int   Nbin = sizeof(var1[i].bin)/sizeof(var1[i].bin[0]);
		for(int j=0; j<Nbin; j++){
			gRandom->SetSeed(0);
			var1[i].bin[j] = gRandom->Uniform(0, 1.0);
			sum += var1[i].bin[j];
		}
		for(int j=0; j<Nbin; j++) var1[i].bin[j] /= sum;
		var1[i].cond = sqrt(pow(var1[i].bin[0]-0.3,2) + pow(var1[i].bin[1]-0.3,2));
		if(var1[i].cond < 0.3){
			var1[i].label = 1;
		}else{
			var1[i].label = 0;
		}
	}


	//blur the boundary region

	int NBlur = (int)(0.25*NObs);

	sort(var1.begin(), var1.end(), sortV1ByCond);
	for(int i=0; i<NBlur; i++){
		gRandom->SetSeed(0);
		float rdm = gRandom->Uniform(0, 1.0);
		if(rdm < 0.1) var1[i].label = 1 - var1[i].label;
		//if(rdm < 0.) var1[i].label = 1 - var1[i].label;
	}


	//write data to a TTree
	TFile *f = new TFile(outName,"recreate");
	TTree *t = new TTree("DataTTree", "DataTTree");

	int tlabel, tvar1Nbin, tvar1label;
	float tvar1[4];
	t->Branch("label"    , &tlabel    , "label/I");
	t->Branch("var1Nbin" , &tvar1Nbin , "var1Nbin/I");
	t->Branch("var1label", &tvar1label, "var1label/I");
	t->Branch("var1"     , tvar1      , "var1[var1Nbin]/F");

	for(int i=0; i<NObs; i++){
		tvar1Nbin = sizeof(var1[i].bin)/sizeof(var1[i].bin[0]);
		tvar1label= var1[i].label;
		for(int j=0; j<tvar1Nbin; j++){ tvar1[j] = var1[i].bin[j]; }
		tlabel = var1[i].label;
		t->Fill();
	}

	t->Write();
	f->Write();
	f->Close();

}


void NonlinearDataGen(void)
{
	NObs = 2048;
	NonlinearPattenHist("TrainNonlinearData.root");

	sleep(1);
	NObs = 512;
	NonlinearPattenHist("PruneNonlinearData.root");

	sleep(1);
	NObs = 2048;
	NonlinearPattenHist("TestNonlinearData.root");

}

