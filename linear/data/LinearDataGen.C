#include "stdio.h"
#include "stdlib.h"
#include <vector>
using namespace std;

//Histogram variable 1 (V1) has 4 bins.
//Histogram variable 2 (V1) has 5 bins.

int NObs;

struct V1
{
	float bin[4];
	float cond; //condition: cond = bin1 + bin2;
	int   label;
};

struct V2
{
	float bin[5];
	float cond; //condition: cond = bin1 + bin2 + bin3;
	int   label;
};


bool sortV1ByCond(const V1 &lhs, const V1 &rhs) { return abs(lhs.cond-0.8) < abs(rhs.cond-0.8); }
bool sortV2ByCond(const V2 &lhs, const V2 &rhs) { return abs(lhs.cond-0.8) < abs(rhs.cond-0.8); }


void LinearPattenHist(const char* outName)
{
	//V1.cond = V1.bin1 + V1.bin2 < 0.8.
	//V2.cond = V2.bin1 + V2.bin2 + V2.bin3 < 0.8
	vector<V1> var1(NObs);
	vector<V2> var2(NObs);

	//randomly assign values
	for(int i=0; i<NObs; i++){
		float sum = 0;
		int   Nbin = sizeof(var1[i].bin)/sizeof(var1[i].bin[0]);
		for(int j=0; j<Nbin; j++){
			var1[i].bin[j] = gRandom->Uniform(0, 1.0);
			sum += var1[i].bin[j];
		}
		for(int j=0; j<Nbin; j++) var1[i].bin[j] /= sum;
		var1[i].cond = var1[i].bin[0] + var1[i].bin[1];
		if(var1[i].cond < 0.8){
			var1[i].label = 1;
		}else{
			var1[i].label = 0;
		}
	}

	for(int i=0; i<NObs; i++){
		float sum = 0;
		int   Nbin = sizeof(var2[i].bin)/sizeof(var2[i].bin[0]);
		for(int j=0; j<Nbin; j++){
			gRandom->SetSeed(0);
			var2[i].bin[j] = gRandom->Uniform(0, 1.0);
			sum += var2[i].bin[j];
		}
		for(int j=0; j<Nbin; j++) var2[i].bin[j] /= sum;
		var2[i].cond = var2[i].bin[0] + var2[i].bin[1] + var2[i].bin[2];
		if(var2[i].cond < 0.8){
			var2[i].label = 1;
		}else{
			var2[i].label = 0;
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

	sort(var2.begin(), var2.end(), sortV2ByCond);
	for(int i=0; i<NBlur; i++){
		gRandom->SetSeed(0);
		float rdm = gRandom->Uniform(0, 1.0);
		if(rdm < 0.1) var2[i].label = 1 - var2[i].label;
		//if(rdm < 0.) var2[i].label = 1 - var2[i].label;
	}


	//write data to a TTree
	TFile *f = new TFile(outName,"recreate");
	TTree *t = new TTree("DataTTree", "DataTTree");

	int tlabel, tvar1Nbin, tvar2Nbin, tvar1label, tvar2label;
	float tvar1[4], tvar2[5];
	t->Branch("label"    , &tlabel    , "label/I");
	t->Branch("var1Nbin" , &tvar1Nbin , "var1Nbin/I");
	t->Branch("var1label", &tvar1label, "var1label/I");
	t->Branch("var1"     , tvar1      , "var1[var1Nbin]/F");
	t->Branch("var2Nbin" , &tvar2Nbin , "var2Nbin/I");
	t->Branch("var2label", &tvar2label, "var2label/I");
	t->Branch("var2"     , tvar2      , "var2[var2Nbin]/F");

	for(int i=0; i<NObs; i++){
		tvar1Nbin = sizeof(var1[i].bin)/sizeof(var1[i].bin[0]);
		tvar2Nbin = sizeof(var2[i].bin)/sizeof(var2[i].bin[0]);
		tvar1label= var1[i].label;
		tvar2label= var2[i].label;
		for(int j=0; j<tvar1Nbin; j++){ tvar1[j] = var1[i].bin[j]; }
		for(int j=0; j<tvar2Nbin; j++){ tvar2[j] = var2[i].bin[j]; }
		if(var1[i].label == 1 && var2[i].label == 1){ tlabel = 1; }else{ tlabel = 0; }

		t->Fill();
	}

	t->Write();
	f->Write();
	f->Close();

}


void LinearDataGen(void)
{
	NObs = 2048;
	LinearPattenHist("TrainLinearData.root");

	sleep(1);
	NObs = 512;
	LinearPattenHist("PruneLinearData.root");

	sleep(1);
	NObs = 2048;
	LinearPattenHist("TestLinearData.root");
}

