#include "stdio.h"
#include "stdlib.h"


char Patterns[2][2][128] = {{"V1_{1}+V1_{2}<0.8", "V1_{1}+V1_{2}>0.8"}, {"(V1_{1}-0.3)^{2}+(V1_{2}-0.3)^{2}<0.3", "(V1_{1}-0.3)^{2}+(V1_{2}-0.3)^{2}>0.3"}};


void V1PlotBin12(const char* fName, int idx, TLegend *l)
{
	TFile *f = TFile::Open(fName);
	TTree *t = (TTree*)f->Get("DataTTree");

	float var1[4];
	int   label;
	int   nentries = t->GetEntries();

	t->SetBranchAddress("var1" , var1);
	t->SetBranchAddress("label", &label);

	TH2F* h0 = new TH2F("var1Pattern0", "Patten of variable 1: bin1 vs. bin2", 20,0,1.0, 20,0,1.0);
	h0->GetXaxis()->SetTitle("bin 1 of variable 1 (V1_{1})");
	h0->GetYaxis()->SetTitle("bin 2 of variable 1 (V1_{2})");
	h0->SetLineColor(kBlue);
	h0->SetMarkerColor(kBlue);
	h0->SetMarkerStyle(20);
	h0->SetMarkerSize(0.2);
	TH2F* h1 = new TH2F("var1Pattern1", "Patten of variable 1: bin1 vs. bin2", 20,0,1.0, 20,0,1.0);
	h1->SetLineColor(kRed);
	h1->SetMarkerColor(kRed);
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(0.2);

	for(int i=0; i<nentries; i++){
		t->GetEntry(i);
		if(label == 0) { h0->Fill(var1[0], var1[1]); }
		if(label == 1) { h1->Fill(var1[0], var1[1]); }
	}

	h0->SetMarkerSize(0.7);
	h1->SetMarkerSize(0.7);

	h1->GetXaxis()->SetTitle("bin 1 of variable 1 (V1_{1})");
	h1->GetYaxis()->SetTitle("bin 2 of variable 1 (V1_{2})");
	h1->Draw("same");
	h0->Draw("same");

	//l->SetBorderSize(0);
	l->SetFillStyle(0);
	l->AddEntry(h1, Patterns[idx][0],"p");
	l->AddEntry(h0, Patterns[idx][1],"p");
	l->AddEntry((TObject*)0, "blurred boundary","");

	l->Draw("same");

}


void DataCheck(void)
{
	TCanvas *clinear = new TCanvas("LinearPattern", "LinearPattern",600,600);
	clinear->SetFillColor(0);

	TLegend* llinear = new TLegend(0.6,0.7,0.95,0.95);
	V1PlotBin12("TrainLinearData.root", 0, llinear);
	TLegend* lsubfigure_a = new TLegend(0.85,0.16,0.96,0.24);
	lsubfigure_a->AddEntry((TObject*)0, "(a)", "");
	lsubfigure_a->SetBorderSize(2);
	lsubfigure_a->Draw();

	clinear->SaveAs("LinearPattern.pdf");

/*
	TCanvas *cnonlinear = new TCanvas("NonlinearPattern", "NonlinearPattern",600,600);
	cnonlinear->SetFillColor(0);

	TLegend* lnonlinear = new TLegend(0.4,0.7,0.95,0.95);
	V1PlotBin12("SyntheticNonlinearData.root", 1, lnonlinear);
	TLegend* lsubfigure_b = new TLegend(0.85,0.16,0.96,0.24);
	lsubfigure_b->AddEntry((TObject*)0, "(b)", "");
	lsubfigure_b->SetBorderSize(2);
	lsubfigure_b->Draw();
*/
}
