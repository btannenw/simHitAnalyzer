#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TH1.h"
#include "TLatex.h"

#include "cmsStyle/tdrStyle.C"
#include "cmsStyle/CMS_lumi.C"

#include <iostream>
#include <iomanip>
#include <fstream>

string topDir;

void addOverflow(TH1D*& histo)
{
  // put overflow bin at end
  int maxBin = histo->GetNbinsX();
  histo->SetBinContent( maxBin, histo->GetBinContent( maxBin ) + histo->GetBinContent( maxBin+1 ) );
  histo->SetBinError  ( maxBin, sqrt( histo->GetBinError(maxBin)*histo->GetBinError(maxBin) + histo->GetBinError(maxBin+1)*histo->GetBinError(maxBin+1) ) );
  histo->SetBinContent( maxBin + 1, 0 );
  histo->SetBinError( maxBin + 1, 0 );

}

void drawFourSampleHist(TCanvas* c0, TFile* fMuon, TFile* fPion, TFile* fKaon, TFile* fProton, string hname, string xname, string modType="")
{
  std::cout << ( "gemSimAnalyzer/" + hname ).c_str() << std::endl;
  TH1D* h_muon   = (TH1D*)fMuon->Get( ( "genSimAnalyzer/" + hname ).c_str() );
  std::cout << ( "gemSimAnalyzer/" + hname ).c_str() << std::endl;
  TH1D* h_pion   = (TH1D*)fPion->Get( ( "genSimAnalyzer/" + hname ).c_str() );
  std::cout << ( "gemSimAnalyzer/" + hname ).c_str() << std::endl;
  TH1D* h_kaon   = (TH1D*)fKaon->Get( ( "genSimAnalyzer/" + hname ).c_str() );
  std::cout << ( "gemSimAnalyzer/" + hname ).c_str() << std::endl;
  TH1D* h_proton = (TH1D*)fProton->Get( ( "genSimAnalyzer/" + hname ).c_str() );
  std::cout << ( "gemSimAnalyzer/" + hname ).c_str() << std::endl;

  addOverflow(h_muon);
  std::cout <<  "mu flow" << std::endl;
  addOverflow(h_pion);
  std::cout <<  "pi flow" << std::endl;
  addOverflow(h_kaon);
  std::cout <<  "ka flow" << std::endl;
  addOverflow(h_proton);
  std::cout <<  "pro flow" << std::endl;
  
  c0->cd();
  c0->SetLogy();

  h_proton->SetMarkerColor(kBlack);
  h_proton->SetLineColor(kBlack);
  h_proton->SetLineWidth(2);
  h_muon->SetLineColor(kRed);
  h_muon->SetMarkerColor(kRed);
  h_muon->SetLineWidth(2);
  h_pion->SetLineColor(kBlue);
  h_pion->SetMarkerColor(kBlue);
  h_pion->SetLineWidth(2);
  h_kaon->SetLineColor(kGreen+2);
  h_kaon->SetMarkerColor(kGreen+2);
  h_kaon->SetLineWidth(2);


  h_proton->GetYaxis()->SetTitle("Entries / Bin");
  h_proton->GetXaxis()->SetTitle( xname.c_str() );

  h_proton->DrawNormalized("e");
  h_proton->GetYaxis()->SetRangeUser(0, 1.1);
  h_muon->DrawNormalized("e same");
  h_kaon->DrawNormalized("e same");
  h_pion->DrawNormalized("e same");

  TLegend* leg = new TLegend();
  leg = new TLegend(0.70, 0.20, .95, .40);
  if ( xname.find("Z-Coordinate")!= string::npos )   leg = new TLegend(0.70, 0.50, .95, .70);

  leg->AddEntry(h_proton, "Proton", "el");
  leg->AddEntry(h_muon, "Muon", "el");
  leg->AddEntry(h_pion, "Pion", "el");
  leg->AddEntry(h_kaon, "Kaon", "el");
  leg->Draw("same");

  TLatex ltx1;
  ltx1.SetTextAlign(9);
  ltx1.SetTextFont(62);
  ltx1.SetTextSize(0.025);
  ltx1.SetNDC();
  string globalRef = "";
  if( xname.find("Endcap") != string::npos && xname.find("X-Coordinate") != string::npos) globalRef = "Oriented Along Global #phi";
  if( xname.find("Endcap") != string::npos && xname.find("Y-Coordinate") != string::npos) globalRef = "Oriented Along Global R";
  if( xname.find("Endcap") != string::npos && xname.find("Z-Coordinate") != string::npos) globalRef = "Oriented Along Global Z";
  if( xname.find("Barrel") != string::npos && xname.find("X-Coordinate") != string::npos) globalRef = "Oriented Along Global #phi";
  if( xname.find("Barrel") != string::npos && xname.find("Y-Coordinate") != string::npos) globalRef = "Oriented Along Global Z";
  if( xname.find("Barrel") != string::npos && xname.find("Z-Coordinate") != string::npos) globalRef = "Oriented Along Global R";
  ltx1.DrawLatex(0.3, 0.92, globalRef.c_str());
  
  if (modType!="")
    ltx1.DrawLatex(0.6, 0.92, modType.c_str());




  CMS_lumi( c0, 0, 33);

  c0->SetLeftMargin(0.15);
  c0->SetRightMargin(0.05);
  c0->SetBottomMargin(0.15);

  c0->Print( (topDir + "/" + hname + ".png").c_str() );
}

void compareSamples(TCanvas* c0, TFile* fMuon, TFile* fPion, TFile* fKaon, TFile* fProton)
{
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "track_pt", "Sim Track p_{T} [GeV]");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_E_loss", "Barrel E_{loss} [MeV]");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localX", "Sim Hit X-Coordinate in Module Frame (Barrel) [mm]");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localY", "Sim Hit Y-Coordinate in Module Frame (Barrel) [mm]");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localZ", "Sim Hit Z-Coordinate in Module Frame (Barrel) [mm]");

  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "endcap_E_loss", "Endcap E_{loss} [MeV]");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "endcap_localX", "Sim Hit X-Coordinate in Module Frame (Endcap) [mm]");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "endcap_localY", "Sim Hit Y-Coordinate in Module Frame (Endcap) [mm]");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "endcap_localZ", "Sim Hit Z-Coordinate in Module Frame (Endcap) [mm]");

  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localX_modType1", "Sim Hit X-Coordinate in Module Frame (Barrel) [mm]", "Module Type 1");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localY_modType1", "Sim Hit Y-Coordinate in Module Frame (Barrel) [mm]", "Module Type 1");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localZ_modType1", "Sim Hit Z-Coordinate in Module Frame (Barrel) [mm]", "Module Type 1");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localX_modType2", "Sim Hit X-Coordinate in Module Frame (Barrel) [mm]", "Module Type 2");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localY_modType2", "Sim Hit Y-Coordinate in Module Frame (Barrel) [mm]", "Module Type 2");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localZ_modType2", "Sim Hit Z-Coordinate in Module Frame (Barrel) [mm]", "Module Type 2");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localX_modType3", "Sim Hit X-Coordinate in Module Frame (Barrel) [mm]", "Module Type 3");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localY_modType3", "Sim Hit Y-Coordinate in Module Frame (Barrel) [mm]", "Module Type 3");
  drawFourSampleHist(c0, fMuon, fPion, fKaon, fProton, "barrel_localZ_modType3", "Sim Hit Z-Coordinate in Module Frame (Barrel) [mm]", "Module Type 3");

}

void drawXYZ(TCanvas* c0, TH1D* h_x, TH1D* h_y, TH1D* h_z, string hname)
{
  h_x->SetYTitle("Entries / 0.5 mm");
  h_y->SetYTitle("Entries / 1.0 mm");
  h_z->SetYTitle("Entries / 0.003 mm");

  c0->Clear();

  c0->Divide(1, 3);
  c0->cd(1);
  h_x->Draw();
  c0->cd(2);
  h_y->Draw();
  c0->cd(3);
  TPad* p1 = (TPad*)(c0->cd(3));
  p1->SetLogy(); 
  h_z->Draw();
  //CMS_lumi( c0, 0, 33);
  c0->cd();

  TLatex ltx1;
  ltx1.SetTextAlign(9);
  ltx1.SetTextFont(62);
  ltx1.SetTextSize(0.025);
  ltx1.SetNDC();
  ltx1.DrawLatex(0.3, 0.82, ("ETL::" + hname).c_str());
  
  c0->SetLeftMargin(0.15);
  c0->SetRightMargin(0.05);
  c0->SetBottomMargin(0.15);

  c0->Print( (topDir + "/endcapLocal_" + hname + ".png").c_str() );
}

void makeXYZPlot(TCanvas* c0, TFile* fMuon)
{
  gDirectory->cd("genSimAnalyzer");
  //gDirectory->ls();
  TTree* t0 = (TTree*)gDirectory->Get("simHitTree");
  std::cout<< "Entries in tree: " << t0->GetEntries() << std::endl;
  TH1D* h_x_entry = new TH1D("h_x_entry", "h_x_entry", 100, -25, 25);
  TH1D* h_y_entry = new TH1D("h_y_entry", "h_y_entry", 100, -50, 50);
  TH1D* h_z_entry = new TH1D("h_z_entry", "h_z_entry", 100, -0.16, 0.16);

  TH1D* h_x_exit = new TH1D("h_x_exit", "h_x_exit", 100, -25, 25);
  TH1D* h_y_exit = new TH1D("h_y_exit", "h_y_exit", 100, -50, 50);
  TH1D* h_z_exit = new TH1D("h_z_exit", "h_z_exit", 100, -0.16, 0.16);

  TH1D* h_x_local = new TH1D("h_x_local", "h_x_local", 100, -25, 25);
  TH1D* h_y_local = new TH1D("h_y_local", "h_y_local", 100, -50, 50);
  TH1D* h_z_local = new TH1D("h_z_local", "h_z_local", 100, -0.16, 0.16);

  t0->Draw("endcap_localX_entry>>+h_x_entry");
  t0->Draw("endcap_localY_entry>>+h_y_entry");
  t0->Draw("endcap_localZ_entry>>+h_z_entry");

  t0->Draw("endcap_localX_exit>>+h_x_exit");
  t0->Draw("endcap_localY_exit>>+h_y_exit");
  t0->Draw("endcap_localZ_exit>>+h_z_exit");

  t0->Draw("endcap_localX>>+h_x_local");
  t0->Draw("endcap_localY>>+h_y_local");
  t0->Draw("endcap_localZ>>+h_z_local");

  h_x_entry->SetXTitle("Local X Entry Coordinates [mm]");
  h_x_exit->SetXTitle("Local X Exit Coordinates [mm]");
  h_x_local->SetXTitle("Local X Coordinates [mm]");

  h_y_entry->SetXTitle("Local Y Entry Coordinates [mm]");
  h_y_exit->SetXTitle("Local Y Exit Coordinates [mm]");
  h_y_local->SetXTitle("Local Y Coordinates [mm]");
  
  h_z_entry->SetXTitle("Local Z Entry Coordinates [mm]");
  h_z_exit->SetXTitle("Local Z Exit Coordinates [mm]");
  h_z_local->SetXTitle("Local Z Coordinates [mm]");

  drawXYZ(c0, h_x_entry, h_y_entry, h_z_entry, "EntryPoint");
  drawXYZ(c0, h_x_exit, h_y_exit, h_z_exit, "ExitPoint");
  drawXYZ(c0, h_x_local, h_y_local, h_z_local, "LocalPosition");

}

void compareZRegion_oneVariable(TCanvas*c0, TTree* t0, double nBins, double xMin, double xMax, string treeVariable, string xname, string yname)
{
  std::cout<< "Entries in tree: " << t0->GetEntries() << std::endl;
  
  TH1D* h_z_sensor = new TH1D( ("h_z_sensor_" + treeVariable).c_str(), ("h_z_sensor_" + treeVariable).c_str(), nBins, xMin, xMax);
  TH1D* h_z_gap = new TH1D( ("h_z_gap_" + treeVariable).c_str(), ("h_z_gap_" + treeVariable).c_str(), nBins, xMin, xMax);

  t0->Draw( ( treeVariable + ">>+h_z_sensor_" + treeVariable).c_str(), "abs(endcap_localZ_exit)>0.15" );
  t0->Draw( ( treeVariable + ">>+h_z_gap_" + treeVariable).c_str(), "abs(endcap_localZ_exit)<0.15" );

  
  c0->cd();
  //if(treeVariable.find("localZ")!=string::npos)   c0->SetLogy();

  h_z_sensor->SetLineColor(kBlack);
  h_z_gap->SetLineColor(kRed);

  h_z_sensor->SetXTitle( xname.c_str() );
  h_z_sensor->SetYTitle( yname.c_str() );

  h_z_sensor->DrawNormalized();
  h_z_sensor->GetYaxis()->SetRangeUser(0, h_z_sensor->GetMaximum()*1.3);
  h_z_sensor->DrawNormalized("e");
  h_z_gap->DrawNormalized("same e");

  TLegend* leg = new TLegend();
  leg = new TLegend(0.70, 0.65, .95, .85);
 
  leg->AddEntry(h_z_sensor, "abs(exitPoint.z()) > 0.15", "el");
  leg->AddEntry(h_z_gap, "abs(exitPoint.z()) < 0.15", "el");
  leg->Draw("same");

  c0->SetLeftMargin(0.15);
  c0->SetRightMargin(0.05);
  c0->SetBottomMargin(0.15);
  c0->Print( (topDir + "/endcapCompareZRegion_" + treeVariable + ".png").c_str() );

}

void compareZRegions(TCanvas* c0, TFile* fMuon)
{
  //gDirectory->cd("genSimAnalyzer");
  //gDirectory->ls();
  TTree* t0 = (TTree*)gDirectory->Get("simHitTree");
  std::cout<< "Entries in tree: " << t0->GetEntries() << std::endl;
 
  compareZRegion_oneVariable(c0, t0, 50, 0, 100, "endcap_absP", "Abs(P) [MeV]", "Entries / 2 MeV");
  compareZRegion_oneVariable(c0, t0, 50, 0, 1, "endcap_E_loss", "E_{Loss} [MeV]", "Entries / 20 KeV");
  compareZRegion_oneVariable(c0, t0, 50, -50, 50, "endcap_localX", "Local X Coordinates [mm]", "Entries / 2.0 mm");
  compareZRegion_oneVariable(c0, t0, 50, -25, 25, "endcap_localY", "Local Y Coordinates [mm]", "Entries / 1.0 mm");
  compareZRegion_oneVariable(c0, t0, 50, -50, 50, "endcap_localX_entry", "Local Entry X Coordinates [mm]", "Entries / 2.0 mm");
  compareZRegion_oneVariable(c0, t0, 50, -25, 25, "endcap_localY_entry", "Local Entry Y Coordinates [mm]", "Entries / 1.0 mm");
  compareZRegion_oneVariable(c0, t0, 50, -50, 50, "endcap_localX_exit", "Local Exit X Coordinates [mm]", "Entries / 2.0 mm");
  compareZRegion_oneVariable(c0, t0, 50, -25, 25, "endcap_localY_exit", "Local Exit Y Coordinates [mm]", "Entries / 1.0 mm");

  compareZRegion_oneVariable(c0, t0, 50, -0.16, 0.16, "endcap_localZ_entry", "Local Entry Z Coordinates [mm]", "Entries / 0.006 mm");
  compareZRegion_oneVariable(c0, t0, 50, -0.16, 0.16, "endcap_localZ", "Local Z Coordinates [mm]", "Entries / 0.006 mm");

  compareZRegion_oneVariable(c0, t0, 50, 8, 15, "endcap_tof", "Time Of Flight [?s]", "Entries / ?s");
}


void produceCombinedPlots()
{
  topDir = "endcapZStudies_12-04-18";

  // *** 0. Drawing options
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  // ** A. CMS Style
  setTDRStyle();
  writeExtraText = true;       // if extra text                                                
  //extraText  = "Internal";  // default extra text is "Preliminary"
  lumi_sqrtS = "#sqrt{s} = 13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)  
  cmsTextSize = 0.3;
  lumiTextSize = 0.3;

  TFile* fMuon = new TFile( (topDir + "/analyzer_singleMuon_zStudies.root").c_str() , "READ");
  //TFile* fKaon = new TFile((topDir + "/analyzer_singleKaon_modType.root").c_str() , "READ");
  //TFile* fPion = new TFile((topDir + "/analyzer_singlePion_modType.root").c_str() , "READ");
  //TFile* fProton = new TFile((topDir + "/analyzer_singleProton_modType.root").c_str() , "READ");

  // *** 1. Compare 4 different sample types
  //compareSamples(c1, fMuon, fPion, fKaon, fProton);

  // *** 2. Compare entry/exit/local points
  makeXYZPlot(c1, fMuon);

  // *** 3. Compare variables in Z sensor area and gap
  compareZRegions(c1, fMuon);

}
