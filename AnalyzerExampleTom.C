//         nicola@cern.ch
//         
//         
//         To use the code:
//         $> root -l
//         root [0] .L AnalyzerExample.C+
//         root [1] AnalyzerExample("InputFile.root")
//         
//         files are available at ( maintained by tommaso.isidori@cern.ch ):
//         https://cernbox.cern.ch/index.php/s/TNWfPTKfFGKBH2A
//         
// 
//           Comparator: starting from vector of time and data, thresh
// 
// 
// 
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <ctime>
#include <map>
#include "TString.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TFitResultPtr.h>
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TProfile.h"
#include "TPave.h"
#include "TLegend.h"




vector<float> Comparator(const std::vector<float>& time, const std::vector<float>& data, const float threshold, const float hysteresis=5e-3, const UInt_t begin=0, const UInt_t end=0) {
  vector<int> i_th_v;
  Int_t i_peakStart=0; 
  float firstCrossing_tmp;
  vector<float> firstCrossing;
  bool above=false;
  bool lockForHysteresis=false;
  TGraph rising;
  UInt_t spline_n=8;
  UInt_t start(spline_n), stop(data.size()-spline_n);
  if ( end != 0 ) {
    start = std::max( begin, spline_n );
    stop = std::min( (UInt_t) data.size()-spline_n, end );
  }
//   cout<< "Inside: From " << start << " to " << stop << endl;

  for (UInt_t i=start; i<stop; ++i) {
    // Look for first edge
    if ( !above && !lockForHysteresis && data.at(i)>threshold ) {
      firstCrossing_tmp = time.at(i);
      i_peakStart = i;
      above = true;
      lockForHysteresis=true;
    }
    // Lock until above threshold+hysteresis
    if ( above && lockForHysteresis && data.at(i)>threshold+hysteresis) {
      lockForHysteresis=false;
    }
    // Look for second edge
    if ( above && !lockForHysteresis && data.at(i)<threshold ) {
      i_th_v.emplace_back( i_peakStart );
      above = false;
    }
    
    if (above && i>=data.size()-spline_n-2) {
        i_th_v.emplace_back( i_peakStart );
        above = false;
        break;
      }
    if ( above && lockForHysteresis && data.at(i)<threshold &&  time.at(i)-firstCrossing_tmp>0.5e-9 ) {
      above = false;
      lockForHysteresis=false;
    }
  }
//   if (above) std::cout<<"Problem: finished above!"<<std::endl;
  //linear interpolation
  for ( const auto& i_th : i_th_v ) {
    for (UInt_t j=0; j<spline_n; ++j) {
      rising.SetPoint(j,data.at(i_th+j-spline_n/2),time.at(i_th+j-spline_n/2));
    }
    firstCrossing.emplace_back( rising.Eval(threshold,NULL,"S") );   
  }

  return firstCrossing;
}




vector<float> CFDMultiple(const std::vector<float>& time, const std::vector<float>& data, const float fraction, const float threshold, const float hysteresis=5e-3) {
  vector<int> i_th_v;
  Int_t i_peakStart=0; 
  float firstCrossing_tmp= .0;
  vector<float> firstCrossing;
  bool above=false;
  bool lockForHysteresis=false;
  TGraph rising;
  UInt_t spline_n=8;
  
  float localMax=0;

  for (UInt_t i=spline_n; i<data.size()-spline_n; ++i) {
    // Look for first edge
    if ( !above && !lockForHysteresis && data.at(i)>threshold ) {
      firstCrossing_tmp = time.at(i);
      i_peakStart = i;
      above = true;
      lockForHysteresis=true;
      localMax = data.at(i);
    }
    // Lock until above threshold+hysteresis
    if ( above && lockForHysteresis && data.at(i)>threshold+hysteresis) {
      lockForHysteresis=false;
      if ( localMax < data.at(i) ) localMax = data.at(i);
    }
    if ( above && data.at(i)>threshold+hysteresis) {      // Search for the maximum
      if ( localMax < data.at(i) ) localMax = data.at(i);
    }
    
    // Look for second edge
    if ( above && !lockForHysteresis && data.at(i)<threshold ) {
      if ( fraction*localMax > threshold ) {
        vector<float> toa_tmp = Comparator(time, data, fraction*localMax, hysteresis, i_peakStart-2, i+10);
        if ( toa_tmp.size() > 0 ) firstCrossing.emplace_back( toa_tmp.at(0) ); 
//         if ( toa_tmp.size() > 1 ) std::cout << "More than one peak in the CFD!" << std::endl;
      }
      above = false;
    }
    
    if (above && i>=data.size()-spline_n-2) {
        if (time.at(i)-firstCrossing_tmp>1e-9) {
          vector<float> toa_tmp = Comparator(time, data, fraction*localMax, hysteresis, i_peakStart-2, i+10);
          if ( toa_tmp.size() > 0 ) firstCrossing.emplace_back( toa_tmp.at(0) );
//           if ( toa_tmp.size() > 1 ) std::cout << "More than one peak in the CFD!" << std::endl;
        }
        above = false;
        break;
      }
    if ( above && lockForHysteresis && data.at(i)<threshold &&  time.at(i)-firstCrossing_tmp>1e-9 ) {
      above = false;
      lockForHysteresis=false;
    }
  }
//   if (above) std::cout<<"Problem: finished above!"<<std::endl;

  return firstCrossing;
}







vector<float> Falltime(const std::vector<float>& time, const std::vector<float>& data, const float fraction_1, const float fraction_2, const float threshold, const float hysteresis=5e-3) {
  vector<int> i_th_v;
  Int_t i_peakStart=0; 
  float firstCrossing_tmp= .0;
  vector<float> firstCrossing;


  bool above=false;
  bool lockForHysteresis=false;
  TGraph rising;
  UInt_t spline_n=8;
  float localMax=0;

  for (UInt_t i=spline_n; i<data.size()-spline_n; ++i) {
    // Look for first edge
    if ( !above && !lockForHysteresis && data.at(i)>threshold ) {
      firstCrossing_tmp = time.at(i);
      i_peakStart = i;
      above = true;
      lockForHysteresis=true;
      localMax = data.at(i);
    }
    // Lock until above threshold+hysteresis
    if ( above && lockForHysteresis && data.at(i)>threshold+hysteresis) {
      lockForHysteresis=false;
      if ( localMax < data.at(i) ) localMax = data.at(i);
    }
    if ( above && data.at(i)>threshold+hysteresis) {      // Search for the maximum
      if ( localMax < data.at(i) ) localMax = data.at(i);
    }
    
    // Look for second edge
    if ( above && !lockForHysteresis && data.at(i)<threshold ) {
      if ( fraction_1*localMax > threshold ) {
        vector<float> toa_tmp_20 = Comparator(time, data, fraction_1*localMax, hysteresis, i_peakStart-2, i+10);
        vector<float> toa_tmp_80 = Comparator(time, data, fraction_2*localMax, hysteresis, i_peakStart-2, i+10);
	if ( toa_tmp_80.size() == 1 && toa_tmp_20.size() == 1 ) firstCrossing.emplace_back( toa_tmp_80.at(0) - toa_tmp_20.at(0) ); 
    }
      above = false;
    }
    
    if (above && i>=data.size()-spline_n-2) {
        if (time.at(i)-firstCrossing_tmp>1e-9) {
          vector<float> toa_tmp = Comparator(time, data, fraction_1*localMax, hysteresis, i_peakStart-2, i+10);
//          if ( toa_tmp.size() > 0 ) firstCrossing.emplace_back( toa_tmp.at(0) );
//           if ( toa_tmp.size() > 1 ) std::cout << "More than one peak in the CFD!" << std::endl;
        }
        above = false;
        break;
      }
    if ( above && lockForHysteresis && data.at(i)<threshold &&  time.at(i)-firstCrossing_tmp>1e-9 ) {
      above = false;
      lockForHysteresis=false;
    }
  }
//   if (above) std::cout<<"Problem: finished above!"<<std::endl;

  return firstCrossing;
}








vector<float> Amplitude(const std::vector<float>& data, const float threshold, const float hysteresis=5e-3, const UInt_t begin=0, const UInt_t end=0) {
  Int_t i_peakStart=0; 
  bool above=false;
  bool lockForHysteresis=false;
  float localMax=0;
  vector<float> Maxampli;
  UInt_t start(0), stop(data.size());
  
  if ( end != 0 ) {
    start = begin;
    stop = std::min( (UInt_t) data.size(), end );
  }

  for (UInt_t i=start; i<stop; ++i) {
    // Look for first edge
    if ( !above && !lockForHysteresis && data.at(i)>threshold ) {
      i_peakStart = i;
      above = true;
      lockForHysteresis=true;
      if ( localMax < data.at(i) ) {
        localMax = data.at(i);
      }
      
    }
    // Lock until above threshold+hysteresis
    if ( above && lockForHysteresis && data.at(i)>threshold+hysteresis) {
      lockForHysteresis=false;
      if ( localMax < data.at(i) ) {
        localMax = data.at(i);
      }
    }
    
    // Look for second edge
    if ( above && !lockForHysteresis && data.at(i)<threshold ) {
      above = false;
       
      Maxampli.emplace_back(localMax);
      localMax = 0;
    }
    
    // Only if the peak is at the end of the vector data
    if (above && i>=data.size()-2) {
      above = false;
      Maxampli.emplace_back(localMax);

      break;
    }
      
    if ( above && lockForHysteresis && data.at(i)<threshold) {
      above = false;
      lockForHysteresis=false;
    }
    

  }

  return Maxampli;
}


            
//  Double_t polinomialFit(vector<Double_t> vec){
//   Double_t p0 =      1660.32;                      
//   Double_t p1 =     -15400.6;                     
//   Double_t p2 = -2.51458e+07;                      
//   Double_t p3 =  2.34461e+09;                      
//   Double_t p4 =  -8.8389e+10;                      
//   Double_t p5 =  1.49665e+12;                      
//   Double_t p6 = -9.02192e+12;   
// 
//      for(int i=0;i<vec.size();i++)
//          function = p0+p1*vec[i]+p2*vec[i]^2+p3*vec[i]^3+p4*vec[i]^4+p5*vec[i]^5+p6*vec[i]^6
//      return function;
// }            
//             
            
Int_t AnalyzerExampleTom(const std::string input_file_name = "Fullscan_6mev_150v_v1.root", const int single_trig=0, const float threshold = 0.03) {
  // Parameters
  float baseline_time = 800;    // the initial 800 ns are used to compute the baseline
  float hysteresis = 5e-3;    
//   float threshold = 0.06;
  float cfd_threshold = 0.5;
  int start_time=clock();
  // Opening input tree
  TFile input_file( input_file_name.c_str() );
  TTree *input_tree = (TTree*) input_file.Get("ttree");
  if (input_tree==NULL) return -1;
  int count=0;
  // Declaration of leaf types
  UInt_t          TrigNumber;
  Float_t         Aaxis;
  Float_t         Baxis;
  Float_t         Caxis;
  Float_t         cax_tmp = 0;
  vector<float>   *Time=NULL;
  vector<float>   thr_vec;
  vector<float>   deltay,derivative;
  vector<float> scansize_vec;
  vector<float>   *Voltage=NULL;
  vector<Double_t> rate_vec;
  vector<Long64_t> vec_tree; 
  vector<Int_t>    Caxis_vec;
  std::map<float, TH1F> RatesMap;
  std::map<float, TGraph> ThrMap;
    Double_t p0 =      1660.32;                      
  Double_t p1 =     -15400.6;                     
  Double_t p2 = -2.51458e+07;                      
  Double_t p3 =  2.34461e+09;                      
  Double_t p4 =  -8.8389e+10;                      
  Double_t p5 =  1.49665e+12;                      
  Double_t p6 = -9.02192e+12;  



  // Pointing branches
  input_tree->SetBranchAddress("TrigNumber", &TrigNumber);
  input_tree->SetBranchAddress("Aaxis", &Aaxis);
  input_tree->SetBranchAddress("Baxis", &Baxis);
  input_tree->SetBranchAddress("Caxis", &Caxis);
  input_tree->SetBranchAddress("Time", &Time);
  input_tree->SetBranchAddress("Voltage", &Voltage);
 
  // Creating output file
  std::string out_file = "out_";
  out_file += "th";
  out_file += std::to_string(threshold);
  out_file += "_";
  out_file += input_file_name;
  TFile output_file(out_file.c_str(),"RECREATE");
  
    TCanvas *c1 = new TCanvas("c1","c1");  
  // Creating histos and graphs
  TGraph CScan_counter_graph;
  CScan_counter_graph.SetName("Counter_vs_Caxis");
  CScan_counter_graph.SetTitle("Counter_vs_Caxis");
      
  TProfile electron_scan_profile;
  electron_scan_profile.SetName("electron_scan_profile");
  electron_scan_profile.SetTitle("Counting's profile on electrons peak");

  TProfile photon_scan_profile;
  photon_scan_profile.SetName("photon_scan_profile");
  photon_scan_profile.SetTitle("Peaks counting profile vs off_line Threshold");
  
  TProfile CScan_counter_profile; 
  CScan_counter_profile.SetName("CScan_counter_profile");
  CScan_counter_profile.SetTitle("Profile of event counter");
                TF1 *fa1 = new TF1("fa1","expo",0,4500);
                
  TGraph CScan_rms_graph;
  CScan_rms_graph.SetName("RMS_vs_Caxis");
  CScan_rms_graph.SetTitle("RMS_vs_Caxis");
  
  TH1F tOfA_h( "timeOfArrivals", "Time Of Arrivals", 5000, 1000, 6000 );
  tOfA_h.GetXaxis()->SetTitle("Time of arrival (ns)");

  TH1F tempo_differenza( "tempo_differenza", "tempo_differenza", 5000,0,4000);
  tempo_differenza.GetXaxis()->SetTitle("Time of arrival differences (ns)");  
  
  TH1F Fall_histo( "Fall time (20-80)", "Fall Time (20-80)", 500, 0, 60 );
  Fall_histo.GetXaxis()->SetTitle("Fall time (ns)");
  
  TH1F Fall_histo_electron( "Electrons peak fall time (20-80)", "Fall Time (20-80)", 500, 0, 50 );
  Fall_histo_electron.GetXaxis()->SetTitle("Fall time (ns)");
  
  TH1F Fall_histo_photon( "Photons peak fall time (20-80)", "Fall Time (20-80)", 500, 0, 50 );
  Fall_histo_photon.GetXaxis()->SetTitle("Fall time (ns)");

  
  TH1F tOfA_cfd_h( "timeOfArrivals_CFD", "Time Of Arrivals using CFD", 5000, 1000, 6000 );
  tOfA_cfd_h.GetXaxis()->SetTitle("Time of arrival (ns)");
 
  TH1F pulseheight( "Pulse Height", "Pulse Height", 300 , 0, 300 );
  pulseheight.GetXaxis()->SetTitle("Amplitude (mV)");
  
  TH2F electron_histo( "electron_histo", "electron_histo", 0.3/0.001, 0, 0.3,1000,0,2000);
  electron_histo.GetXaxis()->SetTitle("threshold [V]");
  electron_histo.GetYaxis()->SetTitle("Peaks");
    
  TGraph hysteresisScan_comparator;
  TGraph hysteresisScan_cfd;
  TGraph hysteresisScan_difference;


  TGraph thresholdScan_comparator;
  TGraph thresholdScan_cfd;
  TGraph thresholdScan_difference;
  
  TProfile tail_fall;
  
  TGraph Voltage_mon;	
  TGraph ratevscaxis;
    ratevscaxis.SetName("Rate vs Caxis");
    ratevscaxis.SetTitle("Rate vs Caxis");
    ratevscaxis.GetXaxis()->SetTitle("Caxis position");
    ratevscaxis.GetYaxis()->SetTitle("Rate [GHz]");
    
  TGraph ThrScanTot;
    ThrScanTot.SetName("Peaks vs threshold");
    ThrScanTot.SetTitle("Peaks vs threshold");
    ThrScanTot.GetXaxis()->SetTitle("threshold [V]");
    ThrScanTot.GetYaxis()->SetTitle("Peaks detected");
    
    TGraph ElThrScan;
    ElThrScan.SetName("El_peaks vs threshold");
    ElThrScan.SetTitle("El_peaks vs threshold");
    ElThrScan.GetXaxis()->SetTitle("threshold [V]");
    ElThrScan.GetYaxis()->SetTitle("Electrons peaks detected");
    
    TGraph PhotThrScan;
    PhotThrScan.SetName("Phot_peaks vs threshold");
    PhotThrScan.SetTitle("Phot_peaks vs threshold");
    PhotThrScan.GetXaxis()->SetTitle("threshold [V]");
    PhotThrScan.GetYaxis()->SetTitle("Photons peaks detected");
        
    
  // Cycle to read the tree
  for (Long64_t tree_i=0; tree_i<input_tree->GetEntries(); ++tree_i) {
    input_tree->GetEntry(tree_i);
    
    vec_tree.push_back(tree_i);
    // Compute the baseline 
    float baseline=.0;
    UInt_t baseline_nOfP=0;
    while (Time->at(baseline_nOfP) - Time->at(0) < baseline_time) {
      baseline += Voltage->at(baseline_nOfP);
      ++baseline_nOfP;
    }
    baseline = baseline/baseline_nOfP;
          
    // Convert in "positive signal" for convenience
    for (UInt_t i=0; i<Voltage->size(); ++i) {
      Voltage->at(i) = baseline-Voltage->at(i);
    }
    //Counting peaks
    vector<float> comparatorTimes = Comparator(*Time, *Voltage, threshold, hysteresis);
//     for(uint h = 0; h <comparatorTimes.size();h++) cout << comparatorTimes.at(h) <<endl; 
    
    //   cout<<TrigNumber<<'\t'<<Caxis<<'\t'<<comparatorTimes.size()<<endl;
    CScan_counter_graph.SetPoint( tree_i, Caxis, comparatorTimes.size() );
    
  
    
    
    
    
    // RMS
    double rms=0;
    for ( auto& v : *Voltage ) rms += std::pow(v,2);
    CScan_rms_graph.SetPoint( tree_i, Caxis, std::sqrt( rms / Voltage->size() ) );
    
    
    for ( const auto& toa : comparatorTimes ) tOfA_h.Fill( toa );
    
    // Ignore trigger with too many hits... TODO: compute baseline just before the hit
    //if ( comparatorTimes.size() > 100 ) continue;        
    
    // CFD
    vector<float> cfdTime = CFDMultiple(*Time, *Voltage, cfd_threshold, threshold, hysteresis);
    for ( const auto& toa : cfdTime ) tOfA_cfd_h.Fill( toa );
    
    
    

//RATEvsCaxis
vector<float> tempo = Comparator(*Time, *Voltage, threshold, hysteresis);
  if ( RatesMap.find( Caxis ) == RatesMap.end() ) { // se non trova Caxis nella mappa significa che non c'e', bisogna creare il TH1F
   std::string title("tempo_differenza_");
   title += std::to_string(Caxis);
   RatesMap[ Caxis ] = TH1F( title.c_str(), title.c_str(), 125,0,4000);
  }
  for(UInt_t i=0;i + 1<tempo.size();i++)
   RatesMap[ Caxis ].Fill( tempo.at(i+1)-tempo.at(i) );

        
        
  
//SINGLE TRIGGER MONITOR AND EVENT DISPLAY
    
    if (tree_i==single_trig){
        int j=0;
        
    for(float hyst=0; hyst<=200e-3;){
        vector<float> cfdTime = CFDMultiple(*Time, *Voltage, cfd_threshold, threshold, hysteresis=hyst);
        vector<float> comparatorTimes = Comparator(*Time, *Voltage, threshold, hysteresis=hyst);
                //cout <<"peaks "<<comparatorTimes.size()<<" - hyst "<<hyst<<endl;
        hysteresisScan_cfd.SetPoint(j,hyst,cfdTime.size());
        hysteresisScan_comparator.SetPoint(j,hyst,comparatorTimes.size());
            hysteresisScan_difference.SetPoint(j,hyst,comparatorTimes.size() - cfdTime.size());
        j++;
        hyst+=1e-3;
    }
    hysteresis=5e-3;
        
        j=0;
         for(float thr=0; thr<=0.3;){
            vector<float> cfdTime = CFDMultiple(*Time, *Voltage, cfd_threshold, thr, hysteresis);
            vector<float> comparatorTimes = Comparator(*Time, *Voltage, thr, hysteresis);
                    // cout <<"peaks "<<comparatorTimes.size()<<" - threshold "<<thr<<endl;
            thresholdScan_cfd.SetPoint(j,thr,cfdTime.size());
            thresholdScan_comparator.SetPoint(j,thr,comparatorTimes.size());
            thresholdScan_difference.SetPoint(j,thr,comparatorTimes.size() - cfdTime.size());
            thr_vec.push_back(thr);
            j++;
            
            thr+=0.001;
        }
    
        for(UInt_t t=0; t<Time->size();t++){
        Voltage_mon.SetPoint(t,Time->at(t),Voltage->at(t));
        }
 
         }
         
         
         

 Double_t thr_max = 0.3; 
 Double_t thr_step = 0.001;

 //Threshold scan for every trigger        
//           if (tree_i==0){
//  int j=0;
//  CScan_counter_profile.SetBins(thr_max/thr_step, 0, thr_max);
//     for(float thr=0; thr<=thr_max;){
//     vector<float> cfdTime_scan_0 = CFDMultiple(*Time, *Voltage, cfd_threshold, thr, hysteresis);
//     vector<float> comparatorTimes_scan_0 = Comparator(*Time, *Voltage, thr, hysteresis);
//     ThrScanTot.SetPoint(j,thr,comparatorTimes_scan_0.size());
//     CScan_counter_profile.Fill(thr,comparatorTimes_scan_0.size(),1);
//     j++;
//     thr+=thr_step;
//                                 }
//     }
//         
// else {
//     for(float thr=0; thr<=thr_max;){
//     vector<float> cfdTime_scan = CFDMultiple(*Time, *Voltage, cfd_threshold, thr, hysteresis);
//     vector<float> comparatorTimes_scan = Comparator(*Time, *Voltage, thr, hysteresis);
//     ThrScanTot.SetPoint(ThrScanTot.GetN(),thr,comparatorTimes_scan.size());
//     CScan_counter_profile.Fill(thr,comparatorTimes_scan.size(),1);
//     thr+=thr_step;
//                                 }
//      }
//      
 
 
 
 
 
 //TOTAL RISE TIME 
vector<float> Fall = Falltime(*Time, *Voltage,0.2,0.8, threshold, hysteresis);   
    for ( uint t=0; t<Fall.size();t++) Fall_histo.Fill(Fall.at(t)); 
 
     
 
 
//Threshold scan for the photons peak
 photon_scan_profile.SetBins(thr_max/thr_step, 0, thr_max);
if(Caxis==5000){      
 if (tree_i==0){
 int j=0;
    for(float thr=0; thr<=thr_max;){
    vector<float> cfdTime_scan_0 = CFDMultiple(*Time, *Voltage, cfd_threshold, thr, hysteresis);
    vector<float> comparatorTimes_scan_0 = Comparator(*Time, *Voltage, thr, hysteresis);
    PhotThrScan.SetPoint(j,thr,comparatorTimes_scan_0.size());
    photon_scan_profile.Fill(thr,comparatorTimes_scan_0.size(),1);
 
    j++;
    thr+=thr_step;}
                                
    vector<float> Fall_0_ph = Falltime(*Time, *Voltage,0.2,0.8, threshold, hysteresis);   
    	for ( uint t=0; t<Fall_0_ph.size();t++) {

    		Fall_histo_photon.Fill(Fall_0_ph.at(t));}
    }
        
else {
    for(float thr=0; thr<=thr_max;){
    vector<float> cfdTime_scan = CFDMultiple(*Time, *Voltage, cfd_threshold, thr, hysteresis);
    vector<float> comparatorTimes_scan = Comparator(*Time, *Voltage, thr, hysteresis);
    PhotThrScan.SetPoint(PhotThrScan.GetN(),thr,comparatorTimes_scan.size());
    photon_scan_profile.Fill(thr,comparatorTimes_scan.size(),1);

    thr+=thr_step;
                                }
    vector<float> Fall = Falltime(*Time, *Voltage,0.2,0.8, threshold, hysteresis);   
    	for ( uint t=0; t<Fall.size();t++) Fall_histo_photon.Fill(Fall.at(t)); 
     }
 }

 


//PULSE HEIGHT ON THE ELCTRON PEAK

if(Caxis > 20000 && Caxis < 30000){
    vector<float> amplitudes = Amplitude(*Voltage, 0.03, hysteresis);
    for (auto& ampl:amplitudes) 
      pulseheight.Fill(1000*ampl);

 }
    
    
    
//Threshold scan for the electron peak
 electron_scan_profile.SetBins(thr_max/thr_step, 0, thr_max);
 if(Caxis==25000){        
//Pulseheight
 if (tree_i==0){
   int j=0;
    for(float thr=0.03; thr<=0.1;){
      vector<float> cfdTime_scan_0 = CFDMultiple(*Time, *Voltage, cfd_threshold, thr, hysteresis);
      vector<float> comparatorTimes_scan_0 = Comparator(*Time, *Voltage, thr, hysteresis);
      electron_scan_profile.Fill(thr,comparatorTimes_scan_0.size(),1);
      ElThrScan.SetPoint(j,thr,comparatorTimes_scan_0.size());
      electron_histo.Fill(thr,comparatorTimes_scan_0.size());

      j++;
      thr+=thr_step;
                                }
        vector<float> Fall_0 = Falltime(*Time, *Voltage,0.2,0.8, threshold, hysteresis);   
    	for ( uint t=0; t<Fall_0.size();t++) Fall_histo_electron.Fill(Fall_0.at(t));
    }
        
else {    
    for(float thr=0; thr<=thr_max;){
      vector<float> cfdTime_scan = CFDMultiple(*Time, *Voltage, cfd_threshold, thr, hysteresis);
      vector<float> comparatorTimes_scan = Comparator(*Time, *Voltage, thr, hysteresis);
      electron_scan_profile.Fill(thr,comparatorTimes_scan.size(),1);
      ElThrScan.SetPoint(ElThrScan.GetN(),thr,comparatorTimes_scan.size());
      electron_histo.Fill(thr,comparatorTimes_scan.size());


//     cout<<"threshold = "<<thr<<"        peaks = "<<comparatorTimes_scan.size()<<endl;
      thr+=thr_step;
                                }
     }
//         tail_fall.SetBins(Time->size(), 0, Time->back());
	vector<float> Fall = Falltime(*Time, *Voltage,0.2,0.8, threshold, hysteresis);   
        
    	for ( uint t=0; t<Fall.size();t++) { Fall_histo_electron.Fill(Fall.at(t));
            
            if(( Fall.at(t) ) >= 10){
                for( uint j = t - 50 ; j < t + 50 ;j++) 
                    tail_fall.Fill(Time->at(j),Voltage->at(j),1); }
                                                                       }   
                                                                            
            
        
	count++;
 	}


    
    
  }
  
  
  
  
  
  
  
  
  
  
  
    
//                                                                           /////// QUI  ///////
    
    int graph_index=0;
  for (auto& histo:RatesMap ) {
   histo.second.Fit("fa1","Q"); //first e second si riferiscono agli elementi della mappa <float,TH1F>
   Double_t slope = fa1->GetParameter(1);
   if ( slope < 0 ) ratevscaxis.SetPoint(graph_index++, histo.first, abs(slope));  
  }

    
      
    
// Double_t fitf = polinomialFit(thr_vec);
//   fitf = p0+p1*vec[i]+p2*vec[i]^2+p3*vec[i]^3+p4*vec[i]^4+p5*vec[i]^5+p6*vec[i]^6
//   TF1 *func = new TF1("fit",fitf,0,0.5,6);
//   func->SetParameter(0,p0);
//   func->SetParameter(1,p1);
//   func->SetParameter(2,p2); 
//   func->SetParameter(3,p3);
//   func->SetParameter(4,p4);
//   func->SetParameter(5,p5);
//   func->SetParameter(6,p6);

  
  
  //NAME AND AXIS HYSTERESIS SCAN
   std::string voltage = "event monitor (trigger # " + std::to_string(single_trig) + ")";
  const char *volt = voltage.c_str();
  Voltage_mon.SetName(volt);
  Voltage_mon.SetTitle("Event monitor");
  Voltage_mon.GetXaxis()->SetTitle("Time");
  Voltage_mon.GetYaxis()->SetTitle("Amplitude (V)");
  
  std::string difference = "hysteresisScan_difference (trigger # " + std::to_string(single_trig) + ")";
  const char *diff = difference.c_str();
  hysteresisScan_difference.SetName(diff);
  hysteresisScan_difference.SetTitle("hysteresisScan_difference");
  hysteresisScan_difference.GetXaxis()->SetTitle("hysteresis (V)");
  hysteresisScan_difference.GetYaxis()->SetTitle("Detected peaks (Comparator - CFD)");
  
  std::string comp = "hysteresisScan_comparator (trigger # " + std::to_string(single_trig) + ")";
  const char *cmp = comp.c_str();
  hysteresisScan_comparator.SetName(cmp);
  hysteresisScan_comparator.SetTitle("hysteresisScan_comparator");
  hysteresisScan_comparator.GetXaxis()->SetTitle("hysteresis (V)");
  hysteresisScan_comparator.GetYaxis()->SetTitle("Detected peaks");

  std::string constfr = "hysteresisScan_cfd (trigger # " + std::to_string(single_trig) + ")";
  const char *cfd = constfr.c_str();
  hysteresisScan_cfd.SetName(cfd);
  hysteresisScan_cfd.SetTitle("hysteresisScan_cfd");
  hysteresisScan_cfd.GetXaxis()->SetTitle("hysteresis (V)");
  hysteresisScan_cfd.GetYaxis()->SetTitle("Detected peaks");  

  //NAME AND AXIS THRESHOLD SCAN
    std::string difference_thr = "thresholdScan_difference (trigger # " + std::to_string(single_trig) + ")";
  const char *diff_thr = difference_thr.c_str();
  thresholdScan_difference.SetName(diff_thr);
  thresholdScan_difference.SetTitle("thresholdScan_difference");
  thresholdScan_difference.GetXaxis()->SetTitle("threshold (V)");
  thresholdScan_difference.GetYaxis()->SetTitle("Detected peaks (Comparator - CFD)");
  std::string comp_thr = "thresholdScan_comparator (trigger # " + std::to_string(single_trig) + ")";
  const char *cmp_thr = comp_thr.c_str();
  thresholdScan_comparator.SetName(cmp_thr);
  thresholdScan_comparator.SetTitle("thresholdScan_comparator");
  thresholdScan_comparator.GetXaxis()->SetTitle("threshold (V)");
  thresholdScan_comparator.GetYaxis()->SetTitle("Detected peaks");

  std::string constfr_thr = "thresholdScan_cfd (trigger # " + std::to_string(single_trig) + ")";
  const char *cfd_thr = constfr_thr.c_str();
  thresholdScan_cfd.SetName(cfd_thr);
  thresholdScan_cfd.SetTitle("thresholdScan_cfd");
  thresholdScan_cfd.GetXaxis()->SetTitle("threshold (V)");
  thresholdScan_cfd.GetYaxis()->SetTitle("Detected peaks");  
  
  
  std::string tail_tmp = "Tail of Fall time distrib.";
  const char *tail = tail_tmp.c_str();
  tail_fall.SetName(tail);
  tail_fall.SetTitle("Tail of Fall time distrib.");
  tail_fall.GetXaxis()->SetTitle("Time (ns)");  
  tail_fall.GetYaxis()->SetTitle("Amplitude (V)");
  
  
  

  Voltage_mon.Write();
  thresholdScan_difference.Write();
  thresholdScan_cfd.Write();
  thresholdScan_comparator.Write();   
  hysteresisScan_difference.Write();
  hysteresisScan_cfd.Write();
  hysteresisScan_comparator.Write();  
  CScan_counter_graph.Write();
  CScan_rms_graph.Write();
  tOfA_h.Write();
  tOfA_cfd_h.Write();
  Fall_histo.Write();
  Fall_histo_electron.Write();
  Fall_histo_photon.Write();
  tempo_differenza.Write();
  tail_fall.Write();
  ratevscaxis.Write();
  pulseheight.Write();
  ElThrScan.Write();
  PhotThrScan.Write();
  electron_scan_profile.Write(); 
  photon_scan_profile.Write();
//   electron_histo.Write();
//   ThrScanTot.Write();
//   CScan_counter_profile.Write();
  
/*  
   auto legend = new TLegend(0.1,0.7,0.48,0.9);
  c1 ->cd(); 
   electron_scan_profile.Rebin(150); 
   photon_scan_profile.Rebin(150);
   photon_scan_profile.SetFillColor(kRed);
   photon_scan_profile.SetFillStyle(3003);
   photon_scan_profile.Draw();
   electron_scan_profile.SetFillColor(kRed);
   electron_scan_profile.SetFillStyle(3003);
   electron_scan_profile.Draw("SAME");   
   legend->SetHeader("Peaks counting profile vs off_line Threshold","C"); // option "C" allows to center the header
   legend->AddEntry("photon_scan_profile","Photons peak","lep");
   legend->AddEntry("electron_scan_profile","Electrons peak","lep");
   legend->Draw();
   c1->Update();*/
        
  
  
  delete input_tree;
  input_file.Close();
  output_file.Close();

  int stop_time=clock();
  cout << "time: " << (stop_time-start_time)/double(CLOCKS_PER_SEC)*1000 << endl;
  return 0;
};
