#include <cmath>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

using namespace std;

template <typename T>
string to_string_with_precision(const T a_value, const int n = 0)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}

template <typename T>
void savePlots(T myPlot, string plotName, bool logY = false, float yMin = 0, float yMax = 1)
{
  TGaxis::SetMaxDigits(3);
  std::string plotPath = "plots/";
  string makeDirectory = "mkdir -p " + plotPath;
  system(makeDirectory.c_str());

  TCanvas *c1  = new TCanvas("myCanvas", "myCanvas",800,800);

  myPlot.GetYaxis()->SetRangeUser(yMin, yMax);

  if (strncmp(typeid(myPlot).name(), "4TH2F", 5) == 0)
  {
    myPlot.Draw("COLZ");
  }
  else
  {
    myPlot.Sumw2();
    if (logY) gPad->SetLogy();
    myPlot.Draw("PE1");
  }

  TPaveText *pt;
  pt = new TPaveText(0.15,0.9,0.95,1., "NDC");
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetTextFont(42);
  TText *pt_LaTex = pt->AddText("#it{#bf{sPHENIX}} Simulation, #sqrt{s} = 200 GeV, pp");
  pt->SetBorderSize(0);
  pt->Draw();
  gPad->Modified();

  string extensions[] = {".C", ".pdf", ".png", ".root"};
  for (auto extension : extensions)
  {
    string output = plotPath + "/" + plotName + extension;
    c1->SaveAs(output.c_str());
  }
}

TH1F makeHisto(int nBins, float min, float max, string xAxisTitle, string unit, int precision, string yAxisTitle = "Geo. Accept.")
{
  TH1F myHisto(xAxisTitle.c_str(), xAxisTitle.c_str(), nBins, min, max);
  
  if (unit != "") xAxisTitle += " [" + unit +  "]";  
  myHisto.GetXaxis()->SetTitle(xAxisTitle.c_str());

  float binWidth = (float) (max - min)/nBins;
  if (unit != "") yAxisTitle += " / " + to_string_with_precision(binWidth, precision) + " " + unit;
  myHisto.GetYaxis()->SetTitle(yAxisTitle.c_str());
  
  return myHisto;
}

TH1F makeHisto(int nBins, float* min, string xAxisTitle, string unit, int precision, string yAxisTitle = "Geo. Accept.")
{
  TH1F myHisto(xAxisTitle.c_str(), xAxisTitle.c_str(), nBins, min);

  if (unit != "") xAxisTitle += " [" + unit +  "]";
  myHisto.GetXaxis()->SetTitle(xAxisTitle.c_str());
  myHisto.GetYaxis()->SetTitle(yAxisTitle.c_str());

  return myHisto;
}

void calcGeomAccept()
{
  string fileName = "files/outputHFTrackEff_Lambda2ppi_CR_2_mode_pTref_3p0_extra_high_pT_points.root";
  TFile *file = new TFile(fileName.c_str());
  TTree* data = (TTree*)file->Get("HFTrackEfficiency");

  string fileName2 = "files/outputHFTrackEff_Kshort2pipi_CR_2_mode_pTref_3p0_extra_high_pT_points.root";
  TFile *file2 = new TFile(fileName2.c_str());
  TTree* data2 = (TTree*)file2->Get("HFTrackEfficiency");

  float Lambda0_mass; data->SetBranchAddress("true_mother_mass",&Lambda0_mass);
  float Lambda0_p; data->SetBranchAddress("true_mother_p",&Lambda0_p);
  float Lambda0_pT; data->SetBranchAddress("true_mother_pT",&Lambda0_pT);
  float Lambda0_eta; data->SetBranchAddress("true_mother_eta",&Lambda0_eta);
  float Lambda0_phi; data->SetBranchAddress("true_mother_phi",&Lambda0_phi);
  int Lambda0_track_1_silicon_seeds; data->SetBranchAddress("reco_track_1_silicon_seeds", &Lambda0_track_1_silicon_seeds);   
  int Lambda0_track_2_silicon_seeds; data->SetBranchAddress("reco_track_2_silicon_seeds", &Lambda0_track_2_silicon_seeds);   

  float K_S0_mass; data2->SetBranchAddress("true_mother_mass",&K_S0_mass);
  float K_S0_p; data2->SetBranchAddress("true_mother_p",&K_S0_p);
  float K_S0_pT; data2->SetBranchAddress("true_mother_pT",&K_S0_pT);
  float K_S0_eta; data2->SetBranchAddress("true_mother_eta",&K_S0_eta);
  float K_S0_phi; data2->SetBranchAddress("true_mother_phi",&K_S0_phi);
  int K_S0_track_1_silicon_seeds; data2->SetBranchAddress("reco_track_1_silicon_seeds", &K_S0_track_1_silicon_seeds);   
  int K_S0_track_2_silicon_seeds; data2->SetBranchAddress("reco_track_2_silicon_seeds", &K_S0_track_2_silicon_seeds);   

  float lower_bin_bounds[] = {0.5, 0.8, 1.1, 1.4, 1.8, 2.2, 3, 4};
  const unsigned int n_variable_bins = sizeof(lower_bin_bounds)/sizeof(lower_bin_bounds[0]) - 1; 

  TH1F Lambda0_all_pT = makeHisto(n_variable_bins, lower_bin_bounds, "pT", "GeV", 1);
  TH1F K_S0_all_pT = makeHisto(n_variable_bins, lower_bin_bounds, "pT", "GeV", 1);

  TH1F Lambda0_in_geometry_pT = makeHisto(n_variable_bins, lower_bin_bounds, "pT", "GeV", 1);
  TH1F K_S0_in_geometry_pT = makeHisto(n_variable_bins, lower_bin_bounds, "pT", "GeV", 1);

  TH1F final_ratio = makeHisto(n_variable_bins, lower_bin_bounds, "pT", "GeV", 1, "#Lambda^{0}/K_{S}^{0} Geo. Accept.");
  TH1F inv_final_ratio = makeHisto(n_variable_bins, lower_bin_bounds, "pT", "GeV", 1, "K_{S}^{0}/#Lambda^{0} Geo. Accept.");

  int n_bins = 15;
  float eta_min_max = 1.1;

  TH1F Lambda0_all_eta = makeHisto(n_bins, -1*eta_min_max, eta_min_max, "#eta", "", 1);
  TH1F K_S0_all_eta = makeHisto(n_bins, -1*eta_min_max, eta_min_max, "#eta", "", 1);

  TH1F Lambda0_in_geometry_eta = makeHisto(n_bins, -1*eta_min_max, eta_min_max, "#eta", "", 1);
  TH1F K_S0_in_geometry_eta = makeHisto(n_bins, -1*eta_min_max, eta_min_max, "#eta", "", 1);

  TH1F final_ratio_eta = makeHisto(n_bins, -1*eta_min_max, eta_min_max, "#eta", "", 1, "#Lambda^{0}/K_{S}^{0} Geo. Accept.");
  TH1F inv_final_ratio_eta = makeHisto(n_bins, -1*eta_min_max, eta_min_max, "#eta", "", 1, "K_{S}^{0}/#Lambda^{0} Geo. Accept.");

  TH1F Lambda0_all_phi = makeHisto(n_bins, -1*M_PI, M_PI, "#phi", "", 1);
  TH1F K_S0_all_phi = makeHisto(n_bins, -1*M_PI, M_PI, "#phi", "", 1);

  TH1F Lambda0_in_geometry_phi = makeHisto(n_bins, -1*M_PI, M_PI, "#phi", "", 1);
  TH1F K_S0_in_geometry_phi = makeHisto(n_bins, -1*M_PI, M_PI, "#phi", "", 1);

  TH1F final_ratio_phi = makeHisto(n_bins, -1*M_PI, M_PI, "#phi", "", 1, "#Lambda^{0}/K_{S}^{0} Geo. Accept.");
  TH1F inv_final_ratio_phi = makeHisto(n_bins, -1*M_PI, M_PI, "#phi", "", 1, "K_{S}^{0}/#Lambda^{0} Geo. Accept.");

  int tmp = 0;
  int barWidth = 50;

  //Lambda0 Geometric Acceptance
  int num_entries = data->GetEntries();
  for (int  l = 0; l < num_entries; ++l)
  {
    if (tmp != (int)100*l/num_entries)
    {
      tmp = (int)100*l/num_entries;
      if ((tmp%1)  == 0)
      {
        cout << "[";
        int pos = barWidth * tmp/100;
        for (int i = 0; i < barWidth; ++i)
        {
          if (i < pos) cout << "=";
          else if (i == pos) cout << ">";
          else cout << " ";
        }
        cout << "] " << tmp << " %\r";
        cout.flush();
      }
    }

    data->GetEntry(l);

    float pz = sqrt(pow(Lambda0_p, 2) - pow(Lambda0_pT, 2));
    float E = sqrt(pow(Lambda0_p, 2) + pow(Lambda0_mass, 2));
    float rapidity = 0.5*log((E + pz)/(E - pz));

    if (abs(rapidity) > 1.) continue;
    
    Lambda0_all_pT.Fill(Lambda0_pT);
    Lambda0_all_eta.Fill(Lambda0_eta);
    Lambda0_all_phi.Fill(Lambda0_phi);

    bool accepted = min(Lambda0_track_1_silicon_seeds, Lambda0_track_2_silicon_seeds) > 0;

    if (accepted)
    {
      Lambda0_in_geometry_pT.Fill(Lambda0_pT);
      Lambda0_in_geometry_eta.Fill(Lambda0_eta);
      Lambda0_in_geometry_phi.Fill(Lambda0_phi);
    }
  }

  tmp = 0; //cout<<"Creating Particles: 0%,"<<flush;
  barWidth = 50;

  //K_S0 Geometric Acceptance
  num_entries = data2->GetEntries();
  for (int  l = 0; l < num_entries; ++l)
  {
    if (tmp != (int)100*l/num_entries)
    {
      tmp = (int)100*l/num_entries;
      if ((tmp%1)  == 0)
      {
        cout << "[";
        int pos = barWidth * tmp/100;
        for (int i = 0; i < barWidth; ++i)
        {
          if (i < pos) cout << "=";
          else if (i == pos) cout << ">";
          else cout << " ";
        }
        cout << "] " << tmp << " %\r";
        cout.flush();
      }
    }

    data2->GetEntry(l);

    float pz = sqrt(pow(K_S0_p, 2) - pow(K_S0_pT, 2));
    float E = sqrt(pow(K_S0_p, 2) + pow(K_S0_mass, 2));
    float rapidity = 0.5*log((E + pz)/(E - pz));

    if (abs(rapidity) > 1.) continue;

    K_S0_all_pT.Fill(K_S0_pT);
    K_S0_all_eta.Fill(K_S0_eta);
    K_S0_all_phi.Fill(K_S0_phi);

    bool accepted = min(K_S0_track_1_silicon_seeds, K_S0_track_2_silicon_seeds) > 0;

    if (accepted)
    {
      K_S0_in_geometry_pT.Fill(K_S0_pT);
      K_S0_in_geometry_eta.Fill(K_S0_eta);
      K_S0_in_geometry_phi.Fill(K_S0_phi);
    }
  }

  cout << "[";
  int pos = barWidth * tmp;
  for (int i = 0; i < barWidth; ++i)
  {
    if (i < pos) cout << "=";
    else if (i == pos) cout << ">";
    else cout << " ";
  }
  cout << "] 100 %\r";
  cout.flush();
  cout<<endl;

  Lambda0_all_pT.Sumw2();
  K_S0_all_pT.Sumw2();

  Lambda0_in_geometry_pT.Sumw2();
  K_S0_in_geometry_pT.Sumw2();

  Lambda0_in_geometry_pT.Divide(&Lambda0_all_pT);
  K_S0_in_geometry_pT.Divide(&K_S0_all_pT);

  savePlots(K_S0_in_geometry_pT, "KS0_geometric_acceptance_ratio", false, 0, 0.035);
  savePlots(Lambda0_in_geometry_pT, "Lambda0_geometric_acceptance_ratio", false, 0, 0.035);

  final_ratio = Lambda0_in_geometry_pT; 
  final_ratio.Divide(&K_S0_in_geometry_pT);

  inv_final_ratio = K_S0_in_geometry_pT; 
  inv_final_ratio.Divide(&Lambda0_in_geometry_pT);

  savePlots(inv_final_ratio, "K_S0_to_Lambda0_geometric_acceptance_ratio", false, 0, 3);
  savePlots(final_ratio, "Lambda0_to_KS0_geometric_acceptance_ratio", false, 0, 1);

  Lambda0_all_eta.Sumw2();
  K_S0_all_eta.Sumw2();

  Lambda0_in_geometry_eta.Sumw2();
  K_S0_in_geometry_eta.Sumw2();

  Lambda0_in_geometry_eta.Divide(&Lambda0_all_eta);
  K_S0_in_geometry_eta.Divide(&K_S0_all_eta);

  savePlots(K_S0_in_geometry_eta, "KS0_geometric_acceptance_ratio_eta", false, 0, 0.035);
  savePlots(Lambda0_in_geometry_eta, "Lambda0_geometric_acceptance_ratio_eta", false, 0, 0.035);

  final_ratio_eta = Lambda0_in_geometry_eta; 
  final_ratio_eta.Divide(&K_S0_in_geometry_eta);

  inv_final_ratio_eta = K_S0_in_geometry_eta; 
  inv_final_ratio_eta.Divide(&Lambda0_in_geometry_eta);

  savePlots(inv_final_ratio_eta, "K_S0_to_Lambda0_geometric_acceptance_ratio_eta", false, 0, 3);
  savePlots(final_ratio_eta, "Lambda0_to_KS0_geometric_acceptance_ratio_eta", false, 0, 1);

  Lambda0_all_phi.Sumw2();
  K_S0_all_phi.Sumw2();

  Lambda0_in_geometry_phi.Sumw2();
  K_S0_in_geometry_phi.Sumw2();

  Lambda0_in_geometry_phi.Divide(&Lambda0_all_phi);
  K_S0_in_geometry_phi.Divide(&K_S0_all_phi);

  savePlots(K_S0_in_geometry_phi, "KS0_geometric_acceptance_ratio_phi", false, 0, 0.035);
  savePlots(Lambda0_in_geometry_phi, "Lambda0_geometric_acceptance_ratio_phi", false, 0, 0.035);

  final_ratio_phi = Lambda0_in_geometry_phi; 
  final_ratio_phi.Divide(&K_S0_in_geometry_phi);

  inv_final_ratio_phi = K_S0_in_geometry_phi; 
  inv_final_ratio_phi.Divide(&Lambda0_in_geometry_phi);

  savePlots(inv_final_ratio_phi, "K_S0_to_Lambda0_geometric_acceptance_ratio_phi", false, 0, 3);
  savePlots(final_ratio_phi, "Lambda0_to_KS0_geometric_acceptance_ratio_phi", false, 0, 1);

  //Geometric + efficiency ratio
  //pT
  std::cout << "| $p_{T}$ [GeV] | $K_S^0$ | $\\Lambda^0$ | $\\Lambda^0 / K_S^0$ | $K_S^0 / \\Lambda^0$ |" << std::endl;
  std::cout << "|:--:|:--:|:--:|:--:|:--:|" << std::endl;
  for (int i = 1; i <= n_variable_bins; ++i)
  {
    std::string low_pT = to_string_with_precision(final_ratio.GetXaxis()->GetBinLowEdge(i), 1);
    std::string high_pT = to_string_with_precision(final_ratio.GetXaxis()->GetBinUpEdge(i), 1);

    std::string Ks0_content = to_string_with_precision(K_S0_in_geometry_pT.GetBinContent(i), 4);
    std::string Ks0_error = to_string_with_precision(K_S0_in_geometry_pT.GetBinError(i), 4);

    std::string Lambda_content = to_string_with_precision(Lambda0_in_geometry_pT.GetBinContent(i), 4);
    std::string Lambda_error = to_string_with_precision(Lambda0_in_geometry_pT.GetBinError(i), 4);

    std::string content = to_string_with_precision(final_ratio.GetBinContent(i), 4);
    std::string error = to_string_with_precision(final_ratio.GetBinError(i), 4);

    std::string inv_content = to_string_with_precision(inv_final_ratio.GetBinContent(i), 4);
    std::string inv_error = to_string_with_precision(inv_final_ratio.GetBinError(i), 4);
    std::cout << "| " << low_pT << " $\\rightarrow$ " << high_pT << " | " << Ks0_content << " $\\pm$ " << Ks0_error << " | " << Lambda_content << " $\\pm$ " << Lambda_error << " | " << content << " $\\pm$ " << error << " | " << inv_content << " $\\pm$ " << inv_error << " |" << std::endl;
  } 

  //eta
  std::cout << "| $\\eta$ | $K_S^0$ | $\\Lambda^0$ | $\\Lambda^0 / K_S^0$ | $K_S^0 / \\Lambda^0$ |" << std::endl;
  std::cout << "|:--:|:--:|:--:|:--:|:--:|" << std::endl;
  for (int i = 1; i <= n_bins; ++i)
  {
    std::string low = to_string_with_precision(final_ratio_eta.GetXaxis()->GetBinLowEdge(i), 2);
    std::string high = to_string_with_precision(final_ratio_eta.GetXaxis()->GetBinUpEdge(i), 2);

    std::string Ks0_content = to_string_with_precision(K_S0_in_geometry_eta.GetBinContent(i), 4);
    std::string Ks0_error = to_string_with_precision(K_S0_in_geometry_eta.GetBinError(i), 4);

    std::string Lambda_content = to_string_with_precision(Lambda0_in_geometry_eta.GetBinContent(i), 4);
    std::string Lambda_error = to_string_with_precision(Lambda0_in_geometry_eta.GetBinError(i), 4);

    std::string content = to_string_with_precision(final_ratio_eta.GetBinContent(i), 2);
    std::string error = to_string_with_precision(final_ratio_eta.GetBinError(i), 2);

    std::string inv_content = to_string_with_precision(inv_final_ratio_eta.GetBinContent(i), 2);
    std::string inv_error = to_string_with_precision(inv_final_ratio_eta.GetBinError(i), 2);
    std::cout << "| " << low << " $\\rightarrow$ " << high << " | " << Ks0_content << " $\\pm$ " << Ks0_error << " | " << Lambda_content << " $\\pm$ " << Lambda_error << " | " << content << " $\\pm$ " << error << " | " << inv_content << " $\\pm$ " << inv_error << " |" << std::endl;
  } 

  //phi
  std::cout << "| $\\phi$ | $K_S^0$ | $\\Lambda^0$ | $\\Lambda^0 / K_S^0$ | $K_S^0 / \\Lambda^0$ |" << std::endl;
  std::cout << "|:--:|:--:|:--:|:--:|:--:|" << std::endl;
  for (int i = 1; i <= n_bins; ++i)
  {
    std::string low = to_string_with_precision(final_ratio_phi.GetXaxis()->GetBinLowEdge(i), 2);
    std::string high = to_string_with_precision(final_ratio_phi.GetXaxis()->GetBinUpEdge(i), 2);

    std::string Ks0_content = to_string_with_precision(K_S0_in_geometry_phi.GetBinContent(i), 4);
    std::string Ks0_error = to_string_with_precision(K_S0_in_geometry_phi.GetBinError(i), 4);

    std::string Lambda_content = to_string_with_precision(Lambda0_in_geometry_phi.GetBinContent(i), 4);
    std::string Lambda_error = to_string_with_precision(Lambda0_in_geometry_phi.GetBinError(i), 4);

    std::string content = to_string_with_precision(final_ratio_phi.GetBinContent(i), 2);
    std::string error = to_string_with_precision(final_ratio_phi.GetBinError(i), 2);

    std::string inv_content = to_string_with_precision(inv_final_ratio_phi.GetBinContent(i), 2);
    std::string inv_error = to_string_with_precision(inv_final_ratio_phi.GetBinError(i), 2);
    std::cout << "| " << low << " $\\rightarrow$ " << high << " | " << Ks0_content << " $\\pm$ " << Ks0_error << " | " << Lambda_content << " $\\pm$ " << Lambda_error << " | " << content << " $\\pm$ " << error << " | " << inv_content << " $\\pm$ " << inv_error << " |" << std::endl;
  } 

  file->Close();
  file2->Close();
}
