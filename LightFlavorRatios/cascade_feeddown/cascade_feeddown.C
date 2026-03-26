#include <iostream>
#include <filesystem>
#include <string>

// Set the axis label max number of digits to 3
TGaxis::SetMaxDigits(3);

template <typename T> T *getObject(const std::string filename, const std::string objname)
{
    TFile *file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return nullptr;
    }
    T *obj = dynamic_cast<T *>(file->Get(objname.c_str()));
    if (!obj)
    {
        std::cerr << "Error getting object: " << objname << " from file: " << filename << std::endl;
        file->Close();
        return nullptr;
    }
    obj->SetDirectory(0);
    file->Close();
    return obj;
}

void draw1Dhistogram(TH1 *hist,                        //
                     bool normalize,                   //
                     bool logx,                        //
                     bool logy,                        //
                     std::string xtitle,               //
                     std::string ytitle,               //
                     std::vector<std::string> addinfo, //
                     std::string plotoption,           //
                     std::string filename              //
)
{
    if (!hist)
    {
        std::cerr << "draw1Dhistogram received a null histogram pointer" << std::endl;
        return;
    }

    if (normalize)
    {
        const double integral = hist->Integral();
        if (integral > 0)
        {
            hist->Scale(1.0 / integral);
        }
    }

    TCanvas *c = new TCanvas("c", "c", 800, 700);
    // gPad->SetRightMargin(0.11);
    gPad->SetLeftMargin(0.18);
    gPad->SetTopMargin(0.07);
    c->cd();
    c->SetLogx(logx);
    c->SetLogy(logy);
    hist->GetXaxis()->SetTitle(xtitle.c_str());
    hist->GetYaxis()->SetTitle(ytitle.c_str());
    hist->GetYaxis()->SetTitleOffset(1.6);
    hist->GetYaxis()->SetRangeUser((logy) ? hist->GetMinimum(0) * 0.5 : 0, (logy) ? hist->GetMaximum() * 25 : hist->GetMaximum() * 1.5);
    // if drawoption contains "PE"
    if (plotoption.find("PE") != std::string::npos)
    {
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.0);
        hist->SetLineWidth(2);
    }
    hist->Draw(plotoption.c_str());
    // Add additional information to the histogram title using TLatex
    TLatex *latex = new TLatex();
    latex->SetTextSize(0.04);
    latex->SetTextAlign(12);
    latex->SetNDC();
    for (size_t i = 0; i < addinfo.size(); ++i)
    {
        latex->DrawLatex(gPad->GetLeftMargin() + 0.04, 0.88 - i * 0.05, addinfo[i].c_str());
    }
    c->SaveAs(Form("%s.png", filename.c_str()));
    c->SaveAs(Form("%s.pdf", filename.c_str()));
    delete c;
}

void cascade_feeddown()
{
    const std::string base_dir = "/sphenix/tg/tg01/hf/hjheng/HF-analysis/simulation/Pythia_ppMinBias/"; // <-- change here to your directory

    std::vector<TH1 *> v_h_lambda_pt_all;
    std::vector<TH1 *> v_h_lambda_pt_from_xi_all;
    std::vector<TH1 *> v_h_lambda_pt_from_xi_charged;
    std::vector<TH1 *> v_h_lambda_pt_from_xi_neutral;
    std::vector<TH1 *> v_h_lambda_rapidity_all;
    std::vector<TH1 *> v_h_lambda_rapidity_from_xi_all;
    std::vector<TH1 *> v_h_lambda_rapidity_from_xi_charged;
    std::vector<TH1 *> v_h_lambda_rapidity_from_xi_neutral;
    std::vector<TH1 *> v_h_lambda_eta_all;
    std::vector<TH1 *> v_h_lambda_eta_from_xi_all;
    std::vector<TH1 *> v_h_lambda_eta_from_xi_charged;
    std::vector<TH1 *> v_h_lambda_eta_from_xi_neutral;
    std::vector<TH1 *> v_h_lambda_phi_all;
    std::vector<TH1 *> v_h_lambda_phi_from_xi_all;
    std::vector<TH1 *> v_h_lambda_phi_from_xi_charged;
    std::vector<TH1 *> v_h_lambda_phi_from_xi_neutral;
    std::vector<TH1 *> v_h_event_counter;

    for (const auto &entry : std::filesystem::directory_iterator(base_dir))
    {
        const std::string filename = entry.path().filename().string();
        if (entry.is_regular_file() && filename.rfind("ppMinBias_Xi_fraction_", 0) == 0 && entry.path().extension() == ".root")
        {
            std::cout << "Found file: " << entry.path().filename() << std::endl;
            TH1D *h_lambda_pt_all = getObject<TH1D>(entry.path().string(), "h_lambda_pt_all");
            TH1D *h_lambda_pt_from_xi_all = getObject<TH1D>(entry.path().string(), "h_lambda_pt_from_xi_all");
            TH1D *h_lambda_pt_from_xi_charged = getObject<TH1D>(entry.path().string(), "h_lambda_pt_from_xi_charged");
            TH1D *h_lambda_pt_from_xi_neutral = getObject<TH1D>(entry.path().string(), "h_lambda_pt_from_xi_neutral");
            TH1D *h_lambda_rapidity_all = getObject<TH1D>(entry.path().string(), "h_lambda_rapidity_all");
            TH1D *h_lambda_rapidity_from_xi_all = getObject<TH1D>(entry.path().string(), "h_lambda_rapidity_from_xi_all");
            TH1D *h_lambda_rapidity_from_xi_charged = getObject<TH1D>(entry.path().string(), "h_lambda_rapidity_from_xi_charged");
            TH1D *h_lambda_rapidity_from_xi_neutral = getObject<TH1D>(entry.path().string(), "h_lambda_rapidity_from_xi_neutral");
            TH1D *h_lambda_eta_all = getObject<TH1D>(entry.path().string(), "h_lambda_eta_all");
            TH1D *h_lambda_eta_from_xi_all = getObject<TH1D>(entry.path().string(), "h_lambda_eta_from_xi_all");
            TH1D *h_lambda_eta_from_xi_charged = getObject<TH1D>(entry.path().string(), "h_lambda_eta_from_xi_charged");
            TH1D *h_lambda_eta_from_xi_neutral = getObject<TH1D>(entry.path().string(), "h_lambda_eta_from_xi_neutral");
            TH1D *h_lambda_phi_all = getObject<TH1D>(entry.path().string(), "h_lambda_phi_all");
            TH1D *h_lambda_phi_from_xi_all = getObject<TH1D>(entry.path().string(), "h_lambda_phi_from_xi_all");
            TH1D *h_lambda_phi_from_xi_charged = getObject<TH1D>(entry.path().string(), "h_lambda_phi_from_xi_charged");
            TH1D *h_lambda_phi_from_xi_neutral = getObject<TH1D>(entry.path().string(), "h_lambda_phi_from_xi_neutral");
            TH1D *h_event_counter = getObject<TH1D>(entry.path().string(), "h_event_counter");

            if (!h_lambda_pt_all || !h_lambda_pt_from_xi_all || !h_lambda_pt_from_xi_charged || !h_lambda_pt_from_xi_neutral ||
                !h_lambda_rapidity_all || !h_lambda_rapidity_from_xi_all || !h_lambda_rapidity_from_xi_charged || !h_lambda_rapidity_from_xi_neutral ||
                !h_lambda_eta_all || !h_lambda_eta_from_xi_all || !h_lambda_eta_from_xi_charged || !h_lambda_eta_from_xi_neutral ||
                !h_lambda_phi_all || !h_lambda_phi_from_xi_all || !h_lambda_phi_from_xi_charged || !h_lambda_phi_from_xi_neutral ||
                !h_event_counter)
            {
                std::cerr << "Skipping file with missing histograms: " << filename << std::endl;
                delete h_lambda_pt_all;
                delete h_lambda_pt_from_xi_all;
                delete h_lambda_pt_from_xi_charged;
                delete h_lambda_pt_from_xi_neutral;
                delete h_lambda_rapidity_all;
                delete h_lambda_rapidity_from_xi_all;
                delete h_lambda_rapidity_from_xi_charged;
                delete h_lambda_rapidity_from_xi_neutral;
                delete h_lambda_eta_all;
                delete h_lambda_eta_from_xi_all;
                delete h_lambda_eta_from_xi_charged;
                delete h_lambda_eta_from_xi_neutral;
                delete h_lambda_phi_all;
                delete h_lambda_phi_from_xi_all;
                delete h_lambda_phi_from_xi_charged;
                delete h_lambda_phi_from_xi_neutral;
                delete h_event_counter;
                continue;
            }

            v_h_lambda_pt_all.push_back(h_lambda_pt_all);
            v_h_lambda_pt_from_xi_all.push_back(h_lambda_pt_from_xi_all);
            v_h_lambda_pt_from_xi_charged.push_back(h_lambda_pt_from_xi_charged);
            v_h_lambda_pt_from_xi_neutral.push_back(h_lambda_pt_from_xi_neutral);
            v_h_lambda_rapidity_all.push_back(h_lambda_rapidity_all);
            v_h_lambda_rapidity_from_xi_all.push_back(h_lambda_rapidity_from_xi_all);
            v_h_lambda_rapidity_from_xi_charged.push_back(h_lambda_rapidity_from_xi_charged);
            v_h_lambda_rapidity_from_xi_neutral.push_back(h_lambda_rapidity_from_xi_neutral);
            v_h_lambda_eta_all.push_back(h_lambda_eta_all);
            v_h_lambda_eta_from_xi_all.push_back(h_lambda_eta_from_xi_all);
            v_h_lambda_eta_from_xi_charged.push_back(h_lambda_eta_from_xi_charged);
            v_h_lambda_eta_from_xi_neutral.push_back(h_lambda_eta_from_xi_neutral);
            v_h_lambda_phi_all.push_back(h_lambda_phi_all);
            v_h_lambda_phi_from_xi_all.push_back(h_lambda_phi_from_xi_all);
            v_h_lambda_phi_from_xi_charged.push_back(h_lambda_phi_from_xi_charged);
            v_h_lambda_phi_from_xi_neutral.push_back(h_lambda_phi_from_xi_neutral);
            v_h_event_counter.push_back(h_event_counter);
        }
    }

    if (v_h_lambda_pt_all.empty())
    {
        std::cerr << "No valid ppMinBias_Xi_fraction_*.root inputs found in " << base_dir << std::endl;
        return;
    }

    // create histograms that are the sums of the corresponding histograms from all input files
    TH1D *h_lambda_pt_all = (TH1D *)v_h_lambda_pt_all[0]->Clone("h_lambda_pt_all");
    TH1D *h_lambda_pt_from_xi_all = (TH1D *)v_h_lambda_pt_from_xi_all[0]->Clone("h_lambda_pt_from_xi_all");
    TH1D *h_lambda_pt_from_xi_charged = (TH1D *)v_h_lambda_pt_from_xi_charged[0]->Clone("h_lambda_pt_from_xi_charged");
    TH1D *h_lambda_pt_from_xi_neutral = (TH1D *)v_h_lambda_pt_from_xi_neutral[0]->Clone("h_lambda_pt_from_xi_neutral");
    TH1D *h_lambda_rapidity_all = (TH1D *)v_h_lambda_rapidity_all[0]->Clone("h_lambda_rapidity_all");
    TH1D *h_lambda_rapidity_from_xi_all = (TH1D *)v_h_lambda_rapidity_from_xi_all[0]->Clone("h_lambda_rapidity_from_xi_all");
    TH1D *h_lambda_rapidity_from_xi_charged = (TH1D *)v_h_lambda_rapidity_from_xi_charged[0]->Clone("h_lambda_rapidity_from_xi_charged");
    TH1D *h_lambda_rapidity_from_xi_neutral = (TH1D *)v_h_lambda_rapidity_from_xi_neutral[0]->Clone("h_lambda_rapidity_from_xi_neutral");
    TH1D *h_lambda_eta_all = (TH1D *)v_h_lambda_eta_all[0]->Clone("h_lambda_eta_all");
    TH1D *h_lambda_eta_from_xi_all = (TH1D *)v_h_lambda_eta_from_xi_all[0]->Clone("h_lambda_eta_from_xi_all");
    TH1D *h_lambda_eta_from_xi_charged = (TH1D *)v_h_lambda_eta_from_xi_charged[0]->Clone("h_lambda_eta_from_xi_charged");
    TH1D *h_lambda_eta_from_xi_neutral = (TH1D *)v_h_lambda_eta_from_xi_neutral[0]->Clone("h_lambda_eta_from_xi_neutral");
    TH1D *h_lambda_phi_all = (TH1D *)v_h_lambda_phi_all[0]->Clone("h_lambda_phi_all");
    TH1D *h_lambda_phi_from_xi_all = (TH1D *)v_h_lambda_phi_from_xi_all[0]->Clone("h_lambda_phi_from_xi_all");
    TH1D *h_lambda_phi_from_xi_charged = (TH1D *)v_h_lambda_phi_from_xi_charged[0]->Clone("h_lambda_phi_from_xi_charged");
    TH1D *h_lambda_phi_from_xi_neutral = (TH1D *)v_h_lambda_phi_from_xi_neutral[0]->Clone("h_lambda_phi_from_xi_neutral");
    TH1D *h_event_counter = (TH1D *)v_h_event_counter[0]->Clone("h_event_counter");

    h_lambda_pt_all->Reset();
    h_lambda_pt_from_xi_all->Reset();
    h_lambda_pt_from_xi_charged->Reset();
    h_lambda_pt_from_xi_neutral->Reset();
    h_lambda_rapidity_all->Reset();
    h_lambda_rapidity_from_xi_all->Reset();
    h_lambda_rapidity_from_xi_charged->Reset();
    h_lambda_rapidity_from_xi_neutral->Reset();
    h_lambda_eta_all->Reset();
    h_lambda_eta_from_xi_all->Reset();
    h_lambda_eta_from_xi_charged->Reset();
    h_lambda_eta_from_xi_neutral->Reset();
    h_lambda_phi_all->Reset();
    h_lambda_phi_from_xi_all->Reset();
    h_lambda_phi_from_xi_charged->Reset();
    h_lambda_phi_from_xi_neutral->Reset();
    h_event_counter->Reset();

    for (size_t i = 0; i < v_h_lambda_pt_all.size(); ++i)
    {
        h_lambda_pt_all->Add(v_h_lambda_pt_all[i]);
        h_lambda_pt_from_xi_all->Add(v_h_lambda_pt_from_xi_all[i]);
        h_lambda_pt_from_xi_charged->Add(v_h_lambda_pt_from_xi_charged[i]);
        h_lambda_pt_from_xi_neutral->Add(v_h_lambda_pt_from_xi_neutral[i]);
        h_lambda_rapidity_all->Add(v_h_lambda_rapidity_all[i]);
        h_lambda_rapidity_from_xi_all->Add(v_h_lambda_rapidity_from_xi_all[i]);
        h_lambda_rapidity_from_xi_charged->Add(v_h_lambda_rapidity_from_xi_charged[i]);
        h_lambda_rapidity_from_xi_neutral->Add(v_h_lambda_rapidity_from_xi_neutral[i]);
        h_lambda_eta_all->Add(v_h_lambda_eta_all[i]);
        h_lambda_eta_from_xi_all->Add(v_h_lambda_eta_from_xi_all[i]);
        h_lambda_eta_from_xi_charged->Add(v_h_lambda_eta_from_xi_charged[i]);
        h_lambda_eta_from_xi_neutral->Add(v_h_lambda_eta_from_xi_neutral[i]);
        h_lambda_phi_all->Add(v_h_lambda_phi_all[i]);
        h_lambda_phi_from_xi_all->Add(v_h_lambda_phi_from_xi_all[i]);
        h_lambda_phi_from_xi_charged->Add(v_h_lambda_phi_from_xi_charged[i]);
        h_lambda_phi_from_xi_neutral->Add(v_h_lambda_phi_from_xi_neutral[i]);
        h_event_counter->Add(v_h_event_counter[i]);
    }

    TH1D *h_feeddown_frac_xi_all = (TH1D *)h_lambda_pt_from_xi_all->Clone("h_feeddown_frac_xi_all");
    h_feeddown_frac_xi_all->Divide(h_lambda_pt_from_xi_all, h_lambda_pt_all, 1.0, 1.0, "B");
    h_feeddown_frac_xi_all->SetTitle(";#Lambda^{0}(+c.c) p_{T} [GeV];f_{#Xi feeddown}");

    TH1D *h_feeddown_frac_xi_rapidity_all = (TH1D *)h_lambda_rapidity_from_xi_all->Clone("h_feeddown_frac_xi_rapidity_all");
    h_feeddown_frac_xi_rapidity_all->Divide(h_lambda_rapidity_from_xi_all, h_lambda_rapidity_all, 1.0, 1.0, "B");
    h_feeddown_frac_xi_rapidity_all->SetTitle(";#Lambda^{0}(+c.c) rapidity;f_{#Xi feeddown}");

    TH1D *h_feeddown_frac_xi_eta_all = (TH1D *)h_lambda_eta_from_xi_all->Clone("h_feeddown_frac_xi_eta_all");
    h_feeddown_frac_xi_eta_all->Divide(h_lambda_eta_from_xi_all, h_lambda_eta_all, 1.0, 1.0, "B");
    h_feeddown_frac_xi_eta_all->SetTitle(";#Lambda^{0}(+c.c) #eta;f_{#Xi feeddown}");

    TH1D *h_feeddown_frac_xi_phi_all = (TH1D *)h_lambda_phi_from_xi_all->Clone("h_feeddown_frac_xi_phi_all");
    h_feeddown_frac_xi_phi_all->Divide(h_lambda_phi_from_xi_all, h_lambda_phi_all, 1.0, 1.0, "B");
    h_feeddown_frac_xi_phi_all->SetTitle(";#Lambda^{0}(+c.c) #phi;f_{#Xi feeddown}");

    std::string plotdir = base_dir + "/cascade_feeddown/figure/"; // <-- change here if you want to save plots in a different directory
    system(("mkdir -p " + plotdir).c_str());
    std::string cutinfo = "#Lambda^{0} |#eta| < 1.3; decay length #leq 2.5 cm";
    draw1Dhistogram(h_feeddown_frac_xi_all, false, false, false, "#Lambda^{0}(+c.c) p_{T} [GeV]", "f_{#Xi feeddown}=N_{#Lambda^{0}(+c.c)}^{#Xi feedown}/N_{#Lambda^{0}(+c.c)}^{Inclusive}", {"PYTHIA8 p+p minimum bias (w. Detroit tune)", cutinfo}, "PE1", plotdir + "feeddown_fraction_xi_all");
    draw1Dhistogram(h_feeddown_frac_xi_rapidity_all, false, false, false, "#Lambda^{0}(+c.c) rapidity", "f_{#Xi feeddown}=N_{#Lambda^{0}(+c.c)}^{#Xi feedown}/N_{#Lambda^{0}(+c.c)}^{Inclusive}", {"PYTHIA8 p+p minimum bias (w. Detroit tune)", cutinfo}, "PE1", plotdir + "feeddown_fraction_xi_rapidity_all");
    draw1Dhistogram(h_feeddown_frac_xi_eta_all, false, false, false, "#Lambda^{0}(+c.c) #eta", "f_{#Xi feeddown}=N_{#Lambda^{0}(+c.c)}^{#Xi feedown}/N_{#Lambda^{0}(+c.c)}^{Inclusive}", {"PYTHIA8 p+p minimum bias (w. Detroit tune)", cutinfo}, "PE1", plotdir + "feeddown_fraction_xi_eta_all");
    draw1Dhistogram(h_feeddown_frac_xi_phi_all, false, false, false, "#Lambda^{0}(+c.c) #phi", "f_{#Xi feeddown}=N_{#Lambda^{0}(+c.c)}^{#Xi feedown}/N_{#Lambda^{0}(+c.c)}^{Inclusive}", {"PYTHIA8 p+p minimum bias (w. Detroit tune)", cutinfo}, "PE1", plotdir + "feeddown_fraction_xi_phi_all");

    TFile *fout = new TFile(Form("%s/cascade_feeddown/Cascade_feeddown_fraction.root", base_dir.c_str()), "RECREATE"); // <-- change here if you want to save the output file in a different directory or with a different name
    h_feeddown_frac_xi_all->Write();
    h_feeddown_frac_xi_rapidity_all->Write();
    h_feeddown_frac_xi_eta_all->Write();
    h_feeddown_frac_xi_phi_all->Write();
    fout->Close();

    const double n_pythia_events = h_event_counter->Integral();
    if (n_pythia_events <= 0)
    {
        std::cerr << "No counted PYTHIA events found." << std::endl;
        return;
    }

    h_lambda_pt_all->Scale(1.0 / n_pythia_events);
    h_lambda_pt_from_xi_all->Scale(1.0 / n_pythia_events);
    h_lambda_pt_from_xi_charged->Scale(1.0 / n_pythia_events);
    h_lambda_pt_from_xi_neutral->Scale(1.0 / n_pythia_events);
    h_lambda_rapidity_all->Scale(1.0 / n_pythia_events);
    h_lambda_rapidity_from_xi_all->Scale(1.0 / n_pythia_events);
    h_lambda_rapidity_from_xi_charged->Scale(1.0 / n_pythia_events);
    h_lambda_rapidity_from_xi_neutral->Scale(1.0 / n_pythia_events);
    h_lambda_eta_all->Scale(1.0 / n_pythia_events);
    h_lambda_eta_from_xi_all->Scale(1.0 / n_pythia_events);
    h_lambda_eta_from_xi_charged->Scale(1.0 / n_pythia_events);
    h_lambda_eta_from_xi_neutral->Scale(1.0 / n_pythia_events);
    h_lambda_phi_all->Scale(1.0 / n_pythia_events);
    h_lambda_phi_from_xi_all->Scale(1.0 / n_pythia_events);
    h_lambda_phi_from_xi_charged->Scale(1.0 / n_pythia_events);
    h_lambda_phi_from_xi_neutral->Scale(1.0 / n_pythia_events);

    draw1Dhistogram(h_lambda_pt_all, false, false, true, "#Lambda^{0}(+c.c) p_{T} [GeV]", "Normalized entries", {"Inclusive #Xi^{-}+#Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_pt_all");
    draw1Dhistogram(h_lambda_pt_from_xi_all, false, false, true, "#Lambda^{0}(+c.c) p_{T} [GeV]", "Normalized entries", {"From #Xi^{-}+#Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_pt_from_xi_all");
    draw1Dhistogram(h_lambda_pt_from_xi_charged, false, false, true, "#Lambda^{0}(+c.c) p_{T} [GeV]", "Normalized entries", {"From #Xi^{-} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_pt_from_xi_charged");
    draw1Dhistogram(h_lambda_pt_from_xi_neutral, false, false, true, "#Lambda^{0}(+c.c) p_{T} [GeV]", "Normalized entries", {"From #Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_pt_from_xi_neutral");
    draw1Dhistogram(h_lambda_rapidity_all, false, false, false, "#Lambda^{0}(+c.c) rapidity", "Normalized entries", {"Inclusive #Xi^{-}+#Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_rapidity_all");
    draw1Dhistogram(h_lambda_rapidity_from_xi_all, false, false, false, "#Lambda^{0}(+c.c) rapidity", "Normalized entries", {"From #Xi^{-}+#Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_rapidity_from_xi_all");
    draw1Dhistogram(h_lambda_rapidity_from_xi_charged, false, false, false, "#Lambda^{0}(+c.c) rapidity", "Normalized entries", {"From #Xi^{-} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_rapidity_from_xi_charged");
    draw1Dhistogram(h_lambda_rapidity_from_xi_neutral, false, false, false, "#Lambda^{0}(+c.c) rapidity", "Normalized entries", {"From #Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_rapidity_from_xi_neutral");
    draw1Dhistogram(h_lambda_eta_all, false, false, false, "#Lambda^{0}(+c.c) #eta", "Normalized entries", {"Inclusive #Xi^{-}+#Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_eta_all");
    draw1Dhistogram(h_lambda_eta_from_xi_all, false, false, false, "#Lambda^{0}(+c.c) #eta", "Normalized entries", {"From #Xi^{-}+#Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_eta_from_xi_all");
    draw1Dhistogram(h_lambda_eta_from_xi_charged, false, false, false, "#Lambda^{0}(+c.c) #eta", "Normalized entries", {"From #Xi^{-} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_eta_from_xi_charged");
    draw1Dhistogram(h_lambda_eta_from_xi_neutral, false, false, false, "#Lambda^{0}(+c.c) #eta", "Normalized entries", {"From #Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_eta_from_xi_neutral");
    draw1Dhistogram(h_lambda_phi_all, false, false, false, "#Lambda^{0}(+c.c) #phi [rad]", "Normalized entries", {"Inclusive #Xi^{-}+#Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_phi_all");
    draw1Dhistogram(h_lambda_phi_from_xi_all, false, false, false, "#Lambda^{0}(+c.c) #phi [rad]", "Normalized entries", {"From #Xi^{-}+#Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_phi_from_xi_all");
    draw1Dhistogram(h_lambda_phi_from_xi_charged, false, false, false, "#Lambda^{0}(+c.c) #phi [rad]", "Normalized entries", {"From #Xi^{-} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_phi_from_xi_charged");
    draw1Dhistogram(h_lambda_phi_from_xi_neutral, false, false, false, "#Lambda^{0}(+c.c) #phi [rad]", "Normalized entries", {"From #Xi^{0} (+c.c)", cutinfo, "(Normalized by N_{PYTHIA events})"}, "hist e1", plotdir + "lambda_phi_from_xi_neutral");
}