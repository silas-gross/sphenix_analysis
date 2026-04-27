#include "analysisHelper.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

// ─────────────────────────────────────────────────────────────────────────────
//  wEEC_doUnfolding.C
//
//  Single-step 3D unfolding of the wEEC, one RooUnfoldBayes call per ΔΦ bin.
//
//  The flat-3D index defined in analysisHelper.h only covers kinematically
//  valid (iL,iS) dijet pT pairs (trueSublPtBins[iS] < trueLeadPtBins[iL+1]),
//  so every truth bin in the response matrix has a nonzero row sum and
//  RooUnfold's Bayesian prior is fully nonzero. No runtime compression needed.
//
//  Dijet-level fakes and misses are encoded directly in the response matrix
//  via RooUnfoldResponse::Fake() and ::Miss() in fillResponse.C. Tower-level
//  fakes/misses within matched dijets are also in the response matrix.
//  RooUnfold handles all of this internally.
//
//  Outputs per ΔΦ bin k:
//    hWEEC3D_meas_{k}      — measured input clone (diagnostic)
//    hWEEC3D_unfolded_{k}  — unfolded flat-3D histogram
//
//  Other outputs:
//    hJetPt_meas, hJetPt_unfolded, hJetPt_truth  — dijet pT unfolding
//
//  Mode::kFull / kHalf: measured from MC response file (hWEEC3D_meas_{k})
//  Mode::kData:         measured from measFile (hWEEC3D_data_{k})
// ─────────────────────────────────────────────────────────────────────────────
void wEEC_doUnfolding(const char* respFile,
                       const char* outFileName,
                       int         nIterWEEC   = 4,
                       Mode        mode        = Mode::kFull,
                       const char* measFile    = nullptr,
                       int         minCounts   = 0)   // bins with fewer raw entries are zeroed
{
    const bool isData    = (mode == Mode::kData);
    const bool isClosure = !isData;
    if (isData && !measFile) {
        std::cerr << "kData mode requires measFile\n"; return;
    }

    TFile* fResp = TFile::Open(respFile, "READ");
    if (!fResp || fResp->IsZombie()) {
        std::cerr << "Cannot open " << respFile << "\n"; return;
    }

    TFile* fMeas = nullptr;
    if (isData) {
        fMeas = TFile::Open(measFile, "READ");
        if (!fMeas || fMeas->IsZombie()) {
            std::cerr << "Cannot open " << measFile << "\n"; return;
        }
    }
    TFile* fMeasSrc = fMeas ? fMeas : fResp;

    TFile* fOut = new TFile(outFileName, "RECREATE");

    auto meas3DName = [&](int k) -> std::string {
        return isData
            ? std::format("hWEEC3D_data_{}", k)
            : std::format("hWEEC3D_meas_{}", k);
    };

    // ── dijet pT unfolding ────────────────────────────────────────────────
    {
        RooUnfoldResponse* respJetPt =
            (RooUnfoldResponse*)fResp->Get("response_jet_pT");
        if (!respJetPt) {
            std::cerr << "Warning: cannot find response_jet_pT — skipping jet pT unfolding\n";
        } else {
            TH1D* hJetPt_input = nullptr;
            if (mode == Mode::kFull) {
                hJetPt_input = (TH1D*)fResp->Get("hJetPt_reco");
                if (hJetPt_input) {
                    hJetPt_input = (TH1D*)hJetPt_input->Clone("hJetPt_meas");
                    hJetPt_input->SetDirectory(0);
                }
            } else {
                hJetPt_input = (TH1D*)fResp->Get("hJetPt_meas");
                if (hJetPt_input) {
                    hJetPt_input = (TH1D*)hJetPt_input->Clone();
                    hJetPt_input->SetDirectory(0);
                }
            }

            if (hJetPt_input) {
                RooUnfoldBayes unfoldJet(respJetPt, hJetPt_input, 4);
                unfoldJet.HandleFakes(true);
                TH1D* hJetPt_unfolded = (TH1D*)unfoldJet.Hunfold();
                hJetPt_unfolded->SetDirectory(0);
                hJetPt_unfolded->SetName("hJetPt_unfolded");

                fOut->cd();
                hJetPt_input->Write("hJetPt_meas");
                hJetPt_unfolded->Write();

                if (isClosure) {
                    TH1D* hJetPt_truth = (TH1D*)fResp->Get("hJetPt_truth");
                    if (hJetPt_truth) {
                        hJetPt_truth->SetDirectory(0);
                        hJetPt_truth->Write();
                    }
                }

                delete hJetPt_input;
                delete hJetPt_unfolded;
            }
        }
    }

    // ── wEEC unfolding — one ΔΦ bin at a time ────────────────────────────
    for (int k = 0; k < nDphi; ++k)
    {
        TH1D* hMeas = (TH1D*)fMeasSrc->Get(meas3DName(k).c_str());
        if (!hMeas || hMeas->Integral() == 0) continue;

        TH1D* hMeasClone = (TH1D*)hMeas->Clone(
            std::format("hWEEC3D_meas_{}", k).c_str());
        hMeasClone->SetDirectory(0);

        RooUnfoldResponse* resp =
            (RooUnfoldResponse*)fResp->Get(
                std::format("response_wEEC3D_{}", k).c_str());
        if (!resp) { delete hMeasClone; continue; }

        // ── Count-based trimming ──────────────────────────────────────────
        // hWEEC3D_counts_{k} is a 2D histogram with the same (reco, truth)
        // binning as the response matrix, filled with unity weight for every
        // matched tower pair that entered a Fill call.  Any cell whose raw
        // count is below minCounts is individually zeroed in the cloned
        // response matrix; cells above threshold are left intact.
        // The measured histogram is not modified — only the response is
        // trimmed, so RooUnfold simply has no path to those truth bins from
        // those reco bins.
        TH2D* hCounts = (TH2D*)fResp->Get(
            std::format("hWEEC3D_counts_{}", k).c_str());
        
        RooUnfoldResponse* respTrimmed = nullptr;
        // hResp2D bin indices: X = reco flat3D (1-based), Y = truth flat3D
        TH2D* hResp2D = (TH2D*)resp->Hresponse();
        if (hCounts) {
            for (int rBin = 1; rBin <= hResp2D->GetNbinsX(); ++rBin)
            for (int tBin = 1; tBin <= hResp2D->GetNbinsY(); ++tBin) {
                double rawCount = hCounts->GetBinContent(rBin, tBin);
                if (rawCount < minCounts)
                    hResp2D->SetBinContent(rBin, tBin, 0.0);
            }
            respTrimmed = new RooUnfoldResponse((TH1D*)resp->Hmeasured()->Clone(), (TH1D*)resp->Htruth()->Clone(), hResp2D);
        }


         RooUnfoldResponse *respToUse = (respTrimmed != nullptr) ? respTrimmed : resp;
    

        std::cout << "Unfolding ΔΦ bin " << k << std::endl;

        // All truth bins in the response matrix are kinematically valid by
        // construction (see analysisHelper.h), so the prior vector is fully
        // nonzero and no runtime compression is needed.
        RooUnfoldBayes unfold(respToUse, hMeasClone, nIterWEEC);
        unfold.HandleFakes(true);
        TH1D* hUnf = (TH1D*)unfold.Hunfold();
        hUnf->SetDirectory(0);
        hUnf->SetName(std::format("hWEEC3D_unfolded_{}", k).c_str());

        fOut->cd();
        hMeasClone->Write();
        hUnf->Write();
        respTrimmed->Write(std::format("response_wEEC3D_trimmed_{}", k).c_str());   // save for drawClosure.C diagnostics

        delete hMeasClone;
        delete hUnf;
        delete respTrimmed;
    }

    fResp->Close(); delete fResp;
    if (fMeas) { fMeas->Close(); delete fMeas; }
    fOut->Close();
    std::cout << "Written: " << outFileName << "\n";
}
