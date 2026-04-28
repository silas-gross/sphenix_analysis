#include "analysisHelper.h"
#include "RooUnfoldResponse.h"

// ─────────────────────────────────────────────────────────────────────────────
//  drawClosure.C
//
//  Draws diagnostic plots from the outputs of:
//    jet_pT_doUnfolding.C  → jetPTFile
//    wEEC_doUnfolding.C    → wEECFile
//    merged response file  → respFile  (truth reference)
//
//  Truth reference for wEEC per-bin plots:
//    kFull — hWEEC3D_truth_{k} from respFile (fullClosure response)
//    kHalf — hWEEC3D_truth_{k} from respFile (halfClosure response, training half)
//    kData — hWEEC3D_truth_{k} from respFile (fullClosure response, full sim)
//
//  Dijet pT plots:
//    jetPT-{label}.png              — 2D truth/unfolded/ratio  (kFull/kHalf)
//    jetPT_projections-{label}.png  — 1D lead/subl pT          (kFull/kHalf)
//    jetPT_recoVsUnfolded-{label}   — reco vs unfolded          (kData)
//
//  wEEC plots:
//    wEEC_inclusive-{label}.png     — inclusive unfolded vs truth
//    wEEC_perBin-{label}.png        — grid: unfolded vs truth per dijet pT bin
//    wEEC_ratio-{label}.png         — grid: unfolded/truth ratio per dijet pT bin
//    wEEC_measVsUnf-{label}.png     — grid: measured vs unfolded (vs truth)
//                                     per dijet pT bin
//
//  All wEEC histograms normalized by NormalizeWEEC (integral then bin width).
// ─────────────────────────────────────────────────────────────────────────────

static const int kColTruth    = kRed   + 1;
static const int kColUnfolded = kBlue  + 1;
static const int kColMeas     = kGreen + 2;
static const int kMarTruth    = 21;
static const int kMarUnfolded = 20;
static const int kMarMeas     = 33;

std::pair<TPad*,TPad*> MakeRatioPads(TVirtualPad* parent)
{
    parent->cd();
    TPad* pTop = new TPad("pTop","",0,0.30,1,1.0);
    TPad* pBot = new TPad("pBot","",0,0.00,1,0.30);
    pTop->SetBottomMargin(0.03); pBot->SetTopMargin(0.03);
    pBot->SetBottomMargin(0.35);
    pTop->Draw(); pBot->Draw();
    return {pTop, pBot};
}

void StyleRatio(TH1D* h, const char* xTitle,
                const char* yTitle = "Unf/Truth",
                double lo = 0.8, double hi = 1.2)
{
    h->SetLineColor(kBlack); h->SetMarkerColor(kBlack); h->SetMarkerStyle(20);
    h->GetYaxis()->SetRangeUser(lo,hi);
    h->GetYaxis()->SetTitle(yTitle);
    h->GetYaxis()->SetTitleSize(0.12); h->GetYaxis()->SetTitleOffset(0.38);
    h->GetYaxis()->SetLabelSize(0.10); h->GetYaxis()->SetNdivisions(504);
    h->GetXaxis()->SetTitle(xTitle);
    h->GetXaxis()->SetTitleSize(0.12); h->GetXaxis()->SetLabelSize(0.10);
    h->SetTitle("");
}

void DrawRefLine(TH1D* h, double y = 1.0)
{
    TLine* l = new TLine(h->GetXaxis()->GetXmin(),y,h->GetXaxis()->GetXmax(),y);
    l->SetLineStyle(2); l->SetLineColor(kGray+1); l->Draw();
}

void StyleLine(TH1D* h, int col, int mar)
{ h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mar); }

TH1D* MakeRatioHist(TH1D* hNum, TH1D* hDen, const char* name)
{
    TH1D* hR = (TH1D*)hNum->Clone(name);
    hR->SetDirectory(0);
    for (int k=1; k<=hR->GetNbinsX(); ++k) {
        double n=hNum->GetBinContent(k), ne=hNum->GetBinError(k);
        double d=hDen->GetBinContent(k), de=hDen->GetBinError(k);
        if (d>0) {
            double r=n/d;
            hR->SetBinContent(k,r);
            hR->SetBinError(k, r*std::sqrt((n>0?(ne/n)*(ne/n):0)+(de/d)*(de/d)));
        } else { hR->SetBinContent(k,0); hR->SetBinError(k,0); }
    }
    return hR;
}

// ─────────────────────────────────────────────────────────────────────────────
void drawClosure(const char* outDir   = ".",
                 const char* respFile = nullptr,
                 Mode        mode     = Mode::kFull)
{
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPaintTextFormat(".2f");

    const bool doTrim    = true;
    const bool usePWMean = true;   // true  → luminosity-weighted mean pair weight per bin
                                   // false → arithmetic bin-center approximation
    const bool isData    = (mode == Mode::kData);
    const bool isClosure = !isData;
    const std::string label   = ModeLabel(mode);
    const std::string plotDir = std::format("{}/Plots", outDir);

    std::string wEECFile  = std::format("{}/wEEC-{}.root",   outDir, label);

    // ════════════════════════════════════════════════════════════════════════
    //  Dijet pT plots
    //  Jet pT histograms are written by wEEC_doUnfolding.C into the same
    //  wEEC output file, so only one input file is needed.
    // ════════════════════════════════════════════════════════════════════════
    TFile* fWEEC = TFile::Open(wEECFile.c_str(),"READ");
    if (!fWEEC || fWEEC->IsZombie()) { std::cerr<<"Cannot open "<<wEECFile<<"\n"; return; }

    // Open response file immediately — needed for truth reference throughout
    TFile* fResp = nullptr;
    if (respFile) {
        fResp = TFile::Open(respFile,"READ");
        if (!fResp || fResp->IsZombie()) {
            std::cerr<<"Cannot open resp file "<<respFile<<"\n"; fResp=nullptr;
        }
    }

    // ── Load pair weight mean histograms (used when usePWMean == true) ────
    // hWEEC3D_pwNum_{k} = sum(evtWeight * pairWeight) per truth flat3D bin
    // hWEEC3D_pwDen_{k} = sum(evtWeight)              per truth flat3D bin
    // Both are nullptr if usePWMean is false or the histograms are absent.
    std::vector<TH1D*> hPWNum(nDphi, nullptr), hPWDen(nDphi, nullptr);
    if (usePWMean && fResp) {
        for (int k = 0; k < nDphi; ++k) {
            hPWNum[k] = (TH1D*)fResp->Get(std::format("hWEEC3D_pwNum_{}", k).c_str());
            hPWDen[k] = (TH1D*)fResp->Get(std::format("hWEEC3D_pwDen_{}", k).c_str());
            if (hPWNum[k]) hPWNum[k]->SetDirectory(0);
            if (hPWDen[k]) hPWDen[k]->SetDirectory(0);
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    //  Project hWEEC3D_unfolded_{k} → per-dijet-pT-bin and inclusive wEEC
    //  This replaces loading pre-projected hWEEC_bin_{ft} / hWEEC_inclusive
    //  which are no longer written by wEEC_doUnfolding.C.
    // ════════════════════════════════════════════════════════════════════════
    std::vector<TH1D*> hUnfBin(nTrueFlat(), nullptr);
    for (int ft=0; ft<nTrueFlat(); ++ft) {
        int iL,iS; TrueFlatToIJ(ft,iL,iS);
        hUnfBin[ft] = new TH1D(
            std::format("hUnfBin_{}",ft).c_str(),
            std::format("Unfolded iL{} iS{};#Delta#phi;wEEC",iL,iS).c_str(),
            nDphi, dPhiBins.data());
        hUnfBin[ft]->Sumw2();
        hUnfBin[ft]->SetDirectory(0);
    }
    TH1D* hInclUnfolded = new TH1D("hWEEC_inclusive",";#Delta#phi;wEEC",nDphi,dPhiBins.data());
    hInclUnfolded->Sumw2();
    hInclUnfolded->SetDirectory(0);

    for (int k=0; k<nDphi; ++k) {
        TH1D* hU3D = (TH1D*)fWEEC->Get(std::format("hWEEC3D_unfolded_{}",k).c_str());
        if (!hU3D) continue;
        // Inclusive projection — IsRecoAccessible gate is inside ProjectWEEC3DSlice
        double ival, ierr;
        if (usePWMean && hPWNum[k] && hPWDen[k])
            ProjectWEEC3DSlice(hU3D, ival, ierr, hPWNum[k], hPWDen[k]);
        else
            ProjectWEEC3DSlice(hU3D, ival, ierr);
        hInclUnfolded->SetBinContent(k+1, ival);
        hInclUnfolded->SetBinError  (k+1, ierr);
        // Per-bin projection — skip bins with no reco support
        for (int ft=0; ft<nTrueFlat(); ++ft) {
            int iL,iS; TrueFlatToIJ(ft,iL,iS);
            if (!IsRecoAccessible(iL,iS)) continue;
            double val,err;
            if (usePWMean && hPWNum[k] && hPWDen[k])
                ProjectWEEC3DSliceForBin(hU3D, iL, iS, val, err, hPWNum[k], hPWDen[k]);
            else
                ProjectWEEC3DSliceForBin(hU3D, iL, iS, val, err);
            hUnfBin[ft]->SetBinContent(k+1, val);
            hUnfBin[ft]->SetBinError  (k+1, err);
        }
    }
    NormalizeWEEC(hInclUnfolded);
    for (int ft=0; ft<nTrueFlat(); ++ft) NormalizeWEEC(hUnfBin[ft]);

    TH1D* hJetPt_unfolded = (TH1D*)fWEEC->Get("hJetPt_unfolded");
    TH1D* hJetPt_meas     = (TH1D*)fWEEC->Get("hJetPt_meas");
    TH1D* hJetPt_truth    = isClosure ? (TH1D*)fWEEC->Get("hJetPt_truth") : nullptr;
    if (!hJetPt_unfolded || !hJetPt_meas) {
        std::cerr<<"Missing jet pT histograms in "<<wEECFile<<"\n"; return;
    }
    hJetPt_unfolded->SetDirectory(0); hJetPt_meas->SetDirectory(0);
    if (hJetPt_truth) hJetPt_truth->SetDirectory(0);

    // Unfolded and truth live in truth-index space; meas lives in reco-index space
    TH2D* h2Unf  = new TH2D("h2Unf", ";Lead p_{T} (GeV);Subl p_{T} (GeV)",
                              nTrueLead,trueLeadPtBins.data(),nTrueSubl,trueSublPtBins.data());
    TH2D* h2Meas = new TH2D("h2Meas",";Lead p_{T} (GeV);Subl p_{T} (GeV)",
                              nRecoLead,recoLeadPtBins.data(),nRecoSubl,recoSublPtBins.data());
    TH2D* h2Truth = isClosure
        ? new TH2D("h2Truth",";Lead p_{T} (GeV);Subl p_{T} (GeV)",
                   nTrueLead,trueLeadPtBins.data(),nTrueSubl,trueSublPtBins.data()) : nullptr;
    TH2D* h2Ratio = isClosure
        ? new TH2D("h2Ratio",";Lead p_{T} (GeV);Subl p_{T} (GeV)",
                   nTrueLead,trueLeadPtBins.data(),nTrueSubl,trueSublPtBins.data()) : nullptr;
    h2Unf ->SetDirectory(0);
    h2Meas->SetDirectory(0);
    if (h2Truth) h2Truth->SetDirectory(0);
    if (h2Ratio) h2Ratio->SetDirectory(0);

    // Fill unfolded and truth (truth-index space)
    for (int f=0; f<nTrueFlat(); ++f) {
        int iL,iS; TrueFlatToIJ(f,iL,iS);
        double u=hJetPt_unfolded->GetBinContent(f+1), ue=hJetPt_unfolded->GetBinError(f+1);
        h2Unf->SetBinContent(iL+1,iS+1,u); h2Unf->SetBinError(iL+1,iS+1,ue);
        if (isClosure && hJetPt_truth) {
            double t=hJetPt_truth->GetBinContent(f+1), te=hJetPt_truth->GetBinError(f+1);
            h2Truth->SetBinContent(iL+1,iS+1,t); h2Truth->SetBinError(iL+1,iS+1,te);
            if (t>0) {
                double r=u/t, re=r*std::sqrt((u>0?(ue/u)*(ue/u):0)+(te/t)*(te/t));
                h2Ratio->SetBinContent(iL+1,iS+1,r); h2Ratio->SetBinError(iL+1,iS+1,re);
            }
        }
    }
    // Fill meas (reco-index space)
    for (int f=0; f<nRecoFlat(); ++f) {
        int iL,iS; RecoFlatToIJ(f,iL,iS);
        double m=hJetPt_meas->GetBinContent(f+1), me=hJetPt_meas->GetBinError(f+1);
        h2Meas->SetBinContent(iL+1,iS+1,m); h2Meas->SetBinError(iL+1,iS+1,me);
    }

    if (isClosure) {
        TCanvas* c1=new TCanvas("c1","Jet pT 2D",1600,550); c1->Divide(3,1);
        c1->cd(1); gPad->SetLogz(); h2Truth->SetTitle("Truth"); h2Truth->Draw("COLZ");
        c1->cd(2); gPad->SetLogz(); h2Unf->SetTitle("Unfolded"); h2Unf->Draw("COLZ");
        c1->cd(3); h2Ratio->SetTitle("Unfolded/Truth");
        h2Ratio->GetZaxis()->SetRangeUser(0.8,1.2); h2Ratio->Draw("COLZ TEXT");
        c1->SaveAs(std::format("{}/jetPT-{}.png",plotDir,label).c_str()); c1->Close(); delete c1;

        TH1D* hLT=(TH1D*)h2Truth->ProjectionX("hLT"), *hLU=(TH1D*)h2Unf->ProjectionX("hLU");
        TH1D* hST=(TH1D*)h2Truth->ProjectionY("hST"), *hSU=(TH1D*)h2Unf->ProjectionY("hSU");
        StyleLine(hLT,kColTruth,kMarTruth); StyleLine(hLU,kColUnfolded,kMarUnfolded);
        StyleLine(hST,kColTruth,kMarTruth); StyleLine(hSU,kColUnfolded,kMarUnfolded);
        TH1D* hLR=MakeRatioHist(hLU,hLT,"hLR"), *hSR=MakeRatioHist(hSU,hST,"hSR");
        TCanvas* c2=new TCanvas("c2","Jet pT Proj",1200,600); c2->Divide(2,1);
        for (int side=0; side<2; ++side) {
            auto [pTop,pBot]=MakeRatioPads(c2->cd(side+1));
            TH1D* hT=(side==0)?hLT:hST, *hU=(side==0)?hLU:hSU, *hR=(side==0)?hLR:hSR;
            const char* xL=(side==0)?"Lead p_{T} (GeV)":"Subl p_{T} (GeV)";
            pTop->cd(); pTop->SetLogy();
            hT->GetYaxis()->SetRangeUser(0.1,std::max(hT->GetMaximum(),hU->GetMaximum())*1.3);
            hT->GetXaxis()->SetLabelSize(0); hT->Draw("E"); hU->Draw("E SAME");
            TLegend* lg=new TLegend(0.55,0.72,0.88,0.88);
            lg->AddEntry(hT,"Truth","lp"); lg->AddEntry(hU,"Unfolded","lp"); lg->Draw();
            pBot->cd(); StyleRatio(hR,xL); hR->Draw("E"); DrawRefLine(hR);
        }
        c2->SaveAs(std::format("{}/jetPT_projections-{}.png",plotDir,label).c_str()); c2->Close(); delete c2;
    } else {
        TH1D* hLM=(TH1D*)h2Meas->ProjectionX("hLM"), *hLU=(TH1D*)h2Unf->ProjectionX("hLU");
        TH1D* hSM=(TH1D*)h2Meas->ProjectionY("hSM"), *hSU=(TH1D*)h2Unf->ProjectionY("hSU");
        StyleLine(hLM,kColMeas,kMarMeas); StyleLine(hLU,kColUnfolded,kMarUnfolded);
        StyleLine(hSM,kColMeas,kMarMeas); StyleLine(hSU,kColUnfolded,kMarUnfolded);
        TH1D* hLR=MakeRatioHist(hLU,hLM,"hLR"), *hSR=MakeRatioHist(hSU,hSM,"hSR");
        TCanvas* c2=new TCanvas("c2","Jet pT Reco vs Unf",1200,600); c2->Divide(2,1);
        for (int side=0; side<2; ++side) {
            auto [pTop,pBot]=MakeRatioPads(c2->cd(side+1));
            TH1D* hM=(side==0)?hLM:hSM, *hU=(side==0)?hLU:hSU, *hR=(side==0)?hLR:hSR;
            const char* xL=(side==0)?"Lead p_{T} (GeV)":"Subl p_{T} (GeV)";
            pTop->cd(); pTop->SetLogy();
            hM->GetYaxis()->SetRangeUser(0.1,std::max(hM->GetMaximum(),hU->GetMaximum())*1.3);
            hM->GetXaxis()->SetLabelSize(0); hM->Draw("E"); hU->Draw("E SAME");
            TLegend* lg=new TLegend(0.55,0.72,0.88,0.88);
            lg->AddEntry(hM,"Measured","lp"); lg->AddEntry(hU,"Unfolded","lp"); lg->Draw();
            pBot->cd(); StyleRatio(hR,xL,"Unf/Meas",0.5,1.5); hR->Draw("E"); DrawRefLine(hR);
        }
        c2->SaveAs(std::format("{}/jetPT_recoVsUnfolded-{}.png",plotDir,label).c_str());
        c2->Close(); delete c2;
    }

    // ════════════════════════════════════════════════════════════════════════
    //  wEEC plots — fWEEC is already open from the jet pT section above
    // ════════════════════════════════════════════════════════════════════════

    // ── Build per-dijet-pT-bin truth reference ────────────────────────────
    // For each ft, project hWEEC3D_truth_{k} for every k using
    // ProjectWEEC3DSliceForBin, assemble into a TH1D vs ΔΦ, then normalize.
    // For kData respFile points to the fullClosure response (full sim truth).
    std::vector<TH1D*> hTruthBin(nTrueFlat(), nullptr);
    if (fResp) {
        for (int ft=0; ft<nTrueFlat(); ++ft) {
            int iL,iS; TrueFlatToIJ(ft,iL,iS);
            if (!IsRecoAccessible(iL,iS)) { hTruthBin[ft] = nullptr; continue; }
            TH1D* hT = new TH1D(
                std::format("hTruthBin_{}",ft).c_str(),
                std::format("Truth iL{} iS{};#Delta#phi;wEEC",iL,iS).c_str(),
                nDphi, dPhiBins.data());
            hT->Sumw2();
            hT->SetDirectory(0);
            for (int k=0; k<nDphi; ++k) {
                TH1D* hT3D=(TH1D*)fResp->Get(
                    std::format("hWEEC3D_truth_{}",k).c_str());
                if (!hT3D) continue;
                double val,err;
                if (usePWMean && hPWNum[k] && hPWDen[k])
                    ProjectWEEC3DSliceForBin(hT3D, iL, iS, val, err, hPWNum[k], hPWDen[k]);
                else
                    ProjectWEEC3DSliceForBin(hT3D, iL, iS, val, err);
                hT->SetBinContent(k+1,val);
                hT->SetBinError  (k+1,err);
            }
            NormalizeWEEC(hT);
            hTruthBin[ft] = hT;
        }
    }

    // ── Build inclusive truth reference ───────────────────────────────────
    TH1D* hInclTruth = nullptr;
    if (fResp) {
        hInclTruth = new TH1D("hInclTruth",";#Delta#phi;wEEC",nDphi,dPhiBins.data());
        hInclTruth->Sumw2();
        hInclTruth->SetDirectory(0);
        for (int k=0; k<nDphi; ++k) {
                TH1D* hT3D=(TH1D*)fResp->Get(
                    std::format("hWEEC3D_truth_{}",k).c_str());
                if (!hT3D) continue;
                double val,err;
                // ProjectWEEC3DSlice already gates on IsRecoAccessible
                if (usePWMean && hPWNum[k] && hPWDen[k])
                    ProjectWEEC3DSlice(hT3D, val, err, hPWNum[k], hPWDen[k]);
                else
                    ProjectWEEC3DSlice(hT3D, val, err);
                hInclTruth->SetBinContent(k+1,val);
                hInclTruth->SetBinError  (k+1,err);
            }
        NormalizeWEEC(hInclTruth);
    }

    // ── Load measured per-bin histograms and normalize ────────────────────
    // For kFull/kHalf: hWEEC3D_meas_{k} from response file, project per bin
    // For kData: hWEEC3D_data_{k} from measFile (not available here, skip)
    std::vector<TH1D*> hMeasBin(nRecoFlat(), nullptr);
    {
        // Use fResp for kFull/kHalf (meas is in response file),
        // skip for kData (measured data file not passed to drawClosure)
        TFile* fMeasSrc = (!isData && fResp) ? fResp : nullptr;
        if (fMeasSrc) {
            for (int ft=0; ft<nRecoFlat(); ++ft) {
                int iL,iS; RecoFlatToIJ(ft,iL,iS);
                TH1D* hM = new TH1D(
                    std::format("hMeasBin_{}",ft).c_str(),
                    std::format("Meas iL{} iS{};#Delta#phi;wEEC",iL,iS).c_str(),
                    nDphi, dPhiBins.data());
                hM->Sumw2();
                hM->SetDirectory(0);
                for (int k=0; k<nDphi; ++k) {
                    TH1D* hM3D=(TH1D*)fMeasSrc->Get(
                        std::format("hWEEC3D_meas_{}",k).c_str());
                    if (!hM3D) continue;
                    double val,err;
                    ProjectWEEC3DSliceForBin(hM3D,iL,iS,val,err);
                    hM->SetBinContent(k+1,val);
                    hM->SetBinError  (k+1,err);
                }
                NormalizeWEEC(hM);
                hMeasBin[ft] = hM;
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    //  Inclusive wEEC
    // ════════════════════════════════════════════════════════════════════════
    if (hInclUnfolded) {
        TCanvas* cI=new TCanvas("cI","Inclusive wEEC",800,600);
        if (hInclTruth) {
            TPad* pTop=new TPad("pTop","",0,0.30,1,1.0);
            TPad* pBot=new TPad("pBot","",0,0.00,1,0.30);
            pTop->SetBottomMargin(0.03); pBot->SetTopMargin(0.03);
            pBot->SetBottomMargin(0.35); pTop->Draw(); pBot->Draw();
            pTop->cd(); pTop->SetLogy();
            StyleLine(hInclTruth,   kColTruth,   kMarTruth);
            StyleLine(hInclUnfolded,kColUnfolded,kMarUnfolded);
            double ymin=1e30,ymax=-1e30;
            for (TH1D* h:{hInclTruth,hInclUnfolded})
                for (int i=1;i<=h->GetNbinsX();++i) {
                    double v=h->GetBinContent(i);
                    if (v>0) ymin=std::min(ymin,v);
                    ymax=std::max(ymax,v);
                }
            hInclTruth->GetYaxis()->SetRangeUser(ymin*0.3,ymax*5.0);
            hInclTruth->GetYaxis()->SetTitle("1/N dN/d#Delta#phi");
            hInclTruth->GetXaxis()->SetLabelSize(0);
            hInclTruth->Draw("E"); hInclUnfolded->Draw("E SAME");
            TLegend* lg=new TLegend(0.55,0.72,0.90,0.88);
            lg->AddEntry(hInclTruth,"Truth (truth towers)","lp");
            lg->AddEntry(hInclUnfolded,"Unfolded","lp"); lg->Draw();
            pBot->cd();
            TH1D* hR=MakeRatioHist(hInclUnfolded,hInclTruth,"hInclRatio");
            StyleRatio(hR,"#Delta#phi"); hR->Draw("E"); DrawRefLine(hR);
        } else {
            cI->cd(); cI->SetLogy();
            StyleLine(hInclUnfolded,kColUnfolded,kMarUnfolded);
            hInclUnfolded->GetYaxis()->SetTitle("1/N dN/d#Delta#phi");
            hInclUnfolded->GetXaxis()->SetTitle("#Delta#phi");
            hInclUnfolded->Draw("E");
        }
        cI->SaveAs(std::format("{}/wEEC_inclusive-{}.png",plotDir,label).c_str());
        cI->Close(); delete cI;
    }

    // ════════════════════════════════════════════════════════════════════════
    //  Per-dijet-pT-bin grid helper
    // ════════════════════════════════════════════════════════════════════════
    const int padW=280, padH=220;

    // Compute shared y range across all populated bins (log scale)
    auto yRange1 = [](const std::vector<TH1D*>& hA,
                      const std::vector<TH1D*>& hB,
                      const std::vector<TH1D*>& hC)
        -> std::pair<double,double>
    {
        double yLo=std::numeric_limits<double>::max();
        double yHi=-std::numeric_limits<double>::max();
        for (int ft=0; ft<nTrueFlat(); ++ft)
            for (TH1D* h : {hA[ft], hB[ft], hC[ft]}) {
                if (!h) continue;
                for (int i=1;i<=h->GetNbinsX();++i) {
                    double v=h->GetBinContent(i);
                    if (v>0) yLo=std::min(yLo,v);
                    yHi=std::max(yHi,v);
                }
            }
        if (yLo>=yHi) return {1e-6,1.0};
        return {yLo*0.3, yHi*5.0};
    };

    // Grid drawing — use a plain nested scope instead of a lambda to avoid
    // Cling interpreter crashes on lambda destruction at function exit.
    // Called three times below via a simple struct-less inline pattern.
    auto DrawGrid = [=](const char* canName, const char* outFile,
                        const std::vector<TH1D*>& hA,
                        const std::vector<TH1D*>& hB,
                        const std::vector<TH1D*>& hC,
                        bool isRatio,
                        double yLo, double yHi,
                        const char* legA, const char* legB, const char* legC)
    {
        TCanvas* cG=new TCanvas(canName,canName,nTrueLead*padW,nTrueSubl*padH);
        cG->Divide(nTrueLead,nTrueSubl,0,0);
        for (int iS=0; iS<nTrueSubl; ++iS)
        for (int iL=0; iL<nTrueLead; ++iL)
        {
            cG->cd(iS*nTrueLead+iL+1);
            gPad->SetLeftMargin(0.18); gPad->SetBottomMargin(0.18);
            gPad->SetRightMargin(0.04); gPad->SetTopMargin(0.10);
            TLatex tex; tex.SetNDC(); tex.SetTextSize(0.09); tex.SetTextAlign(13);
            std::string bl=std::format("{:.0f}-{:.0f}/{:.0f}-{:.0f} GeV",
                trueLeadPtBins[iL],trueLeadPtBins[iL+1],trueSublPtBins[iS],trueSublPtBins[iS+1]);

            int ft = TrueFlatIndex(iL,iS);

            // Gray out kinematically forbidden bins (subl lower edge >= lead upper edge).
            // All other bins are drawn even if outside the primary reporting region.
            if (ft < 0 || !hA[ft]) {
                gPad->SetFillColor(kGray);
                tex.DrawLatex(0.22,0.88,bl.c_str()); continue;
            }

            if (isRatio) {
                TH1D* hH=hA[ft];
                hH->SetLineColor(kBlack); hH->SetMarkerColor(kBlack);
                hH->SetMarkerStyle(20);
                hH->GetYaxis()->SetRangeUser(yLo,yHi);
                hH->GetYaxis()->SetTitle("Unf/Truth");
                hH->GetYaxis()->SetTitleSize(0.09); hH->GetYaxis()->SetLabelSize(0.08);
                hH->GetYaxis()->SetNdivisions(504);
                hH->GetXaxis()->SetTitle("#Delta#phi");
                hH->GetXaxis()->SetTitleSize(0.09); hH->GetXaxis()->SetLabelSize(0.08);
                hH->SetTitle("");
                hH->Draw("E"); DrawRefLine(hH);
            } else {
                gPad->SetLogy();
                TH1D* hBase = hB[ft] ? hB[ft] : (hC[ft] ? hC[ft] : hA[ft]);
                hBase->GetYaxis()->SetRangeUser(yLo,yHi);
                hBase->GetYaxis()->SetTitle("1/N dN/d#Delta#phi");
                hBase->GetYaxis()->SetTitleSize(0.09); hBase->GetYaxis()->SetTitleOffset(0.9);
                hBase->GetYaxis()->SetLabelSize(0.08);
                hBase->GetXaxis()->SetTitle("#Delta#phi");
                hBase->GetXaxis()->SetTitleSize(0.09); hBase->GetXaxis()->SetLabelSize(0.08);
                hBase->SetTitle("");
                StyleLine(hA[ft],kColUnfolded,kMarUnfolded);
                if (hB[ft]) StyleLine(hB[ft],kColTruth,kMarTruth);
                if (hC[ft]) StyleLine(hC[ft],kColMeas,  kMarMeas);
                hBase->Draw("E");
                hA[ft]->Draw("E SAME");
                if (hB[ft] && hB[ft]!=hBase) hB[ft]->Draw("E SAME");
                if (hC[ft] && hC[ft]!=hBase) hC[ft]->Draw("E SAME");
                TLegend* lg=new TLegend(0.40,0.72,0.96,0.88);
                lg->SetTextSize(0.07);
                lg->AddEntry(hA[ft],legA,"lp");
                if (hB[ft]) lg->AddEntry(hB[ft],legB,"lp");
                if (hC[ft]) lg->AddEntry(hC[ft],legC,"lp");
                lg->Draw();
            }
            tex.DrawLatex(0.22,0.88,bl.c_str());
        }
        cG->cd();
        TLatex ct; ct.SetNDC(); ct.SetTextSize(0.022); ct.SetTextAlign(22);
        ct.DrawLatex(0.5,0.003,"Lead jet p_{T} bin #rightarrow");
        ct.SetTextAngle(90); ct.DrawLatex(0.006,0.5,"Subl jet p_{T} bin #rightarrow");
        cG->SaveAs(outFile); cG->Close(); delete cG;
    };

    // ════════════════════════════════════════════════════════════════════════
    //  Plot 1: unfolded vs truth  (wEEC_perBin)
    // ════════════════════════════════════════════════════════════════════════
    {
        std::vector<TH1D*> hNull(nTrueFlat(), nullptr);
        auto [yLo,yHi] = yRange1(hUnfBin, hTruthBin, hNull);
        DrawGrid("cPerBin",
                 std::format("{}/wEEC_perBin-{}.png",plotDir,label).c_str(),
                 hUnfBin, hTruthBin, hNull,
                 false, yLo, yHi,
                 "Unfolded", "Truth (truth towers)", "");
    }

    // ════════════════════════════════════════════════════════════════════════
    //  Plot 2: ratio unfolded/truth  (wEEC_ratio)  — only if truth available
    // ════════════════════════════════════════════════════════════════════════
    if (fResp) {
        std::vector<TH1D*> hRatioBin(nTrueFlat(), nullptr);
        for (int ft=0; ft<nTrueFlat(); ++ft) {
            if (hUnfBin[ft] && hTruthBin[ft])
                hRatioBin[ft] = MakeRatioHist(hUnfBin[ft], hTruthBin[ft],
                                    std::format("hRatioBin_{}",ft).c_str());
        }
        std::vector<TH1D*> hNull(nTrueFlat(), nullptr);
        DrawGrid("cRatio",
                 std::format("{}/wEEC_ratio-{}.png",plotDir,label).c_str(),
                 hRatioBin, hNull, hNull,
                 true, 0.5, 1.5,
                 "Unf/Truth", "", "");
        for (int ft=0; ft<nTrueFlat(); ++ft) delete hRatioBin[ft];
    }

    // ════════════════════════════════════════════════════════════════════════
    //  Plot 3: measured vs unfolded vs truth  (wEEC_measVsUnf)
    // ════════════════════════════════════════════════════════════════════════
    {
        // Expand reco-indexed meas into truth-indexed space for DrawGrid
        std::vector<TH1D*> hMeasBinTrue(nTrueFlat(), nullptr);
        for (int ftR=0; ftR<nRecoFlat(); ++ftR) {
            if (!hMeasBin[ftR]) continue;
            int iLr,iSr; RecoFlatToIJ(ftR,iLr,iSr);
            int iLt=FindBin(recoLeadPtBins[iLr],trueLeadPtBins);
            int iSt=FindBin(recoSublPtBins[iSr], trueSublPtBins);
            if (iLt>=0 && iSt>=0) hMeasBinTrue[TrueFlatIndex(iLt,iSt)]=hMeasBin[ftR];
        }
        std::vector<TH1D*> hNull(nTrueFlat(), nullptr);
        auto [yLo,yHi] = yRange1(hUnfBin, hTruthBin, !isData ? hMeasBinTrue : hNull);
        DrawGrid("cMeasVsUnf",
                 std::format("{}/wEEC_measVsUnf-{}.png",plotDir,label).c_str(),
                 hUnfBin,
                 hTruthBin,
                 !isData ? hMeasBinTrue : hNull,
                 false, yLo, yHi,
                 "Unfolded", "Truth (truth towers)",
                 !isData ? "Measured" : "");
    }

    // ════════════════════════════════════════════════════════════════════════
    //  1D closure projections: lead pT, subl pT, pair weight
    //  wEEC_projections-{label}.png
    // ════════════════════════════════════════════════════════════════════════
    if (isClosure && fResp) {
        TH1D* hUnfLead = new TH1D("hUnfLead","",nTrueLead,trueLeadPtBins.data());
        TH1D* hUnfSubl = new TH1D("hUnfSubl","",nTrueSubl,trueSublPtBins.data());
        TH1D* hUnfPW   = new TH1D("hUnfPW",  "",nPairWeight,pairWeightBins.data());
        TH1D* hTruLead = new TH1D("hTruLead","",nTrueLead,trueLeadPtBins.data());
        TH1D* hTruSubl = new TH1D("hTruSubl","",nTrueSubl,trueSublPtBins.data());
        TH1D* hTruPW   = new TH1D("hTruPW",  "",nPairWeight,pairWeightBins.data());
        for (TH1D* h:{hUnfLead,hUnfSubl,hUnfPW,hTruLead,hTruSubl,hTruPW})
            { h->SetDirectory(0); h->Sumw2(); }

        // Project from the already-normalized per-bin wEEC histograms so that
        // truth and unfolded use exactly the same normalization convention.
        // hUnfBin[ft] and hTruthBin[ft] are normalized wEEC vs ΔΦ.
        // Summing over ΔΦ bins gives the ΔΦ-integrated shape weight, and the
        // errors from RooUnfold's covariance are properly carried through.
        // Lead and subl pT projections come from summing over ΔΦ bins.
        for (int ft=0; ft<nTrueFlat(); ++ft) {
            int iL, iS; TrueFlatToIJ(ft, iL, iS);
            if (!IsRecoAccessible(iL, iS)) continue;
            double lCen = 0.5*(trueLeadPtBins[iL]+trueLeadPtBins[iL+1]);
            double sCen = 0.5*(trueSublPtBins[iS]+trueSublPtBins[iS+1]);
            int bL = hUnfLead->FindBin(lCen);
            int bS = hUnfSubl->FindBin(sCen);
            if (hUnfBin[ft]) {
                double sumV=0, sumE2=0;
                for (int k=1; k<=hUnfBin[ft]->GetNbinsX(); ++k) {
                    sumV  += hUnfBin[ft]->GetBinContent(k);
                    sumE2 += std::pow(hUnfBin[ft]->GetBinError(k), 2);
                }
                hUnfLead->AddBinContent(bL, sumV);
                hUnfSubl->AddBinContent(bS, sumV);
                hUnfLead->SetBinError(bL, std::hypot(hUnfLead->GetBinError(bL), std::sqrt(sumE2)));
                hUnfSubl->SetBinError(bS, std::hypot(hUnfSubl->GetBinError(bS), std::sqrt(sumE2)));
            }
            if (hTruthBin[ft]) {
                double sumV=0, sumE2=0;
                for (int k=1; k<=hTruthBin[ft]->GetNbinsX(); ++k) {
                    sumV  += hTruthBin[ft]->GetBinContent(k);
                    sumE2 += std::pow(hTruthBin[ft]->GetBinError(k), 2);
                }
                hTruLead->AddBinContent(bL, sumV);
                hTruSubl->AddBinContent(bS, sumV);
                hTruLead->SetBinError(bL, std::hypot(hTruLead->GetBinError(bL), std::sqrt(sumE2)));
                hTruSubl->SetBinError(bS, std::hypot(hTruSubl->GetBinError(bS), std::sqrt(sumE2)));
            }
        }

        // Pair weight projection: read from raw 3D histograms since pair weight
        // is the innermost flat3D dimension and hUnfBin/hTruthBin integrate over it.
        // Read errors before AddBinContent to avoid Sumw2 corruption.
        for (int k=0; k<nDphi; ++k) {
            TH1D* hU3D = (TH1D*)fWEEC->Get(std::format("hWEEC3D_unfolded_{}",k).c_str());
            TH1D* hT3D = (TH1D*)fResp->Get(std::format("hWEEC3D_truth_{}",k).c_str());
            for (int ft=0; ft<nTrueFlat(); ++ft)
            for (int iPw=0; iPw<nPairWeight; ++iPw) {
                int iL, iS; TrueFlatToIJ(ft, iL, iS);
                if (!IsRecoAccessible(iL, iS)) continue;
                int iF = ft * nPairWeight + iPw;
                double den   = (usePWMean && hPWDen[k]) ? hPWDen[k]->GetBinContent(iF+1) : 0.0;
                double wMean = (den > 0)
                    ? hPWNum[k]->GetBinContent(iF+1) / den
                    : 0.5*(pairWeightBins[iPw]+pairWeightBins[iPw+1]);
                int bP = hUnfPW->FindBin(wMean);
                if (hU3D) {
                    double v=hU3D->GetBinContent(iF+1), e=hU3D->GetBinError(iF+1);
                    double ePprev = hUnfPW->GetBinError(bP);
                    hUnfPW->AddBinContent(bP, v);
                    hUnfPW->SetBinError(bP, std::hypot(ePprev, e));
                }
                if (hT3D) {
                    double v=hT3D->GetBinContent(iF+1), e=hT3D->GetBinError(iF+1);
                    double ePprev = hTruPW->GetBinError(bP);
                    hTruPW->AddBinContent(bP, v);
                    hTruPW->SetBinError(bP, std::hypot(ePprev, e));
                }
            }
        }

        struct ProjSpec { TH1D* hU; TH1D* hT; const char* xL; bool logx; };
        std::vector<ProjSpec> projs = {
            {hUnfLead,hTruLead,"Lead p_{T} (GeV)",false},
            {hUnfSubl,hTruSubl,"Subl p_{T} (GeV)",false},
            {hUnfPW,  hTruPW,  "Pair weight",      true }
        };
        TCanvas* cProj=new TCanvas("cProj","wEEC projections",1800,600);
        cProj->Divide(3,1);
        for (int p=0; p<3; ++p) {
            auto [hU,hT,xL,logx]=projs[p];
            StyleLine(hT,kColTruth,kMarTruth);
            StyleLine(hU,kColUnfolded,kMarUnfolded);
            auto [pTop,pBot]=MakeRatioPads(cProj->cd(p+1));
            pTop->cd(); pTop->SetLogy(); if (logx) pTop->SetLogx();
            double ymax=std::max(hT->GetMaximum(),hU->GetMaximum());
            double ymin=1e30;
            for (int b=1;b<=hT->GetNbinsX();++b)
                if (hT->GetBinContent(b)>0) ymin=std::min(ymin,hT->GetBinContent(b));
            if (ymin>ymax) ymin=ymax*1e-4;
            hT->GetYaxis()->SetRangeUser(ymin*0.3,ymax*3.0);
            hT->GetXaxis()->SetLabelSize(0);
            hT->GetYaxis()->SetTitle("Pairs (arb.)");
            hT->Draw("E"); hU->Draw("E SAME");
            TLegend* lg=new TLegend(0.55,0.72,0.88,0.88);
            lg->AddEntry(hT,"Truth","lp"); lg->AddEntry(hU,"Unfolded","lp"); lg->Draw();
            pBot->cd(); if (logx) pBot->SetLogx();
            TH1D* hR=MakeRatioHist(hU,hT,std::format("hProjR_{}",p).c_str());
            StyleRatio(hR,xL); hR->Draw("E"); DrawRefLine(hR);
        }
        cProj->SaveAs(std::format("{}/wEEC_projections-{}.png",plotDir,label).c_str());
        cProj->Close(); delete cProj;
        delete hUnfLead; delete hUnfSubl; delete hUnfPW;
        delete hTruLead; delete hTruSubl; delete hTruPW;
    }

    // ════════════════════════════════════════════════════════════════════════
    //  Pair weight closure: per ΔΦ bin grid (pairWeight_closure_dPhi)
    //  and per (ΔΦ, dijet pT) bin PDF (pairWeight_closure_perBin)
    // ════════════════════════════════════════════════════════════════════════
    if (isClosure && fResp) {
        // Build pair weight projections: hPWUnf[k] and hPWTru[k] (summed over dijet pT)
        // and hPWUnfBin[k*nTrueFlat()+ft] / hPWTruBin[k*nTrueFlat()+ft] (per dijet pT bin)
        std::vector<TH1D*> hPWUnf(nDphi,nullptr), hPWTru(nDphi,nullptr);
        std::vector<TH1D*> hPWUnfBin(nDphi*nTrueFlat(),nullptr);
        std::vector<TH1D*> hPWTruBin(nDphi*nTrueFlat(),nullptr);
        std::cout << "nPairWeight: " << nPairWeight << std::endl;
        for(int k=0; k<=nPairWeight; k++)
        {
            std::cout << "   pw " << k << ": " << pairWeightBins[k] << std::endl;
        }
        for (int k=0; k<nDphi; ++k) {
            hPWUnf[k]=new TH1D(std::format("hPWUnf_{}",k).c_str(),"",nPairWeight,pairWeightBins.data());
            hPWTru[k]=new TH1D(std::format("hPWTru_{}",k).c_str(),"",nPairWeight,pairWeightBins.data());
            hPWUnf[k]->SetDirectory(0); hPWTru[k]->SetDirectory(0);
            hPWUnf[k]->Sumw2();        hPWTru[k]->Sumw2();
            for (int ft=0; ft<nTrueFlat(); ++ft) {
                hPWUnfBin[k*nTrueFlat()+ft]=new TH1D(
                    std::format("hPWUnfBin_{}_{}",k,ft).c_str(),"",nPairWeight,pairWeightBins.data());
                hPWTruBin[k*nTrueFlat()+ft]=new TH1D(
                    std::format("hPWTruBin_{}_{}",k,ft).c_str(),"",nPairWeight,pairWeightBins.data());
                hPWUnfBin[k*nTrueFlat()+ft]->SetDirectory(0);
                hPWTruBin[k*nTrueFlat()+ft]->SetDirectory(0);
                hPWUnfBin[k*nTrueFlat()+ft]->Sumw2();
                hPWTruBin[k*nTrueFlat()+ft]->Sumw2();
            }
        }

        // Read unfolded 3D histograms directly from already-open fWEEC
        for (int k=0; k<nDphi; ++k) {
                TH1D* hU=(TH1D*)fWEEC->Get(std::format("hWEEC3D_unfolded_{}",k).c_str());
                TH1D* hT=(TH1D*)fResp->Get(std::format("hWEEC3D_truth_{}",k).c_str());
                for (int ft=0; ft<nTrueFlat(); ++ft)
                for (int iPw=0; iPw<nPairWeight; ++iPw) {
                    int iL, iS; TrueFlatToIJ(ft, iL, iS);
                    if (!IsRecoAccessible(iL, iS)) continue;
                    int iF = ft * nPairWeight + iPw;
                    double den   = (usePWMean && hPWDen[k]) ? hPWDen[k]->GetBinContent(iF+1) : 0.0;
                    double wMean = (den > 0)
                        ? hPWNum[k]->GetBinContent(iF+1) / den
                        : 0.5*(pairWeightBins[iPw]+pairWeightBins[iPw+1]);
                    if (hU) {
                            double v=hU->GetBinContent(iF+1), e=hU->GetBinError(iF+1);
                            int bP=hPWUnf[k]->FindBin(wMean);
                            int bPb=hPWUnfBin[k*nTrueFlat()+ft]->FindBin(wMean);
                            double bPe = hPWUnf[k]->GetBinError(bP);
                            double bPbe = hPWUnfBin[k*nTrueFlat()+ft]->GetBinError(bPbe);
                            hPWUnf[k]->Fill(wMean,v);
                            hPWUnfBin[k*nTrueFlat()+ft]->Fill(wMean,v);
                            hPWUnf[k]->SetBinError(bP,std::hypot(bPe,e));
                            hPWUnfBin[k*nTrueFlat()+ft]->SetBinError(bPb,std::hypot(bPbe,e));
                    }
                    if (hT) {
                            double v=hT->GetBinContent(iF+1);
                            double e=hT->GetBinError(iF+1);
                            int bP = hPWTru[k]->FindBin(wMean);
                            int bPb = hPWTruBin[k*nTrueFlat()+ft]->FindBin(wMean);
                            double bPe = hPWTru[k]->GetBinError(bP);
                            double bPbe = hPWTruBin[k*nTrueFlat()+ft]->GetBinError(bPb);
                            hPWTru[k]->Fill(wMean,v);
                            hPWTruBin[k*nTrueFlat()+ft]->Fill(wMean,v);
                            hPWTru[k]->SetBinError(bP,std::hypot(bPe,e));
                            hPWTruBin[k*nTrueFlat()+ft]->SetBinError(bPb,std::hypot(bPbe,e));
                    }
                }
        }

        // ── Plot 1: 4×8 ΔΦ grid ──────────────────────────────────────────
        const int nColsPW=4, nRowsPW=(nDphi+nColsPW-1)/nColsPW;
        TCanvas* cPWGrid=new TCanvas("cPWGrid","PW closure dPhi",nColsPW*250,nRowsPW*500);
        cPWGrid->Divide(nColsPW,nRowsPW,0,0);
        for (int k=0; k<nDphi; ++k) {
            StyleLine(hPWTru[k],kColTruth,kMarTruth);
            StyleLine(hPWUnf[k],kColUnfolded,kMarUnfolded);
            auto [pTop,pBot]=MakeRatioPads(cPWGrid->cd(k+1));
            pTop->cd(); pTop->SetLogy(); pTop->SetLogx();
            double ymax=std::max(hPWTru[k]->GetMaximum(),hPWUnf[k]->GetMaximum());
            double ymin=1e30;
            for (int b=1;b<=hPWTru[k]->GetNbinsX();++b)
                if (hPWTru[k]->GetBinContent(b)>0) ymin=std::min(ymin,hPWTru[k]->GetBinContent(b));
            if (ymin>ymax) ymin=ymax*1e-4;
            hPWTru[k]->GetYaxis()->SetRangeUser(ymin*0.3,ymax*3.0);
            hPWTru[k]->GetXaxis()->SetLabelSize(0);
            hPWTru[k]->GetYaxis()->SetTitle("Pairs");
            hPWTru[k]->GetYaxis()->SetTitleSize(0.07);
            hPWTru[k]->SetTitle(std::format("{:.2f}-{:.2f}",dPhiBins[k],dPhiBins[k+1]).c_str());
            hPWTru[k]->Draw("E"); hPWUnf[k]->Draw("E SAME");
            TLegend* lg=new TLegend(0.55,0.72,0.92,0.88);
            lg->SetTextSize(0.07);
            lg->AddEntry(hPWTru[k],"Truth","lp"); lg->AddEntry(hPWUnf[k],"Unfolded","lp"); lg->Draw();
            pBot->cd(); pBot->SetLogx();
            TH1D* hR=MakeRatioHist(hPWUnf[k],hPWTru[k],std::format("hPWRatio_{}",k).c_str());
            StyleRatio(hR,"Pair weight"); hR->Draw("E"); DrawRefLine(hR);
        }
        cPWGrid->SaveAs(std::format("{}/pairWeight_closure_dPhi-{}.png",plotDir,label).c_str());
        cPWGrid->Close(); delete cPWGrid;

        // ── Plot 2: PDF — one page per ΔΦ bin, nTrueLead×nTrueSubl grid ──
        std::string pdfName=std::format("{}/pairWeight_closure_perBin-{}.pdf",plotDir,label);
        for (int k=0; k<nDphi; ++k) {
            TCanvas* cPage=new TCanvas(
                std::format("cPWPage_{}",k).c_str(),"",nTrueLead*220,nTrueSubl*220);
            cPage->Divide(nTrueLead,nTrueSubl,0,0);
            for (int iS=0; iS<nTrueSubl; ++iS)
            for (int iL=0; iL<nTrueLead; ++iL) {
                int ft=TrueFlatIndex(iL,iS);
                auto [pTop,pBot]=MakeRatioPads(cPage->cd(iS*nTrueLead+iL+1));
                std::string bl=std::format("{:.0f}-{:.0f}/{:.0f}-{:.0f} GeV",
                    trueLeadPtBins[iL],trueLeadPtBins[iL+1],
                    trueSublPtBins[iS],trueSublPtBins[iS+1]);
                // Gray out kinematically forbidden bins
                if (ft < 0) {
                    gPad->SetFillColor(kGray);
                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.08);
                    tx.DrawLatex(0.15,0.85,bl.c_str()); continue;
                }
                TH1D* hU=hPWUnfBin[k*nTrueFlat()+ft];
                TH1D* hT=hPWTruBin[k*nTrueFlat()+ft];
                StyleLine(hT,kColTruth,kMarTruth);
                StyleLine(hU,kColUnfolded,kMarUnfolded);
                // Draw all bins, even those outside the reporting region —
                // gray only kinematically forbidden bins (ft < 0, handled above)
                // and bins with no truth entries (empty histogram).
                if (hT->Integral()==0) {
                    gPad->SetFillColor(kGray+1);  // darker gray = truth-only / empty
                    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.08);
                    tx.DrawLatex(0.15,0.85,bl.c_str()); continue;
                }
                pTop->cd(); pTop->SetLogy(); pTop->SetLogx();
                double ymax=std::max(hT->GetMaximum(),hU->GetMaximum());
                double ymin=1e30;
                for (int b=1;b<=hT->GetNbinsX();++b)
                    if (hT->GetBinContent(b)>0) ymin=std::min(ymin,hT->GetBinContent(b));
                if (ymin>ymax) ymin=ymax*1e-4;
                hT->GetYaxis()->SetRangeUser(ymin*0.3,ymax*3.0);
                hT->GetXaxis()->SetLabelSize(0);
                hT->GetYaxis()->SetTitleSize(0.08); hT->GetYaxis()->SetTitle("Pairs");
                hT->SetTitle(""); hT->Draw("E"); hU->Draw("E SAME");
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.08); tx.SetTextAlign(13);
                tx.DrawLatex(0.15,0.88,bl.c_str());
                pBot->cd(); pBot->SetLogx();
                TH1D* hR=MakeRatioHist(hU,hT,std::format("hPWBinR_{}_{}",k,ft).c_str());
                StyleRatio(hR,"Pair weight"); hR->Draw("E"); DrawRefLine(hR);
            }
            cPage->cd();
            TLatex ct; ct.SetNDC(); ct.SetTextSize(0.015); ct.SetTextAlign(22);
            ct.DrawLatex(0.5,0.002,std::format("#Delta#phi: {:.2f}-{:.2f} (bin {})",
                dPhiBins[k],dPhiBins[k+1],k).c_str());
            std::string pageOpt=(k==0)?(pdfName+"("):((k==nDphi-1)?(pdfName+")"):pdfName);
            cPage->SaveAs(pageOpt.c_str());
            cPage->Close(); delete cPage;
        }

        // Cleanup pair weight histograms
        for (int k=0; k<nDphi; ++k) {
            delete hPWUnf[k]; delete hPWTru[k];
            for (int ft=0; ft<nTrueFlat(); ++ft) {
                delete hPWUnfBin[k*nTrueFlat()+ft];
                delete hPWTruBin[k*nTrueFlat()+ft];
            }
        }
    }


    // ════════════════════════════════════════════════════════════════════════
    //  Response matrix visualization
    //  Two plots, both read from fResp:
    //
    //  respMatrix_pairWeight-{label}.png
    //    Grid of 8 representative ΔΦ bins showing the pair weight response:
    //    the nPairWeight × nPairWeight matrix obtained by marginalizing the
    //    3D response over all (iL, iS) dijet pT bins.  Reco pair weight on X,
    //    truth pair weight on Y.  Column-normalized so each reco bin sums to 1,
    //    making the migration probability visible.
    //
    //  respMatrix_dijetPt-{label}.png
    //    The nFlat × nFlat = 49 × 49 dijet pT response matrix (lead pT × subl pT
    //    in flat index) summed over all ΔΦ bins and pair weight bins.
    //    Column-normalized.  Shows the jet energy scale migration.
    // ════════════════════════════════════════════════════════════════════════
    if (fResp) {
        // ── Full 3D response matrix grid: all ΔΦ bins ────────────────────
        // Each pad shows the full nFlat3D×nFlat3D response matrix for one
        // ΔΦ bin, with reco flat3D index on X and truth flat3D index on Y.
        // Saved as PDF for lossless high-resolution rendering.
        const int nCols3D = 4, nRows3D = (nDphi + nCols3D - 1) / nCols3D;
        const int pad3D = 400;
        TCanvas* cResp = new TCanvas("cResp", "Full 3D response matrix",
                                     nCols3D*pad3D, nRows3D*pad3D);
        cResp->Divide(nCols3D, nRows3D, 0, 0);

        //vector<TH2D*> hDraw(nDphi, nullptr);

        for (int k = 0; k < nDphi; ++k) {
            cResp->cd(k+1);
            gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
            gPad->SetRightMargin(0.15); gPad->SetTopMargin(0.10);
            gPad->SetLogz();

            RooUnfoldResponse* resp =
                (RooUnfoldResponse*)fResp->Get(
                    std::format("response_wEEC3D_{}", k).c_str());
            if (!resp) continue;

            TH2D* hFull = (TH2D*)resp->Hresponse();
            if (!hFull || hFull->Integral() == 0) continue;

            // Clone so we can set title/style without modifying the stored object
            TH2D *hDraw = (TH2D*)hFull->Clone(
                std::format("hResp3D_{}",k).c_str());
            hDraw->SetDirectory(0);
            hDraw->SetTitle(
                std::format("{:.2f}<#Delta#phi<{:.2f};Reco flat3D;Truth flat3D",
                    dPhiBins[k], dPhiBins[k+1]).c_str());
            hDraw->GetXaxis()->SetTitleSize(0.05);
            hDraw->GetYaxis()->SetTitleSize(0.05);
            hDraw->GetXaxis()->SetLabelSize(0.04);
            hDraw->GetYaxis()->SetLabelSize(0.04);
            hDraw->GetZaxis()->SetLabelSize(0.04);
            //double zMax = hDraw->GetMaximum();
            //if (zMax > 0) hDraw->GetZaxis()->SetRangeUser(zMax*1e-4, zMax);
            hDraw->Draw("COLZ");

            TLatex tex; tex.SetNDC(); tex.SetTextSize(0.06); tex.SetTextAlign(13);
            tex.DrawLatex(0.14, 0.92,
                std::format("{:.2f}-{:.2f}",
                    dPhiBins[k], dPhiBins[k+1]).c_str());
            //delete hDraw;
        }

        cResp->SaveAs(std::format("{}/respMatrix_full3D-{}.pdf",
                                  plotDir, label).c_str());

        cResp->Close(); delete cResp;

        if(doTrim)
        {
            TCanvas* cRespTrimmed = new TCanvas("cRespTrimmed", "Full 3D response matrix Trimmed",
                                        nCols3D*pad3D, nRows3D*pad3D);
            cRespTrimmed->Divide(nCols3D, nRows3D, 0, 0);

            //vector<TH2D*> hDraw(nDphi, nullptr);

            for (int k = 0; k < nDphi; ++k) {
                cRespTrimmed->cd(k+1);
                gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
                gPad->SetRightMargin(0.15); gPad->SetTopMargin(0.10);
                gPad->SetLogz();

                RooUnfoldResponse* resp =
                    (RooUnfoldResponse*)fWEEC->Get(
                        std::format("response_wEEC3D_trimmed_{}", k).c_str());
                if (!resp) continue;

                TH2D* hFull = (TH2D*)resp->Hresponse();
                if (!hFull || hFull->Integral() == 0) continue;

                // Clone so we can set title/style without modifying the stored object
                TH2D *hDraw = (TH2D*)hFull->Clone(
                    std::format("hResp3D_trimmed_{}",k).c_str());
                hDraw->SetDirectory(0);
                hDraw->SetTitle(
                    std::format("{:.2f}<#Delta#phi<{:.2f};Reco flat3D;Truth flat3D",
                        dPhiBins[k], dPhiBins[k+1]).c_str());
                hDraw->GetXaxis()->SetTitleSize(0.05);
                hDraw->GetYaxis()->SetTitleSize(0.05);
                hDraw->GetXaxis()->SetLabelSize(0.04);
                hDraw->GetYaxis()->SetLabelSize(0.04);
                hDraw->GetZaxis()->SetLabelSize(0.04);
                //double zMax = hDraw->GetMaximum();
                //if (zMax > 0) hDraw->GetZaxis()->SetRangeUser(zMax*1e-4, zMax);
                hDraw->Draw("COLZ");

                TLatex tex; tex.SetNDC(); tex.SetTextSize(0.06); tex.SetTextAlign(13);
                tex.DrawLatex(0.14, 0.92,
                    std::format("{:.2f}-{:.2f}",
                        dPhiBins[k], dPhiBins[k+1]).c_str());
                //delete hDraw;
            }

            cRespTrimmed->SaveAs(std::format("{}/respMatrixTrimmed_full3D-{}.pdf",
                                    plotDir, label).c_str());
        }
        
        // ── Cleanup ───────────────────────────────────────────────────────────
        fWEEC->Close(); delete fWEEC;

        // ── dijet pT response: sum over all ΔΦ and pair weight bins ──────
        TH2D* hJetResp = new TH2D("hJetResp",
            "Dijet pT response (from wEEC resp, marginalized);"
            "Reco flat dijet pT bin;Truth flat dijet pT bin",
            nRecoFlat(), 0, nRecoFlat(), nTrueFlat(), 0, nTrueFlat());
        hJetResp->SetDirectory(0);

        for (int k=0; k<nDphi; ++k) {
            RooUnfoldResponse* resp =
                (RooUnfoldResponse*)fResp->Get(
                    std::format("response_wEEC3D_{}", k).c_str());
            if (!resp) continue;
            TH2D* hFull = (TH2D*)resp->Hresponse();
            if (!hFull || hFull->Integral() == 0) continue;

            for (int ftR=0; ftR<nRecoFlat(); ++ftR)
            for (int ftT=0; ftT<nTrueFlat(); ++ftT) {
                double sum = 0;
                for (int rIPw=0; rIPw<nPairWeight; ++rIPw)
                for (int tIPw=0; tIPw<nPairWeight; ++tIPw) {
                    int rFlat = ftR * nPairWeight + rIPw;
                    int tFlat = ftT * nPairWeight + tIPw;
                    sum += hFull->GetBinContent(rFlat+1, tFlat+1);
                }
                if (sum > 0)
                    hJetResp->AddBinContent(
                        hJetResp->GetBin(ftR+1, ftT+1), sum);
            }
        }

        // Column-normalize
        for (int rBin=1; rBin<=nRecoFlat(); ++rBin) {
            double colSum=0;
            for (int tBin=1; tBin<=nTrueFlat(); ++tBin)
                colSum += hJetResp->GetBinContent(rBin,tBin);
            if (colSum>0)
                for (int tBin=1; tBin<=nTrueFlat(); ++tBin)
                    hJetResp->SetBinContent(rBin,tBin,
                        hJetResp->GetBinContent(rBin,tBin)/colSum);
        }

        TCanvas* cJR = new TCanvas("cJR","Dijet pT response (marginalized)",800,700);
        cJR->cd(); gPad->SetLogz();
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);
        gPad->SetRightMargin(0.15); gPad->SetTopMargin(0.08);
        double jrMax = hJetResp->GetMaximum();
        if (jrMax > 0) hJetResp->GetZaxis()->SetRangeUser(jrMax*1e-3, jrMax);
        hJetResp->GetXaxis()->SetTitle("Reco flat dijet pT bin");
        hJetResp->GetYaxis()->SetTitle("Truth flat dijet pT bin");
        hJetResp->Draw("COLZ");
        cJR->SaveAs(std::format("{}/respMatrix_dijetPt-{}.png",
                                plotDir,label).c_str());
        cJR->Close(); delete cJR;
        delete hJetResp;

        fResp->Close(); delete fResp;
    }

    for (int ft=0; ft<nTrueFlat(); ++ft) {
        delete hTruthBin[ft];
        delete hUnfBin[ft];
    }
    for (int ft=0; ft<nRecoFlat(); ++ft)
        delete hMeasBin[ft];

    std::cout << "Done.\n";
    gSystem->Exit(0);
}