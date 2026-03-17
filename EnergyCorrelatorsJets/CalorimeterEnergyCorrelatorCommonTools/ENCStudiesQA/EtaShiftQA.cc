EtaShiftQA::EtaShiftQA()
{
	EMCALQA = new PerCaloQAPlots("EMCAL");
	IHCALQA	= new PerCaloQAPlots("IHCAL");
	OHCALQA = new PerCaloQAPlots("OHCAL");
	MetaEQA = new PerCaloQAPlots("Meta Calo, EMCAL radius");
	MetaIQA = new PerCaloQAPlots("Meta Calo, IHCAL radius");
	MetaOQA = new PerCaloQAPlots("Meta Calo, OHCAL radius");
	
	caluclatedJetpt 	= new TH1F(
				"h_calc_Jet_pt", 
				"Jet p_{T} calculated from unshifted consituents; p_{T}^{jet}[GeV]; N_{jet}",
				50, -1, 99
			);
	caluclatedShiftedJetpt 	= new TH1F(
				"h_shift_calc_Jet_pt", 
				"Jet p_{T} calculated from shifted consituents; p_{T}^{jet} [GeV]; N_{jet}",
				50, -1, 99
			);
	hzVTX			= new TH1F(
				"h_zvtx", 
				"Vertex position; z_{vtx} [cm]; N_{events}"
				300, -150.5, 149.5
			);
				
	for(int i = 0; i<6; i++)
	{
		int z_range = (i+1)*10;
		EMCAL_Z_QA->at(i)	= new PerCaloQAPlots(std::format("EMCAL, |z|_{vtx} < {}", z_range));
		IHCAL_Z_QA->at(i) 	= new PerCaloQAPlots(std::format("IHCAL, |z|_{vtx} < {}", z_range));
		OHCAL_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("OHCAL, |z|_{vtx} < {}", z_range));
		MetaE_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, EMCAL radius, |z|_{vtx} < {}", z_range));
		MetaI_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, IHCAL radius, |z|_{vtx} < {}", z_range));
		MetaO_Z_QA->at(i)  	= new PerCaloQAPlots(std::format("Meta Calo, OHCAL radius, |z|_{vtx} < {}", z_range));
	}
}
