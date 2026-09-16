#include "TopoClusterMatching.h"

TopoClusterMatching::TopoClusterMatching(float mpT, const std::string name)
{
	minpT=mpt;
	getBins(mpt);
	h_ClusterPairEtAll=new TH1F(
			"h_ClusterPairEtAll", "Cluster Pair E_{T}; E_{T, i} #times E_{T,j} / < E_{T, dijet} >^{2}; N_{pair}",
			clusterPtBin.size(), clusterPtBin.data());
	h_TruthPairEtAll=new TH1F(
			"h_TruthPairEtAll", "Truth Pair E_{T}; E_{T, i} #times E_{T,j} / < E_{T, dijet} >^{2}; N_{pair}",
			clusterPtBin.size(), clusterPtBin.data());
	for(int i=0; i<(int)sub_bucket.size()-1; i++)
	{
		for(int j=0; j<(int)sub_bucket.at(i).size()-1; j++)
		{
			h_ClusterPairAllEt_Div[i][j]=new TH1F(
					std::format("h_ClPairETAll_{}_{}", lead_bucket[i], sub_bucket[i][j]).c_str();
					Form("Cluster Pair E_{T}, %d< p_{T, lead}^{jet} < %d, %d< p_{T, sub}^{jet}; E_{T, i} #times E_{T,j} / < E_{T, dijet} >^{2}; N_{pair}", lead_bucket[i], lead_bucket[i+1], sub_bucket[i][j], sub_bucket[i][j]).c_str(),
					clusterPtBin.size(), clusterPtBin.data());
			h_TruthPairEtAll_Div[i][j]=new TH1F(
					std::format("h_TrPairETAll_{}_{}", lead_bucket[i], sub_bucket[i][j]).c_str();
					Form("Truth Pair E_{T}, %d< p_{T, lead}^{jet} < %d, %d< p_{T, sub}^{jet}; E_{T, i} #times E_{T,j} / < E_{T, dijet} >^{2}; N_{pair}", lead_bucket[i], lead_bucket[i+1], sub_bucket[i][j], sub_bucket[i][j]).c_str(),
					clusterPtBin.size(), clusterPtBin.data());
		}
	}
	h_ClusterPairEt=new TH1F(
			"h_ClusterPairEt", "Cluster Pair E_{T}; E_{T, i} #times E_{T,j} / < E_{T, dijet} >^{2}; N_{pair}",
			clusterPtBin.size(), clusterPtBin.data());
	h_TruthPairEt=new TH1F(
			"h_TruthPairEt", "Truth Pair E_{T}; E_{T, i} #times E_{T,j} / < E_{T, dijet} >^{2}; N_{pair}",
			clusterPtBin.size(), clusterPtBin.data());
	for(int i=0; i<(int)sub_bucket.size()-1; i++)
	{
		for(int j=0; j<(int)sub_bucket.at(i).size()-1; j++)
		{
			h_ClusterPairEt_Div[i][j]=new TH1F(
					std::format("h_ClPairET_{}_{}", lead_bucket[i], sub_bucket[i][j]).c_str();
					Form("Cluster Pair E_{T}, %d< p_{T, lead}^{jet} < %d, %d< p_{T, sub}^{jet}; E_{T, i} #times E_{T,j} / < E_{T, dijet} >^{2}; N_{pair}", lead_bucket[i], lead_bucket[i+1], sub_bucket[i][j], sub_bucket[i][j]).c_str(),
					clusterPtBin.size(), clusterPtBin.data());
			h_TruthPairEt_Div[i][j]=new TH1F(
					std::format("h_TrPairET_{}_{}", lead_bucket[i], sub_bucket[i][j]).c_str();
					Form("Truth Pair E_{T}, %d< p_{T, lead}^{jet} < %d, %d< p_{T, sub}^{jet}; E_{T, i} #times E_{T,j} / < E_{T, dijet} >^{2}; N_{pair}", lead_bucket[i], lead_bucket[i+1], sub_bucket[i][j], sub_bucket[i][j]).c_str(),
					clusterPtBin.size(), clusterPtBin.data());
		}
	}
					
	event_cut = new DijetEventCuts(); //require a leading jet of 12 GeV sublead 7 GeV, keep it in |eta|<0.7, set dPhi > 3 pi/4
	
}

void TopoCluterMatching::findMatchingCluster(std::vector<CandidateObj*> cls, CandidateObj* trpt ) 
{
	for(int i = 0; i<(int)cls.size(); i++)
	{
		CandidateObj cl = cls[i];
		bool goodpos = cl->isIndR(trpt);
		bool goodE = cl->isE(trpt);
		if(goodE && goodpos)
		{
			trpt->SetMatch(i);
			break;
		}
		else continue;
	}
	return;

}

void TopoClusterMatching::getBins(float minpt)
{
	//building the logarithmic bins 
	int nbins 	= 100;
	float min 	= std::pow(minpt/30., 2); //use 30 GeV as Q approximation for this purpose
	float max 	= 0.5;
       	float binwidth 	= std::log(max) - std::log(min);	
	binwidth	= binwidth/(float)nbins; //linear spacing in log
	clusterPtBin.push_back(1e-8);
	cluterPtBin.push_back(min);
	for(int i=0; i<nbins; i++)
	{
		float logbinlow	= std::log(clusterPtBin[i]);
		logbinlow 	= logbinlow + binwidth;
		
		clusterPtBin.push_back(std::pow(10, logbinlow));
	}
	clusterPtBin.push_back(1.);

	//create the buckets for which subleading and leading jets we have 
	int nbuckets 	= 10; 
	float minLJet	= 12.;
	float minSLRat	= minLJet/7.;
	float maxJet	= 100.;
	float bucketW	= (maxJet - minLJet)/((float)nbuckets);
	
	for(int i=0; i<nbuckets; i++)
	{
		float jetbucket = minLJet + bucketW*i;
		lead_bucket.push_back(jetbucket);
		std::vector<float> tempbins;
		float subw = bucketW*((i+1) - minSLRat*i);
		for(int j=0; j<nbuckets; j++)
		{
			tempbins.push_back(jetbucket * minSLRat + subw * j);
		}
		sub_bucket.push_back(tempbins);
	}
			
	return;

}

float TopoClusterMatching::getR(CandidateObj* o1, CandidateObj* o2)
{
	float deta = o1->eta - o2->eta;
	float dphi = o1->phi - o2->phi;
	if(dphi > M_PI) dphi += - 2*M_PI; 
	dphi= std::abs(dphi);
	float R2 = std::pow(dphi,2) + std::pow(deta, 2);
	return std::sqrt(R2);
}

bool TopoClusterMatching::runDijetCut(PHCompositeNode* topNode)
{
	auto jets 	= findNode::getClass<JetContainer>(topNode, jet_node_name);
	bool isGood 	= event_cut->passesTheCut(jets);	
	return isGood;
}
std::pair<int, int> TopoClusterMatching::findBucket(float lpt, float slpt)
{
	//this just figures out which bin of the lead-sublead pt matrix to fill in 
	int l_i = 0; 
	int s_i	= 0;
 
	for(int i = 0; i<(int)lead_buckets.size(); i++)
		if ( lpt < lead_buckes.at(i)) l_i = i-1;
	for(int i = 0; i<(int)sublead_buckets.at(l_i).size(); i++)
		if ( slpt < sublead_buckets.at(i)) s_i = i-1;
	return std::make_pair(l_i, s_i);
}
bool TopoClusterMatching::isAMatch(CandidateObj* t, CandidateObj* c)
{
	bool matched = false;
	float erat = t->E / c->E;
	float dR = getR(t, c);
	if(erat > eratmin && dR < dRmax) matched = true;
	t->isMatch = matched;
	c->isMatch = matched;
	return matched;
}
void TopoClusterMatching::MatchAllTruthToClusters(PHCompositeNode* topNode) 
{
	auto truthMap = getClass::findNode<PHG4TruthInfoContainer*>(topNode, "G4TruthInfo");
	std::vector<CandidateObj*> valid_truth {};
	if(!truthinfo) return;
	for(
		auto iter = truthinfo->GetSPHENIXPrimaryParticleRange().first; 
		iter != truthinfo->GetSPHENIXPrimaryParticleRange().second; 
		++iter
	   )
	{
		if(!iter) continue;
		PHG4Particle* p = iter->second;
		if(!p) continue;
		bool goodKin = isGoodKin(p);
		if(goodKin){
			CandidateObj* pt = new CandidateObj(p);
			valid_truth.push_back(pt);
			continue;
		}
		else continue;
	}
	std::vector<CandidateObj*> valid_topo {};
	auto clusters = findNode::getClass<RawClusterContainer>(topNode, "TOPOCLUSTER_ALLCALO");
	if(!clusters)return;
	for(
		auto iter=clusters->getClusters().begin;
		iter !=clusters->getClusters().end;
		++iter
	   )
	{
		if(!iter) continue;
		bool goodKin = isGoodKin(*iter);
		if(goodKin)
		{
			CandidateObj* cl = new CandidateObj(*iter);
			valid_topo.push_back(pt);
			continue;
		}
		else continue;
	}
	for(auto t:valid_truth)
		for(auto c:valid_topo)
		{
			bool matched = isAMatch(t, c);
			if(isAMatch)
			{
				h_MatchedTruth->Fill(t->pt);
			      	h_RealCluster->Fill(c->pt);	
				break;
			}
			else continue;
		}
	return;

}
void TopoClusterMatching::getEEC(std::vector<CandidateObj*> objs, bool isTruth, std::pair<int, int> jet_bin_index)
{
	for(int i=0; i<(int)objs.size()-1; i++)
	{
	
//		if(objs.at(i)->isMatched == false) continue;
		for(int j=i+1; j<(int)objs.size()-1; i++)
		{
//			if(objs.at(j)->isMatched == false) continue;
			float pairEt = objs[i]->Et * ojbs[j]->Et;
			if(isTruth){
				if(objs[j]->isMatched && objs[i]->isMatched){
			       		h_TruthPairEt->Fill(pairEt/Q2);
					h_TruthPairEt_Div[jet_bin_index.first][jet_bin_index.second]->Fill(pairEt/Q2);
				}
			       	h_TruthPairEtAll->Fill(pairEt/Q2);
				h_TruthPairEtAll_Div[jet_bin_index.first][jet_bin_index.second]->Fill(pairEt/Q2);
			}
			else
			{
				if(objs[j]->isMatched && objs[i]->isMatched){
					h_ClusterPairEt->Fill(pairEt/Q2);
					h_ClusterPairEt_Div[jet_bin_index.first][jet_bin_index.second]->Fill(pairEt/Q2);
				}
				h_ClusterPairEtAll->Fill(pairEt/Q2);
				h_ClusterPairEtAll_Div[jet_bin_index.first][jet_bin_index.second]->Fill(pairEt/Q2);
			}
		}
	}
	return;
}
bool TopoClusterMatching::KinCuts(PHG4Particle* p)
{
	bool kingood {false};
	if(!p) return kingood;
	float px {p->get_px()};
	float py {p->get_py()};
	float pz {p->get_pz()};
	float e	 {p->get_e()};
	float eta { std::atanh(pz / e)}; 
	bool isEM {false};
	if( std::abs(eta) <= 1.1)
	{
		int pid = std::abs(p->get_pid());
		if( pid == 11 || pid == 13 || pid == 22) 
			isEM = true;
		else if( pid < 11 || pid > 16 ) 
			isEM = false;
		else return kingood;
		float threhold = isEM ? 0.2 : 0.5;
		if(e > threshold) kingood = true;
	}
	return kingood;
}
bool TopoClusterMatching::KinCuts(RawCluster* p)
{
	bool kingood {false};
	if(!p) return kingood;
	float phi {p->get_phi()};
	float e	 {p->get_energy()};
	float eta {p->get_eta}; 
	if( std::abs(eta) <= 1.1)
		if(e > minpt) kingood = true;
	return kingood;
}
	 
int TopoClusterMatching::process_event(PHCompositeNode* topNode)
{
	if(verbose > 1) std::cout<<"event number: " <<n_evt<<std::endl;
	n_evt++;

	bool isDijet = runDijetCut(topNode);	
	if(!isDijet) return Fun4AllReturnCodes::EVENT_OK;

	jet_bin_index	= findBucket(event_cut->getLeadPt(), event_cut->getSubleadPt());
	Q2 = 0.5* ( std::pow(event_cut->getLeadPt(), 2) + std::pow(event_cut->getLeadPt(), 2));
	MatchAllTruthToClusters(topNode, jet_bin_index);	
	return Fun4AllReturnCodes::EVENT_OK;
}

