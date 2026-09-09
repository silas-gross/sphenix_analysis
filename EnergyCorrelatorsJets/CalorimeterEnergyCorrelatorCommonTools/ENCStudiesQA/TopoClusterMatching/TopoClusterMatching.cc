#include "TopoClusterMatching.h"

TopoClusterMatching::TopoClusterMatching(float mpT, const std::string name)
{
	minpT=mpt;
	getBins(mpt);
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
int TopoClusterMatching::process_event(PHCompositeNode* topNode)
{
	if(verbose > 1) std::cout<<"event number: " <<n_evt<<std::endl;
	n_evt++;

	bool isDijet = runDijetCut(topNode);	
	if(!isDijet) return Fun4AllReturnCodes::EVENT_OK;

	std::pair<int, int> jet_bin_index	= findBucket(event_cut->getLeadPt(), event_cut->getSubleadPt());
	
	return Fun4AllReturnCodes::EVENT_OK;
}
