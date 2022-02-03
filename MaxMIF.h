#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include <algorithm>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<list>
#include <math.h>
#include <string.h>
 
using namespace std; 

#define STR_SEPTABLE "\t"
#define SMAllESTNUM 1e-15

/* Create Class named MaxMIF */
class MaxMIF
{
	typedef struct PPI
	{
		int Id;
		double Idw;
		PPI(int _id,double _idw)
		{
			this->Id=_id;
		    this->Idw=_idw;
		}
	}PPI;
	
public:
	void ReadCancerGene();
	void ReadGeneID();
	bool ReadMutData(string ifspathstr);
	void ComputeMut();
	bool ReadNetwork(string ifspathstr);
	void ComputeRank();
    void OutputRank(string rankName);
	void MutDataClear();
	void NetworkClear();
	
private:
	void StringSplit(string longstr,string str1, vector<string>& vectstr);
	
private:
	int m_intgid, m_intgid1,m_intSampNum;
	double m_doubMinWeight,m_doubMinMutScore;
	string m_filepath,m_strtemp,m_genename;
	vector<string> m_vecTemp,m_vectUNMutGenes;
	map<int,string> m_mapIdGene;
	map<string,int> m_mapGeneId;
	map<string,vector<int> > m_geneMutData;
	map<int,list<PPI> > m_netWork;
	map<int,list<PPI> >::iterator iterPPI;
	map<string,double> m_mapGeneMut;
	multimap<double,string> m_GeneScore;
	set<int> m_CancerGene[5];
	ifstream m_ifs;
	ofstream m_ofs;	
};

/* the functions of MaxMIF */

// Read reference cancer gene sets
void MaxMIF::ReadCancerGene()
{
	string path[5]={"Background//CancerGeneSet//CGC","Background//CancerGeneSet//CGCpointMut",
	"Background//CancerGeneSet//Rule2020","Background//CancerGeneSet//HCD",
	"Background//CancerGeneSet//MouseMut"};
	for(int k = 0; k < 5;++k)
	{
		m_ifs.open(path[k].c_str());
		while(getline(m_ifs,m_strtemp))
		{
			m_intgid=atoi(m_strtemp.c_str());
			m_CancerGene[k].insert(m_intgid);
		}
		m_ifs.close();
	}
}

// Read gene Entrez ID file
void MaxMIF::ReadGeneID()
{
	string idpath="Background//geneID";
	m_ifs.open(idpath.c_str());
	while(getline(m_ifs, m_strtemp))
	{
		StringSplit(m_strtemp, STR_SEPTABLE,  m_vecTemp);
		m_intgid=atoi(m_vecTemp[1].c_str());
		m_mapGeneId[m_vecTemp[0]]=m_intgid;
		m_mapIdGene[m_intgid]=m_vecTemp[0];
	}
	m_ifs.close();
}

void MaxMIF::StringSplit(string longstr,string sep, vector<string>& vectstr)
{
	vectstr.clear();
	int nend=0;
    int nbegin=0;
    while(nend != string::npos)
    {
        nend = longstr.find(sep, nbegin);
        if(nend == string::npos)
            vectstr.push_back(longstr.substr(nbegin, longstr.length()-nbegin));
        else
            vectstr.push_back(longstr.substr(nbegin, nend-nbegin));
        nbegin = nend + sep.length();
    }
}

// Read the gene-level somatic mutation data
bool MaxMIF::ReadMutData(string ifspathstr)
{
	m_ifs.open(ifspathstr.c_str());
	if(m_ifs == NULL)
	{
		cout<<"\n\n"<<"Fail to open the Mutation data: "<<ifspathstr<<"! \n"
		<<"Please check if your input is correct and try again later!"<<"\n";
		return false;
	}
	getline(m_ifs,m_strtemp);
	StringSplit(m_strtemp, STR_SEPTABLE,  m_vecTemp);
	m_intSampNum=m_vecTemp.size()-1;
	// if the input data is  mutation frequency
	if( m_intSampNum == 1)
	{
		cout<<"\n"<<"Type: Mutation frequency data ..."<<"\n";
		while(getline(m_ifs,m_strtemp))
		{
			StringSplit(m_strtemp, STR_SEPTABLE,  m_vecTemp);
			m_genename=m_vecTemp[0].c_str();
			float gmut = atof(m_vecTemp[1].c_str());
			if(gmut != 0.0f)
			{
				m_mapGeneMut[m_genename]=gmut;
			}	
		}
		m_ifs.close();
		return true;
	}
	// if the input data is  mutation matrix 
	cout<<"\n"<<"Type: Mutation matrix data; Sample size: "<<m_intSampNum<<"\n";
	vector<int> vectinttemp;
	string s1="1",s2="1.0";
	while(getline(m_ifs,m_strtemp))
	{
		StringSplit(m_strtemp, STR_SEPTABLE,  m_vecTemp);
		m_genename=m_vecTemp[0].c_str();
		for (int j=0;j< m_intSampNum;++j)
		{
			if(m_vecTemp[j+1].c_str()==s1 || m_vecTemp[j+1].c_str()==s2)
			{
				vectinttemp.push_back(j);
			}
		}
		// symbol genes
		if(!vectinttemp.empty())
		{
			m_geneMutData[m_genename]=vectinttemp;
			vectinttemp.clear();
		}
		else
		{
			// deal with special genes
			map<string,int>::iterator itergid=m_mapGeneId.find(m_genename);
			if(itergid != m_mapGeneId.end())
				m_vectUNMutGenes.push_back(m_genename);
		}
	}
	m_ifs.close();
	return true;
}

// Compute the mutation-score of every mutated genes
void MaxMIF::ComputeMut()
{
	vector<double> vectSampMutNum;
	map<string,vector<int> >::iterator iter1;
	int maxMutNum=1;
	for(int j=0;j<m_intSampNum;++j)
	{
		double sMutNum=0;
		for(iter1=m_geneMutData.begin();iter1 !=m_geneMutData.end(); iter1++)
		{
			m_genename=iter1->first;
			vector<int> vecSampleID= iter1-> second;
			vector<int>::iterator iter2;
			iter2=find(vecSampleID.begin(),vecSampleID.end(),j);
			if (iter2!=vecSampleID.end())
			{
				sMutNum++;
			}
		}
		if(sMutNum !=0)
		{
			vectSampMutNum.push_back((double)1/sMutNum);
			if(maxMutNum < sMutNum)
				maxMutNum = sMutNum;
		}
		else
			vectSampMutNum.push_back(0);	
	}
	
	for(iter1=m_geneMutData.begin();iter1 !=m_geneMutData.end(); iter1++)
	{
		double gidMut=0;
		m_genename=iter1->first;
		vector<int> vecSampleID= iter1-> second;
		vector<int>::iterator iter3=vecSampleID.begin();
		for(; iter3 !=vecSampleID.end();iter3++)
		{
			int samplabel=*iter3;
			gidMut=gidMut+vectSampMutNum[samplabel];	
		}
		m_mapGeneMut[m_genename]=gidMut;
	}
	vectSampMutNum.clear();
	m_geneMutData.clear();
	// genes without mutation 
	m_doubMinMutScore = (double)1/(maxMutNum);
	////  Deal with genes who are not mutated
	for(int k=0;k < m_vectUNMutGenes.size(); k++ )
	{
		m_genename = m_vectUNMutGenes[k];
		m_mapGeneMut[m_genename]= m_doubMinMutScore;
	}
	m_vectUNMutGenes.clear();
}

// Read functional network
bool MaxMIF::ReadNetwork(string ifspathstr)
{
	m_doubMinWeight=0.01;  // set a small weight;
	double ggweight;
	m_ifs.open(ifspathstr.c_str());
	if(m_ifs == NULL)
	{
		cout<<"\n"<<"Fail to open the Network data: "<<ifspathstr<<"! \n"
		<<"Please check if your input is correct and try again later!"<<"\n";
		return false;
	}
	while(getline(m_ifs,m_strtemp))
	{
		StringSplit(m_strtemp, STR_SEPTABLE,  m_vecTemp);
		m_intgid=atoi(m_vecTemp[0].c_str());
		m_intgid1=atoi(m_vecTemp[1].c_str());
		map<int,string>::iterator iteridg,itergid1g;
		iteridg=m_mapIdGene.find(m_intgid);
		itergid1g=m_mapIdGene.find(m_intgid1);
		if((iteridg==m_mapIdGene.end() ) | (itergid1g==m_mapIdGene.end() ) )
			continue;
		ggweight=atof(m_vecTemp[2].c_str());
		if(ggweight < m_doubMinWeight)
			m_doubMinWeight = ggweight;
		iterPPI = m_netWork.find(m_intgid);
		if( iterPPI != m_netWork.end() )
		{
			list<PPI>& temp1 = m_netWork[m_intgid];
			temp1.push_back(PPI(m_intgid1,ggweight));
		}
		else
		{
			list<PPI > pid_w11;
			pid_w11.push_back(PPI(m_intgid1,ggweight));
			m_netWork[m_intgid] = pid_w11;
		}
		iterPPI = m_netWork.find(m_intgid1);
		if( iterPPI != m_netWork.end() )
		{
			list<PPI>& temp2 = m_netWork[m_intgid1];
			temp2.push_back(PPI(m_intgid,ggweight));
		}
		else
		{
			list<PPI > pid_w12; 
			pid_w12.push_back(PPI(m_intgid,ggweight));
			m_netWork[m_intgid1] = pid_w12;
		}
	}
	m_ifs.close();
	return true;	
}

// Compute the scores of genes by MaxMIF and rank in descending order
void MaxMIF::ComputeRank()
{
	map<string,double>::iterator iterMut=m_mapGeneMut.begin();
	for(;iterMut != m_mapGeneMut.end(); ++iterMut)
	{
		m_genename=iterMut-> first;
		double MaxMIFScore= pow(m_mapGeneMut[m_genename]*m_doubMinWeight,2);
		if(MaxMIFScore <= 0)
			continue;
		
		// deal with special genes
		map<string,int>::iterator itergid=m_mapGeneId.find(m_genename);
		if(itergid == m_mapGeneId.end())
		{
			pair<double,string> p1(1e-3*MaxMIFScore,m_genename);
			m_GeneScore.insert(p1);
			continue;
		}
		
		// deal with normal genes
		m_intgid=m_mapGeneId[m_genename];
		iterPPI = m_netWork.find(m_intgid);
		if( iterPPI == m_netWork.end())
		{
			if( m_mapGeneMut[m_genename] == m_doubMinMutScore)
			{
				pair<double,string> p1(0.1*MaxMIFScore,m_genename);
				m_GeneScore.insert(p1);
				continue;
			}
			else
			{
				// when the gene is mutated 
				pair<double,string> p1(MaxMIFScore+SMAllESTNUM*m_mapGeneMut[m_genename],m_genename);
				m_GeneScore.insert(p1);
				continue;
			}
		}
		list<PPI> pid_w1=m_netWork[m_intgid];
		double  ggmaxforce=MaxMIFScore,forceTemp,rij;
		string genename1;
		list<PPI>::iterator iter4=pid_w1.begin();
		for(;iter4!=pid_w1.end(); ++iter4)
		{
			PPI gid_weight=*iter4;	
			genename1=m_mapIdGene[gid_weight.Id];
			map<string,double>::iterator iterMut1=m_mapGeneMut.find(genename1);
			if(iterMut1 !=m_mapGeneMut.end())
			{
				rij=1/gid_weight.Idw;
				forceTemp=m_mapGeneMut[m_genename]*m_mapGeneMut[genename1]/(rij*rij);
				if(ggmaxforce<forceTemp)
					ggmaxforce=forceTemp;
			}
		}
		if(ggmaxforce >0)
		{
			pair<double,string> p1(  (ggmaxforce + (SMAllESTNUM*m_mapGeneMut[m_genename]) ),m_genename);
			m_GeneScore.insert(p1);
		}	
	}
}

// Output the result of ranking by MaxMIF
void MaxMIF::OutputRank(string rankName)
	{
		m_ofs.open(rankName.c_str());
		m_ofs<<"Rank"<<"\t"<<"GeneSymbol"<<"\t"<<"EntrezID"<<"\t"<<
		"CGC"<<"\t"<<"CGCpointMut"<<"\t"<<"Rule2020"<<"\t"<<"HCD"<<"\t"<<"MouseMut"<<
		"\t"<<"Mutation-score"<<"\t"<<"MaxMIF-Score"<<"\n";
		int grank=0,countAll[5] = {0};
		multimap<double,string>::reverse_iterator iter5 = m_GeneScore.rbegin();
		for(;iter5 != m_GeneScore.rend();iter5++)
		{
			double MaxMIFScore=iter5 ->first;
			m_genename=iter5-> second;
			map<string,int>::iterator itergid=m_mapGeneId.find(m_genename);
			if(itergid !=m_mapGeneId.end())
				m_intgid=m_mapGeneId[m_genename];
			else
				m_intgid= 0;
		    grank++;
		    m_ofs<<grank<<"\t"<<m_genename<<"\t"<<m_intgid;
		    for(int k = 0; k < 5;++k)
			{
				if(m_CancerGene[k].count(m_intgid))
				{
					m_ofs<<"\t"<<"Y";
				    ++countAll[k];
				}
				else 
			    m_ofs<<"\t"<<"N";
			}
		    m_ofs<<"\t"<<m_mapGeneMut[m_genename]<<"\t"<<MaxMIFScore<<"\n";
		}
		m_ofs.close();
		m_GeneScore.clear();
	}
	
void MaxMIF::MutDataClear()
{
	m_mapGeneMut.clear();	
}

void MaxMIF::NetworkClear()
{
	m_netWork.clear();
}
