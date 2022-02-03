#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <cmath>
#include <cstring>
//#include <filesystem> // zmi TODO replace all file path operations to this
#include <experimental/filesystem> // zmi TODO replace all file path operations to this
#include "ProgressBar.hpp" // for custom progress
#include "json.hpp" // for json config parser
#include <cassert>

using namespace std;
namespace fs = std::experimental::filesystem;
//namespace fs = std::filesystem;
void CheckPathExist(const string &sPath);

/* Create Class MaxMIF */
class MaxMIF {
public:
    static inline bool bVerbose{true};
    static inline bool bUseCancerReferenceData{true};

    typedef struct PPI {
        int Id;
        double Idw;

        PPI(int _id, double _idw) {
            this->Id = _id;
            this->Idw = _idw;
        }
    } PPI;

public:
    MaxMIF(); // constructor
    bool ReadMaxMIFConfig();

    void ReadCancerGene();

    void ReadGeneID();

    bool ReadMutData(const string &ifspathstr, const bool &bIsMutationFrequency);

    void ComputeMut();

    bool ReadNetwork(const string &sInputNetworkPath, const string &sDefaultNetworkID);

    void ComputeRank();

    void OutputRank(const string &sOutputFile);

    void MutDataClear();

    void NetworkClear();

private:
    static void StringSplit(string longstr, const string &str1, vector<string> &vectstr);

    static inline const double EPSILON{1e-15};
    static inline const string SEPARATOR{"\t"};
    static inline const string CONFIGFILENAME{"MaxMIF.json"};
    static inline const unsigned long m_barSize{80};

private:
    int m_intgid, m_intgid1, m_intSampNum;
    double m_doubMinWeight, m_doubMinMutScore;
    string m_strtemp, m_genename;
    vector<string> m_vecTemp, m_vectUNMutGenes;
    map<int, string> m_mapIdGene;
    map<string, int> m_mapGeneId;
    map<string, vector<int> > m_geneMutData;
    map<int, list<PPI> > m_netWork;
    map<int, list<PPI> >::iterator iterPPI;
    map<string, double> m_mapGeneMut;
    multimap<double, string> m_GeneScore;
    set<int> m_CancerGene[6];
    ifstream m_ifs;
    ofstream m_ofs;

    string m_CGC_CancerGeneSet_path;
    string m_CGCpointMut_CancerGeneSet_path;
    string m_Rule2020_CancerGeneSet_path;
    string m_HCD_CancerGeneSet_path;
    string m_MouseMut_CancerGeneSet_path;
    string m_CGC2019_CancerGeneSet_path;
    string m_geneID_path;
    string m_HumanNet_path;
    string m_STRING_path;

};



/* the functions of MaxMIF */

// TODO make it optional for real problem
// Read reference cancer gene sets
void MaxMIF::ReadCancerGene() {
    cout << "Reading reference cancer gene sets..." << endl;
    string path[6] = {
            m_CGC_CancerGeneSet_path,
            m_CGCpointMut_CancerGeneSet_path,
            m_Rule2020_CancerGeneSet_path,
            m_HCD_CancerGeneSet_path,
            m_MouseMut_CancerGeneSet_path,
            m_CGC2019_CancerGeneSet_path};
    for (int k = 0; k < 6; ++k) {
        CheckPathExist(path[k]);
        m_ifs.open(path[k].c_str(), ifstream::in);
        if (!m_ifs.is_open()) {
            cout << endl << "Failed to open the CancerGeneSet data: " << path[k].c_str() << "!" << endl
                 << "Please check if your input is correct and try again later!" << endl;
            exit(-1);
        }
        while (getline(m_ifs, m_strtemp)) {
            //cout << "m_strtemp: "<< m_strtemp << endl;
            m_intgid = std::stoi(m_strtemp); // atoi(m_strtemp.c_str());
            //cout <<"k: "<<k<<" current: "<< m_strtemp << " m_intgid: "<< m_intgid << endl;
            assert(m_intgid != 0);
            m_CancerGene[k].insert(m_intgid);
        }
        m_ifs.close();
    }
}

// Read gene Entrez ID file
void MaxMIF::ReadGeneID() {
    cout << "Reading gene Entrez ID set..." << endl;
    string idpath = m_geneID_path;
    CheckPathExist(idpath);
    m_ifs.open(idpath.c_str(), ifstream::in);
    if (!m_ifs.is_open()) {
        cout << endl << "Failed to open the geneID data: " << idpath.c_str() << "!" << endl
             << "Please check if your input is correct and try again later!" << endl;
        exit(-1);
    }
    while (getline(m_ifs, m_strtemp)) {
        StringSplit(m_strtemp, SEPARATOR, m_vecTemp);
        m_intgid = stoi(m_vecTemp[1]); //atoi(m_vecTemp[1].c_str());
        m_mapGeneId[m_vecTemp[0]] = m_intgid;
        m_mapIdGene[m_intgid] = m_vecTemp[0];
    }
    m_ifs.close();
}

void MaxMIF::StringSplit(string longstr, const string &sep, vector<string> &vectstr) {
    vectstr.clear();
    unsigned long nend{0L};
    unsigned long nbegin{0L};
    while (nend != string::npos) {
        nend = longstr.find(sep, nbegin);
        if (nend == string::npos)
            vectstr.push_back(longstr.substr(nbegin, longstr.length() - nbegin));
        else
            vectstr.push_back(longstr.substr(nbegin, nend - nbegin));
        nbegin = nend + sep.length();
    }
}

// Read the gene-level somatic mutation data
bool MaxMIF::ReadMutData(const string &ifspathstr, const bool &bIsMutationFrequency) {
    CheckPathExist(ifspathstr);
    m_ifs.open(ifspathstr.c_str(), ifstream::in);
    if (!m_ifs.is_open()) {
        cout << endl << "Failed to open the Mutation data: " << ifspathstr << "!" << endl
             << "Please check if your input is correct and try again later!" << endl;
        return false;
    }

    // get number of lines in file
    unsigned long numLines{};
    cout << endl;
    while (getline(m_ifs, m_strtemp)) {
        ++numLines;
        //cout<<"\rLines: "<<numLines;
    }
    // reset stream
    m_ifs.clear();
    m_ifs.seekg(0, std::ios::beg);

    getline(m_ifs, m_strtemp);
    StringSplit(m_strtemp, SEPARATOR, m_vecTemp);
    m_intSampNum = m_vecTemp.size() - 1;
    const unsigned long limit = numLines;

    if (true == bIsMutationFrequency) {
        if (!m_intSampNum == 1) {
            cout << endl << "Failed to identify the Mutation frequency data: " << ifspathstr << "!" << endl
                 << "Please check if your input is correct (2 columns only) and try again later!" << endl;
            return false;
        } else {
            // initialize the bar
            ProgressBar progressBar(limit, m_barSize);
            // if the input data is  mutation frequency
            cout << endl << "Reading Mutation frequency data ..." << endl;
            while (getline(m_ifs, m_strtemp)) {
                ++progressBar;
                StringSplit(m_strtemp, SEPARATOR, m_vecTemp);
                //m_genename = m_vecTemp[0].c_str();
                int geneID_tmp = stoi(m_vecTemp[0]);
                m_genename = m_mapIdGene[geneID_tmp]; // tempo for entrezID input
                if(m_genename.empty()){ m_genename =  m_vecTemp[0];}
                //cout << "m_genename "<< m_genename << endl; // tempo
                double gmut = stod(m_vecTemp[1]);
                if (gmut > std::numeric_limits<double>::epsilon()) {
                    m_mapGeneMut[m_genename] = gmut;
                }
                // display the bar
                if (MaxMIF::bVerbose) {
                    progressBar.display();
                }
            }
            m_ifs.close();
            // tell the bar to finish
            progressBar.done();
            return true;
        }
    }

    if (false == bIsMutationFrequency) {
        if (m_intSampNum == 1) {
            cout << endl << "Failed to identify the Mutation matrix data: " << ifspathstr << "!" << endl
                 << "Please check if your input is correct and try again later!" << endl;
            return false;
        } else {
            // if the input data is  mutation matrix
            cout << endl << "Reading Mutation matrix data for " << m_intSampNum << " samples..." << endl;
            vector<int> vectinttemp;
            string s1 = "1", s2 = "1.0";
            // initialize the bar
            ProgressBar progressBar(limit, m_barSize);
            while (getline(m_ifs, m_strtemp)) {
                ++progressBar;
                StringSplit(m_strtemp, SEPARATOR, m_vecTemp);
                //m_genename = m_vecTemp[0].c_str();
                int geneID_tmp = stoi(m_vecTemp[0]);
                m_genename = m_mapIdGene[geneID_tmp]; // tempo for entrezID input
                if(m_genename.empty()){ m_genename =  m_vecTemp[0];}
                //cout << "m_genename "<< m_genename << endl; // tempo

                for (int j = 0; j < m_intSampNum; ++j) {
                    if (m_vecTemp[j + 1].c_str() == s1 || m_vecTemp[j + 1].c_str() == s2) {
                        vectinttemp.push_back(j);
                    }
                }
                // symbol genes
                if (!vectinttemp.empty()) {
                    m_geneMutData[m_genename] = vectinttemp;
                    vectinttemp.clear();
                } else {
                    // deal with special genes
                    map<string, int>::iterator itergid = m_mapGeneId.find(m_genename);
                    if (itergid != m_mapGeneId.end())
                        m_vectUNMutGenes.push_back(m_genename);
                }
                // display the bar
                if (MaxMIF::bVerbose) {
                    progressBar.display();
                }
            }
            m_ifs.close();
            // tell the bar to finish
            progressBar.done();
            return true;
        }
    }
    return false;
}

// Compute the mutation-score of every mutated genes
void MaxMIF::ComputeMut() {
    vector<double> vectSampMutNum;
    map<string, vector<int> >::iterator iter1;
    int maxMutNum = 1;
    // initialize the bar
    const unsigned long limit1{static_cast<unsigned long>(m_intSampNum)};
    ProgressBar progressBar1(limit1, m_barSize);
    cout << "Computing mutations for " << limit1 << " samples and " << m_geneMutData.size() << " genes... " << endl;
    for (int j = 0; j < m_intSampNum; ++j) {
        ++progressBar1;
        int nMutNum = 0;
        for (iter1 = m_geneMutData.begin(); iter1 != m_geneMutData.end(); iter1++) {
            m_genename = iter1->first;
            vector<int> vecSampleID = iter1->second;
            vector<int>::iterator iter2;
            iter2 = find(vecSampleID.begin(), vecSampleID.end(), j);
            if (iter2 != vecSampleID.end()) {
                nMutNum++;
            }
        }
        if (nMutNum != 0) {
            vectSampMutNum.push_back((double) 1.0 / (double) nMutNum);
            if (maxMutNum < nMutNum)
                maxMutNum = nMutNum;
        } else {
            vectSampMutNum.push_back(0);
        }

        if (MaxMIF::bVerbose) {
            progressBar1.display();
        }
    }
    progressBar1.done();

    // initialize the bar
    const unsigned long limit2{m_geneMutData.size()};
    ProgressBar progressBar2(limit2, m_barSize);

    for (iter1 = m_geneMutData.begin(); iter1 != m_geneMutData.end(); iter1++) {
        ++progressBar2;
        double gidMut = 0;
        m_genename = iter1->first;
        vector<int> vecSampleID = iter1->second;
        vector<int>::iterator iter3 = vecSampleID.begin();
        for (; iter3 != vecSampleID.end(); iter3++) {
            int samplabel = *iter3;
            gidMut = gidMut + vectSampMutNum[samplabel];
        }
        m_mapGeneMut[m_genename] = gidMut;
        if (MaxMIF::bVerbose) {
            progressBar2.display();
        }
    }
    progressBar2.done();

    vectSampMutNum.clear();
    m_geneMutData.clear();
    // genes without mutation
    m_doubMinMutScore = (double) 1.0 / ((double) maxMutNum);

    // initialize the bar
    const unsigned long limit3{m_vectUNMutGenes.size()};
    ProgressBar progressBar3(limit3, m_barSize);

    ////  Deal with genes who are not mutated
    for (int k = 0; k < m_vectUNMutGenes.size(); k++) {
        ++progressBar3;
        m_genename = m_vectUNMutGenes[k];
        m_mapGeneMut[m_genename] = m_doubMinMutScore;
        if (MaxMIF::bVerbose) {
            progressBar3.display();
        }
    }
    progressBar3.done();
    m_vectUNMutGenes.clear();
}

// Read functional network
bool MaxMIF::ReadNetwork(const string &sInputNetworkPath, const string &sDefaultNetworkID) {
    m_doubMinWeight = 0.01;  // set a small weight;
    double ggweight;
    string sNetworkPath;
    if (sDefaultNetworkID == "STRING10") {
        sNetworkPath = m_STRING_path;
    } else if (sDefaultNetworkID == "HumanNet1") {
        sNetworkPath = m_HumanNet_path;
    } else {
        sNetworkPath = sInputNetworkPath;
    }

    CheckPathExist(sNetworkPath);
    m_ifs.open(sNetworkPath.c_str(), ifstream::in);
    if (!m_ifs.is_open()) {
        cout << endl << "Failed to open the Network data: " << sNetworkPath << "!" << endl
             << "Please check if your input is correct and try again later!" << endl;
        return false;
    }

    // get number of lines in file
    unsigned long numLines{};
    cout << endl;
    while (getline(m_ifs, m_strtemp)) {
        ++numLines;
        //cout<<"\rLines: "<<numLines;
    }
    // reset stream
    m_ifs.clear();
    m_ifs.seekg(0, std::ios::beg);

    const unsigned long limit{numLines};
    ProgressBar progressBar(limit, m_barSize);
    cout << "Importing " << numLines << " lines ..." << endl;
    while (getline(m_ifs, m_strtemp)) {
        ++progressBar;
        StringSplit(m_strtemp, SEPARATOR, m_vecTemp);
        m_intgid = stoi(m_vecTemp[0]);// atoi(m_vecTemp[0].c_str());
        m_intgid1 = stoi(m_vecTemp[1]); // atoi(m_vecTemp[1].c_str());
        map<int, string>::iterator iteridg, itergid1g;
        iteridg = m_mapIdGene.find(m_intgid);
        itergid1g = m_mapIdGene.find(m_intgid1);
        if ((iteridg == m_mapIdGene.end()) | (itergid1g == m_mapIdGene.end()))
            continue;
        ggweight =stod(m_vecTemp[2]); // atof(m_vecTemp[2].c_str());
        if (ggweight < m_doubMinWeight)
            m_doubMinWeight = ggweight;
        iterPPI = m_netWork.find(m_intgid);
        if (iterPPI != m_netWork.end()) {
            list<PPI> &temp1 = m_netWork[m_intgid];
            temp1.push_back(PPI(m_intgid1, ggweight));
        } else {
            list<PPI> pid_w11;
            pid_w11.push_back(PPI(m_intgid1, ggweight));
            m_netWork[m_intgid] = pid_w11;
        }
        iterPPI = m_netWork.find(m_intgid1);
        if (iterPPI != m_netWork.end()) {
            list<PPI> &temp2 = m_netWork[m_intgid1];
            temp2.push_back(PPI(m_intgid, ggweight));
        } else {
            list<PPI> pid_w12;
            pid_w12.push_back(PPI(m_intgid, ggweight));
            m_netWork[m_intgid1] = pid_w12;
        }
        if (MaxMIF::bVerbose) {
            progressBar.display();
        }
    }
    m_ifs.close();
    progressBar.done();
    return true;
}

// Compute the scores of genes by MaxMIF and rank in descending order
void MaxMIF::ComputeRank() {
    const unsigned long limit{m_mapGeneMut.size()};
    ProgressBar progressBar(limit, m_barSize);
    cout << "Calculating " << limit << " ranks ..." << endl;

    map<string, double>::iterator iterMut = m_mapGeneMut.begin();
    for (; iterMut != m_mapGeneMut.end(); ++iterMut) {
        ++progressBar;
        m_genename = iterMut->first;
        double MaxMIFScore = pow(m_mapGeneMut[m_genename] * m_doubMinWeight, 2);
        if (MaxMIFScore <= 0)
            continue;

        // deal with special genes
        map<string, int>::iterator itergid = m_mapGeneId.find(m_genename);
        if (itergid == m_mapGeneId.end()) {
            pair<double, string> p1(1e-3 * MaxMIFScore, m_genename);
            m_GeneScore.insert(p1);
            continue;
        }

        // deal with normal genes
        m_intgid = m_mapGeneId[m_genename];
        iterPPI = m_netWork.find(m_intgid);
        if (iterPPI == m_netWork.end()) {
            if (m_mapGeneMut[m_genename] == m_doubMinMutScore) {
                pair<double, string> p1(0.1 * MaxMIFScore, m_genename);
                m_GeneScore.insert(p1);
                continue;
            } else {
                // when the gene is mutated
                pair<double, string> p1(MaxMIFScore + EPSILON * m_mapGeneMut[m_genename], m_genename);
                m_GeneScore.insert(p1);
                continue;
            }
        }
        list<PPI> pid_w1 = m_netWork[m_intgid];
        double ggmaxforce = MaxMIFScore, forceTemp, rij;
        string genename1;
        list<PPI>::iterator iter4 = pid_w1.begin();
        for (; iter4 != pid_w1.end(); ++iter4) {
            PPI gid_weight = *iter4;
            genename1 = m_mapIdGene[gid_weight.Id];
            map<string, double>::iterator iterMut1 = m_mapGeneMut.find(genename1);
            if (iterMut1 != m_mapGeneMut.end()) {
                rij = 1 / gid_weight.Idw;
                forceTemp = m_mapGeneMut[m_genename] * m_mapGeneMut[genename1] / (rij * rij);
                if (ggmaxforce < forceTemp)
                    ggmaxforce = forceTemp;
            }
        }
        if (ggmaxforce > 0) {
            pair<double, string> p1((ggmaxforce + (EPSILON * m_mapGeneMut[m_genename])), m_genename);
            m_GeneScore.insert(p1);
        }
        if (MaxMIF::bVerbose) {
            progressBar.display();
        }
    }
    progressBar.done();
}

// Output the result of ranking by MaxMIF
void MaxMIF::OutputRank(const string &sOutputFile) {
    // check that output not exists
    fs::path path = sOutputFile;
    if (fs::exists(path)) {
        cout << "ERROR: Output file " << sOutputFile << " already exist" << endl;
        exit(-1);
    };

    m_ofs.open(sOutputFile.c_str(), ifstream::out);
    if (!m_ofs.is_open()) {
        cout << endl << "Failed to open the Output path: " << sOutputFile << "!" << endl
             << "Please check if your input is correct and try again later!" << endl;
        exit(-1);
    }

    if (MaxMIF::bUseCancerReferenceData) {
        m_ofs << "Rank" << "\t" << "GeneSymbol" << "\t" << "EntrezID" << "\t" <<
              "CGC" << "\t" << "CGCpointMut" << "\t" << "Rule2020" << "\t" << "HCD" << "\t" << "MouseMut" << "\t" << "CGC2019" <<
              "\t" << "Mutation-score" << "\t" << "MaxMIF-Score" << endl;
    } else {
        m_ofs << "Rank" << "\t" << "GeneSymbol" << "\t" << "EntrezID" << "\t" <<
              "Mutation-score" << "\t" << "MaxMIF-Score" << endl;
    }

    int grank = 0;
    int countAll[6] = {0};
    multimap<double, string>::reverse_iterator iter5 = m_GeneScore.rbegin();
    for (; iter5 != m_GeneScore.rend(); iter5++) {
        double MaxMIFScore = iter5->first;
        m_genename = iter5->second;
        map<string, int>::iterator itergid = m_mapGeneId.find(m_genename);
        string sGeneID;
        if (itergid != m_mapGeneId.end()){
            m_intgid = m_mapGeneId[m_genename];
            sGeneID = to_string(m_intgid);
        } else{
            m_intgid = 0;
            sGeneID = m_genename; // if not found, copy gene name (it is replaced by id already)
        }
        grank++;

        m_ofs << grank << "\t" << m_genename << "\t" << sGeneID;
        if (MaxMIF::bUseCancerReferenceData) {
            for (int k = 0; k < 6; ++k) {
                if (m_CancerGene[k].count(m_intgid)) {
                    m_ofs << "\t" << "Y";
                    ++countAll[k];
                } else
                    m_ofs << "\t" << "N";
            }
        }
        m_ofs << "\t" << std::scientific << std::setprecision(16) << m_mapGeneMut[m_genename] << "\t" << MaxMIFScore << endl;
    }
    m_ofs.close();
    m_GeneScore.clear();
}

void MaxMIF::MutDataClear() {
    m_mapGeneMut.clear();
}

void MaxMIF::NetworkClear() {
    m_netWork.clear();
}

bool MaxMIF::ReadMaxMIFConfig() {
    cout << string(120, '-') << endl;
    fs::path sCurrentExePath = fs::current_path();
    fs::path sConfigPath = sCurrentExePath / CONFIGFILENAME;
    CheckPathExist(sConfigPath.string());

    ifstream ifs(sConfigPath.string());
    nlohmann::json m_jConfig; // remove m_jConfig from class,  add just parsed and checked data fields with paths instead
    ifs >> m_jConfig; // read json
    cout << "Successfully read config file " << sConfigPath.string() << endl;
    // DONE read all paths into public data fields  and validate if files exist.
    if (m_jConfig["parameters"]["general"]["verbose"] == "false") {
        bVerbose = false;
    }
    if (m_jConfig["parameters"]["general"]["useCancerSet"] == "false") {
        bUseCancerReferenceData = false;
    }
    m_CGC_CancerGeneSet_path = m_jConfig["parameters"]["CancerGeneSets"]["CGC"]["path"];
    CheckPathExist(m_CGC_CancerGeneSet_path);
    m_CGCpointMut_CancerGeneSet_path = m_jConfig["parameters"]["CancerGeneSets"]["CGCpointMut"]["path"];
    CheckPathExist(m_CGCpointMut_CancerGeneSet_path);
    m_Rule2020_CancerGeneSet_path = m_jConfig["parameters"]["CancerGeneSets"]["Rule2020"]["path"];
    CheckPathExist(m_Rule2020_CancerGeneSet_path);
    m_HCD_CancerGeneSet_path = m_jConfig["parameters"]["CancerGeneSets"]["HCD"]["path"];
    CheckPathExist(m_HCD_CancerGeneSet_path);
    m_MouseMut_CancerGeneSet_path = m_jConfig["parameters"]["CancerGeneSets"]["MouseMut"]["path"];
    CheckPathExist(m_MouseMut_CancerGeneSet_path);
    m_CGC2019_CancerGeneSet_path = m_jConfig["parameters"]["CancerGeneSets"]["CGC2019"]["path"];
    CheckPathExist(m_CGC2019_CancerGeneSet_path);
    m_geneID_path = m_jConfig["parameters"]["geneID"]["path"];
    CheckPathExist(m_geneID_path);
    m_HumanNet_path = m_jConfig["parameters"]["networks"]["HUMANNET"]["path"];
    CheckPathExist(m_HumanNet_path);
    m_STRING_path = m_jConfig["parameters"]["networks"]["STRING"]["path"];
    CheckPathExist(m_STRING_path);
    cout << "Successfully checked for existence of all input files from config" << endl;
    cout << string(120, '-') << endl;
    return true;
}

// constructor
MaxMIF::MaxMIF() {
    if (!ReadMaxMIFConfig()) {
        exit(1);
    }

}

void CheckPathExist(const string &sPath) {
    fs::path path = sPath;
    if (!fs::exists(path)) {
        cout << "ERROR: Expected file " << sPath << " does not exist" << endl;
        exit(1);
    };

}
