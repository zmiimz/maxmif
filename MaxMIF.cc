/* Author:  Yingnan Hou<houyingnan3@gmail.com> 2017
 * Usage: This is MaxMIF package. Use, redistribution, modify without limitations
 */

/***********************************************************************/

#include "MaxMIF.h"

/***********************************************************************/

int main (int argc, char** argv)
{
	MaxMIF maxmif1;
	
	string Makers="MaxMIF_",outPath="Output//",
	defnetworkPath="Background//Network//",newNetName="HumanNet",networkPath, MutDataPath, MutName;
	if(argc == 2)
	{
		MutDataPath= argv[1],MutName= "new";
	}
	else if(argc == 3  )
	{
		MutDataPath= argv[1],MutName= argv[2];
	}
	else if(argc == 4  )
	{
		MutDataPath= argv[1],MutName= argv[2], networkPath=argv[3];
	}
	else
	{
		MutDataPath= argv[1],MutName= argv[2], networkPath=argv[3], newNetName=argv[4];
	}
	
	/* Read in the reference gene sets */
	maxmif1.ReadCancerGene();
	maxmif1.ReadGeneID();
	
	/* Read in the mutation data */
	string outputName,outputPath;
	cout<<"\n"<<"Read Mutation data: "<<MutDataPath<<"   ...";
	bool readTF=maxmif1.ReadMutData(MutDataPath);
	if(readTF== false) return 0;
	
	/* choose network  */ 
	int nettype=2; // use the default network: HumanNet and STRINGv10
	maxmif1.ComputeMut();
	for(int nnet =1; nnet <= nettype; nnet ++)
	{
		if(nnet == 2) newNetName= "STRINGv10";
		if(argc >= 4) nettype=1;
		else
		{
			networkPath = defnetworkPath + newNetName;
		}
		cout<<"\n"<<"Read Network data:  "<< newNetName <<"   ...";
		
		/* Read in the Network data */
		bool readNetTF=maxmif1.ReadNetwork(networkPath);
		if(readNetTF== false) return 0;
		
		/* Compute MaxMIF Score */
		maxmif1.ComputeRank();
		outputName=Makers+MutName+"_"+newNetName+".txt";
		outputPath=outPath+outputName;
		maxmif1.OutputRank(outputPath);
		cout<<"\n\n"<<"Output file:"<<"\n"<<"\t\t"<< outputName<<"  ..."<<"\n";
		maxmif1.NetworkClear();
	}
	maxmif1.MutDataClear();
	return 0;	
}