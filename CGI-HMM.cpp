#include <fstream>
#include <math.h>
#include <vector>
#include <iostream>
#include <sstream>	//for stringstream class
#include <sarray.h>
#include <darray.h>
#include <code.h>
#include <getdate.h>
#include <CGI-HMM.h>
using namespace std;

int main(int argc,char* argv[])
{
	char o;
	const char* parampath="prob.hmm";
	char* seqpath=NULL;
	const char* outpath="CGIres";
	int seqlen=0;
	double GCcont=0.5;
	double OEratio=0.6;
	
	bool sta1=true;
	while((o=getopt(argc,argv,"p:s:l:c:r:o:Sh"))!=-1)
	{	switch (o)
		{	case 'p':parampath=optarg;
				break;
			case 's':seqpath=optarg;
				break;
			case 'l':seqlen=stn<int>(optarg);
				break;
			case 'c':GCcont=stn<double>(optarg);
				break;
			case 'r':OEratio=stn<double>(optarg);
				break;
			case 'o':outpath=optarg;
				break;
			case 'S':sta1=false;
				break;
			case 'h':system("less Instruction_Files/CGI-HMM");
				exit(0);
			case '?':exit(0);
		}
	}

	if(!seqpath)
	{	cout<<"Need to specify a sequence file with -s option\n";
		exit(0);
	}

	vector<string> states,symbols;
	sarray<double> iprob;
	darray<double> tprob,eprob;
	cout<<"Loading HMM transition probability\n";
	readHMM(parampath,states,symbols,iprob,tprob,eprob);

	cout<<"Loading sequence file\n";
	readseq(seqpath,outpath,iprob,tprob,eprob,seqlen,GCcont,OEratio,sta1);

}
