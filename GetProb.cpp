#include <fstream>
#include <math.h>
#include <vector>
#include <iostream>
#include <sstream>	//for stringstream class
#include <sarray.h>
#include <darray.h>
#include <GetProb.h>
using namespace std;

int main(int argc,char* argv[])
{
	char o;
	char* seqpath=NULL;
	char* pospath=NULL;
	const char* outpath="res.hmm";
	while((o=getopt(argc,argv,"s:p:o:h"))!=-1)
	{	switch (o)
		{	case 's':seqpath=optarg;
				break;
			case 'p':pospath=optarg;
				break;
			case 'o':outpath=optarg;
				break;
			case 'h':system("less Instruction_Files/GetProb");
				exit(0);
			case '?':exit(0);
		}
	}
	if(!seqpath || !pospath)
	{	cout<<"Need sequence and position file. Use -h for help\n";
		exit(0);
	}
	darray<int> pos;
	cout<<"Loading position file...\n";
	readpos(pospath,pos);		
	cout<<"Loading sequence file...\n";
	sarray<int> seq;
	readseq(seqpath,seq);

	darray<double> tprob;
	sarray<double> iprob;
	cout<<"Compute transition probability...\n";
	transprob(seq,pos,tprob,iprob);
	cout<<"writing output to "<<outpath<<endl;
	output(iprob,tprob,outpath);
}
