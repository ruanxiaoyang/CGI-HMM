#include <float.h>

void readHMM(const char* _path,vector<string> states,vector<string> symbols,sarray<double> & iprob,darray<double> & tprob,darray<double> & eprob)
{
	ifstream IN(_path);
	if(!IN)
	{	cout<<"can\'t find file:"<<_path<<endl;
		exit(0);
	}
	char a;
	string str,word;
	double prob;
	int iter,stn,emn,row=0,col=0;
	while(IN.get(a))
	{	if(a=='>')
		{	getline(IN,str);
			while(getline(IN,str) && str.find(">")==string::npos)
			{	stringstream ss(str);
				while(ss>>word)
				{	states.push_back(word);
				}
			}
			stn=states.size();
			tprob.fast_resize(stn,stn);

			while(getline(IN,str) && str.find(">")==string::npos)
			{	stringstream ss(str);
				while(ss>>word)
				{	symbols.push_back(word);
				}
			}
			emn=symbols.size();
			eprob.fast_resize(stn,emn);
			
			iter=0;
			iprob.resize(stn,(double)1/((double)stn));	
			while(getline(IN,str) && str.find(">")==string::npos)
			{	stringstream ss(str);
				while(ss>>prob && iter<iprob.size())
				{	iprob[iter]=prob;
					iter+=1;
				}
			}
			if(iter!=iprob.size())
			{	cout<<"States count/probability count mismatch\n";
				exit(0);
			}

			
			while(getline(IN,str) && str.find(">")==string::npos)
			{	stringstream ss(str);	
				while(ss>>prob)
				{	if(row>=tprob.getrnum() || col>=tprob.getcnum())
					{	cout<<"Transition probability matrix incorrect row/column number\n";
						exit(0);
					}
					tprob(row,col)=prob;
					col+=1;
					if(col==stn)
					{	row+=1;
						col=0;
					}
				}
			}		
	
			row=0;col=0;
			while(getline(IN,str) && str.find(">")==string::npos)
			{	stringstream ss(str);	
				while(ss>>prob)
				{	if(row>=eprob.getrnum() || col>=eprob.getcnum())
					{	cout<<"Emission probability matrix incorrect row/column number\n";
						exit(0);
					}
					eprob(row,col)=prob;
					col+=1;
					if(col==emn)
					{	row+=1;
						col=0;
					}
				}
			}	
		}		
	}	
}

sarray<double> logtrans(sarray<double> & prob)
{	
	sarray<double> res(prob.size());
	for(int i=0;i<prob.size();++i)
	{	if(prob[i]!=0)
			res[i]=-1*log(prob[i]);
		else
			res[i]=-1*log(0.0001);
	}
	return res;
}

darray<double> logtrans(darray<double> & prob)
{
	darray<double> res(prob.getrnum(),prob.getcnum());
	for(int i=0;i<prob.getrnum();++i)
	{	for(int j=0;j<prob.getcnum();++j)
		{	if(prob(i,j)!=0)
				res(i,j)=-1*log(prob(i,j));
			else
				res(i,j)=-1*log(0.0001);
		}
	}
	return res;
}

void probchain(sarray<double> & _iprob,darray<double> & _tprob,darray<double> & _eprob,sarray<short> & _seq,sarray<double> & finalprob,darray<short> & wp)
{
	sarray<double> liprob;
	liprob=logtrans(_iprob);
	darray<double> ltprob,leprob;
	ltprob=logtrans(_tprob);
	leprob=logtrans(_eprob);
	
sarray<double> prevp(liprob.size(),0);
sarray<double> currp(liprob.size(),0);
sarray<double>* pp=&prevp;
sarray<double>* pc=&currp;
sarray<double>* swap;

	wp.fast_resize(2,_seq.size());
	for(int i=0;i<liprob.size();++i)
	{	prevp[i]=liprob[i]+leprob(i,_seq[0]);
	}
	int row=liprob.size(),prevmaxI,ITER;
	double max,tmp;
	for(int j=1;j<_seq.size();++j)
	{	ITER=0;
	 	for(int i=0;i<row;++i)
		{	max=DBL_MAX;
			if(_seq[j]!=i && _seq[j]!=i-4)
			{	pc->operator[](i)=DBL_MAX;
				continue;
			}
			for(int I=0;I<row;++I)
			{	if(_seq[j-1]!=I && _seq[j-1]!=I-4)
					continue;
				tmp=pp->operator[](I)+ltprob(I,i)+leprob(i,_seq[j]);
				if(tmp<max)
				{	max=tmp;
					prevmaxI=I;
				}
			}
			pc->operator[](i)=max;
			wp(ITER,j)=(prevmaxI>3 ? 1 : 0);
			ITER+=1;
		}
		swap=pp;
		pp=pc;
		pc=swap;
	}	
	finalprob=*pp;
}

void hiddenstate(sarray<double> & _iprob,darray<double> & _tprob,darray<double> & _eprob,sarray<short> & _seq,sarray<short> & state)
{	
	sarray<double> finalprob;
	darray<short> wp;
	state.resize(_seq.size());
	sarray<int> indexarr;

	clock_t sta=clock();
	probchain(_iprob,_tprob,_eprob,_seq,finalprob,wp);
	double maxprob=finalprob.smin(indexarr);
	cout<<"(-log)maxprob:"<<maxprob<<endl;

	state[state.size()-1]=(indexarr[0]>3 ? 1 : 0);
	for(int j=state.size()-2;j>=0;--j)
		state[j]=wp(state[j+1],j+1);
	clock_t end=clock();

	ofstream LOG("log",ios_base::app);
	LOG<<"Find hidden states for sequence of length "<<_seq.size()<<"bp\n"
		<<"Time lapse "<<(end-sta)/1e6<< " sec\n";
}


void state2CGI(sarray<short> & _state,sarray<short> & _seq,darray<double> & CGI,bool _sta1=true)
{
	CGI.clear();
	int iter=0,obs=0,start=0,end=0,ccnt=0,gcnt=0,cgcnt=0,offset=1;
	if(_sta1==false)
		offset=0;

	clock_t sta=clock();	
	while(iter<_state.size())
	{	while(_state[iter]==1 && iter<_state.size())
		{	iter+=1;
			continue;
		}
		start=iter;
		ccnt=0;
		gcnt=0;
		cgcnt=0;
		while(_state[iter]==0 && iter<_state.size())
		{
			if(_seq[iter]==2)
				ccnt+=1;
			if(_seq[iter]==3)
			{	gcnt+=1;
				if(_seq[iter-1]==2)
					cgcnt+=1;
			}
			iter+=1;
			continue;
		}
		end=iter-1;
		if(end>start)
		{
			CGI.fill(obs,0,obs+1);
			CGI.fill(obs,1,start+offset);
			CGI.fill(obs,2,end+offset);
			CGI.fill(obs,3,end-start+1);
			CGI.fill(obs,4,(ccnt+gcnt)/((double)(end-start+1)));
			CGI.fill(obs,5,(cgcnt)/((double)(ccnt*gcnt))*(end-start+1));
			CGI.fill(obs,6,cgcnt);
			CGI.fill(obs,7,(end-start+1)/(double)cgcnt);
			obs+=1;
		}
	}
	CGI.resize(obs,8);
	clock_t tend=clock();

	ofstream LOG("log",ios_base::app);
        LOG<<"Locate CGI\n"
                <<"Time lapse "<<(tend-sta)/1e6<< " sec\n\n";
}
				
void reportCGI(string & _path,darray<double> & _CGI,const int & _seqlen,const double & _GCcont,const double & _OEratio)
{
	int iter=1;
	ofstream OUT(_path.c_str());
	OUT<<"ID\tSta\tEnd\tLen\tGCprop\tOEratio\tCpG_Num\tMeanDist\n";
	for(int i=0;i<_CGI.getrnum();++i)
	{	if(_CGI(i,3)>=_seqlen && _CGI(i,4)>=_GCcont && _CGI(i,5)>=_OEratio)
		{	OUT<<iter<<"\t"
			<<(int)_CGI(i,1)<<"\t"
			<<(int)_CGI(i,2)<<"\t"
			<<(int)_CGI(i,3)<<"\t";
			OUT.precision(3);
			OUT<<_CGI(i,4)<<"\t"
			<<_CGI(i,5)<<"\t"
			<<(int)_CGI(i,6)<<"\t"
			<<_CGI(i,7)<<"\n";
			iter+=1;
		}
	}
}

//read in fasta seq one by one and process immediately. Discard unused seq to save space  
void readseq(const char* _path,const char* _outpath,sarray<double> & _iprob,darray<double> & _tprob, darray<double> & _eprob,const int & _seqlen,const double & _GCcont,const double & _OEratio,bool _sta1) 
{ 
        ifstream IN(_path); 
        if(! IN) 
        {       cout<<_path<<" not found\n"; 
                exit(0); 
        } 
        char a; 
        int iter,seqcount=0; 
        string str; 
	sarray<short> seq;
	sarray<short> state;
	darray<double> CGI;
	
	clock_t sta=clock();
        while(IN.get(a)) 
        {       if(a=='>') 
                {       getline(IN,str);
			stringstream ss(str);
			ss>>str;
			cout<<str<<endl;
			iter=0;
			seqcount+=1;	
                        while(IN.get(a) && a!='>') 
                        {       if(a=='\n' || a==' ') 
                                        continue; 
                                seq.fill(iter,na2id(a),10000); 
                                iter+=1; 
                        }
			seq.fill_trim(); 
			cout<<"Find hidden state for "<<str<<"\n"
                        <<"Length:"<<seq.size()<<"\n";
                	hiddenstate(_iprob,_tprob,_eprob,seq,state);

                	cout<<"Generate output for "<<str<<"\n";
                	state2CGI(state,seq,CGI,_sta1);

                	string outputname=_outpath+string(".")+str+".cgi";
                	reportCGI(outputname,CGI,_seqlen,_GCcont,_OEratio);
                	seq.clear();
                } 
		IN.unget();
        }
	clock_t end=clock();

	ofstream LOG("log",ios_base::app);
	LOG<<getdate()<<endl;
	LOG<<"Read in and processed "<<seqcount<<" sequence(s)\n"
		<<"Time lapse "<<(end-sta)/1e6<< " sec\n"
		<<"Filter threshold\n"
		<<"CGI length:"<<_seqlen<<"\n"
		<<"GC content:"<<_GCcont<<"\n"
		<<"OE ratio:  "<<_OEratio<<"\n\n";
}


