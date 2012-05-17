#include <code.h>
void readpos(const char* _path,darray<int> & pos)
{
	ifstream IN(_path);
	if(! IN)
	{	cout<<_path<<" not found\n";
		exit(0);
	}
	int tmp;
	int iter=0,row=1;
	while(IN>>tmp)
	{	row=(row==1 ? 0 : 1);
		pos.fill(row,iter,tmp,100);
		if(row==1)
			iter+=1;
	}
	pos.resize(2,iter);
}



void readseq(const char* _path,sarray<int> & seq)
{
	ifstream IN(_path);
	if(! IN)
	{	cout<<_path<<" not found\n";
		exit(0);
	}
	char a;
	int iter=0;
	string str;
	while(IN.get(a))
	{	if(a=='>')
		{	getline(IN,str);
			while(IN.get(a))
			{	if(a=='\n' || a==' ')
					continue;
				seq.fill(iter,na2id(a),10000);
				iter+=1;
			}
		}
	}
	seq.fill_trim();
}

void transprob(sarray<int> & _seq,const int & _sta,const int & _end,const int & _base,darray<int> & tcount)
{	
	for(int j=_sta+1;j<=_end;++j)
	{	if(_seq[j-1]<0 || _seq[j]<0)
			continue;
		tcount(_seq[j-1]+_base,_seq[j]+_base)+=1;
	}
}

void count2prob(darray<int> & _tcount,darray<double> & tprob,sarray<double> & iprob)
{
	tprob.fast_resize(_tcount.getrnum(),_tcount.getcnum());	
	sarray<int> rowsum;
	rowsum=_tcount.rowsum();
	for(int i=0;i<_tcount.getrnum();++i)
	{	for(int j=0;j<_tcount.getcnum();++j)
		{	tprob(i,j)=(double)_tcount(i,j)/rowsum[i];
		}
	}
	iprob.resize(rowsum.size());
	int total=rowsum.sum();
	for(int i=0;i<rowsum.size();++i)
		iprob[i]=rowsum[i]/(double)total;
	return;
}

void transprob(sarray<int> & _seq,darray<int> & _pos,darray<double> & tprob,sarray<double> & iprob)
{
	darray<int> tcount(8,8,0);
	int nstart,nend;
	for(int p=0;p<_pos.getcnum();++p)
	{	if(p==0)
			nstart=0;
		else
			nstart=_pos(1,p-1)+1;
		nend=_pos(0,p)-1;
		transprob(_seq,nstart,nend,4,tcount);
		transprob(_seq,_pos(0,p),_pos(1,p),0,tcount);
	}
	transprob(_seq,_pos(1,_pos.getcnum()-1)+1,_seq.size()-1,4,tcount);
	for(int i=0;i<4;++i)
	{	for(int j=4;j<8;++j)
		{	tcount(i,j)=tcount(j,i)=(int)(_pos.getcnum()/16.0);
		}
	}
	count2prob(tcount,tprob,iprob);
	return;
}	

void output(sarray<double> & _iprob,darray<double> & _tprob,const char * outpath)
{
	darray<int> eprob(8,4,0);
	for(int i=0;i<4;++i)
		eprob(i,i)=eprob(i+4,i)=1;	
	ofstream OUT(outpath);
	OUT<<">states\n"
		<<"A+\tT+\tC+\tG+\tA-\tT-\tC-\tG-\n";

	OUT<<">symbol\n"
		<<"A\tT\tC\tG\n";

	OUT<<">initial\n";
	_iprob.record(OUT);
	OUT<<"\n";

	OUT<<">transition\n";
	_tprob.record(OUT);

	OUT<<">emission\n";
	eprob.record(OUT);	
	return;
}
	

