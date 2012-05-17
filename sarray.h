#include <iostream>
#include <stdlib.h>
using namespace std;

template<class type>		//forward declaration for friend function
class sarray;

template<class type>
ostream & operator << (ostream & os, const sarray<type> & temp);


template<class type>
class sarray
{
public:		int _size;
		type * sparray;
		type _defau;
		void addheap();
		int _sizelim;
	sarray();
	sarray(const int & size);
	sarray(const int & size,const type & defau);
	sarray(const sarray & temp);
	sarray(const std::string & temp);
	void init();
	void init(const int & size);
	void init(const int & size, const type & defau);
	type * getarr();
	void operator=(const sarray<type> & temp);
	sarray<type> operator+(const sarray<type> & temp);
	void operator=(const string & temp);
	bool operator==(const string & temp);
	bool operator==(const sarray & temp);
	type & operator[](const int & index){return sparray[index];}
 	friend ostream & operator << <> (ostream & os,const sarray<type> & temp);
	int size(){return _size;}
	void resize(const int & size);
	void resize(const int &size,const type & defau);
	type smax();
	type smax(sarray<int> & indexarr);
	type smin();
	type smin(sarray<int> & indexarr);
	void pushback(const type & temp);
	void merge(const type & temp);
	bool find(const type & findval);
	bool empty();
	sarray<type> unique();		//must sort first
	sarray<type> subseq(const int & start,const int & end);
	sarray<type> c(sarray<type> & b);
	string tos();
	void freem(){free(sparray);sparray=NULL;}
	void clear(){sparray=(type*)realloc(sparray,sizeof(type)*0);_size=0;}	
	~sarray();
	void fill(const int & index,const type & val,const int & buffer);
	void fill_trim();
	type sum();
	void record(ostream & os);	
};

template<class type>
void sarray<type>::record(ostream &os)
{
        if(_size==0)
        {	os<<"Empty array!\n";
                return;
        }
	for(int i=0;i<_size;++i)
       		os<<sparray[i]<<"\t";
}

template<class type>
type sarray<type>::sum()
{
	type total=(type)0;
	for(int i=0;i<_size;++i)
		total+=sparray[i];
	return total;
}

//self expand size to hold new element
//primarily used for importing large amount of data from file
//under most circumstances should be followed by fill_trim() function to discard unused space. 
template<class type>
void sarray<type>::fill(const int & index,const type & val,const int & buffer=1000)
{
	if(index<_size)
	{	_sizelim=index+1;
		sparray[index]=val;
		return;
	}
	if(index>=_size)
	{	this->resize(_size+buffer);
	}
	this->fill(index,val,buffer);
	return;
}	

template<class type>
void sarray<type>::fill_trim()
{
	if(_sizelim<_size)
		this->resize(_sizelim);
	return;
}


//constructors
/***********************************/
template<class type>
sarray<type>::sarray()
{
	init();
}

template<class type>
sarray<type>::sarray(const int & size)
{
	init(size);
}
template<class type>
sarray<type>::sarray(const int &size, const type &defau)
{
	init(size,defau);
}

template<class type>
void sarray<type>::init()
{
	_size=0;
	_sizelim=0;
	_defau=(type)0;
	sparray=(type*)malloc(sizeof(type)*0);
}
template<class type>
void sarray<type>::init(const int & size)
{
	_size=size;
	_defau=(type)0;
	sparray=(type*)malloc(sizeof(type)*_size);
}
template<class type>
void sarray<type>::init(const int & size, const type & defau)
{
	_size=size;
	_defau=defau;
	sparray=NULL;
	sparray=(type*)malloc(sizeof(type)*_size);
	for(int i=0;i<_size;++i)
		sparray[i]=_defau;
}
template<class type>
sarray<type>::sarray(const sarray<type> &temp)
{
	if(this==&temp)
		return;
	_size=temp._size;
	_defau=temp._defau;
	sparray=(type*)malloc(sizeof(type)*_size);
	for(int i=0;i<_size;++i)
		sparray[i]=temp.sparray[i];
	return;
}
template<class type>
sarray<type>::sarray(const std::string & temp)
{
	this->operator=(temp);
	return;
}
/*********************************************/
//destructor
/****************************************/
template<class type>
sarray<type>::~sarray()
{
	this->freem();
}
/****************************************/

//operator
/****************************************/
template<class type>
sarray<type> sarray<type>::operator+(const sarray<type> & temp)
{
	sarray<type> result(_size+temp._size);
    for(int i=0;i<_size;++i)
		result[i]=sparray[i];
	for(int i=0;i<temp._size;++i)
		result[_size+i]=temp.sparray[i];
	return result;
}

template<class type>
void sarray<type>::operator =(const std::string & temp)
{
	_size=temp.size();
	_defau='a';
	sparray=(char*)malloc(sizeof(char)*_size);
	for(int i=0;i<_size;++i)
		sparray[i]=temp[i];
	return;
}
template<class type>
void sarray<type>::operator =(const sarray<type> &temp)
{
	if(this==&temp)
		return;
	_size=temp._size;
	_defau=temp._defau;
	sparray=(type*)realloc(sparray,sizeof(type)*_size);
	for(int i=0;i<_size;++i)
		sparray[i]=temp.sparray[i];
	return;
}

template<class type>
bool sarray<type>::operator ==(const sarray<type> & temp)
{
	if(_size!=temp._size)
		return false;
	if(_size==0 && temp._size==0)
		return true;
	for(int i=0;i<_size;++i)
	{
		if(sparray[i]!=temp.sparray[i])
			return false;
	}
	return true;
}
template<class type>
bool sarray<type>::operator ==(const string & temp)
{
	if(_size!=temp.size())
		return false;
	for(int i=0;i<_size;++i)
	{
		if(sparray[i]!=temp[i])
			return false;
	}
	return true;
}

template<class type>
ostream & operator<<(ostream & os,const sarray<type> & temp)
{
	if(temp._size==0)
	{
		os<<"Empty array!";
		return os;
	}
	if(sizeof(type)>4)
	{
		os.flags(ios::fixed);os.precision(3);os.setf(ios::left);
	    for(int i=0;i<temp._size;++i)
	    {
		   os.width(9);
		   os<<temp.sparray[i]<<" ";
	    }
	}
	else if(sizeof(type)==4 || sizeof(type)==2)
	{	for(int i=0;i<temp._size;++i)
		   os<<temp.sparray[i]<<" ";
	}
	else
	{	for(int i=0;i<temp._size;++i)
		   os<<temp.sparray[i];
	}
	return os;
}
/************************************************/

template<class type>
type* sarray<type>::getarr()
{
	return this->sparray;
}
template<class type>
void sarray<type>::addheap()
{
    _size+=1;
    sparray=(type*)realloc(sparray,sizeof(type)*_size);
}
template<class type>
void sarray<type>::pushback(const type & temp)
{    
    this->addheap();
    sparray[_size-1]=temp;
}
template<class type>
sarray<type> sarray<type>::subseq(const int &start, const int &end)
{
	int size=end-start+1;
	sarray<type> temp(size);
	for(int i=0;i<size;++i)
		temp[i]=sparray[start+i];
	return temp;
}
template<class type>
type sarray<type>::smax()
{
	type max;
	max=sparray[0];
	for(int i=1;i<_size;++i)
	{
		if(max<sparray[i])
			max=sparray[i];
	}
	return max;
}
template<class type>
type sarray<type>::smax(sarray<int> & indexarr)
{
	type max;
	int index=0;
	max=sparray[0];
	for(int i=1;i<_size;++i)
	{
		if(max<sparray[i])
		{
			max=sparray[i];
			index=i;
		}
	}
	for(int i=index;i<_size;++i)
	{
		if(max==sparray[i])
			indexarr.pushback(i);
	}
	return max;
}
template<class type>
type sarray<type>::smin()
{
	type min;
	min=sparray[0];
	for(int i=1;i<_size;++i)
	{
		if(min>sparray[i])
			min=sparray[i];
	}
	return min;
}

template <class type>
type sarray<type>::smin(sarray<int> & indexarr)
{
	type min=this->smin();
	for(int i=0;i<_size;++i)
	{
		if(sparray[i]==min)
			indexarr.pushback(i);
	}
	return min;
}
template<class type>
void sarray<type>::merge(const type & temp)
{
	for(int i=0;i<_size;++i)
	{
		if(sparray[i]==temp)
			return;
	}
	this->pushback(temp);
}


template<class type>
void sarray<type>::resize(const int &size)
{
	if(_size==size || size<0)
		return;
	sparray=(type*)realloc(sparray,sizeof(type)*size);
	_size=size;
	return;
}
template<class type>
void sarray<type>::resize(const int &size,const type & defau)
{
	if(_size==size || size<0)
		return;
	_defau=defau;
	sparray=(type*)realloc(sparray,sizeof(type)*size);
	if(_size<size)
	{	for(int i=_size;i<size;++i)
			sparray[i]=_defau;
	}
	_size=size;
	return;
}
template<class type>
bool sarray<type>::find(const type & findval)
{
	for(int i=0;i<_size;++i)
	{
		if(sparray[i]==findval)
			return true;
	}
	return false;
}
template<class type>
sarray<type> sarray<type>::unique()
{
	sarray<type> res;
	res.pushback(sparray[0]);
	int i=0;
	while(i<_size)
	{
		if(res[res.size()-1]!=sparray[i])
			res.pushback(sparray[i]);
		i+=1;
	}
	return res;
}
template<class type>
string sarray<type>::tos()
{
	string str;
	stringstream oss;
	for(int i=0;i<_size;++i)
	{
		oss<<sparray[i];
		str+=oss.str()+" ";
		oss.clear();oss.str("");
	}
	return str;
}

template<class type>
bool sarray<type>::empty()
{
	return (_size==0 ? true : false);
}

template<class type>
sarray<type> sarray<type>::c(sarray<type> & b)
{
	int size=_size+b.size();
	sarray<type> res(size);
	for(int i=0;i<_size;++i)
		res[i]=sparray[i];
	for(int i=0;i<b.size();++i)
		res[_size+i]=b[i];
	return res;
}
