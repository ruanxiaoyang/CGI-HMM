int na2id(const char & _na)
{
	if(_na=='A' || _na=='a')
		return 0;
	if(_na=='T' || _na=='t')
		return 1;
	if(_na=='C' || _na=='c')
		return 2;
	if(_na=='G' || _na=='g')
		return 3;
	return (int)rand()%2;
}

template <class T>
string ntos(const T& num)
{
	std::ostringstream oss;
	oss<<num;
	return oss.str();
}

template <class T>
T stn(const std::string& s, std::ios_base& (*f)(std::ios_base&)=std::dec)	//*f is pointer to function which takes ios_base and return ios_base
{
	T t;
	std::istringstream iss(s);
	if(!(iss >> f >> t).fail())
	{
		return t;
	}
	else
	{
		cout<<"error"<<endl;
		return 0;
	}
}
