#include <stdio.h>
#include <time.h>

using namespace std;
string getdate()
{
	time_t t=time(NULL);
	tm* now=localtime(&t);
	char buffer[20];
	sprintf(buffer,"%d/%d/%d %d:%d",now->tm_year+1900,now->tm_mon+1,now->tm_mday,now->tm_hour,now->tm_min);
	return string(buffer);
}
