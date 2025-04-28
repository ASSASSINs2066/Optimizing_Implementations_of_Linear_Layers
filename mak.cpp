#include<bits/stdc++.h>
using namespace std;

#define int long long
const int M=1e7+7;
int re()
{
	int x=0,w=1;
	char y=getchar();
	while(y<'0'||y>'9')
	{
		if(y=='-')w=-1;
		y=getchar();
	}
	while(y<='9'&&y>='0')
	{
		x=x*10+y-'0';
		y=getchar();
	}return x*w;
}
int n,m;


signed main()
{
    n=32;
    freopen("A.txt","r",stdin);
	freopen("S.txt","w",stdout);
    for(int i=1;i<=n;i++)
    {
        int y=0;
        cout<<"0x";
        for(int j=1;j<=n;j++)
        {
            int x=re();
            y=y*2+x;
            if(j%4==0)
            {
                printf("%x",y);
                y=0;
            }
        }cout<<","<<endl;
    }
	return 0;
}