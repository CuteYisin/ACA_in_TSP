#include<bits/stdc++.h>
#define ll long long
#define random(x) rand()%(x)
using namespace std;

int main()
{
    freopen("in.txt","w",stdout);
    srand((int)time(0));
    for(int i=1;i<=30;i++)
    {
        for(int j=i+1;j<=30;j++)
        {
            cout<<i<<" "<<j<<" "<<random(50)+1<<endl;
        }
    }
    return 0;
}
