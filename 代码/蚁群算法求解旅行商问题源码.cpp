#include<bits/stdc++.h>
#include <ctime>
#define ll long long
#define U unsigned
#define sqr(x) ((x)*(x))
#define rep0(x,n) for(int x=0;x<n;x++)
#define endl "\n"
using namespace std;

const int N=30,M=40;
const int Max=500;
double alpha=2,beta=5,rou=0.1,alpha1=0.1;
//信息启发因子,期望启发式因子,全局信息素挥发系数,局部信息素挥发系数

double dis[N][N];
double Lnn;

//蚁群系统
class AntColonySystem
{
private:
	double info[N][N], visible[N][N];//节点之间的信息素量,节点之间的启发式信息量
public:
	AntColonySystem(){}

	//计算当前节点到下一节点转移的概率
	double Transition(int i, int j)
	{
	    if(i!=j)  return (pow(info[i][j],alpha)*pow(visible[i][j],beta));
        else  return 0.0;
	}

	//根据式4局部更新
	void UpdateLocalPathRule(int i, int j)
	{
	    info[i][j]=info[j][i]=(1.0-alpha1)*info[i][j]+alpha1*(1.0/(N*Lnn));
	}

	//初始化
	void InitParameter(double value)
	{
	    rep0(i,N)
        {
            rep0(j,N)
            {
                info[i][j]=info[j][i]=value;
                if(i!=j) visible[i][j]=visible[j][i]=1.0/dis[i][j];
            }
        }
	}

	//全局信息素更新
	void UpdateGlobalPathRule(int* bestTour, int globalBestLength)
	{
	    rep0(i,N)
        {
            int row=*(bestTour+2*i);
            int col=*(bestTour+2*i+1);
            info[row][col]=info[col][row]=(1.0-rou)*info[row][col]+rou*(1.0/globalBestLength);
        }
	}
};

//蚂蚁个体
class ACSAnt
{
private:
	AntColonySystem* antColony;
protected:
	int startCity, cururentCity;//初始城市编号，当前城市编号
	int allowed[N];             //禁忌表
	int Tour[N][2];             //一个个路径段序列组成，即（currentcity，nextcity）
	int currentTourIndex;       //当前路径索引，从0开始，存储蚂蚁经过城市的编号
public:
	ACSAnt(AntColonySystem* acs, int start)
	{
		antColony=acs;
		startCity=start;
	}

	//选择下一节点
	int Choose()
	{
	    int nextCity=-1;
        double q=rand()/(double)RAND_MAX;    //产生一个0~1之间的随机数q

        if(q<=0.1)	                         //随机采样和比较
        {
            double probability=-1.0;
            rep0(i,N)
            {
                if(allowed[i])
                {
                    double prob=antColony->Transition(cururentCity,i);  //计算当前节点转移到下一节点的概率
                    if(prob>probability)
                    {
                        nextCity=i;
                        probability=prob;
                    }
                }
            }
        }
        else
        {
            //按概率转移
            double p=rand()/(double)RAND_MAX;	//生成一个随机数,用来判断落在哪个区间段
            double sum=0.0;
            double probability=0.0;	            //概率的区间点，p落在哪个区间段，则该点是转移的方向

            rep0(i,N)
            {
                if(allowed[i])  sum+=antColony->Transition(cururentCity,i);
            }

            rep0(j,N)
            {
                if(allowed[j]&&sum>0)
                {
                    probability+=antColony->Transition(cururentCity,j)/sum; //往城市j转移的概率
                    if(probability>=p||(p>0.9999&&probability>0.9999))
                    {
                        nextCity=j;
                        break;
                    }
                }
            }
        }
        return nextCity;
	}

	//开始搜索
	int* Search()
	{
	    cururentCity=startCity;
        int toCity;
        currentTourIndex=0;
        rep0(i,N)
        {
            allowed[i]=1;
        }
        allowed[cururentCity]=0;
        int endCity;
        int count=0;
        do
        {
            count++;
            endCity=cururentCity;
            toCity=Choose();	                                  //选择下一个节点
            if(toCity>=0)
            {
                MoveToNextCity(toCity);                           //移动到下一个节点
                antColony->UpdateLocalPathRule(endCity, toCity);  //进行局部更新
                cururentCity=toCity;
            }
        }while(toCity>=0);
        MoveToNextCity(startCity);
        antColony->UpdateLocalPathRule(endCity,startCity);

        return *Tour;
	}

	//移动到下一节点
	void MoveToNextCity(int nextCity)
	{
	    allowed[nextCity]=0;
        Tour[currentTourIndex][0]=cururentCity;
        Tour[currentTourIndex][1]=nextCity;
        currentTourIndex++;
        cururentCity=nextCity;
	}
};

//选择下一个最邻近节点
int ChooseNextNode(int currentNode, int visitedNode[])
{
    int nextNode=-1;
	double shortDistance=0.0;
	rep0(i,N)
	{
		if(visitedNode[i])
		{
			if(shortDistance==0.0)
			{
				shortDistance=dis[currentNode][i];
				nextNode=i;
			}
			if (shortDistance<dis[currentNode][i])
			{
				nextNode=i;
			}
		}
	}
	return nextNode;
}

//给一个节点由最近邻距离方法计算长度Lnn
double CalAdjacentDistance(int node)
{
    double sum=0.0;
	int visitedNode[N];
	rep0(j,N)
	{
		visitedNode[j]=1;
	}
	visitedNode[node]=0;
	int currentNode=node;
	int nextNode;
	do
	{
		nextNode=ChooseNextNode(currentNode, visitedNode);
		if(nextNode>=0)
		{
			sum+=dis[currentNode][nextNode];
			currentNode=nextNode;
			visitedNode[currentNode]=0;
		}
	}while(nextNode>=0);
	sum+=dis[currentNode][node];
	return sum;
}

//由矩阵表示两两城市之间的距离
void calculateAllDistance()
{
    freopen("in.txt","r",stdin);
	rep0(i,N)
	{
		for(int j=i+1;j<N;j++)
		{
			int a,b,c;cin>>a>>b>>c;
			dis[i][j]=dis[j][i]=c;
		}
	}
}

//获得经过n个城市的路径长度
double calculateSumOfDistance(int* tour)
{
    double sum=0;
	rep0(i,N)
	{
		int row= *(tour+2*i);
		int col= *(tour+2*i+1);
		sum+=dis[row][col];
	}
	return sum;
}

int main()
{
	calculateAllDistance();

	AntColonySystem* acs = new AntColonySystem();
	ACSAnt* ants[M];
	for (int k = 0; k < M; k++)  ants[k] = new ACSAnt(acs, (int)(k%N));

	time_t timer;time(&timer);
	unsigned long seed = timer;
	seed%=56000;
	srand((unsigned int)seed);

	int node=rand()%N;	                //随机选择一个节点计算由最近邻方法得到Lnn
	Lnn = CalAdjacentDistance(node);
	double initInfo = 1 / (N * Lnn);
	acs->InitParameter(initInfo);	    //根据式1初始化蚁群信息素强度

	int globalTour[N][2];               //全局最优路径序列
	double globalBestLength=0.0;	    //全局最优路径长度

	static clock_t Start,Finish,s,f;
	s=clock();

	rep0(i,Max)
	{
        Start=clock();

		int localTour[N][2];	        //局部最优路径序列
		double localBestLength=0.0;	    //局部最优路径长度
		double tourLength;	            //当前路径长度

		rep0(j,M)
		{
			int* tourPath=ants[j]->Search();
			tourLength=calculateSumOfDistance(tourPath);

			//局部比较，并记录路径和长度
			if(tourLength<localBestLength||abs(localBestLength-0.0)<0.000001)
			{
				rep0(m,N)
				{
					int row=*(tourPath+2*m);
					int col=*(tourPath+2*m+1);
					localTour[m][0]=row;
					localTour[m][1]=col;
				}
				localBestLength=tourLength;
			}
		}

		//全局比较，并记录路径和长度
		if(localBestLength<globalBestLength||abs(globalBestLength-0.0)<0.000001)
		{
			rep0(m,N)
			{
				globalTour[m][0]=localTour[m][0];
				globalTour[m][1]=localTour[m][1];
			}
			globalBestLength=localBestLength;
		}

		acs->UpdateGlobalPathRule(*globalTour,globalBestLength);

		Finish = clock();
        double time_second=double(Finish-Start)/CLOCKS_PER_SEC;

		cout<<"第"<<i+1<<"迭代最优路径:"<<localBestLength<<" "<<endl;
		printf("此次迭代运行时间:%fs\n",time_second);
		rep0(m,N)  cout<<localTour[m][0]+1<<" ";
		cout<<endl;
		cout<<endl;
	}

    f=clock();
    double t=double(f-s)/CLOCKS_PER_SEC;

	cout<<"全局最优路径长度:"<<globalBestLength<<endl;
	printf("总运行时间:%fs\n",t);
	cout<<"全局最优路径:";
	rep0(m,N)  cout<<globalTour[m][0]+1<<" ";
	cout<<endl;
	system("pause");
	return 0;
}

