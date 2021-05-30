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
//��Ϣ��������,��������ʽ����,ȫ����Ϣ�ػӷ�ϵ��,�ֲ���Ϣ�ػӷ�ϵ��

double dis[N][N];
double Lnn;

//��Ⱥϵͳ
class AntColonySystem
{
private:
	double info[N][N], visible[N][N];//�ڵ�֮�����Ϣ����,�ڵ�֮�������ʽ��Ϣ��
public:
	AntColonySystem(){}

	//���㵱ǰ�ڵ㵽��һ�ڵ�ת�Ƶĸ���
	double Transition(int i, int j)
	{
	    if(i!=j)  return (pow(info[i][j],alpha)*pow(visible[i][j],beta));
        else  return 0.0;
	}

	//����ʽ4�ֲ�����
	void UpdateLocalPathRule(int i, int j)
	{
	    info[i][j]=info[j][i]=(1.0-alpha1)*info[i][j]+alpha1*(1.0/(N*Lnn));
	}

	//��ʼ��
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

	//ȫ����Ϣ�ظ���
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

//���ϸ���
class ACSAnt
{
private:
	AntColonySystem* antColony;
protected:
	int startCity, cururentCity;//��ʼ���б�ţ���ǰ���б��
	int allowed[N];             //���ɱ�
	int Tour[N][2];             //һ����·����������ɣ�����currentcity��nextcity��
	int currentTourIndex;       //��ǰ·����������0��ʼ���洢���Ͼ������еı��
public:
	ACSAnt(AntColonySystem* acs, int start)
	{
		antColony=acs;
		startCity=start;
	}

	//ѡ����һ�ڵ�
	int Choose()
	{
	    int nextCity=-1;
        double q=rand()/(double)RAND_MAX;    //����һ��0~1֮��������q

        if(q<=0.1)	                         //��������ͱȽ�
        {
            double probability=-1.0;
            rep0(i,N)
            {
                if(allowed[i])
                {
                    double prob=antColony->Transition(cururentCity,i);  //���㵱ǰ�ڵ�ת�Ƶ���һ�ڵ�ĸ���
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
            //������ת��
            double p=rand()/(double)RAND_MAX;	//����һ�������,�����ж������ĸ������
            double sum=0.0;
            double probability=0.0;	            //���ʵ�����㣬p�����ĸ�����Σ���õ���ת�Ƶķ���

            rep0(i,N)
            {
                if(allowed[i])  sum+=antColony->Transition(cururentCity,i);
            }

            rep0(j,N)
            {
                if(allowed[j]&&sum>0)
                {
                    probability+=antColony->Transition(cururentCity,j)/sum; //������jת�Ƶĸ���
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

	//��ʼ����
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
            toCity=Choose();	                                  //ѡ����һ���ڵ�
            if(toCity>=0)
            {
                MoveToNextCity(toCity);                           //�ƶ�����һ���ڵ�
                antColony->UpdateLocalPathRule(endCity, toCity);  //���оֲ�����
                cururentCity=toCity;
            }
        }while(toCity>=0);
        MoveToNextCity(startCity);
        antColony->UpdateLocalPathRule(endCity,startCity);

        return *Tour;
	}

	//�ƶ�����һ�ڵ�
	void MoveToNextCity(int nextCity)
	{
	    allowed[nextCity]=0;
        Tour[currentTourIndex][0]=cururentCity;
        Tour[currentTourIndex][1]=nextCity;
        currentTourIndex++;
        cururentCity=nextCity;
	}
};

//ѡ����һ�����ڽ��ڵ�
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

//��һ���ڵ�������ھ��뷽�����㳤��Lnn
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

//�ɾ����ʾ��������֮��ľ���
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

//��þ���n�����е�·������
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

	int node=rand()%N;	                //���ѡ��һ���ڵ����������ڷ����õ�Lnn
	Lnn = CalAdjacentDistance(node);
	double initInfo = 1 / (N * Lnn);
	acs->InitParameter(initInfo);	    //����ʽ1��ʼ����Ⱥ��Ϣ��ǿ��

	int globalTour[N][2];               //ȫ������·������
	double globalBestLength=0.0;	    //ȫ������·������

	static clock_t Start,Finish,s,f;
	s=clock();

	rep0(i,Max)
	{
        Start=clock();

		int localTour[N][2];	        //�ֲ�����·������
		double localBestLength=0.0;	    //�ֲ�����·������
		double tourLength;	            //��ǰ·������

		rep0(j,M)
		{
			int* tourPath=ants[j]->Search();
			tourLength=calculateSumOfDistance(tourPath);

			//�ֲ��Ƚϣ�����¼·���ͳ���
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

		//ȫ�ֱȽϣ�����¼·���ͳ���
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

		cout<<"��"<<i+1<<"��������·��:"<<localBestLength<<" "<<endl;
		printf("�˴ε�������ʱ��:%fs\n",time_second);
		rep0(m,N)  cout<<localTour[m][0]+1<<" ";
		cout<<endl;
		cout<<endl;
	}

    f=clock();
    double t=double(f-s)/CLOCKS_PER_SEC;

	cout<<"ȫ������·������:"<<globalBestLength<<endl;
	printf("������ʱ��:%fs\n",t);
	cout<<"ȫ������·��:";
	rep0(m,N)  cout<<globalTour[m][0]+1<<" ";
	cout<<endl;
	system("pause");
	return 0;
}

