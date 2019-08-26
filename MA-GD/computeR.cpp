#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<vector>
#include"GD.h"
using namespace std;
int max_cluster_num=0;
int max_cluster_size=0;
double compute(vector<node> &rnode, vector<vector<int>> &G)//计算网络G的R值
{
	int N = rnode.size();
	int Maxdegreep;
	double largest = 0;
	int mid;
	double res = 0;
	vector<vector<int>> G2;//创建用于被攻击的邻接图
	G2.resize(N, vector<int>(N, 0));
	G2 = G;
	for (int i = 0; i < N - 1; i++)
	{
		Maxdegreep = MaxDegreePosition(rnode);
		DeleteNode(rnode, Maxdegreep, G2);
		rnode[Maxdegreep].attacked_rank = i;
		mid = MaxConnectedCluster(rnode, G2);
		largest = largest + (double)mid / N;
	}
	res = (double)largest / N;
	for (int i = 0; i < N; i++)//将点的状态全部还原
	{
		rnode[i].activation = 1;
	}
	calculate_degree(rnode, G);
	return res;
}

int MaxDegreePosition(vector<node> &rnode)//计算网络中degree最大的节点的位置
{
	int N = rnode.size();
	int maxnump = 0;
	int maxdegree = rnode[0].degree;
	for (int i = 1; i < N; i++)
	{
		if (rnode[i].degree > maxdegree)
		{
			maxnump = i;
			maxdegree = rnode[i].degree;
		}
	}
	if (maxdegree == 0)
	{
		for (int i = 0; i < rnode.size(); i++)
		{
			if (rnode[i].activation == 1)
			{
				maxnump = i;
				return maxnump;
			}
		}
	}
	return maxnump;
}

void calculate_degree(vector<node> &rnode, vector<vector<int>> &G)//是计算整个网络所有节点的degree
{
	int N = rnode.size();
	for (int p = 0; p < N; p++)
	{
		rnode[p].degree = 0;
		for (int i = 0; i < N; i++)
		{
			if (G[p][i] == 1)
			{
				rnode[p].degree++;
			}
		}
	}
}

void DeleteNode(vector<node> &rnode, int p, vector<vector<int>> &G)//p是删除的节点位置
{
	int N = rnode.size();
	rnode[p].activation = 0;//冻结节点
	for (int i = 0; i < N; i++)
	{
		if (G[p][i] == 1)
		{
			G[p][i] = 0;
			G[i][p] = 0;
			rnode[p].degree--;
			rnode[i].degree--;
		}
	}
}

int MaxConnectedCluster(vector<node> &rnode, vector<vector<int>> &G)//返回最大连接簇的节点数量
{
	int N = rnode.size();
	vector<int> search(N);//存储节点是否被检索过的信息
	int maxconnectnum = 1;
	int p = 0;
	int IN = 0;
	int maxscmf = 0;
	for (int i = 0; i < N; i++)
	{
		search[i] = i;
		rnode[i].max_cluster_flag = 0;
	}
	max_cluster_num = 1;
	while (whether_search_end(search, rnode))
	{
		while (search[p] == -1 || rnode[p].activation == 0)
		{
			p++;
		}
		IN = Include_Nodes(p, G, search, rnode);
		if (maxconnectnum <= IN)
		{
			maxconnectnum = IN;
			maxscmf = max_cluster_num;
			set_max_cluster_flag_0(rnode, 0);
		}
		else
		{
			set_max_cluster_flag_0(rnode,1);
		}
		max_cluster_num++;
	}
	max_cluster_num=maxscmf;
	return maxconnectnum;
}

int MaxConnectedCluster_nosignal(vector<node> &rnode, vector<vector<int>> &G)//返回最大连接簇的节点数量,不带标记
{
	int N = rnode.size();
	vector<int> search(N);//存储节点是否被检索过的信息
	int maxconnectnum = 1;
	int p = 0;
	int IN = 0;
	for (int i = 0; i < N; i++)
	{
		search[i] = i;
	}
	while (whether_search_end(search, rnode))
	{
		while (search[p] == -1 || rnode[p].activation == 0)
		{
			p++;
		}
		IN = Include_Nodes(p, G, search, rnode);
		if (maxconnectnum <= IN)
		{
			maxconnectnum = IN;
			
		}

	}
	return maxconnectnum;
}

int Include_Nodes_nosingal(int p, vector<vector<int>> &G, vector<int> &search, vector<node> &rnode)//输入节点位置p，检索该节点连接的未被检索过的节点数，返回其值.不带标记
{
	int N = rnode.size();
	int res;
	if (search[p] != -1 && rnode[p].activation != 0)
	{
		search[p] = -1;
		int maxconnectnum = 1;
		for (int i = 0; i < rnode[p].neibor.size(); i++)
		{
			if (G[p][rnode[p].neibor[i]] == 1 && search[rnode[p].neibor[i]] != -1)
			{
				maxconnectnum = maxconnectnum + Include_Nodes_nosingal(rnode[p].neibor[i], G, search, rnode);
			}
		}
		res = maxconnectnum;
	}
	else
	{
		res = 0;
	}
	return res;
}

bool whether_search_end(vector<int> &search, vector<node> &rnode)//返回false说明全部点都被检索过了，返回true则还有点未被检索
{
	int N = rnode.size();
	bool u = false;
	for (int i = 0; i < N; i++)
	{
		if (search[i] != -1 && rnode[i].activation != 0)
		{
			u = true;
		}
	}
	return u;
}

int Include_Nodes(int p, vector<vector<int>> &G, vector<int> &search, vector<node> &rnode)//输入节点位置p，检索该节点连接的未被检索过的节点数，返回其值
{
	int N = rnode.size();
	int res;
	if (search[p] != -1 && rnode[p].activation != 0)
	{
		search[p] = -1;
		rnode[p].max_cluster_flag = max_cluster_num;
		int maxconnectnum = 1;
		for (int i = 0; i < rnode[p].neibor.size(); i++)
		{
			if (G[p][rnode[p].neibor[i]] == 1 && search[rnode[p].neibor[i]] != -1)
			{
				maxconnectnum = maxconnectnum + Include_Nodes(rnode[p].neibor[i], G, search, rnode);
			}
		}
		res = maxconnectnum;
	}
	else
	{
		res = 0;
	}
	return res;
}




void backupdegree(vector<node> &rnode, vector<vector<int>> &G)
{
	for (int i = 0; i < rnode.size(); i++)
	{
		for (int j = 0; j < rnode.size(); j++) {
			if (G[i][j] == 1) {
				rnode[i].checkdegree++;
			}
		}
	}
}

bool checkdegree(vector<node> &rnode, vector<vector<int>> &G)
{
	for (int i = 0; i < rnode.size(); i++)
	{
		rnode[i].degree = 0;
		for (int j = 0; j < rnode.size(); j++) {
			if (G[i][j] == 1) {
				rnode[i].degree++;
			}
		}
	}
	for (int i = 0; i < rnode.size(); i++)
	{
		if (rnode[i].degree != rnode[i].checkdegree)
		{
			return false;
		}
	}
	return true;

}


void translate_G_to_neibor(vector<node> &rnode, vector<vector<int>> &G)
{
	for (int i = 0; i < rnode.size(); i++) {
		rnode[i].degree = 0;
		vector<int>().swap(rnode[i].neibor);
		for (int j = 0; j < rnode.size(); j++) {
			if (G[i][j] == 1 && rnode[j].activation != 0) {
				rnode[i].neibor.push_back(j);
				rnode[i].degree++;
			}
		}
	}
}


void recovernode(vector<node> &rnode)
{
	for (int i = 0; i < rnode.size(); i++) {
		rnode[i].search_flag = 1;
	}
}

void set_max_cluster_flag_0(vector<node> &rnode, int pattern)//pattern:1是当前max_cluster_num的删除；pattern:0不是当前max_cluster_num的删除
{
	for (int i = 0; i < rnode.size(); i++) {
		if (pattern == 1) {
			if (rnode[i].max_cluster_flag == max_cluster_num) {
				rnode[i].max_cluster_flag = 0;
			}
		}
		else {
			if (rnode[i].max_cluster_flag != max_cluster_num) {
				rnode[i].max_cluster_flag = 0;
			}
		}
	}
}

void calculateGD(vector<node> &rnode, vector<vector<int>> &G)//计算eij--最优选择概率
{
	int Maxdegreep;
	vector<vector<int>> G2;//创建用于被攻击的邻接图
	G2.resize(rnode.size(), vector<int>(rnode.size(), 0));
	G2 = G;
	for (int i = 0; i < rnode.size(); i++) {
		for (int j = 0; j < rnode.size(); j++) {
			rnode[i].delta_cp_current[j] = 0;
		}
	}
	for (int i = 0; i < rnode.size() - 1; i++)
	{
		//rnode的邻居节点不会删除，只参与节点度计算
		Maxdegreep = MaxDegreePosition(rnode);
		DeleteNode(rnode, Maxdegreep, G2);
		rnode[Maxdegreep].attacked_rank = i;
		max_cluster_size = MaxConnectedCluster(rnode, G2);
		calculate_cp(rnode, G2);
	}
	update_cp(rnode);
	for (int i = 0; i < rnode.size(); i++)//将点的状态全部还原
	{
		rnode[i].activation = 1;
	}
	calculate_degree(rnode, G);
}

void calculate_cp(vector<node> &rnode, vector<vector<int>> &G)//计算边连接概率
{

	double mccn = 0,loop1;
	for (int i = 0; i < rnode.size(); i++) {
		if (rnode[i].activation == 1) {
			for (int j = 0; j < i; j++) {
				if (rnode[j].activation == 1&&i!=j) {

					if (rnode[i].max_cluster_flag != max_cluster_num) {
						if (rnode[j].max_cluster_flag == max_cluster_num) {
							rnode[i].delta_cp_current[j] += gradient;
							//decrease_other_cp(rnode[i], j, (double)gradient / (rnode.size() - 1));
							rnode[j].delta_cp_current[i] += gradient;
							//decrease_other_cp(rnode[j], i, (double)gradient / (rnode.size() - 1));
						}
					}
					else
					{
						if (rnode[j].max_cluster_flag != max_cluster_num) {
							rnode[i].delta_cp_current[j] += gradient;
							//decrease_other_cp(rnode[i], j, (double)gradient / (rnode.size() - 1));
							rnode[j].delta_cp_current[i] += gradient;
							//decrease_other_cp(rnode[j], i, (double)gradient / (rnode.size() - 1));
						}
						else
						{
							if (G[i][j] == 1) {
								G[i][j] = 0;
								G[j][i] = 0;
								rnode[i].neibor.erase(find(rnode[i].neibor.begin(), rnode[i].neibor.end(), j));
								rnode[j].neibor.erase(find(rnode[j].neibor.begin(), rnode[j].neibor.end(), i));
								rnode[i].degree--;
								rnode[j].degree--;
								mccn = MaxConnectedCluster_nosignal(rnode, G);
								if (max_cluster_size > mccn) {
									loop1 = max_cluster_size - mccn;
									rnode[i].delta_cp_current[j] += gradient *loop1*g_rate;
									double middd;
									middd = gradient *loop1* g_rate;
									//decrease_other_cp(rnode[i], j, middd / (double)(rnode.size() - 1));
									rnode[j].delta_cp_current[i] += gradient *loop1*g_rate;
									//decrease_other_cp(rnode[j], i, middd / (double)(rnode.size() - 1));
								}
								else
								{
									rnode[i].delta_cp_current[j] -= gradient * g_rate;
									double middd;
									middd = gradient * g_rate;
									//decrease_other_cp(rnode[i], j, -1.0*middd / (double)(rnode.size() - 1));
									rnode[j].delta_cp_current[i] -= gradient * g_rate;
									//decrease_other_cp(rnode[j], i, -1.0*middd / (double)(rnode.size() - 1));
								}
								G[i][j] = 1;
								G[j][i] = 1;
								rnode[i].neibor.push_back(j);
								rnode[j].neibor.push_back(i);
								rnode[i].degree++;
								rnode[j].degree++;
							}
						}
					}
				}
			}
			for (int j = i + 1; j < rnode.size(); j++) {
				if (rnode[j].activation == 1&&i!=j) {

					if (rnode[i].max_cluster_flag != max_cluster_num) {
						if (rnode[j].max_cluster_flag == max_cluster_num) {
							rnode[i].delta_cp_current[j] += gradient;
							//decrease_other_cp(rnode[i], j, (double)gradient / (rnode.size() - 1));
							rnode[j].delta_cp_current[i] += gradient;
							//decrease_other_cp(rnode[j], i, (double)gradient / (rnode.size() - 1));
						}
					}
					else
					{
						if (rnode[j].max_cluster_flag != max_cluster_num) {
							rnode[i].delta_cp_current[j] += gradient;
							//decrease_other_cp(rnode[i], j, (double)gradient / (rnode.size() - 1));
							rnode[j].delta_cp_current[i] += gradient;
							//decrease_other_cp(rnode[j], i, (double)gradient / (rnode.size() - 1));
						}
						else
						{
							if (G[i][j] == 1) {
								G[i][j] = 0;
								G[j][i] = 0;
								rnode[i].neibor.erase(find(rnode[i].neibor.begin(), rnode[i].neibor.end(), j));
								rnode[j].neibor.erase(find(rnode[j].neibor.begin(), rnode[j].neibor.end(), i));
								rnode[i].degree--;
								rnode[j].degree--;
								mccn = MaxConnectedCluster_nosignal(rnode, G);
								if (max_cluster_size > mccn) {
									loop1 = max_cluster_size - mccn;
									rnode[i].delta_cp_current[j] += gradient * loop1*g_rate;
									double middd;
									middd = gradient * loop1* g_rate;
									//decrease_other_cp(rnode[i], j, middd / (double)(rnode.size() - 1));
									rnode[j].delta_cp_current[i] += gradient * loop1*g_rate;
									//decrease_other_cp(rnode[j], i, middd / (double)(rnode.size() - 1));
								}
								else
								{
									rnode[i].delta_cp_current[j] -= gradient * g_rate;
									double middd;
									middd = gradient * g_rate;
									//decrease_other_cp(rnode[i], j, -1.0*middd / (double)(rnode.size() - 1));
									rnode[j].delta_cp_current[i] -= gradient * g_rate;
									//decrease_other_cp(rnode[j], i, -1.0*middd / (double)(rnode.size() - 1));
								}
								G[i][j] = 1;
								G[j][i] = 1;
								rnode[i].neibor.push_back(j);
								rnode[j].neibor.push_back(i);
								rnode[i].degree++;
								rnode[j].degree++;
							}
						}
					}
				}
			}
		}
		
	}
}

void decrease_other_cp(node &rnode, int j, double d_value)//j为增加的那条边,其他边减少值
{
	vector<double> random_value;
	int random_p;
	if (d_value == 0) return;
	for (int i = 0; i < rnode.cp.size()-1; i++) {
		random_value.push_back(i*0.4/(rnode.cp.size()-1));
	}
	random_p = (double)(rand() / (RAND_MAX + 1.0))*random_value.size();
	for (int i = 0; i < j; i++) {
		
		rnode.delta_cp_current[i] -= d_value*(random_value[(random_p+i)%(rnode.cp.size()-1)]+0.8);


	//	rnode.delta_cp_current[i] -= d_value/(rnode.cp.size()-1.0);
	}
	for (int i = j + 1; i < rnode.delta_cp_current.size(); i++) {
		rnode.delta_cp_current[i] -= d_value * (random_value[(random_p + i) % (rnode.cp.size() - 1)]+0.8);

	//	rnode.delta_cp_current[i] -= d_value/ (rnode.cp.size() - 1.0);
	}
}

void update_cp(vector<node> &rnode)
{
	double lr = 0;
	lr = learning_rate;
	for (int i = 0; i < rnode.size(); i++) {
		for (int j = 0; j < rnode.size(); j++) {
			decrease_other_cp(rnode[i], j, rnode[i].delta_cp_current[j]);
		}
	}
	for (int i = 0; i < rnode.size(); i++) {
		for (int j = 0; j < rnode.size(); j++) {
			rnode[i].cp[j] = rnode[i].cp[j] + lr * (rnode[i].delta_cp_current[j] - belta*rnode[i].delta_cp_before[j]);
			rnode[i].delta_cp_before[j] = rnode[i].delta_cp_current[j];
		}
	}
}