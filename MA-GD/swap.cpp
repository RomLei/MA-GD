#include"GD.h"
#include"vector"
#include<algorithm>
using namespace std;

int update_network(vector<node> &rnode, vector<vector<int>> &G)
{
	int k=-1, l=-1;
	int sucflag = 0;
	int ranki_in_j, rankj_in_i;
	bool suc = true;
	for (int i = 0,search_bound=0; i < rnode.size(); i++) {
		search_bound = rnode[i].degree;
		for (int jnumber = 0,j,start_num1=0,start_num2=0; jnumber < search_bound&&jnumber<rnode.size(); jnumber++) {
			k = -1;
			l = -1;
			suc = true;
			j = rnode[i].prefer_edge_rank[jnumber];
			if (jnumber != 0 && rnode[i].cp[j] == rnode[i].cp[rnode[i].prefer_edge_rank[jnumber - 1]]) {
				search_bound++;
			}
			ranki_in_j = find_rank(rnode, j, i);
			if (G[i][j] == 0&&i!=j) {
				if (ranki_in_j < rnode[j].degree)
				{
						k = find_change_neibor(rnode, i, jnumber, start_num1);
						if (k != -1) {
							l = find_change_neibor(rnode, j, ranki_in_j, start_num2);
						}
						if (l != -1 &&k!=l) {
							if (G[k][l] == 0) {
								rnode[i].neibor.erase(find(rnode[i].neibor.begin(), rnode[i].neibor.end(), k));
								rnode[k].neibor.erase(find(rnode[k].neibor.begin(), rnode[k].neibor.end(), i));
								rnode[i].neibor.push_back(j);
								rnode[j].neibor.push_back(i);
								rnode[j].neibor.erase(find(rnode[j].neibor.begin(), rnode[j].neibor.end(), l));
								rnode[l].neibor.erase(find(rnode[l].neibor.begin(), rnode[l].neibor.end(), j));
								rnode[l].neibor.push_back(k);
								rnode[k].neibor.push_back(l);
								G[i][j] = 1;
								G[j][i] = 1;
								G[k][l] = 1;
								G[l][k] = 1;
								G[i][k] = 0;
								G[k][i] = 0;
								G[j][l] = 0;
								G[l][j] = 0;
								sucflag++;
							}
						}
				}
			}
		}
	}
	return sucflag;
}

void update_perfer_edge_rank(vector<node> &rnode)//将边连接概率按索引排序
{
	int i = 0;
	for (; i < rnode.size(); i++) {
		sort_part(rnode[i]);
	}
}

void sort_part(node &rnode)
{
	vector<vector<double>> ranklist(2,vector<double>(rnode.cp.size(),0));
	for (int i = 0; i < rnode.cp.size(); i++) {
		ranklist[0][i] = rnode.cp[i];
		ranklist[1][i] = i;
	}
	quick_sort(ranklist, 0, rnode.cp.size() - 1);
	for (int j = 0; j < rnode.cp.size(); j++) {
		rnode.prefer_edge_rank[j] = ranklist[1][j];
	}
}

void quick_sort(vector<vector<double>> &ranklist, int start, int end)
{
	if (start > end) {
		return;
	}
	int i = start, j = end;
	double mid = ranklist[0][start];
	int midrank = ranklist[1][start];
	while (i < j) {
		while (i<j&&ranklist[0][j]<=mid) {
			j--;
		}
		ranklist[0][i] = ranklist[0][j];
		ranklist[1][i] = ranklist[1][j];
		while (i<j&&ranklist[0][i]>mid) {
			i++;
		}
		ranklist[0][j] = ranklist[0][i];
		ranklist[1][j] = ranklist[1][i];
	}
	ranklist[0][j] = mid;
	ranklist[1][j] = midrank;
	quick_sort(ranklist, start, i - 1);
	quick_sort(ranklist, i + 1, end);
}

int find_rank(vector<node> &rnode, int node1, int node2)//node2在node1中的排位
{
	for (int i = 0,t=0; i < rnode.size(); i++,t++) {
		if (i != 0 && rnode[node1].cp[rnode[node1].prefer_edge_rank[i]] == rnode[node1].cp[rnode[node1].prefer_edge_rank[i - 1]]) {
			t--;
		}
		if (rnode[node1].prefer_edge_rank[i] == node2) {
			return t;
		}
	}
	return -1;
}

int find_change_neibor(vector<node> &rnode, int node1, int jnumber,int start_num)
{
	vector<int> k;
	for (int t = start_num, kkk = 0; t < rnode[node1].neibor.size(); t++) {
		kkk = rnode[node1].neibor[t];
		if (find(rnode[node1].prefer_edge_rank.begin(), rnode[node1].prefer_edge_rank.end(), kkk) > rnode[node1].prefer_edge_rank.begin() + jnumber) {
			k.push_back(kkk);
		}
	}
	if (!k.size()) {
		return -1;
	}
	return k[(rand()/(RAND_MAX+1.0))*k.size()];
}

bool swap_edge(vector<node> &rnode, vector<vector<int>> &G, int &falsenum)//随机交换一层中的两条边,返回false说明交换失败，true成功
{
	int node1, node1_neibor, node2, node2_neibor;
	node1 = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
	node1_neibor = ((double)rand() / RAND_MAX)*(rnode[node1].neibor.size() - 1);
	node1_neibor = rnode[node1].neibor[node1_neibor];
	node2 = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
	node2_neibor = ((double)rand() / RAND_MAX)*(rnode[node2].neibor.size() - 1);
	node2_neibor = rnode[node2].neibor[node2_neibor];
	if (find(rnode[node1].neibor.begin(), rnode[node1].neibor.end(), node2_neibor) != rnode[node1].neibor.end())
	{
		falsenum++;
		return false;
	}
	else if (find(rnode[node2].neibor.begin(), rnode[node2].neibor.end(), node1_neibor) != rnode[node2].neibor.end())
	{
		falsenum++;
		return false;
	}
	else if (node1 == node2 || node1_neibor == node2 || node1 == node2_neibor || node2 == node1_neibor)
	{
		falsenum++;
		return false;
	}
	else
	{
		rnode[node1].neibor.erase(find(rnode[node1].neibor.begin(), rnode[node1].neibor.end(), node1_neibor));
		rnode[node1_neibor].neibor.erase(find(rnode[node1_neibor].neibor.begin(), rnode[node1_neibor].neibor.end(), node1));
		rnode[node2].neibor.erase(find(rnode[node2].neibor.begin(), rnode[node2].neibor.end(), node2_neibor));
		rnode[node2_neibor].neibor.erase(find(rnode[node2_neibor].neibor.begin(), rnode[node2_neibor].neibor.end(), node2));
		rnode[node1].neibor.push_back(node2_neibor);
		rnode[node2_neibor].neibor.push_back(node1);
		rnode[node2].neibor.push_back(node1_neibor);
		rnode[node1_neibor].neibor.push_back(node2);
		G[node1][node1_neibor] = 0;
		G[node1_neibor][node1] = 0;
		G[node2][node2_neibor] = 0;
		G[node2_neibor][node2] = 0;
		G[node1][node2_neibor] = 1;
		G[node2_neibor][node1] = 1;
		G[node2][node1_neibor] = 1;
		G[node1_neibor][node2] = 1;
	}
	return true;
}

bool swap_2(vector<node> &rnode, vector<vector<int>> &G, int &falsenum,double rate1)//新的两条边之间的最优交换概率之和要大于之前的的两条边之和
{
	int node1, node1_neibor, node2, node2_neibor;
	node1 = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
	node1_neibor = ((double)rand() / RAND_MAX)*(rnode[node1].neibor.size() - 1);
	node1_neibor = rnode[node1].neibor[node1_neibor];
	node2 = ((double)rand() / RAND_MAX)*(rnode.size() - 1);
	node2_neibor = ((double)rand() / RAND_MAX)*(rnode[node2].neibor.size() - 1);
	node2_neibor = rnode[node2].neibor[node2_neibor];
	if (find(rnode[node1].neibor.begin(), rnode[node1].neibor.end(), node2_neibor) != rnode[node1].neibor.end())
	{
		falsenum++;
		return false;
	}
	else if (find(rnode[node2].neibor.begin(), rnode[node2].neibor.end(), node1_neibor) != rnode[node2].neibor.end())
	{
		falsenum++;
		return false;
	}
	else if (node1 == node2 || node1_neibor == node2 || node1 == node2_neibor || node2 == node1_neibor)
	{
		falsenum++;
		return false;
	}
	else if ((rnode[node1].cp[node1_neibor] + rnode[node1_neibor].cp[node1] + rnode[node2].cp[node2_neibor] + rnode[node2_neibor].cp[node2])+rate1 > rnode[node1].cp[node2_neibor] + rnode[node2_neibor].cp[node1] + rnode[node2].cp[node1_neibor] + rnode[node1_neibor].cp[node2]) {
		falsenum++;
		return false;
		
	}
	else
	{
		rnode[node1].neibor.erase(find(rnode[node1].neibor.begin(), rnode[node1].neibor.end(), node1_neibor));
		rnode[node1_neibor].neibor.erase(find(rnode[node1_neibor].neibor.begin(), rnode[node1_neibor].neibor.end(), node1));
		rnode[node2].neibor.erase(find(rnode[node2].neibor.begin(), rnode[node2].neibor.end(), node2_neibor));
		rnode[node2_neibor].neibor.erase(find(rnode[node2_neibor].neibor.begin(), rnode[node2_neibor].neibor.end(), node2));
		rnode[node1].neibor.push_back(node2_neibor);
		rnode[node2_neibor].neibor.push_back(node1);
		rnode[node2].neibor.push_back(node1_neibor);
		rnode[node1_neibor].neibor.push_back(node2);
		G[node1][node1_neibor] = 0;
		G[node1_neibor][node1] = 0;
		G[node2][node2_neibor] = 0;
		G[node2_neibor][node2] = 0;
		G[node1][node2_neibor] = 1;
		G[node2_neibor][node1] = 1;
		G[node2][node1_neibor] = 1;
		G[node1_neibor][node2] = 1;
	}
	return true;
}