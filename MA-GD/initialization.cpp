#include"vector"
#include"GD.h"

using namespace std;

void initilizationBA(vector<node> &rnode, vector<vector<int>> &G)
{
	int degree_ALL = 0;
	double pp = 0, ppi = 0;
	G[0][1] = 1; G[1][0] = 1;
	G[2][0] = 1; G[0][2] = 1;
	G[2][1] = 1; G[1][2] = 1;
	rnode[0].degree = 2;
	rnode[1].degree = 2;
	rnode[2].degree = 2;
	for (int i = 3, k = 0, j = 0; i < rnode.size(); i++)
	{
		k = 0; j = 0;
		ppi = 0;
		degree_ALL = 0;
		for (int u = 0; u < i; u++)
		{
			degree_ALL += rnode[u].degree;
		}
		pp = (double)rand() / RAND_MAX;
		while (rnode[i].degree < 1)
		{
			if ((pp - ppi) <= (double)rnode[k].degree / degree_ALL)
			{
				G[i][k] = 1; G[k][i] = 1;
				rnode[i].degree++;
				rnode[k].degree++;
			}
			else
			{
				ppi += (double)rnode[k].degree / degree_ALL;
				k++;
			}
		}
		degree_ALL = 0;
		for (int u = 0; u < i; u++)
		{
			if (u != k)
			{
				degree_ALL += rnode[u].degree;
			}
		}
		pp = (double)rand() / RAND_MAX;
		ppi = 0;
		if (k == 0)
		{
			j++;
		}
		while (rnode[i].degree < 2)
		{
			if ((pp - ppi) <= (double)rnode[j].degree / degree_ALL)
			{
				if (G[i][j] == 1)
				{
					pp = (double)rand() / RAND_MAX;
					ppi = 0;
					j = 0;
				}
				else
				{
					G[i][j] = 1; G[j][i] = 1;
					rnode[i].degree++;
					rnode[j].degree++;
				}
			}
			else
			{
				ppi += (double)rnode[j].degree / degree_ALL;
				j++;
				if (j == k)
				{
					j++;
				}
			}
		}
	}
}

void initialization_SF_r(vector<node> &rnode)//初始化完全随机SF网络
{
	int degree_ALL = 0, degree_MID = 0, addNodeNum;
	double rn, rn2;
	//创建初始的三个点
	int n1, n2, n3;
	for (int i = 0; i < rnode.size(); i++) {
		rnode[i].degree = 0;
	}
	n1 = ((double)rand() / (RAND_MAX + 1))*rnode.size();
	do {
		n2 = ((double)rand() / (RAND_MAX + 1))*rnode.size();
	} while (n1 == n2);
	do {
		n3 = ((double)rand() / (RAND_MAX + 1))*rnode.size();
	} while (n1 == n3 || n2 == n3);
	rnode[n1].neibor.push_back(n2);
	rnode[n1].neibor.push_back(n3);
	rnode[n2].neibor.push_back(n1);
	rnode[n2].neibor.push_back(n3);
	rnode[n3].neibor.push_back(n1);
	rnode[n3].neibor.push_back(n2);
	rnode[n1].degree += 2;
	rnode[n2].degree += 2;
	rnode[n3].degree += 2;
	rnode[n1].activation = 1;
	rnode[n2].activation = 1;
	rnode[n3].activation = 1;
	vector<int> a(rnode.size() - 3, -1);
	for (int i = 0, j = 0; i < rnode.size() - 3; j++, i++) {
		if (j != n1 && j != n2 && j != n3) {
			a[i] = j;
		}
		else {
			i--;
		}
	}
	for (int i = 0; i <rnode.size() - 3; i++) {
		swap(a[i], a[((double)rand() / (RAND_MAX + 1))*(rnode.size() - 3)]);
	}
	//加入其他点
	for (int i = 0; i < rnode.size() - 3; i++)
	{
		degree_ALL = calculate_degree2(rnode);
		rn = (double)rand() / RAND_MAX;
		degree_MID = 0;
		for (int j = 0; j < rnode.size(); j++)
		{
			if (a[i] != j && rnode[j].activation == 1 && (find(rnode[j].neibor.begin(), rnode[j].neibor.end(), a[i]) == rnode[j].neibor.end())) {
				degree_MID += rnode[j].neibor.size();
				if ((rn <= ((double)degree_MID / degree_ALL)))
				{
					rnode[a[i]].neibor.push_back(j);
					rnode[j].neibor.push_back(a[i]);
					rnode[a[i]].degree++;
					rnode[j].degree++;
					rn = (double)rand() / RAND_MAX;
					degree_ALL -= rnode[j].degree + 1;
					j = -1;
					degree_MID = 0;
				}
				if (rnode[a[i]].degree == 2)
				{
					rnode[a[i]].activation = 1;
					break;
				}
			}
		}
	}
}


void translate_rnode_G(vector<node> &rnode, vector<vector<int>> &G)
{
	for (int i = 0; i < rnode.size(); i++) {
		for (int j = 0; j < rnode[i].neibor.size(); j++) {
			G[i][rnode[i].neibor[j]] = 1;
			G[rnode[i].neibor[j]][i] = 1;
		}
	}
}

int calculate_degree2(vector<node> &rnode)//计算总的度
{
	int degree_ALL = 0;
	for (int i = 0; i < rnode.size(); i++)
	{
		degree_ALL += rnode[i].degree;
	}
	return degree_ALL;
}

void initialize_cp(vector<node> &rnode)
{
	for (int i = 0; i < rnode.size(); i++) {
		rnode[i].cp.resize(rnode.size(), 0);
		for (int j = 0; j < rnode[i].neibor.size(); j++) {
			rnode[i].cp[rnode[i].neibor[j]] = 1.0;
		}
		rnode[i].delta_cp_before.resize(rnode.size(), 0);
		rnode[i].delta_cp_current.resize(rnode.size(), 0);
		rnode[i].prefer_edge_rank.resize(rnode.size(), 0);
	}
}