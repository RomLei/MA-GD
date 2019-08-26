#pragma once
#include"vector"
using namespace std;
#define g_rate 2;
#define learning_rate 1.2;
#define gradient 0.01
#define belta 0

class node
{
public:
	int degree;//�洢��
	int activation;//�����
	int search_flag;//�Ƿ����
	int checkdegree;//����
	int max_cluster_flag;//�Ƿ���������Ӵ�
	int attacked_rank;//������˳��
	vector<double> cp;//�����Ӹ���
	vector<double> delta_cp_current;//��ǰ����ֵcp
	vector<double> delta_cp_before;//��������ֵcp
	vector<int> prefer_edge_rank;//���ӱ߸������У��Ӹߵ���
	vector<int> neibor;
	node()
	{
		activation = 1;
		search_flag = 1;
		max_cluster_flag = 1;
	}
};
double mainloop();//Դ.cpp

void initilizationBA(vector<node> &rnode, vector<vector<int>> &G);//initialization.cpp
void initialize_cp(vector<node> &rnode);
void initialization_SF_r(vector<node> &rnode);
int calculate_degree2(vector<node> &rnode);
void translate_rnode_G(vector<node> &rnode, vector<vector<int>> &G);

int update_network(vector<node> &rnode, vector<vector<int>> &G);//swap.cpp
void sort_part(node &rnode);
void update_perfer_edge_rank(vector<node> &rnode);
int find_rank(vector<node> &rnode, int node1, int node2);
int find_change_neibor(vector<node> &rnode, int node1, int jnumber, int start_num);
void quick_sort(vector<vector<double>> &ranklist, int start, int end);
bool swap_edge(vector<node> &rnode, vector<vector<int>> &G, int &falsenum);
bool swap_2(vector<node> &rnode, vector<vector<int>> &G, int &falsenum, double rate1);

double compute(vector<node> &rnode, vector<vector<int>> &G);//computeR.cpp
int MaxDegreePosition(vector<node> &rnode);
void calculate_degree(vector<node> &rnode, vector<vector<int>> &G);
void DeleteNode(vector<node> &rnode, int p, vector<vector<int>> &G);
int MaxConnectedCluster(vector<node> &rnode, vector<vector<int>> &G);
bool whether_search_end(vector<int> &search, vector<node> &rnode);
int Include_Nodes(int p, vector<vector<int>> &G, vector<int> &search, vector<node> &rnode);
void set_max_cluster_flag_0(vector<node> &rnode, int pattern);
void recovernode(vector<node> &rnode);
void translate_G_to_neibor(vector<node> &rnode, vector<vector<int>> &G);
int MaxConnectedCluster_nosignal(vector<node> &rnode, vector<vector<int>> &G);
int Include_Nodes_nosingal(int p, vector<vector<int>> &G, vector<int> &search, vector<node> &rnode);

void calculateGD(vector<node> &rnode, vector<vector<int>> &G);//calculateGD.cpp
void decrease_other_cp(node &rnode, int j, double d_value);
void update_cp(vector<node> &rnode);
void calculate_cp(vector<node> &rnode, vector<vector<int>> &G);

void print_network(vector<node> &rnode, int nettype);//printnetwork.cpp