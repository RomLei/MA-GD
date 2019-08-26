#include"vector"
#include<iostream>
#include<fstream>
#include<string>
#include"GD.h"

using namespace std;
int cycnum = 1;
int main()
{
	ofstream result;
	result.open("result.txt");
	double R = 0, Rmax = 0, Rmin = 100, Rmid = 0;
	double k = 1;
	for (int i = 0; i < k; i++) {
		Rmid = mainloop();
		if (Rmid > Rmax) {
			Rmax = Rmid;
		}
		if (Rmid < Rmin)
		{
			Rmin = Rmid;
		}
		R += Rmid;
	}
	R = R / k;
	result << "the final R is" << R << endl;
	result << "the max R is " << Rmax << "; the min R is " << Rmin << endl;
	result.close();
	return 0;
}

double mainloop()
{
	srand(cycnum);
	cycnum++;
	double R = 0, Rmax = 0, R_before = 0;
	int gen = 20000, falsenum = 0;
	int countnum = 0;
	//
	//int N = 10;
	ifstream translationresult("translationresult.txt");//读取真实网络文件
	char realworldgraph;
	int jack = 0;
	int N = 0;
	realworldgraph = translationresult.get();
	while (realworldgraph != 10)
	{
		N++;
		realworldgraph = translationresult.get();
		printf("reading data %d\n", jack); jack++;
	}
	translationresult.seekg(0);
	string realworldline;
	vector<node> rnode(N);
	vector<node> rnodeb(N);
	vector<node> rnodeb2(N);
	vector<node> rnodemax(N);
	vector<vector<int>> G(N, vector<int>(N, 0));
	vector<vector<int>> Gb(N, vector<int>(N, 0));
	vector<vector<int>> Gb2(N, vector<int>(N, 0));
	vector<vector<int>> Gmax(N, vector<int>(N, 0));
	//创建真实网络
	for (int i = 0, abama = 0, abama2 = 0; i < N; i++)
	{
		abama2 = 0;
		getline(translationresult, realworldline);
		for (int j = 0; j < N; j++)
		{
			G[i][j] = realworldline[abama2] - 48;
			printf("creating network %d\n", abama); abama++;
			abama2++;
		}
		abama++;
	}
	translate_G_to_neibor(rnode, G);


	//initialization_SF_r(rnode);
	//translate_rnode_G(rnode, G);
	//calculate_degree(rnode,G);

	//初始化网络
//	initilizationBA(rnode, G);
	initialize_cp(rnode);
	//预训练
	print_network(rnode, 0);
	//开始循环
	int succ = false;
	double rate1 = 0.5;
	for (int i = gen; i > 0; i--) {

		rnodeb = rnode;
		Gb = G;

		

		if ((rand() / (RAND_MAX + 1.0))<exp(-i * 4.0 / gen)) {
		calculateGD(rnode, G);
		//	update_perfer_edge_rank(rnode);
		while (!swap_2(rnode, G, falsenum,rate1));
		rate1 = rate1 + 0.0007;
		//	succ = update_network(rnode, G);
		}
		else {
			while (!swap_edge(rnode, G, falsenum));
		}

		//

		R = compute(rnode, G);


		if (R > Rmax)
		{
			rnodemax = rnode;
			Gmax = G;
			Rmax = R;
		}
		else {
			rnode = rnodeb;
			G = Gb;
		}
		//对比R
		//succ = 0;

		R_before = R;


		printf("gen:%d, cycnum:%d , R: %f, Rmax: %f，falsenum:%d rnode[26].cp[134]:%f, status'10-5':%d, suc:%d\n", i, cycnum, R, Rmax, falsenum, rnode[26].cp[134], G[26][134], succ);
	}
	print_network(rnode, 1);
	return Rmax;
}