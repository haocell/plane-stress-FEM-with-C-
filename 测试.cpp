#include<iostream>
#include<math.h>  
#include<fstream>
#include<iomanip>
#include"matrix.h"  
using namespace std;

int main()
{	
	char fname[20];
	cout << "请输入数据文件名称: ";	//读取文件
	cin >> fname;
	ifstream infile;
	infile.open(fname, ios_base::in);
	if (!infile)
	{
		cout << fname << ": Can not open this file" << endl;
		return 0;
	}
	//节点数量，单元数量，x固定节点数量，y固定节点数量，x荷载数量，y荷载数量
	int numNode, numEle, xFixNum, yFixNum,xLoadNum, yLoadNum;     //定义辅助变量  
	double Emodule, poisson;	//弹性模量，泊松比
	infile >> numNode >> numEle >> xFixNum>> yFixNum>> xLoadNum>> yLoadNum;
	if (numNode > 700) {
		cout << "节点数量过多！";
		return 0;
	}
	if (numEle > 500) {
		cout << "单元数量过多！";
		return 0;
	}
	infile >> Emodule >> poisson;
	Matrix nodeCoor(numNode,2);//节点坐标
	Matrix nodes(numEle, 3);//单元节点编号
	for (int i = 0; i < numNode; i++) {//读取节点坐标
		infile >> nodeCoor.p[i][0];
		infile >> nodeCoor.p[i][1];
	}
	for (int i = 0; i < numEle; i++) {//读取单元
		infile >> nodes.p[i][0];
		infile >> nodes.p[i][1];
		infile >> nodes.p[i][2];
	}
	int nn = xFixNum + yFixNum;
	Matrix fixNum(nn,1);
	for (int i = 0; i < nn; i++) {//读取固定位移
		infile >> fixNum.p[i][0];
	}
	int nf = xLoadNum + yLoadNum;
	Matrix fload(nf, 2);
	for (int i = 0; i < nf; i++) {//读取荷载
		infile >> fload.p[i][0];
		infile >> fload.p[i][1];
	}
	/*输出坐标
	for (int i = 0; i < numNode;  i++) {
		cout<< nodeCoor.p[i][0]<<"  ";
		cout << nodeCoor.p[i][1] << endl;;
	}
	*/
	infile.close();
	
	////////////////////////////////////
	/// 输入控制参数
	///////////////////////////////////
	const int lx = 16;	//x段数
	const int ly = 8;	//y段数
	int sysDof = 2 * numNode;	//系统总自由度
	const int eleDof = 6;	//每个单元自由度
	////////////////////////////////////
	/// 矩阵和向量初始化
	///////////////////////////////////
	Matrix ff(sysDof,1);	//系统力向量
	double k[eleDof][eleDof] = { 0 };	//单元刚度矩阵
	Matrix  kk(sysDof,sysDof);	//系统总刚度矩阵
	//double disp[sysDof];	//系统位移
	double eleDisp[eleDof] = { 0 };	//单元位移
	Matrix stress(numEle,3);	//应力矩阵
	Matrix strain(numEle,3);	//应变矩阵
	int index[eleDof];	//索引向量
	Matrix D(3, 3);
	double B[3][6] = { 0 };
	//定义D矩阵
	double ev = Emodule / (1 - pow(poisson, 2));
	D.p[0][0] = D.p[1][1] = ev;
	D.p[0][1] = D.p[1][0] = ev * poisson;
	D.p[2][2] = ev * (1 - poisson) / 2;
	////////////////////////////////////
	/// 引入边界条件
	///////////////////////////////////
	Matrix bcDofNum(nn, 1);
	for (int i = 0; i < nn; i++) {
		if (i < xFixNum) {
			bcDofNum.p[i][0] = 2*(fixNum.p[i][0]-1);
		}
		else {
			bcDofNum.p[i][0] = 2 * (fixNum.p[i][0]-1)+1;
		}
	}
	////////////////////////////////////
	/// 计算单元刚度矩阵
	///////////////////////////////////
	int nd[3];
	double x[3] = { 0 }, y[3] = { 0 };
	for (int iel = 0; iel < numEle; iel++) {
		nd[0] = nodes.p[iel][0];
		nd[1] = nodes.p[iel][1];
		nd[2] = nodes.p[iel][2];
		x[0] = nodeCoor.p[nd[0] - 1][0];
		x[1] = nodeCoor.p[nd[1] - 1][0];
		x[2] = nodeCoor.p[nd[2] - 1][0];
		y[0] = nodeCoor.p[nd[0] - 1][1];
		y[1] = nodeCoor.p[nd[1] - 1][1];
		y[2] = nodeCoor.p[nd[2] - 1][1];
		//节点位移对应矩阵编号
		int n = 0;
		for (int i = 0; i < 3; i++) {
			int start = (nd[i] - 1) * 2;
			for (int j = 0; j < 2; j++) {
				index[n] = start + j;
				n = n + 1;
			}
		}

		double area = 0.5 * (x[0] * y[1] + x[1] * y[2] + x[2] * y[0] - x[0] * y[2] - x[1] * y[0] - x[2] * y[1]);
		double bi[3] = { (y[1] - y[2]),(y[2] - y[0]),(y[0] - y[1]) };
		double ci[3] = { (x[2] - x[1]),(x[0] - x[2]),(x[1] - x[0]) };
		for (int i = 0; i < 3; i++) {
			int i1 = 2 * i;
			int i2 = i1 + 1;
			B[0][i1] = bi[i];
			B[1][i2] = ci[i];
			B[2][i1] = ci[i];
			B[2][i2] = bi[i];
		}

		Matrix Bi(3, 6);
		Matrix BiT(6, 3);
		Matrix Di(3, 3);
		Bi = (*B);
		Di = D;
		BiT = BiT.T(Bi);	//Bi转置
		Di *= Bi;
		BiT *= Di;
		Matrix k(eleDof, eleDof);
		k = BiT;
		//cout << iel << endl;
		//k.Show();//输出单元刚度矩阵
		k *= (1 / (4 * area));

		//刚度组装
		int ii, jj;
		for (int i = 0; i < 6; i++) {
			ii = index[i];
			for (int j = 0; j < 6; j++) {
				jj = index[j];
				kk.p[ii][jj] = kk.p[ii][jj] + k.p[i][j];
			}
		}
	}
	//引入荷载阵
	for (int i = 0; i < nf; i++) {
		if (i < xLoadNum) {
			int j = 2 * (fload.p[i][0]-1);
			ff.p[j][0] = fload.p[i][1];
		}
		else {
			int j = 2 * (fload.p[i][0] - 1) + 1;
			ff.p[j][0]= fload.p[i][1];
		}
	}
	//刚度矩阵置一法
	//int nn = 18;	//固定边界自由度
	for (int i = 0; i < nn; i++) {
		int c = bcDofNum.p[i][0];
		for (int j = 0; j < sysDof; j++) {
			kk.p[c][j] = 0;
		}
		kk.p[c][c] = 1;
		ff.p[c][0] = 0;
	}
	Matrix kk1(sysDof, sysDof);
	Matrix ff1(sysDof, 1);
	kk1 = kk;
	ff1 = ff;

	////////////////////////////////////
	///计算求解矩阵
	///////////////////////////////////
	//节点位移求解
	Matrix disp(sysDof, 1);
	disp.Solve(kk1, ff1);
	//单元应力计算
	Matrix stress1(numEle, 3);	//应力矩阵
	Matrix strain1(numEle, 3);	//应变矩阵
	int nd1[3];
	double x1[3] = { 0 }, y1[3] = { 0 };
	for (int ielp = 0; ielp < numEle; ielp++) {
		nd1[0] = nodes.p[ielp][0];
		nd1[1] = nodes.p[ielp][1];
		nd1[2] = nodes.p[ielp][2];
		x1[0] = nodeCoor.p[nd1[0] - 1][0];
		x1[1] = nodeCoor.p[nd1[1] - 1][0];
		x1[2] = nodeCoor.p[nd1[2] - 1][0];
		y1[0] = nodeCoor.p[nd1[0] - 1][1];
		y1[1] = nodeCoor.p[nd1[1] - 1][1];
		y1[2] = nodeCoor.p[nd1[2] - 1][1];
		//节点位移对应编号
		int n = 0;
		//cout << "单元" << ielp;
		for (int i = 0; i < 3; i++) {
			int start = (nd1[i] - 1) * 2;
			for (int j = 0; j < 2; j++) {
				index[n] = start + j;
				//cout << "  " << index[n] + 1;
				n = n + 1;
			}
		}
		//cout << endl;

		double area = 0.5 * (x1[0] * y1[1] + x1[1] * y1[2] + x1[2] * y1[0] - x1[0] * y1[2] - x1[1] * y1[0] - x1[2] * y1[1]);
		double bi[3] = { (y1[1] - y1[2]),(y1[2] - y1[0]),(y1[0] - y1[1]) };
		double ci[3] = { (x1[2] - x1[1]),(x1[0] - x1[2]),(x1[1] - x1[0]) };
		for (int i = 0; i < 3; i++) {
			int i1 = 2 * i;
			int i2 = i1 + 1;
			B[0][i1] = bi[i];
			B[1][i2] = ci[i];
			B[2][i1] = ci[i];
			B[2][i2] = bi[i];
		}
		//单元节点位移
		Matrix eleDisp(6, 1);

		for (int i = 0; i < 6; i++) {
			eleDisp.p[i][0] = disp.p[index[i]][0];
		}

		//cout << endl;
		Matrix Bi(3, 6);
		Matrix BiT(6, 3);
		Matrix Di(3, 3);
		Bi = (*B);
		Bi *= (1 / (2 * area));

		Matrix estrain(3, 1); //单元应变
		Matrix estress(3, 1);//单元应力
		estrain = Bi;
		estrain *= eleDisp;
		estress = D;
		estress *= estrain;
		//将单元应力汇总
		for (int i = 0; i < 3; i++) {
			strain1.p[ielp][i] = estrain.p[i][0];
			stress1.p[ielp][i] = estress.p[i][0];
		}


	}
	////////////////////////////////////
	///数据输出
	///////////////////////////////////


	ofstream file;
	file.open("04.text", ios::out);
	//file << setiosflags(ios::fixed);	//设置输出位数
	file << setiosflags(ios::right);	//设置右对齐
	file << "x,y,u,v" << endl;
	for (int i = 0; i < numNode; i++) {
		file << setprecision(2) << setw(6) << nodeCoor.p[i][0] << setw(8) << nodeCoor.p[i][1];
		file << setw(10) << disp.p[2 * i][0] << setw(10) << disp.p[2 * i + 1][0] << endl;
	}
	//输出每个单元应力
	file << "每个单元sigmaX，sigmaY，sigmaXY" << endl;
	for (int i = 0; i < numEle; i++) {
		file << "单元" << i + 1 << setw(10) << stress1.p[i][0] << setw(10) << stress1.p[i][1] << setw(10) << stress1.p[i][2] << endl;
	}
	/*适用于该题的输出
	//按矩形图像分别输出节点位移
	file << endl << "位移u" << endl;
	for (int i = ly + 1; i > 0; i--) {
		for (int j = 0; j < lx + 1; j++) {
			file << setw(10) << disp.p[2 * ((ly + 1) * j + i) - 2][0];
		}
		file << endl;
	}
	file << endl << "位移v" << endl;
	for (int i = ly + 1; i > 0; i--) {
		for (int j = 0; j < lx + 1; j++) {
			file << setw(10) << disp.p[2 * ((ly + 1) * j + i) - 1][0];
		}
		file << endl;
	}
	//按矩形图像分别输出单元应力
	file << endl << "sigmaX" << endl;
	for (int i = 2 * ly; i > 0; i--) {
		for (int j = 0; j < lx; j++) {
			file << setw(10) << stress1.p[2 * ly * j + i - 1][0];
		}
		file << endl;
	}
	file << endl << "igmaY" << endl;;
	for (int i = 2 * ly; i > 0; i--) {
		for (int j = 0; j < lx; j++) {
			file << setw(10) << stress1.p[2 * ly * j + i - 1][1];
		}
		file << endl;
	}
	file << "sigmaXY" << endl;
	for (int i = 2 * ly; i > 0; i--) {
		for (int j = 0; j < lx; j++) {
			file << setw(10) << stress1.p[2 * ly * j + i - 1][2];
		}
		file << endl;
	}
	*/
	file.close();	//关闭文件
	cout << "程序运行成功，打开文件\"01.text\"查看结果";
	return 0;
	
}

/*
int feelDof(int* nd) {
	int n = 0;
	int index[6];
	for (int i = 0; i < 3; i++) {
		int start = (nd[i] - 1) * 2;
		for (int j = 0; j < 2; j++) {
			index[n] = start + j;
			n = n + 1;
		}
	}
	return *index;
}
*/
