#include<iostream>
#include<math.h>  
#include<fstream>
#include<iomanip>
#include"matrix.h"  
using namespace std;

int main()
{	
	char fname[20];
	cout << "�����������ļ�����: ";	//��ȡ�ļ�
	cin >> fname;
	ifstream infile;
	infile.open(fname, ios_base::in);
	if (!infile)
	{
		cout << fname << ": Can not open this file" << endl;
		return 0;
	}
	//�ڵ���������Ԫ������x�̶��ڵ�������y�̶��ڵ�������x����������y��������
	int numNode, numEle, xFixNum, yFixNum,xLoadNum, yLoadNum;     //���帨������  
	double Emodule, poisson;	//����ģ�������ɱ�
	infile >> numNode >> numEle >> xFixNum>> yFixNum>> xLoadNum>> yLoadNum;
	if (numNode > 700) {
		cout << "�ڵ��������࣡";
		return 0;
	}
	if (numEle > 500) {
		cout << "��Ԫ�������࣡";
		return 0;
	}
	infile >> Emodule >> poisson;
	Matrix nodeCoor(numNode,2);//�ڵ�����
	Matrix nodes(numEle, 3);//��Ԫ�ڵ���
	for (int i = 0; i < numNode; i++) {//��ȡ�ڵ�����
		infile >> nodeCoor.p[i][0];
		infile >> nodeCoor.p[i][1];
	}
	for (int i = 0; i < numEle; i++) {//��ȡ��Ԫ
		infile >> nodes.p[i][0];
		infile >> nodes.p[i][1];
		infile >> nodes.p[i][2];
	}
	int nn = xFixNum + yFixNum;
	Matrix fixNum(nn,1);
	for (int i = 0; i < nn; i++) {//��ȡ�̶�λ��
		infile >> fixNum.p[i][0];
	}
	int nf = xLoadNum + yLoadNum;
	Matrix fload(nf, 2);
	for (int i = 0; i < nf; i++) {//��ȡ����
		infile >> fload.p[i][0];
		infile >> fload.p[i][1];
	}
	/*�������
	for (int i = 0; i < numNode;  i++) {
		cout<< nodeCoor.p[i][0]<<"  ";
		cout << nodeCoor.p[i][1] << endl;;
	}
	*/
	infile.close();
	
	////////////////////////////////////
	/// ������Ʋ���
	///////////////////////////////////
	const int lx = 16;	//x����
	const int ly = 8;	//y����
	int sysDof = 2 * numNode;	//ϵͳ�����ɶ�
	const int eleDof = 6;	//ÿ����Ԫ���ɶ�
	////////////////////////////////////
	/// �����������ʼ��
	///////////////////////////////////
	Matrix ff(sysDof,1);	//ϵͳ������
	double k[eleDof][eleDof] = { 0 };	//��Ԫ�նȾ���
	Matrix  kk(sysDof,sysDof);	//ϵͳ�ܸնȾ���
	//double disp[sysDof];	//ϵͳλ��
	double eleDisp[eleDof] = { 0 };	//��Ԫλ��
	Matrix stress(numEle,3);	//Ӧ������
	Matrix strain(numEle,3);	//Ӧ�����
	int index[eleDof];	//��������
	Matrix D(3, 3);
	double B[3][6] = { 0 };
	//����D����
	double ev = Emodule / (1 - pow(poisson, 2));
	D.p[0][0] = D.p[1][1] = ev;
	D.p[0][1] = D.p[1][0] = ev * poisson;
	D.p[2][2] = ev * (1 - poisson) / 2;
	////////////////////////////////////
	/// ����߽�����
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
	/// ���㵥Ԫ�նȾ���
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
		//�ڵ�λ�ƶ�Ӧ������
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
		BiT = BiT.T(Bi);	//Biת��
		Di *= Bi;
		BiT *= Di;
		Matrix k(eleDof, eleDof);
		k = BiT;
		//cout << iel << endl;
		//k.Show();//�����Ԫ�նȾ���
		k *= (1 / (4 * area));

		//�ն���װ
		int ii, jj;
		for (int i = 0; i < 6; i++) {
			ii = index[i];
			for (int j = 0; j < 6; j++) {
				jj = index[j];
				kk.p[ii][jj] = kk.p[ii][jj] + k.p[i][j];
			}
		}
	}
	//���������
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
	//�նȾ�����һ��
	//int nn = 18;	//�̶��߽����ɶ�
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
	///����������
	///////////////////////////////////
	//�ڵ�λ�����
	Matrix disp(sysDof, 1);
	disp.Solve(kk1, ff1);
	//��ԪӦ������
	Matrix stress1(numEle, 3);	//Ӧ������
	Matrix strain1(numEle, 3);	//Ӧ�����
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
		//�ڵ�λ�ƶ�Ӧ���
		int n = 0;
		//cout << "��Ԫ" << ielp;
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
		//��Ԫ�ڵ�λ��
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

		Matrix estrain(3, 1); //��ԪӦ��
		Matrix estress(3, 1);//��ԪӦ��
		estrain = Bi;
		estrain *= eleDisp;
		estress = D;
		estress *= estrain;
		//����ԪӦ������
		for (int i = 0; i < 3; i++) {
			strain1.p[ielp][i] = estrain.p[i][0];
			stress1.p[ielp][i] = estress.p[i][0];
		}


	}
	////////////////////////////////////
	///�������
	///////////////////////////////////


	ofstream file;
	file.open("04.text", ios::out);
	//file << setiosflags(ios::fixed);	//�������λ��
	file << setiosflags(ios::right);	//�����Ҷ���
	file << "x,y,u,v" << endl;
	for (int i = 0; i < numNode; i++) {
		file << setprecision(2) << setw(6) << nodeCoor.p[i][0] << setw(8) << nodeCoor.p[i][1];
		file << setw(10) << disp.p[2 * i][0] << setw(10) << disp.p[2 * i + 1][0] << endl;
	}
	//���ÿ����ԪӦ��
	file << "ÿ����ԪsigmaX��sigmaY��sigmaXY" << endl;
	for (int i = 0; i < numEle; i++) {
		file << "��Ԫ" << i + 1 << setw(10) << stress1.p[i][0] << setw(10) << stress1.p[i][1] << setw(10) << stress1.p[i][2] << endl;
	}
	/*�����ڸ�������
	//������ͼ��ֱ�����ڵ�λ��
	file << endl << "λ��u" << endl;
	for (int i = ly + 1; i > 0; i--) {
		for (int j = 0; j < lx + 1; j++) {
			file << setw(10) << disp.p[2 * ((ly + 1) * j + i) - 2][0];
		}
		file << endl;
	}
	file << endl << "λ��v" << endl;
	for (int i = ly + 1; i > 0; i--) {
		for (int j = 0; j < lx + 1; j++) {
			file << setw(10) << disp.p[2 * ((ly + 1) * j + i) - 1][0];
		}
		file << endl;
	}
	//������ͼ��ֱ������ԪӦ��
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
	file.close();	//�ر��ļ�
	cout << "�������гɹ������ļ�\"01.text\"�鿴���";
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
