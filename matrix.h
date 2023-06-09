
#ifndef __MATRIX_CLL_H__
#define __MATRIX_CCL_H__
#include "pch.h"

class Matrix {
public:
	int rows_num, cols_num;
	double** p;
	void initialize();//��ʼ������

public:
	Matrix(int, int);
	Matrix(int, int, double);//Ԥ��ֿռ�
	virtual ~Matrix();//��������Ӧ�����麯�������Ǵ��಻��������
	Matrix& operator=(const Matrix&);//����ĸ���
	Matrix& operator=(double*);//�������ֵ��������
	Matrix& operator+=(const Matrix&);//�����+=����
	Matrix& operator-=(const Matrix&);//-=
	Matrix& operator*=(const Matrix&);//*=
	Matrix& operator*=(double);//*=
	Matrix operator*(const Matrix& m)const;
	Matrix& Solve(const Matrix&, const Matrix&);//������Է�����Ax=b
	void Show() const;//������ʾ
	void swapRows(int, int);
	double det();//����������ʽ
	double Point(int i, int j) const;
	static Matrix inv(Matrix);//�����������
	static Matrix eye(int);//����һ����λ����
	int row() const;
	int col() const;
	Matrix& T(const Matrix& m);//����ת�õ�ʵ��
	Matrix gaussianEliminate();//��˹��Ԫ��
	friend std::istream& operator>>(std::istream&, Matrix&);//ʵ�־��������
};

#endif