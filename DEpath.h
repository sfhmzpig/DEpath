//
// Created by MBA on 18/1/21.
//

#ifndef DEPATH_DEPATH_H
#define DEPATH_DEPATH_H

#include <vector>
#include <array>
#include <string>
#include <memory>
#include <fstream>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <ctime>

const double precession = 0.0000000001;

int X = 600; //长
int Y = 600; //宽
int Z = 210; //高
int L = 1; //每个网格边长
const int Pop = 256; //种群个数
const int NH = 1; //特殊种群个数
const double F = 0.7; //变异力
const double Fmin = 0.3;
const double Fmax = 0.9;
double CRl = 0.1;
double CRu = 0.6;
double CR = 0.6; //交叉时选择变异体的概率
const int maxIte = 20000; //最大迭代次数
const int reIni = 3000; //判断种群最好个体是否搜寻到更好的结果的迭代次数
const int mut1 = 60000;
const int minNum = 99999; //最小数常数
int w1 = 10;
const double w2 = 0.1;
int no = 2;
const int T = 5;

std::array<int,3> start_point = {13,502,149}; //起点
std::array<int,3> end_point = {503,108,200}; //终点
int threat_point[T][3] = {{99,307,47},{246,140,60},{287,246,59},{325,357,67}};
int D; //个体染色体个数
const double PI = 3.14159265358979323846;
/*
using XYZ = std::array<std::array<std::array<int,Z+1>,Y+1>,X+1>; //三维数组，模拟环境上每个网格处是否有障碍物
std::shared_ptr<XYZ> mapInfo; //环境太大，栈空间不足时用堆空间储存
XYZ mapInf; //栈空间储存环境信息
 */
int *** map_info;

double lineDis;

std::vector<std::array<int,3>> W; //采用几何方法计算出得一个特殊个体
std::vector<std::array<int,2>> P; //更新W时为更新函数提供面（plane），二维数组，第一维是切面（plane），第二维是切面的障碍物个数

int *** popu_now; //种群信息
int *** popu_next; //种群信息（变异后，由popu_now计算）
int *** popu_cross; //种群信息（交叉后）
double fit_now,fit_next; //暂时无用
int signal1;
const double misPre = 0.000001;

void MakeMap(); //构建地图信息函数

std::vector<std::array<int,3>> IniW(std::vector<std::array<int,3>>); //初始化W

std::vector<std::array<int,2>> SepPlane(std::vector<std::array<int,2>> P,std::array<int,3> sp,std::array<int,3> ep); //构建P函数

std::vector<std::array<int,3>> SelectPAN(std::vector<std::array<int,3>> W,int plane); //选择前置和后置点函数

std::array<int,3> SelectPoint(std::array<int,3> pre, std::array<int,3> next); //选择被切面上的点

std::vector<std::array<int,3>> ReW(std::vector<std::array<int,3>> W,std::array<int,3> addP); //更新W

double CalFes(std::array<int,3> CanP); //计算适应度

double distance(double x1,double y1,double z1,double x2,double y2,double z2); //计算两点间距离

double Collosion(int **); //计算个体的碰撞度，取值［0，1］

double NorDis(int **); //计算个体的距离度，取值［0，1］

double Throught(int **);

std::array<double,Pop> PopFits(int ***,int); //计算种群每个个体适应度

//std::array<int,3> Mut(int ***popu_now,std::array<int,3> pop,std::array<std::array<int,3>,D>);

double PopFit(int **pop_now,int); //计算一个个体适应度
#endif //DEPATH_DEPATH_H
