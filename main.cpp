#include <iostream>
#include "DEpath.h"

void MakeMap(){ //构建地图，mapinf三维数组储存地图信息，0为可行，1为障碍
    std::ifstream inZ{"mapZ.txt"}; //储存 xy面的z值（高度）信息，储存入mapZ，方便导入matlab
    std::istream_iterator<double > iiZ{inZ};
    std::array<std::array<double,X>,Y> high;
    for (auto &i : high){
        for (auto &j : i){
            ++iiZ;
            ++iiZ;
            j = *iiZ++;
        }
    }
    for (int i = 0; i != X+1; ++i){
        for (int j = 0; j != Y+1; ++j){
            for (int k = 0; k != Z+1; ++k) {
                if (high[j][i] - k >= precession){
                    mapInf[i][j][k] = 1;
                }
            }
        }
    }
}

void MakeMap1(){
    double tmpDis = 0;
    map_info = new int **[X+1];
    for (int i = 0; i!= X+1; ++i){
        map_info[i] = new int *[Y+1];
        for (int j = 0; j!= Y+1; ++j){
            map_info[i][j] = new int[Z+1];
        }
    }
    for (int i = 0; i != T; ++i){
        for (int j = 0; j != X+1; ++j){
            for (int k = 0; k != Y+1; ++k){
                tmpDis = distance((double)threat_point[i][0],(double)threat_point[i][1],0,(double)j,(double)k,0);
                if (tmpDis - threat_point[i][2] <= precession){
                    for (int ii = 0; ii != Z+1; ++ii){
                        map_info[j][k][ii] = 1;
                    }
                }
            }
        }
    }
}

void RandTIme(){
    clock_t time_rand = clock();
    srand((unsigned int)time_rand);
}

std::vector<std::array<int,3>> IniW(std::vector<std::array<int,3>> W,std::array<int,3> sp,std::array<int,3> ep){ //初始化W
    W.push_back(sp);
    W.push_back(ep);
    return W;
}

bool compare1(std::array<int,2> a,std::array<int,2> b){ //用于P的排序功能函数
    return a[1] > b[1];
}
std::vector<std::array<int,2>> SepPlane(std::vector<std::array<int,2>> P,std::array<int,3> sp,std::array<int,3> ep){ //创建P，储存起点到终点之间每个面的障碍物点数目
    int planeNum = ep[0] - sp[0] - 1; //起点和终点间面的个数
    std::array<int,2> tmp;
    int inFeaPoint; //面的障碍物个数

    for (int i = 0; i != planeNum; ++i){
        tmp[0] = i + sp[0] + 1; //P的第一维储存是第几个面
        inFeaPoint = 0;
        for (int j = 0; j != Y+1; ++j){
            for (int k = 0; k != Z+1; ++k){
                //if (mapInf[tmp[0]][j][k] == 1)
                if (map_info[tmp[0]][j][k])
                    ++inFeaPoint;
            }
        }
        tmp[1] = inFeaPoint; //P的第二维储存面有多少个障碍物
        P.push_back(tmp);
    }

    std::sort(P.begin(),P.end(),compare1);

    return P;
}



std::vector<std::array<int,3>> SelectPAN(std::vector<std::array<int,3>> W,int plane){ //输入哪个面，和W，来为面选择前置点和后置点
    std::vector<std::array<int,3>> selectPoint; //第一个为后置点，第二个为前置点
    std::array<int,3> tmp;

    for (int i = 0; i != W.size(); ++i){
        if (plane < W[i][0]){ //第一个大于面的W，这个W的点为后置点，它的前一个点为前置点
            tmp = W[i];
            selectPoint.push_back(tmp);
            tmp = W[i-1];
            selectPoint.push_back(tmp);
            break;
            /*
            if (plane == W[i-1][0]){
                tmp = W[i-2];
                selectPoint.push_back(tmp);
            } else {
                tmp = W[i-1];
                selectPoint.push_back(tmp);
            }
             */
        }
    }
//    if (plane == 45)
//        std::cout << "pAN ok" <<std::endl;
    return selectPoint;
};

std::array<int,3> ThrPoint(std::array<int,3> pre, std::array<int,3> next, int plane){ //计算两个点的连线穿越的面的交点
    std::array<int,3> thrPoint;

    thrPoint[0] = plane;
    thrPoint[1] = (int)((plane - pre[0]) * (next[1] - pre[1]) / (next[0] - pre[0]) + pre[1]);
    thrPoint[2] = (int)((plane - pre[0]) * (next[2] - pre[2]) / (next[0] - pre[0]) + pre[2]);
//    if (plane == 45)
//        std::cout << "tp ok" <<std::endl;
    return thrPoint;
};

std::array<int,3> SelecPoint(std::array<int,3> thrP){ //通过交点选择新的点
    std::array<int,3> selecP;
    std::array<int,3> tmp_thrP;
    double minFes = minNum,tmpFes;
    int fesNum = 0,radius = L; //判断候选点中是否有可行点，如果没有radius增加一个L

    while (fesNum == 0){
        for (int i = thrP[1]-radius; i != thrP[1]+radius+1; ++i){
            if (i >= 0 && i <= Y){
                for (int j = thrP[2]-radius; j != thrP[2]+radius+1; ++j){
                    if (j >= 0 && j <= Z){
                        if (map_info[thrP[0]][i][j] == 0){
                            ++fesNum;
                            tmp_thrP[0] = thrP[0];
                            tmp_thrP[1] = i;
                            tmp_thrP[2] = j;
                            tmpFes = CalFes(tmp_thrP);
                            if (minFes - tmpFes > precession){
                                minFes = tmpFes;
                                selecP = tmp_thrP;
                            }
                        }
                    }
                }
            }

        }
        radius += L;
    }
//    if (thrP[0] == 45)
//        std::cout << "selecP ok" <<std::endl;
    return selecP;
};

double CalFes(std::array<int,3> CanP){ //计算候选点的适应度
    double f_p = minNum,d_p;
    double tmp_fp;

    for (int i = 0; i != X+1; ++i){
        for (int j = 0; j != Y+1; ++j){
            for (int k = 0; k != Z+1; ++k){
                if (map_info[i][j][k] == 1){
                    tmp_fp = distance(i,j,k,CanP[0],CanP[1],CanP[2]);
                    if (f_p - tmp_fp > precession){
                        f_p = tmp_fp;
                    }
                }
            }
        }
    }

    d_p = distance(start_point[0],start_point[1],start_point[2],CanP[0],CanP[1],CanP[2]) + distance(end_point[0],end_point[1],end_point[2],CanP[0],CanP[1],CanP[2]);
    return f_p / d_p;
}

bool compare2(std::array<int,3> a,std::array<int,3> b){ //用于排序W的功能函数，将W按内含的点的x坐标从小到大排序
    return a[0] < b[0];
}
std::vector<std::array<int,3>> ReW(std::vector<std::array<int,3>> W,std::array<int,3> addP){
    W.push_back(addP);
    std::sort(W.begin(),W.end(),compare2);

//    std::sort(W.begin(),W.end(),[](std::array a,std::array b)->bool{return a[0] < b[0];});
    return W;
};

std::vector<std::array<int,3>> SetW(std::vector<std::array<int,3>> W, std::vector<std::array<int,2>> P){ //更新W
    int selecNum = end_point[0] - start_point[0] - 1;
    for (int i = 0; i != selecNum; ++i){
        auto preANext = SelectPAN(W,P[i][0]); //选择前，后置点
        auto thrP = ThrPoint(preANext[1],preANext[0],P[i][0]); //选择被切点
        auto selecP = SelecPoint(thrP); //通过被切点选择新点
        W = ReW(W,selecP); //更新W
        std::cout << "plane:" << P[i][0] << "finished." << std::endl;
    }
    return W;
};

double distance(double x1,double y1,double z1,double x2,double y2,double z2){
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

int *** IniPop(std::vector<std::array<int,3>> W,int *** popu_now, int signal1){ //初始化种群函数
    /*
    popu_now = new int **[Pop];
    for (int i = 0; i != Pop; ++i){
        popu_now[i] = new int*[D];
        for (int j = 0; j != D; ++j){
            popu_now[i][j] = new int[3];
        }
    }
     */
    signal1 = 0;
    for (int i = 0; i != D; ++i){ //初始化时第一个个体为特殊个体，它的值按照集合方法计算（W）
        popu_now[0][i][0] = W[i][0];
        popu_now[0][i][1] = W[i][1];
        popu_now[0][i][2] = W[i][2];
    }
    for (int i = 1; i != Pop; ++i){ //其余个体，x坐标不变，y，z坐标在［0，Y］和［0，Z］中随机
        for (int j = 1; j!= D-1; ++j){
            popu_now[i][j][0] = popu_now[0][j][0];
            popu_now[i][j][1] = (int)(Y * (rand() / (RAND_MAX+1.0)));
            popu_now[i][j][2] = (int)(Z * (rand() / (RAND_MAX+1.0)));
        }
        popu_now[i][0][0] = popu_now[0][0][0]; //每个个体第一个和最后一个染色体不变，为起点和终点
        popu_now[i][0][1] = popu_now[0][0][1];
        popu_now[i][0][2] = popu_now[0][0][2];
        popu_now[i][D-1][0] = popu_now[0][D-1][0];
        popu_now[i][D-1][1] = popu_now[0][D-1][1];
        popu_now[i][D-1][2] = popu_now[0][D-1][2];

    }
//    std::cout << "pop_ini ok" << std::endl;
    return popu_now;
};

std::array<int,3> Mut(int ***popu_now,std::array<int,3> pop,std::array<int,3> best_xyz,double F_tmp){ //变异的小函数，为每个个体的每个染色体提供变异
    int best_y = 0,best_z = 0,r1,r2,x,y,z,s_num = 0; //best为当前种群最好个体
    best_y = best_xyz[1];
    best_z = best_xyz[2];
    y = z = -1;
    x = pop[0];
    while (y < 0 || y > Y || z < 0 || z > Z){ //当变异超越环境，重新变异
        r1 = (int)(rand()/(RAND_MAX + 1.0) * Pop);
        r2 = (int)(rand()/(RAND_MAX + 1.0) * Pop);
        y = (int)(pop[1] + F_tmp * (best_y - pop[1]) + F_tmp * (popu_now[r1][x-start_point[0]][1] - popu_now[r2][x-start_point[0]][1]));
        z = (int)(pop[2] + F_tmp * (best_z - pop[2]) + F_tmp * (popu_now[r1][x-start_point[0]][2] - popu_now[r2][x-start_point[0]][2]));
//        ++s_num;
//        if (s_num > 5)
//            RandTIme();
    }
//    std::cout << x << '\t' << y << '\t' << z << std::endl;
    return std::array<int,3>{x,y,z};
}
std::array <int,3> Mut1(std::array<int,3> pop){
    int x,y,z;
    double r1,r2;
    double r3 = 0.5;
    if (map_info[pop[0]][pop[1]][pop[2]] == 1){
        y = z = -1;
        x = pop[0];
        while (y < 0 || y > Y || z < 0 || z > Z || (map_info[x][y][z] == 1)){
            r1 = rand()/(RAND_MAX + 1.0);
            if (r3 - r1 >= precession){
                r2 = rand()/(RAND_MAX + 1.0);
            } else {
                r2 =  - (rand()/(RAND_MAX + 1.0));
            }
            y = (int)(pop[1] + r2 * (Y / 3));
            z = (int)(pop[2] + r2 * (Z / 3));
        }
    } else {
        x = pop[0];
        y = pop[1];
        z = pop[2];
    }
    return std::array<int,3>{x,y,z};
}


int *** Mutation(int *** popu_now, int *** popu_next, int signal1, int iteTh){ //变异函数，调用Mut完成种群变异
    /*
    popu_next = new int **[Pop];
    for (int i = 0; i != Pop; ++i){
        popu_next[i] = new int *[D];
        for (int j = 0; j!= D; ++j){
            popu_next[i][j] = new int[3];
        }
    }
     */

    double F_tmp = Fmin + (Fmax - Fmin) * (maxIte - iteTh)/maxIte;
    auto popu_fit = PopFits(popu_now,iteTh);
    std::array<int,3> pop,pop_next;
    double minFit = minNum;
    std::vector<std::array<int,3>> best;
    std::array<int,3> tmp;
    std::array<int,3> best_xyz;
    int index = 0;

    if (signal1 == 0){
        for (int i = 0; i!= Pop; ++i){ //寻找当前种群最好个体
            if (minFit - popu_fit[i] > precession){
                index = i;
                minFit = popu_fit[i];
                best.clear();
                for (int j = 0; j != D; ++j){
                    tmp[0] = popu_now[i][j][0];
                    tmp[1] = popu_now[i][j][1];
                    tmp[2] = popu_now[i][j][2];
                    best.push_back(tmp);
                }
            }
        }
//    std::cout << "best:" << index << std::endl;
        for (int i = 0; i != Pop ; ++i){ //每个个体的每条染色体输入Mut函数分别变异
            for (int j = 1; j != D-1; ++j){
                pop[0] = popu_now[i][j][0];
                pop[1] = popu_now[i][j][1];
                pop[2] = popu_now[i][j][2];
                best_xyz = best[j];
                pop_next = Mut(popu_now,pop,best_xyz,F_tmp);
//            std:: cout << j << std::endl;
                popu_next[i][j][0] = pop_next[0];
                popu_next[i][j][1] = pop_next[1];
                popu_next[i][j][2] = pop_next[2];
            }
            popu_next[i][0][0] = popu_now[i][0][0]; //第一和最后一条染色体不变异
            popu_next[i][0][1] = popu_now[i][0][1];
            popu_next[i][0][2] = popu_now[i][0][2];
            popu_next[i][D-1][0] = popu_now[i][D-1][0];
            popu_next[i][D-1][1] = popu_now[i][D-1][1];
            popu_next[i][D-1][2] = popu_now[i][D-1][2];
        }
    } else {
        for (int i = 0; i != Pop ; ++i) { //每个个体的每条染色体输入Mut函数分别变异
            for (int j = 1; j != D - 1; ++j) {
                pop[0] = popu_now[i][j][0];
                pop[1] = popu_now[i][j][1];
                pop[2] = popu_now[i][j][2];
                pop_next = Mut1(pop);
//            std:: cout << j << std::endl;
                popu_next[i][j][0] = pop_next[0];
                popu_next[i][j][1] = pop_next[1];
                popu_next[i][j][2] = pop_next[2];
            }
            popu_next[i][0][0] = popu_now[i][0][0]; //第一和最后一条染色体不变异
            popu_next[i][0][1] = popu_now[i][0][1];
            popu_next[i][0][2] = popu_now[i][0][2];
            popu_next[i][D - 1][0] = popu_now[i][D - 1][0];
            popu_next[i][D - 1][1] = popu_now[i][D - 1][1];
            popu_next[i][D - 1][2] = popu_now[i][D - 1][2];
        }
    }

//    std::cout << "Mutation ok" << std::endl;
    return popu_next;
}

double SetCR(double ave_fit,double min_fit,double max_fit,double me_fit){
    double CR;
    if (me_fit - ave_fit >= precession){ //revise1
        CR = CRl + (CRu - CRl)*(me_fit - min_fit)/(max_fit - min_fit);
    } else {
        CR = CRl;
    }
    return CR;
}

double AveFit(int *** popu_now, int iteTh){
    double sum_fit = 0;
    auto pops_fit = PopFits(popu_now,iteTh);
    for(auto &i : pops_fit){
        sum_fit += i;
    }
    return sum_fit/Pop;
}

std::array<double,2> MaxAMinFit(int *** popu_now, int iteTh){
    double min_fit = minNum;
    double max_fit = 0;
    auto pops_fit = PopFits(popu_now,iteTh);
    for (auto &i : pops_fit){
        if (min_fit - i >= precession)
            min_fit = i;
        if (i - max_fit >= precession)
            max_fit = i;
    }
    return std::array<double,2>{min_fit,max_fit};
}

int *** Crossover(int *** popu_now, int *** popu_next, int *** popu_cross, int iteTh){ //交叉函数
    int jrand;
    double j_rand;
    double min_fit,max_fit,ave_fit,me_fit;
    /*
    int *** popu_cross = new int **[Pop];
    for (int i = 0; i != Pop; ++i){
        popu_cross[i] = new int*[D];
        for (int j = 0; j != D; ++j){
            popu_cross[i][j] = new int[3];
        }
    }
     */
    double CR;
    ave_fit = AveFit(popu_now,iteTh);
    auto tmp_fit = MaxAMinFit(popu_now,iteTh);
    min_fit = tmp_fit[0];
    max_fit = tmp_fit[1];
    for (int i = 0; i != Pop; ++i){
        jrand = (int)((D-1) * (rand()/(RAND_MAX+1.0)));
        me_fit = PopFit(popu_now[i], iteTh);
        CR = SetCR(ave_fit,min_fit,max_fit,me_fit);
        for (int j = 0; j != D; ++j){
            j_rand = rand()/(RAND_MAX+1.0);
//            std::cout << popu_next[i][j][0] << '\t' << popu_next[i][j][1] << '\t' << popu_next[i][j][2] << '\t' << popu_now[i][j][0] << '\t' << popu_now[i][j][1] << '\t' << popu_now[i][j][2] << std::endl;
            if (CR - j_rand >= precession || j == jrand){ //满足条件，选择变异的染色体作为交叉个体的染色体，否则选择原来的个体作为染色体
                popu_cross[i][j][0] = popu_next[i][j][0];
                popu_cross[i][j][1] = popu_next[i][j][1];
                popu_cross[i][j][2] = popu_next[i][j][2];
            } else {
                popu_cross[i][j][0] = popu_now[i][j][0];
                popu_cross[i][j][1] = popu_now[i][j][1];
                popu_cross[i][j][2] = popu_now[i][j][2];
            }
//            std::cout << popu_next[i][j][0] << '\t' << popu_next[i][j][1] << '\t' << popu_next[i][j][2] << '\t' << popu_now[i][j][0] << '\t' << popu_now[i][j][1] << '\t' << popu_now[i][j][2] << std::endl;

        }
    }
//    std::cout << "Crossover ok" << std::endl;
    return popu_cross;
}

int SetSignal(int ini, int signal1){
    if (ini >= mut1){
        signal1 = 1;
    } else {
        signal1 = 0;
    }
    return signal1;
}

int *** Selection(int *** popu_now, int *** popu_cross, int signal1, int iteTh){ //选择函数，从原来种群和交叉种群的每个个体中选择表现更好的个体作为新个体
    double now_fit,cross_fit;
    double now_max = 0,cross_max = 0;
    /*
    for (int i = 0; i != Pop; ++i){
        for (int j = 0; j != D; ++j){
            std::cout << popu_now[i][j][0] << '\t' << popu_now[i][j][1] << '\t' << popu_now[i][j][2] << '\t' << popu_cross[i][j][0] << '\t' << popu_cross[i][j][1] << '\t' << popu_cross[i][j][2] << std::endl;

        }
    }
     */
    for (int i = 0; i != Pop; ++i){
        now_fit = PopFit(popu_now[i],iteTh); //原来种群的适应度
        cross_fit = PopFit(popu_cross[i],iteTh); //交叉种群的适应度
        if (now_fit - now_max >= precession) //modify1
            now_max = now_fit;
        if (cross_fit - cross_max >= precession)
            cross_max = cross_fit;
        if (cross_fit - now_fit >= precession){ //选择
            for (int j = 0; j != D; ++j){
                popu_now[i][j][0] = popu_cross[i][j][0];
                popu_now[i][j][1] = popu_cross[i][j][1];
                popu_now[i][j][2] = popu_cross[i][j][2];
            }
        }
    }
    std::cout << now_max << '\t' << cross_max << std::endl;

//    std::cout << "Selection ok" << std::endl;
    return popu_now;
}

int ** GetOne(int *** popu_now, int iteTh){ //从种群中选择表现最好的个体
    std::array<double,Pop> pop_fits = PopFits(popu_now,iteTh);
    double maxFit = 0;
    int pop_index = 0;
    for (int i = 0; i!= Pop; ++i){ //modify1
        if (pop_fits[i] - maxFit >= precession){
            maxFit = pop_fits[i];
            pop_index = i;
        }
    }
//    std::cout << "GetOne ok" << std::endl;
    return popu_now[pop_index];
}

std::string WriteOne(int ** get_one){ //将结果写入'result.txt'
    double dis = 0;
    std::ofstream os{"result"};
    std::ofstream ods{"distance"};
    for (int i = 0; i != D; ++i){
        os << get_one[i][0] << '\t' << get_one[i][1] << '\t' << get_one[i][2] << std::endl;
        if (i != D - 1){
            dis += distance(get_one[i][0],get_one[i][1],get_one[i][2],get_one[i+1][0],get_one[i+1][1],get_one[i+1][2]);
        }
    }
    ods << dis << std::endl;
    return "Compelete!";
}

double PopFit(int **pop_now, int iteTh){ //个体适应度
//    std::cout << Collosion(pop_now) << '\t' << NorDis(pop_now) << std::endl;
    double w;
    w = 100 * iteTh/maxIte + 1;
    return (Collosion(pop_now) + NorDis(pop_now));
    //return NorDis(pop_now);
}

std::array<double,Pop> PopFits(int ***popu_now, int iteTh){ //种群适应度
    std::array<double,Pop> popu_fit;

    double w;
    w = 100 * iteTh/maxIte + 1;
    for (int i = 0; i != Pop; ++i){
        popu_fit[i] = Collosion(popu_now[i]) + NorDis(popu_now[i]);
        //popu_fit[i] = NorDis(popu_now[i]);
    }
    return popu_fit;
};

double Collosion(int **pop){ //个体包含障碍物适应度
    int pop_num = 0;
    int x,y,z;
    for (int i = 1; i!= D-1; ++i){
        x = pop[i][0];
        y = pop[i][1];
        z = pop[i][2];
        if (map_info[x][y][z] == 0) //modify1
            ++pop_num;
    }
//    std::cout << ((double)pop_num)/(D-2) << std::endl;
    return ((double)pop_num)/(D-2);
}

double NorDis(int **pop){ //个体距离适应度
    double pop_dis = 0;
    double dis;
    int x,y,z,x1,y1,z1;
    for (int i = 0; i!= D-1; ++i){
        x = pop[i][0];
        y = pop[i][1];
        z = pop[i][2];
        x1 = pop[i+1][0];
        y1 = pop[i+1][1];
        z1 = pop[i+1][2];
        pop_dis += distance(x,y,z,x1,y1,z1);
    }
    dis = distance((double)start_point[0],(double)start_point[1],(double)start_point[2],end_point[0],end_point[1],end_point[2]);
//    std::cout << pop_dis << std::endl;

    return dis/pop_dis; //modify1
/*
    if (pop_dis - dis <= 100){
        return atan(pop_dis)/(PI/2);
    } else {
        return atan(pop_dis/10)/(PI/2);
    }
    */
//    std::cout << pop_dis/dis << std::endl;
//    return pop_dis/dis;
}

std::ofstream WriteResult(std::vector<std::array<int,3>> W){ //将W写入'W.txt'
    std::string W_result{"W.txt"};
    std::ofstream os{W_result};
    for (auto &i : W){
        os << i[0] << '\t' << i[1] << '\t' << i[2] << '\n';
    }
    return os;
}

int ** DEPhase(std::vector<std::array<int,3>> W, int ***popu_now, int ***popu_next ,int ***popu_cross){ //DE
    double tmpFit = 0;
    double tmpFit1 = 0;
    std::ofstream fit_stream{"fit.txt"};
    popu_now = IniPop(W,popu_now,signal1); //初始化种群
    auto get_one = GetOne(popu_now,0); //选择最好的个体
    tmpFit = PopFit(get_one,0);
    for (int i = 0; i != maxIte; ++i){
        popu_next = Mutation(popu_now,popu_next,signal1,i); //种群变异
        popu_cross = Crossover(popu_now,popu_next,popu_cross,i); //交叉
        popu_now = Selection(popu_now,popu_cross,signal1,i); //选择
        signal1 = SetSignal(i,signal1);
        if (i%500 == 0)
            RandTIme();
        if ((i+1)%reIni == 0){ //判断是否更好
            auto get_one1 = GetOne(popu_now,i);
            tmpFit1 = PopFit(get_one1,i);
            fit_stream << tmpFit << '\t' << tmpFit1 << std::endl;
            std::cout << tmpFit << '\t' << tmpFit1 << std::endl;
            if (tmpFit - tmpFit1 >= precession) { // 没有更好就重新初始化
                popu_now = IniPop(W,popu_now,signal1);
                tmpFit = PopFit(GetOne(popu_now,i),i);
                i = 0;
            } else {
                tmpFit = tmpFit1;
            }
            RandTIme(); //随机数种子
        }
        std::cout << i << std::endl;
    }
    get_one = GetOne(popu_now,maxIte);
    return get_one;
}

int *** NewSpa(int *** popu){ //为种群分配空间
    popu = new int **[Pop];
    for (int i = 0; i != Pop; ++i){
        popu[i] = new int *[D];
        for (int j = 0; j != D; ++j){
            popu[i][j] = new int[3];
        }
    }
    return popu;
}

void FreePopu(int *** popu){ //释放空间
    for (int  i = 0; i != Pop; ++i){
        for (int j = 0; j != D; ++j){
            free(popu[i][j]);
        }
    }
    for (int i = 0; i!= Pop; ++i){
        free(popu[i]);
    }
    free(popu);
}

std::vector<std::array<int,3>> MakeW(std::vector<std::array<int,3>> W){
    std::ifstream is{"W.txt"};
    std::istream_iterator<int> ii{is};
    std::array<int,3> tmp_p;
    for (int i = 0; i != D; ++i){
        for (auto &j : tmp_p){
            j = *ii++;
        }
        W.push_back(tmp_p);
    }

    return W;
};

int main() {
    MakeMap1();
    RandTIme(); //随机数种子
//    W = IniW(W,start_point,end_point); //初始化W
//    P = SepPlane(P,start_point,end_point); //创建P
//    W = SetW(W,P); //创建W
//    WriteResult(W);
    W = MakeW(W);
    popu_now = NewSpa(popu_now); //分配空间
    popu_next = NewSpa(popu_next);
    popu_cross = NewSpa(popu_cross);
    auto get_one = DEPhase(W,popu_now,popu_next,popu_cross); //DE
    std::cout << WriteOne(get_one) << std::endl;
    FreePopu(popu_now); //释放空间
    FreePopu(popu_next);
    FreePopu(popu_cross);
    return 0;
}