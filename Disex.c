#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define N 100//200
#define a_max 1.0 // 4.767
#define a_min 0.0
#define b_max 3.0 // 4.767
#define b_min 0.0
#define mu 0.5
#define v_muu 0.0//2.0
#define d 0.01
#define dt 1.0/(b_max+d+1.0)
#define M  (int)b_max+1000 //500
#define T  500
#define T_f 2000
#define Q 1    //平均を取るための試行回数


int main(){
    double con_a[N][N],P_list[2]={0.0},con_b[N][N],F_list[7]={0.5};//{0.1,0.5,1.0,2.0,3.0,10.0,100.0}
    double R_list[15]={3.0};
    int i,j,k,l,m,n,h,q,con[N][N],Rsize,count;//{2,3,4,5,7,9,13,17,20,30,50,60,80,100,140,180,200,300};
    double eps,z,ramda,a,b,pr,def,c,q_ss,q_is,q_os,q_si,q_ii,q_oi,q_so,q_io,q_oo,q_sj,q_ij,q_jj,q_oj;
    int ns,ni,ne,no,nj,i_s,i_i,i_o,s_s,s_i,s_o,o_s,o_i,o_o,j_o,j_s,j_i,j_j;
    ns=0;ni=0;ne=0;i_s=0;i_i=0;i_o=0;s_s=0;s_i=0;s_o=0;o_s=0;o_i=0;o_o=0;

    
    FILE *gp,*data1,*data2,*data3,*data4,*data5,*data6;
    char *data_file1,*data_file2,*data_file3,*data_file4,*data_file5,*data_file6;

    int t[15000],I[15000],S[15000],E[15000],O[15000];
    double x[T],y1[T],y2[T],y3[T],y4[T],y5[T],y6[T],y7[T],y8[T],y9[T],y10[T],y11[T],a_low[T],a_hi[T];

    srand((unsigned int)time(NULL));

    printf("%d",Rsize);

    eps=0.8093;z=4.0;
    for(i=0;i<70;i++){
        x[i]=(double)i*0.01;
        a=(double)i*0.01;
        y1[i]=(double)((z+1)*a*a+a*(-2*(1-eps)+2*z*(1-eps)+pow(z,2))+(-4*eps*(z-1)+2*z*(1-eps)-pow(z,2)*(1-2*eps))+pow(pow(a+(z-2),2)*((z+1)*(z+1)*a*a+a*(4*eps-2*z+2*z*z*(3-2*eps))+(4*eps*eps+4*eps*z*(1-2*eps)+z*z*(1-2*eps)*(1-2*eps))),1/2))/(2*(-a+(1-eps)*(z-1)*(z-2)));
        y2[i]=(double)a*z*(a*z*(1-2*eps)-eps*(z-1))/(a*a*eps*z*(z-1)+a*(z*z*(1-eps)-z*(1-2*eps))-eps*(z-1)*(z-1));
    }
    
    
    //図の描画
        data_file1="EXline.dat";
        // data_file1="outadde_f.dat";
        data1=fopen(data_file1,"w");
        for(l=0;l<70;l++){
            fprintf(data1,"%f\t%f\t%f\n", x[l],y1[l],y2[l]);
        }
        fclose(data1);
        gp=popen("gnuplot -persist","w");
        fprintf(gp,"set terminal png\n");
        //fprintf(gp,"set logscale\n");
        // fprintf(gp,"set output 'Addefunction_f_%2f_P_%2f.png'\n",F_list[0],P_list[0]);
        fprintf(gp,"set output 'EXline.png'\n");

        
        fprintf(gp,"set xrange [0:%f]\n",0.61);
        fprintf(gp,"set xlabel '1/m_+'\n");
        fprintf(gp,"set yrange [0:%d]\n",30);//5.0
        fprintf(gp,"set ylabel 'm__/m_+'\n");
        // fprintf(gp,"f(x)=2*x/(x-1)\n");

        fprintf(gp,"plot \'%s\'using 1:2 with lines linetype 1 linecolor rgb 'red' title \"Pathegen Driven Extinction \",\'%s\' using 1:3 with lines linetype 1 linecolor rgb 'black' title \"Disease free\"\n",data_file1,data_file1);
        // fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"f=%f \"",data_file1,F_list[0]);
        // for(l=1;l<6;l++){
        //     fprintf(gp,",\'%s\' using 1:%d with lines linetype %d title \"f=%f \"",data_file1,l+2,l+1,F_list[l]);
        // }
        // fprintf(gp,",\'%s\' using 1:%d with lines linetype %d title \"f=%f \"\n",data_file1,9,8,F_list[6]);
        
        
        
        pclose(gp);



    }







