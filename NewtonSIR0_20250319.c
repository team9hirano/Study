#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define N 200//200
#define a_max 4.767 // 4.767
#define a_min 0.2
#define b_max 200.0 // 4.767
#define b_min 5.0
#define mu 0.5
#define ep 0.8093
#define v_muu 2.0//2.0
#define d 1.0
#define dt 1.0/(b_max+d+1.0)
#define M  (int)b_max+d+800+3000 //b_max+1000
#define T  500
#define T_f 2000
#define z 4.0
#define Q 1    //平均を取るための試行回数

double infection(int con[][N],double con_a[][N],double con_b[][N],double P_list[],double F_list[],double c,double r,int P_size,int f_size,int x[],double y[],int S[],int I[],int O[],int g);


int main(){
    double J[N][N];//{0.1,0.5,1.0,2.0,3.0,10.0,100.0}
    double F_list[9],r[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},A_list[T],B_list[T],R_list[T],D_list[T],A2_list[T],B2_list[T];//5.0,7.0,8.0,9.0,10.0,15.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0
    int i,j,k,l,m,h,q,con[N][N],Rsize,n_max,n_min,n,number,count,count_max;
    double rhos,rhoi,rhor,A,B,C,D,E,F,G,H,I,beta,alpha,gamma,omega,m1,m2,m3,R,befrhos,befrhoi,befrhor,befA,befB,befC,befD,befE,befF,befG;
    
    FILE *gp,*data1,*data2,*data3,*data4,*data5,*data6;
    char *data_file1,*data_file2,*data_file3,*data_file4,*data_file5,*data_file6;

    // int x[T],t[15000],I[15000],S[15000],E[15000],O[15000];
    // double y1[T],y2[T],y3[T],y4[T],y5[T],y6[T],y7[T],y8[T],y9[T],y10[T],y11[T];

    srand((unsigned int)time(NULL));
    //絶滅境界線
    count=0;
    count_max=10;
    for(m=1;m<T;m++){
        count=0;
        if(m<=10){
            R=(double)m/10+1.0;
       }else{
            R=(double)m-8.0;
       }
        // printf("alpha=%f\n",alpha);
        A_list[m]=R;
        B_list[m]=(double)b_max;
        
        for(h=(int)b_max;h>0;h--){
            beta=(double)h;
            if(count > count_max)break;
            A=0.25;
            B=0.25;
            C=0.75;
            D=0.25;
            E=0.25;
            F=0.55;
            G=0.25;
            rhos=0.25;
            rhoi=0.25;
            rhor=0.25;
            //R=5.0;
            //beta=100.0;
            befA=A;befB=B;befC=C;befD=D;befF=F;befG=G;befrhos=rhos;befrhoi=rhoi;befrhor=rhor;
            alpha=1.0;
            gamma=1.0;
            omega=10.0;
        
            for(l=0;l<100000;l++){
                J[0][0]=R*A-beta*B-d;
                J[0][1]=0;
                J[0][2]=0;
                J[0][3]=R*rhos;
                J[0][4]=-beta*rhos;
                J[0][5]=0;
                J[0][6]=0;
                J[0][7]=0;
                J[0][8]=0;
                J[0][9]=0;
                J[1][0]=beta*B;
                J[1][1]=-(d+alpha+gamma);
                J[1][2]=0;
                J[1][3]=0;
                J[1][4]=beta*rhos;
                J[1][5]=0;
                J[1][6]=0;
                J[1][7]=0;
                J[1][8]=0;
                J[1][9]=0;
                J[2][0]=0;
                J[2][1]=gamma;
                J[2][2]=-d;
                J[2][3]=0;
                J[2][4]=0;
                J[2][5]=0;
                J[2][6]=0;
                J[2][7]=0;
                J[2][8]=0;
                J[2][9]=0;
                J[3][0]=(-R*(1-1/z)*A*A*(1-rhoi-rhor))/((1-rhos-rhoi-rhor)*(1-rhos-rhoi-rhor));
                J[3][1]=(-R*(1-1/z)*A*A*rhos)/((1-rhos-rhoi-rhor)*(1-rhos-rhoi-rhor));
                J[3][2]=(-R*(1-1/z)*A*A*rhos)/((1-rhos-rhoi-rhor)*(1-rhos-rhoi-rhor));
                J[3][3]=R*(1-1/z)*D+beta*B/z-R*(1/z+2*(1-1/z)*rhos*A/(1-(rhos+rhoi+rhor)))-2*R*A;
                J[3][4]=alpha+beta*A/z;
                J[3][5]=0;
                J[3][6]=R*(1-1/z)*A;
                J[3][7]=0;
                J[3][8]=0;
                J[3][9]=0;
                J[4][0]=0;
                J[4][1]=0;
                J[4][2]=0;
                J[4][3]=R*(1-1/z)*E-R*B;
                J[4][4]=beta*(1-1/z)*C-(2*d+alpha+gamma)-beta*(1-2*B)/z-R*A+d;
                J[4][5]=beta*(1-1/z)*B;
                J[4][6]=0;
                J[4][7]=R*(1-1/z)*A;
                J[4][8]=0;
                J[4][9]=0;
                J[5][0]=2*R*A*A*(1-(rhoi+rhor))/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)));
                J[5][1]=2*R*A*A*rhos/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)));
                J[5][2]=2*R*A*A*rhos/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)));
                J[5][3]=2*R*(1/z+2*(1-1/z)*rhos*A/(1-(rhos+rhoi+rhor)))-R*C;
                J[5][4]=-beta*(1-2/z)*C;
                J[5][5]=-d-beta*(1-2/z)*B;
                J[5][6]=0;
                J[5][7]=0;
                J[5][8]=0;
                J[5][9]=0;
                J[6][0]=-(R*(1-2/z)*A*D*rhos+d*D*(rhos+rhoi+rhor)+alpha*D*rhoi)/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)))-(R*(1-2/z)*A*D+d*D)/(1-(rhos+rhoi+rhor));
                J[6][1]=-(R*(1-2/z)*A*D*rhos+d*D*(rhos+rhoi+rhor)+alpha*D*rhoi)/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)))-(d+alpha)*D/(1-(rhos+rhoi+rhor));
                J[6][2]=-(R*(1-2/z)*A*D*rhos+d*D*(rhos+rhoi+rhor)+alpha*D*rhoi)/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)))-d*D/(1-(rhos+rhoi+rhor));
                J[6][3]=-R*(1-2/z)*D*rhos/(1-(rhos+rhoi+rhor));
                J[6][4]=0;
                J[6][5]=0;
                J[6][6]=-2*d-(R*(1-2/z)*A*rhos+d*(rhos+rhoi+rhor)+alpha*rhoi)/(1-(rhos+rhoi+rhor));
                J[6][7]=-2*d;
                J[6][8]=0;
                J[6][9]=0;
                J[7][0]=(beta*(1-1/z)*A*B*rhos+d*B*rhos+(d+alpha)*F*rhoi+d*G*rhoi-d*(rhos+rhoi+rhor)*E-alpha*E*rhoi+R*A*E*rhos/z)/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)))+(beta*(1-1/z)*A*B+d*B-d*E+R*A*E/z)/(1-(rhos+rhoi+rhor));
                J[7][1]=(beta*(1-1/z)*A*B*rhos+d*B*rhos+(d+alpha)*F*rhoi+d*G*rhoi-d*(rhos+rhoi+rhor)*E-alpha*E*rhoi+R*A*E*rhos/z)/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)))+((d+alpha)*F+d*G-d*E-alpha*E)/(1-(rhos+rhoi+rhor));
                J[7][2]=(beta*(1-1/z)*A*B*rhos+d*B*rhos+(d+alpha)*F*rhoi+d*G*rhoi-d*(rhos+rhoi+rhor)*E-alpha*E*rhoi+R*A*E*rhos/z)/((1-(rhos+rhoi+rhor))*(1-(rhos+rhoi+rhor)))+d*E/(1-(rhos+rhoi+rhor));
                J[7][3]=(beta*(1-1/z)*B*rhos-R*E*rhos/z)/(1-(rhos+rhoi+rhor));
                J[7][4]=(beta*(1-1/z)*A*rhos+d*rhos)/(1-(rhos+rhoi+rhor));
                J[7][5]=0;
                J[7][6]=0;
                J[7][7]=(-d*(rhos+rhoi+rhor)-alpha*rhoi+R*A*rhos/z)/(1-(rhos+rhoi+rhor))-(d+alpha+gamma);
                J[7][8]=(d+alpha)*rhoi/(1-(rhos+rhoi+rhor));
                J[7][9]=d*rhoi/(1-(rhos+rhoi+rhor));
                J[8][0]=2*beta*(1/z+(1-1/z)*B)*B/rhoi-B*F/rhoi;
                J[8][1]=-2*beta*(1/z+(1-1/z)*B)*B*rhos/(rhoi*rhoi)+rhos*B*F/(rhoi*rhoi);
                J[8][2]=0;
                J[8][3]=0;
                J[8][4]=2*beta*(1/z+2*(1-1/z)*B)-rhos*F/rhoi;
                J[8][5]=0;
                J[8][6]=0;
                J[8][7]=0;
                J[8][8]=-(d+alpha+gamma)-rhos*B/rhoi;
                J[8][9]=0;
                J[9][0]=R*A-beta*B-d;
                J[9][1]=0;
                J[9][2]=0;
                J[9][3]=R*rhos;
                J[9][4]=-beta*rhos;
                J[9][5]=0;
                J[9][6]=0;
                J[9][7]=0;
                J[9][8]=0;
                J[9][9]=0;
                F_list[0]=-(beta*(1-1/z)*(C-A-B)*x-d*x-beta*(1/z+(1-1/z)*A)*x-beta*x*x);
                F_list[1]=-(beta*(1-1/z)*A*(C-A-B)-(d+alpha+gamma)*A-beta*(1/z+(1-1/z)*A)*A-gamma*A*B+beta*A*A);
                F_list[2]=-(R*(1-1/z)*B-d*B+alpha*A+d+beta*A*B/z-R*B/z-gamma*B*B);
                F_list[3]=-(-gamma*A+(1-C)*(d+R*B-beta*A+beta*(1-1/z)*A));
                F_list[4]=-(beta*(1-1/z)*x*B+d*x+(d+alpha)*(s-x-y)+d*(1-s)-beta*x*y);
                F_list[5]=-(-beta*(1-1/z)*x*(1-C)-gamma*(s-x-y)+(1-s)*(d+beta*x));
                F_list[6]=;
                F_list[7]=;
                F_list[8]=;
                if(h==b_max&&m==1&&l==0){
                    printf("J11=%f,J12=%f,J13=%f,J14=%f,J15=%f,J16=%f\n",J[0][0],J[0][1],J[0][2],J[0][3],J[0][4],J[0][5]);
                    printf("J21=%f,J22=%f,J23=%f,J24=%f,J25=%f,J26=%f\n",J[1][0],J[1][1],J[1][2],J[1][3],J[1][4],J[1][5]);
                    printf("J31=%f,J32=%f,J33=%f,J34=%f,J35=%f,J36=%f\n",J[2][0],J[2][1],J[2][2],J[2][3],J[2][4],J[2][5]);
                    printf("J41=%f,J42=%f,J43=%f,J44=%f,J45=%f,J46=%f\n",J[3][0],J[3][1],J[3][2],J[3][3],J[3][4],J[3][5]);
                    printf("J51=%f,J52=%f,J53=%f,J54=%f,J55=%f,J56=%f\n",J[4][0],J[4][1],J[4][2],J[4][3],J[4][4],J[4][5]);
                    printf("J61=%f,J62=%f,J63=%f,J64=%f,J65=%f,J66=%f\n",J[5][0],J[5][1],J[5][2],J[5][3],J[5][4],J[5][5]);
                    printf("F0=%f,F1=%f,F2=%f,F3=%f,F4=%f,F5=%f\n",F_list[0],F_list[1],F_list[2],F_list[3],F_list[4],F_list[5]);
                }
                //ガウスの消去法
                number=6;
                for(i=0;i<number;i++){
                    for(j=i+1;j<number;j++){
                    m1=J[j][i]/J[i][i];
                    for(k=i+1;k<number;k++){
                        J[j][k]=J[j][k]-m1*J[i][k];
                        
                    } 
                    F_list[j]=F_list[j]-m1*F_list[i];
                    }
                }
                for(i=number-1;i>=0;i--){
                    m2=0.0;
                    if(i==number-1){
                        m2=0.0;
                    }
                    else{
                        for(j=i+1;j<number;j++){
                            m2=m2+J[i][j]*r[j];
                        }
                    }
                    
                    r[i]=(1/J[i][i])*(F_list[i]-m2);
                }

                if(fabs(r[0])<0.001&&fabs(r[1])<0.001&&fabs(r[2])<0.001&&fabs(r[3])<0.001&&fabs(r[4])<0.001&&fabs(r[5])<0.001){
                    
                    if(R*B-d<0&&beta*x-(d+alpha+gamma)<0){   //r*B-d,beta*x-(d+alpha+gamma) or R*B-beta*A-d
                
                        B_list[m]=beta;
                        
                        
                
                    }else{
                        count=count+1;
                    }
                    //printf("でたよ");
                    break;
                    
                }
                befx=x;befy=y;befs=s;befA=A;befB=B;befC=C;
                    x=x+r[0];
                    A=A+r[1];
                    B=B+r[2];
                    C=C+r[3];
                    y=y+r[4];
                    s=s+r[5];
                

            }
        
        }
     
        printf("r=%f,x=%f,A=%f,B=%f,C=%f,y=%f,s=%f\n",R,x,A,B,C,y,s);
    }

    //diseasefree境界線
    count=0;
    count_max=10;
    for(m=1;m<T;m++){
       count=0;
       if(m<=10){
            R=(double)m/10+1.0;
       }else{
            R=(double)m-8.0;
       }
        //R=(double)m;
        // printf("alpha=%f\n",alpha);
        A2_list[m]=R;
        B2_list[m]=0.0;
        for(h=1;h<=(int)b_max;h++){
            beta=(double)h;
            if(count > count_max)break;
            x=0.3;
            y=0.3;
            s=0.75;
            A=0.0;
            B=0.5;
            C=1.0;
            
            rhos=0.5;
            //beta=100.0;
            gamma=1.0;
            alpha=1.0;
            omega=10.0;
        
            for(l=0;l<100000;l++){
                J[0][0]=beta*(1-1/z)*(C-A-B)-d-beta*(1/z+(1-1/z)*A)-2*beta*x;
                J[0][1]=R*(1-1/z)*(rhos/(1-rhos))*B;
                J[0][2]=0;
                J[0][3]=R*(1-1/z)*(rhos/(1-rhos))*y-beta*(1-1/z)*x;
                J[0][4]=beta*(1-1/z)*x;
                J[0][5]=R*(1-1/z)*B*y*pow((1-rhos),-2.0);
                J[1][0]=beta*(1-1/z)*B+d-(d+alpha)-beta*y;
                J[1][1]=-(d+alpha)-R*(1-1/z)*(rhos/(1-rhos))*B-beta*x;
                J[1][2]=alpha;
                J[1][3]=beta*(1-1/z)*x-R*(1-1/z)*(rhos/(1-rhos))*y;
                J[1][4]=0;
                J[1][5]=R*B*(1-1/z)*y*pow(1-rhos,-2);
                J[2][0]=-beta*(1-1/z)*(1-C)+gamma+beta*(1-s);
                J[2][1]=gamma;
                J[2][2]=-gamma-beta*x-d;
                J[2][3]=0;
                J[2][4]=beta*(1-1/z)*x;
                J[2][5]=0;
                J[3][0]=0;
                J[3][1]=0;
                J[3][2]=0;
                J[3][3]=R*(1-1/z)*(1-2*B*(rhos/(1-rhos)))-d-R*(1/z+2*B*(1-1/z)*(rhos/(1-rhos)))-2*gamma*B;
                J[3][4]=0;
                J[3][5]=-2*R*(1-1/z)*pow(1-rhos,-2);
                J[4][0]=0;
                J[4][1]=0;
                J[4][2]=0;
                J[4][3]=R*(1-C);
                J[4][4]=-(d+R*B);
                J[4][5]=0;
                J[5][0]=0;
                J[5][1]=0;
                J[5][2]=0;
                J[5][3]=R*rhos;
                J[5][4]=0;
                J[5][5]=R*B-d;
                F[0]=-(R*(1-1/z)*(rhos/(1-rhos))*B*y+beta*(1-1/z)*(C-A-B)*x-d*x-beta*(1/z+(1-1/z)*A)*x-beta*x*x);
                F[1]=-(beta*(1-1/z)*x*B+d*x+(d+alpha)*(s-x-y)+d*(1-s)-R*(1-1/z)*y*(rhos/(1-rhos))*B-beta*x*y);
                F[2]=-(-beta*(1-1/z)*x*(1-C)-gamma*(s-x-y)+(1-s)*(d+beta*x));
                F[3]=-(R*(1-1/z)*(1-(rhos/(1-rhos))*B)*B-d*B+alpha*A+d+beta*A*B/z-R*B*(1/z+(1-1/z)*(rhos/(1-rhos))*B)-gamma*B*B);//R*(1-1/z)*(1-(rhos/(1-rhos))*B)*B-d*B+alpha*A+d+beta*A*B/z-R*B*(1/z+(1-1/z)*(rhos/(1-rhos))*B)-gamma*B*B
                F[4]=-(-gamma*A+(1-C)*(d+R*B-beta*A+beta*(1-1/z)*A));//-gamma*A+(1-C)*(d+R*B-beta*A+beta*(1-1/z)*A);
                F[5]=-(rhos*(R*B-d));
                //ガウスの消去法
                number=6;
                for(i=0;i<number;i++){
                    for(j=i+1;j<number;j++){
                        m1=J[j][i]/J[i][i];
                        for(k=i+1;k<number;k++){
                            J[j][k]=J[j][k]-m1*J[i][k];
                            
                        } 
                        F_list[j]=F_list[j]-m1*F_list[i];
                    }
                }
                for(i=number-1;i>=0;i--){
                    m2=0.0;
                    if(i==number-1){
                        m2=0.0;
                    }
                    else{
                        for(j=i+1;j<number;j++){
                            m2=m2+J[i][j]*r[j];
                        }
                    }
                    
                    r[i]=(1/J[i][i])*(F_list[i]-m2);
                }

                if(fabs(r[0])<0.001&&fabs(r[1])<0.001&&fabs(r[2])<0.001&&fabs(r[3])<0.001&&fabs(r[4])<0.001&&fabs(r[5])<0.001){
                    
                    if(R*B-d<0&&beta*x-(d+alpha+gamma)<0){   //r*B-d,beta*x-(d+alpha+gamma)
                

                        B2_list[m]=beta;
                        
                        
                
                    }else{
                        count=count+1;
                    }
                    break;
                    
                }
                
                    x=x+r[0];
                    //A=A-r[1];
                    B=B+r[3];
                    C=C+r[4];
                    rhos=rhos+r[5];
                    y=y+r[1];
                    s=s+r[2];
                

            }
        
        }
    
    }



    //printf("B_list[0]=%f\n",B_list[0]);
    //S,Iなどの時間変遷図(after)
                 data_file3="Newtonline.dat";
                 data_file2="SIRDot.dat";
                 data3=fopen(data_file3,"w");
                 for(l=0;l<T;l++){
                     fprintf(data3,"%f\t%f\t%f\t%f\n",A_list[l],B_list[l],A2_list[l],B2_list[l]);
                 }
                 fclose(data3);
                 gp=popen("gnuplot -persist","w");
                 fprintf(gp,"set logscale\n");
                 fprintf(gp,"set terminal png\n");
                
                 fprintf(gp,"set output 'NewtonSIR0line.png'\n");
            

                
                 fprintf(gp,"set xrange [0:%d]\n",T+500);
                 fprintf(gp,"set xlabel 'r'\n");
                 fprintf(gp,"set yrange [0:%d]\n",T+500);
                 fprintf(gp,"set ylabel 'beta'\n");
                  //5.0
                
                 //fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"a= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,A_list[h],data_file3,data_file3);
                 fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"exline\",\'%s\' using 3:4 with lines linetype 2 title \"diseasefree \",\'%s\' using 1:2 with points pointtype 1 title \"point \"\n",data_file3,data_file3,data_file2);
                 pclose(gp);


}







// void infection(int con[][N],double con_a[][N],double con_b[][N],double P_list[],double F_list[],int P_size,int f_size,int x[],double y[],int S[],int E[],int I[],int O[],int a);
double infection(int con[][N],double con_a[][N],double con_b[][N],double P_list[],double F_list[],double c,double r,int P_size,int f_size,int x[],double y[],int S[],int I[],int O[],int g){
    int i,j,k,m,l,MM;
    double v_mu,pr,a;
    //a=1.0;
    FILE *gp,*data1,*data2,*data3,*data4,*data5,*data6;
    char *data_file1,*data_file2,*data_file3,*data_file4,*data_file5,*data_file6;
    srand((unsigned int)time(NULL));

    if(g==0){v_mu=0;MM=30;}else{v_mu=v_muu;MM=M;} //100,M
    for(j=0;j<T;j++){ //M
        double vir;
        int ni,ne;
        int ss,ee,ii,oo;
        a=r;
        for(k=0;k<MM*N*N;k++){
            int x,y,n;
            x=(int)rand()%N,y=(int)rand()%N;
            n=(int)rand()%4;pr=(double)rand()/RAND_MAX;
            if(con[x][y]==1){
                            
                if(pr < r*dt){
                    if(n==0){
                                
                        if(con[(x-1)%N][y]==0){ //(x-1)%Nについて確認
                                    
                            con[(x-1)%N][y]=1;
                        }
                                
                    }
                    else if(n==1){
                        if(con[(x+1)%N][y]==0){
                            con[(x+1)%N][y]=1;
                        }
                                    
                    }
                    else if(n==2){
                        if(con[x][(y-1)%N]==0){
                            con[x][(y-1)%N]=1;
                        }
                                    
                    }
                    else{
                        if(con[x][(y+1)%N]==0){
                            con[x][(y+1)%N]=1;
                        }
                                    
                                    
                    }
                }
                
            }
            
            else if(con[x][y]==2){
                int xx,yy;
                xx=(int)rand()%N,yy=(int)rand()%N;
                pr=(double)rand()/RAND_MAX;
                            
                if(pr < con_b[x][y]*dt){
                    if( (double)rand()/RAND_MAX < P_list[P_size]){
                        if(con[xx][yy]==1){
                            con[xx][yy]=2;
                            if((double)rand()/RAND_MAX<mu*dt){
                                if((double)rand()/RAND_MAX<0.5){
                                    con_b[xx][yy]=con_b[x][y]+v_mu;
                                    if(con_b[xx][yy]>b_max){
                                        con_b[xx][yy]=b_max;
                                    }
                                }else{
                                    con_b[xx][yy]=con_b[x][y]-v_mu;
                                    if(con_b[xx][yy]<b_min){
                                        con_b[xx][yy]=b_min;
                                    }
                                                        
                                }
                                con_a[xx][yy]=con_a[x][y];
                            }
                            else{
                                con_a[xx][yy]=con_a[x][y];
                                con_b[xx][yy]=con_b[x][y];
                            }
                                        
                        }
                    }
                    else{
                        if(n==0){
                            if(con[(x-1)%N][y]==1){
                                con[(x-1)%N][y]=2;con_a[(x-1)%N][y]=con_a[x][y];
                                if((double)rand()/RAND_MAX<mu*dt){
                                    if((double)rand()/RAND_MAX<0.5){
                                        con_b[(x-1)%N][y]=con_b[x][y]+v_mu;
                                        if(con_b[(x-1)%N][y]>b_max){
                                            con_b[(x-1)%N][y]=b_max;
                                        }
                                    }else{
                                        con_b[(x-1)%N][y]=con_b[x][y]-v_mu;
                                        if(con_b[(x-1)%N][y]<b_min){
                                            con_b[(x-1)%N][y]=b_min;
                                        }
                                                            
                                    }
                                    // con_b[(x-1)%N][y]=3*con_a[(x-1)%N][y];
                                }
                                else{
                                    con_b[(x-1)%N][y]=con_b[x][y];
                                        // con_b[(x-1)%N][y]=3*con_a[(x-1)%N][y];
                                }
                                
                            }
                        }
                        else if(n==1){
                            if(con[(x+1)%N][y]==1){
                                con[(x+1)%N][y]=2;con_a[(x+1)%N][y]=con_a[x][y];
                                if((double)rand()/RAND_MAX<mu*dt){
                                    if((double)rand()/RAND_MAX<0.5){
                                        con_b[(x+1)%N][y]=con_b[x][y]+v_mu;
                                        if(con_b[(x+1)%N][y]>b_max){
                                            con_b[(x+1)%N][y]=b_max;
                                        }
                                    }else{
                                        con_b[(x+1)%N][y]=con_b[x][y]-v_mu;
                                        if(con_b[(x+1)%N][y]<b_min){
                                            con_b[(x+1)%N][y]=b_min;
                                        }
                                                            
                                    }
                                    // con_b[(x+1)%N][y]=3*con_a[(x-1)%N][y];
                                }
                                else{
                                    con_b[(x+1)%N][y]=con_b[x][y];
                                        // con_b[(x+1)%N][y]=3*con_a[(x-1)%N][y];
                                }
                            }
                            
                        }
                            
                        else if(n==2){
                            if(con[x][(y-1)%N]==1){
                                con[x][(y-1)%N]=2;con_a[x][(y-1)%N]=con_a[x][y];
                                if((double)rand()/RAND_MAX<mu*dt){
                                    if((double)rand()/RAND_MAX<0.5){
                                        con_b[x][(y-1)%N]=con_b[x][y]+v_mu;
                                        if(con_b[x][(y-1)%N]>b_max){
                                            con_b[x][(y-1)%N]=b_max;
                                            }
                                    }else{
                                        con_b[x][(y-1)%N]=con_b[x][y]-v_mu;
                                        if(con_b[x][(y-1)%N] < b_min){
                                            con_b[x][(y-1)%N]=b_min;
                                        }
                                                            
                                    }
                                    // con_b[x][(y-1)%N]=3*con_a[x][(y-1)%N];
                                }
                                else{
                                    con_b[x][(y-1)%N]=con_b[x][y];
                                    // con_b[x][(y-1)%N]=3*con_a[x][(y-1)%N];
                                }
                            }
                            
                        }
                        else{
                            if(con[x][(y+1)%N]==1){
                                con[x][(y+1)%N]=2;con_a[x][(y+1)%N]=con_a[x][y];
                                if((double)rand()/RAND_MAX<mu*dt){
                                    if((double)rand()/RAND_MAX<0.5){
                                        con_b[x][(y+1)%N]=con_b[x][y]+v_mu;
                                        if(con_b[x][(y+1)%N]>b_max){
                                            con_b[x][(y+1)%N]=b_max;
                                            }
                                    }else{
                                        con_b[x][(y+1)%N]=con_b[x][y]-v_mu;
                                        if(con_b[x][(y+1)%N] < b_min){
                                            con_b[x][(y+1)%N]=b_min;
                                        }
                                                            
                                    }
                                    // con_b[x][(y-1)%N]=3*con_a[x][(y-1)%N];
                                }
                                else{
                                    con_b[x][(y+1)%N]=con_b[x][y];
                                    // con_b[x][(y-1)%N]=3*con_a[x][(y-1)%N];
                                }
                            }
                            
                        }
                    }
                }
                else if(pr < con_b[x][y]*dt+con_a[x][y]*dt){
                    con[x][y]=1;
                    con_a[x][y]=0;
                    con_b[x][y]=0;
                }
                            
                // else{
                //     if((double)rand()/RAND_MAX<mu*dt){
                //             if((double)rand()/RAND_MAX<0.5){
                //                 con_a[x][y]=con_a[x][y]+v_mu;
                //                 if(con_b[x][y]>b_max){
                //                     con_b[x][y]=b_max;
                //                 }
                                            
                //             }
                //             else{
                //                 con_b[x][y]=con_b[x][y]-v_mu;
                //                 if(con_b[x][y]<b_min){
                //                     con_b[x][y]=b_min;
                //                 }
                                            
                //             }
                //         // con_b[x][y]=3*con_a[x][y];
                //     }
                                    
                                
                // }    
                            
                            
                            
                            
                            
            }

                        



        }
        ss=0,ii=0,ee=0,oo=0;
        for(l=0;l<N;l++){
            for(m=0;m<N;m++){
                if(con[l][m]==2){
                    ii=ii+1;
                }
                else if(con[l][m]==1){
                    ss=ss+1;
                }
                else if(con[l][m]==0){
                    oo=oo+1;
                }
                

            }

        }
        S[j]=ss;I[j]=ii;O[j]=oo;
                    




        vir=0.0,ni=0,ne=0;
        for(l=0;l<N;l++){
            for(m=0;m<N;m++){
                if(con[l][m]==2){
                    ni=ni+1;
                    vir=vir+con_b[l][m];
                                
                }
                            

            }

        }
                    
                    
        vir=(double)vir/ni;
    x[j]=j;//printf("%d\n",ni);
    y[j]=vir;

    
    // if(f_size==0){
    //     y1[j]=y1[j]+vir;
                    
    // }
    // else if(f_size==1) {y2[j]=y2[j]+vir;}
    // else if(f_size==2){y3[j]=y3[j]+vir;}else if(f_size==3){y4[j]=y4[j]+vir;}else if(f_size==4){y5[j]=y5[j]+vir;}else if(f_size==5){y6[j]=y6[j]+vir;}else if(f_size==6){y7[j]=y7[j]+vir;}//else if(i==7){y8[j]=ns;}
                

    }
    
    //SIS_betaの図
    data_file1="SIS_beta.dat";
    //     // data_file1="outadde_f.dat";
         data1=fopen(data_file1,"w");
         for(l=0;l<T;l++){
             fprintf(data1,"%d\t%f\n",l,y[l]);
         }
         fclose(data1);
         gp=popen("gnuplot -persist","w");
         fprintf(gp,"set terminal png\n");
         // fprintf(gp,"set logscale\n");
         // fprintf(gp,"set output 'Addefunction_f_%2f_P_%2f.png'\n",F_list[0],P_list[0]);
         fprintf(gp,"set output 'SIS_beta_%2f.png'\n",a);

        
         fprintf(gp,"set xrange [0:%d]\n",T);
         fprintf(gp,"set xlabel 'T'\n");
         fprintf(gp,"set yrange [0:%f]\n",b_max+10.0);//5.0
         fprintf(gp,"set ylabel 'beta'\n");
        

         fprintf(gp,"plot \'%s\'using 1:2 with lines linetype 1 linecolor rgb 'red'\n",data_file1);
         // fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"f=%f \"",data_file1,F_list[0]);
         // for(l=1;l<6;l++){
         //     fprintf(gp,",\'%s\' using 1:%d with lines linetype %d title \"f=%f \"",data_file1,l+2,l+1,F_list[l]);
         // }
         // fprintf(gp,",\'%s\' using 1:%d with lines linetype %d title \"f=%f \"\n",data_file1,9,8,F_list[6]);
        
        
        
         pclose(gp);



    pr=0.0;
    for(l=T-50;l<T;l++){
        pr=pr+y[l];
    }
    pr=pr/50.0;
    
    return pr;
    
}




