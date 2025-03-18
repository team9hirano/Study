#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define N 200//200
#define a_max 4.767 // 4.767
#define a_min 0.2
#define b_max 500.0 // 4.767
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
#define w 1.0
#define Q 1    //平均を取るための試行回数
//回復率

int main(){
    double J[N][N];//{0.1,0.5,1.0,2.0,3.0,10.0,100.0}
    double F[6],r[6]={0.0,0.0,0.0,0.0,0.0,0.0},A_list[T],B_list[T],R_list[T],D_list[T],A2_list[T],B2_list[T];//5.0,7.0,8.0,9.0,10.0,15.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0
    int i,j,k,l,m,h,q,con[N][N],Rsize,n_max,n_min,n,number,count,count_max;
    double x,y,s,A,B,C,f,G,H,I,beta,alpha,gamma,omega,m1,m2,m3,R,rhos,U,rhoi,V;
    
    FILE *gp,*data1,*data2,*data3,*data4,*data5,*data6;
    char *data_file1,*data_file2,*data_file3,*data_file4,*data_file5,*data_file6;

    //diseasefree境界線
    count=0;
    count_max=10;
    for(m=1;m<T;m++){
        for(i=0;i<number;i++){
            r[i]=0.0;
        }
        count=0;
        gamma=(double)m;
        // printf("alpha=%f\n",alpha);
        A2_list[m]=gamma;
        B2_list[m]=(double)b_max;
        for(h=1;h<=(int)b_max;h++){
            beta=(double)h;
            if(count > count_max){
                   
                break;
            }
            x=gamma/beta;
            U=0.5;
            alpha=1.0;
            R=1.0;
            rhos=0.3;
            rhoi=0.3;
            V=0.3;
            //beta=100.0;
            
            //gamma=5.0;
            omega=1.0;
        
            for(l=0;l<100000;l++){
                J[0][0]=-w;
                J[0][1]=-w-beta*x;
                J[0][2]=-beta*rhoi;
                J[0][3]=0;
                J[0][4]=0;

                J[1][0]=0;
                J[1][1]=-gamma+beta*x;
                J[1][2]=beta*rhoi;
                J[1][3]=0;
                J[1][4]=0;

                J[2][0]=beta*(1-1/z)*x*rhoi/(rhos*rhos);
                J[2][1]=-beta*(1-1/z)*x/rhos;
                J[2][2]=beta*(1-1/z)*V-beta*(1/z+(1-1/z)*(rhoi/rhos))-2*beta*x;
                J[2][3]=-w;
                J[2][4]=beta*(1-1/z)*x;

                J[3][0]=-beta*(1-1/z)*x*rhoi/(rhos*rhos);
                J[3][1]=beta*(1-1/z)*x/rhos;
                J[3][2]=-beta*(1-1/z)*(1-2*rhoi*x/rhos-V)+gamma+beta*(1-U);
                J[3][3]=-gamma-w-beta*x;
                J[3][4]=beta*(1-1/z)*x;

                J[4][0]=2*w*x*rhoi/(rhos*rhos);
                J[4][1]=-2*w*x/rhos-beta*(1-2/z)*x*x*V/rhos+V*w/rhos;
                J[4][2]=-2*w*rhoi/rhos-2*x*V*rhoi*beta*(1-2/z)/rhos;
                J[4][3]=0;
                J[4][4]=-2*w-rhoi*x*x*beta*(1-2/z)/rhos-w*(1/rhos-1-rhoi/rhos);






                F[0]=w*(1-rhos-rhoi)-beta*x*rhoi;
                F[1]=beta*rhoi*x-gamma*rhoi;
                F[2]=w*(1-U)+beta*(1-1/z)*x*V-beta*(1/z+(1-1/z)*rhoi/rhos)*x-beta*x*x;
                F[3]=-beta*(1-1/z)*(1-rhoi*x/rhos-V)*x-gamma*(U-x)+w*(1-U)+beta*x*(1-U);
                F[4]=2*w*(1-rhoi*x/rhos-V)-V*rhoi*x*x*beta*(1-2/z)/rhos-V*w*(1/rhos-1-rhoi/rhos);
                //if(h==0&&m==0&&l==0)printf("J11=%f,J12=%f,J21=%f,J22=%f,F[0]=%f,F[1]=%f\n",J[0][0],J[0][1],J[1][0],J[1][1],F[0],F[1]);
                //ガウスの消去法
                number=5;
                for(i=0;i<number;i++){
                    for(j=i+1;j<number;j++){
                        m1=J[j][i]/J[i][i];
                        for(k=i+1;k<number;k++){
                            J[j][k]=J[j][k]-m1*J[i][k];
                        
                        } 
                        F[j]=F[j]-m1*F[i];
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
                    
                    r[i]=(1/J[i][i])*(F[i]-m2);
                }
                
                if(fabs(r[0])<0.001&&fabs(r[1])<0.001&&fabs(r[2])<0.001&&fabs(r[3])<0.001&&fabs(r[4])<0.001){
                    //printf("l=%d\n",l);


                //     if(beta*x-gamma<0){   //r*B-d,beta*x-(d+alpha+gamma)
                // //printf("安定,beta=%f",beta);
                //         B2_list[m]=beta;
                //         //printf("beta=%f\n",B_list[m]);
                //         //if(m==0)printf("beta=%f\n",B_list[m]);
                        
                
                //     }else{
                         if(count<=count_max){
                             count=count+1;
                         }            
                //         }else{
                //             printf("beta=%f\n",B2_list[m]);
                //         }
                //     }
                    break;
                    
                }
                rhos=rhos-r[0];
                rhoi=rhoi-r[1];
                x=x-r[2];
                U=U-r[3];
                V=V-r[4];
                //if(h==0&&m==0&&l==0)printf("x=%f,U=%f,r[0]=%f,r[2]=%f,F[0]=%f,F[1]=%f\n",x,U,r[0],r[1],F[0],F[1]);
                // B=B-r[2];
                // C=C-r[3];
                // y=y-r[4];
                // s=s-r[5];

            }
        //printf("x=%f,A=%f,B=%f,C=%f\n",x,A,B,C);
        
            //printf("不安定");
            
            
        
        // printf("beta=%f\n",beta);
        // printf("x=%f,A=%f,B=%f,C=%f\n",x,A,B,C);
        }
        printf("gamma=%f,beta=%f,rhos=%f,rhoi=%f,x=%f,U=%f,V=%f\n",gamma,beta,rhos,rhoi,x,U,V);
    //    printf("beta=%f\n",B_list[m]);
    //    printf("x=%f,A=%f,B=%f,C=%f,y=%f,s=%f\n",x,A,B,C,y,s); 
    }



    //printf("B_list[0]=%f\n",B_list[0]);
    //S,Iなどの時間変遷図(after)
                //  data_file3="NewtonSIRSline.dat";
                //  data3=fopen(data_file3,"w");
                //  for(l=0;l<T;l++){
                //      fprintf(data3,"%f\t%f\n",A2_list[l],B2_list[l]);
                //  }
                //  fclose(data3);
                //  gp=popen("gnuplot -persist","w");
                //  //fprintf(gp,"set logscale\n");
                //  fprintf(gp,"set terminal png\n");
                
                //  fprintf(gp,"set output 'NewtonSIRSline_beta_gamma.png'\n");
            

                
                //  fprintf(gp,"set xrange [0:%d]\n",T);
                //  fprintf(gp,"set xlabel 'gamma'\n");
                //  fprintf(gp,"set yrange [0:%d]\n",T);
                //  fprintf(gp,"set ylabel 'beta'\n");
                //   //5.0
                
                //  //fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"a= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,A_list[h],data_file3,data_file3);
                // fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"disease free\"\n",data_file3);
                //  pclose(gp);






}