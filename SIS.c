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
#define M  (int)b_max+d+800 //b_max+1000
#define T  500
#define T_f 2000
#define z 4.0
#define Q 1    //平均を取るための試行回数

double infection(int con[][N],double con_a[][N],double con_b[][N],double P_list[],double F_list[],double c,double r,int P_size,int f_size,int x[],double y[],int S[],int I[],int O[],int g);


int main(){
    double con_a[N][N],P_list[2]={0.0},con_b[N][N],F_list[7]={0.5};//{0.1,0.5,1.0,2.0,3.0,10.0,100.0}
    double A_list[15]={0.1,0.2,0.3,0.4,0.5,0.7,0.9,1.0,2.0,3.0,4.0,6.0,7.0,8.0,10.0};//5.0,7.0,8.0,9.0,10.0,15.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0
    int i,j,k,l,m,h,q,con[N][N],Rsize,n_max,n_min,n;
    double a,b,c,pr,def,r;
    a=1.0;
    b=10.0;Rsize=15;c=0.0;
    FILE *gp,*data1,*data2,*data3,*data4,*data5,*data6;
    char *data_file1,*data_file2,*data_file3,*data_file4,*data_file5,*data_file6;

    int x[T],t[15000],I[15000],S[15000],E[15000],O[15000];
    double y1[T],y2[T],y3[T],y4[T],y5[T],y6[T],y7[T],y8[T],y9[T],y10[T],y11[T];

    srand((unsigned int)time(NULL));



    for(q=0;q<Q;q++){
        for(h=0;h<Rsize;h++){ 
            for(i=0;i<1;i++){
                //c=0.05+i*0.05;
                printf("r=%f\n",A_list[h]);
                a=A_list[h];
                /*conの初期化(part1)*/
                for(j=0;j<N;j++){
                    for(k=0;k<N;k++){
                        con[j][k]=0;
                        con_a[j][k]=0;
                        con_b[j][k]=0;
                    }
                }
                /*conの初期化(part2)*/
                for(j=0;j<N;j++){
                    for(k=0;k<N;k++){
                        n_max=2,n_min=1;
                        n =(int)rand()%(n_max-n_min+1)+n_min;
                        if(n==1){
                            con[j][k]=1;
                            con_a[j][k]=0;
                            con_b[j][k]=0;
                        }
                        else{
                            con[j][k]=2;
                            con_a[j][k]=a;
                            con_b[j][k]=b;

                        }
                    }
                }
                // for(j=0;j<N;j++){
                //     for(k=0;k<N;k++){
                //        if((double)rand()/RAND_MAX < 1.0/100){
                //             con[j][k]=2;
                //             con_a[j][k]=a;
                //             con_b[j][k]=b;
                //             // printf("侵入完了");
                //        } 
                //     }
                // }

            def=infection(con,con_a,con_b,P_list,F_list,c,A_list[h],0,0,x,y1,S,I,O,0);
               
            printf("でたよー");    
                /*変異株導入*/
                for(j=0;j<N;j++){
                    for(k=0;k<N;k++){
                        if(con[j][k]==2){
                            if((double)rand()/RAND_MAX < 0.001){
                                if((double)rand()/RAND_MAX < 0.5){
                                    // con_a[j][k]=a+v_muu;
                                    con_b[j][k]=con_b[j][k]+v_muu;
                                }
                                else{
                                    con_b[j][k]=con_b[j][k]-v_muu;
                                    if(con_b[j][k]<b_min){
                                        // con_a[j][k]=a;
                                        con_b[j][k]=b_min;
                                    }
                                }

                            }
                            else{
                                con_a[j][k]=1.0;
                                con_b[j][k]=b;
                            }
                        }
                    }
                }
                
                int ns,ni,ne,no;
                ns=0,ni=0,ne=0;
                for(l=0;l<N;l++){
                    for(m=0;m<N;m++){
                        if(con[l][m]==2){
                            ni=ni+1;
                        }
                        else if(con[l][m]==1){
                            ns=ns+1;
                        }
                        else{
                            no=no+1;
                        }
                        
                    }

                }
                data_file2="SISbef.dat";
                data2=fopen(data_file2,"w");
                for(l=0;l<T;l++){
                     fprintf(data2,"%d\t%d\t%d\t%d\n", l, S[l],I[l],O[l]);
                }
                
                fclose(data2);
                gp=popen("gnuplot -persist","w");
                fprintf(gp,"set terminal png\n");
                
                fprintf(gp,"set output 'SISbef_%.3f.png'\n",c);
            

                
                fprintf(gp,"set xrange [0:%d]\n",T);
                fprintf(gp,"set yrange [0:%d]\n",N*N); //5.0
                
                fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"P= %f  S \",\'%s\' using 1:3 with lines linetype 2 title \"I \",\'%s\' using 1:4 with lines linetype 3 title \"0\"\n",data_file2,P_list[0],data_file2,data_file2);
                
                pclose(gp);









                /*メインの平衡状態を作る*/
                printf("OK\n");
                y2[h]=infection(con,con_a,con_b,P_list,F_list,c,A_list[h],0,0,x,y1,S,I,O,1);

                data_file3="SISafter.dat";
                data3=fopen(data_file3,"w");
                for(l=0;l<T;l++){
                    fprintf(data3,"%d\t%d\t%d\t%d\n", l, S[l],I[l],O[l]);
                }
                fclose(data3);
                gp=popen("gnuplot -persist","w");
                fprintf(gp,"set terminal png\n");
                
                fprintf(gp,"set output 'SISafter_%.3f.png'\n",r);
            

                
                fprintf(gp,"set xrange [0:%d]\n",T);
                fprintf(gp,"set yrange [0:%d]\n",N*N); //5.0
                
                fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"P= %f  S \",\'%s\' using 1:3 with lines linetype 3 title \"I \",\'%s\' using 1:4 with lines linetype 4 title \"0\"\n",data_file3,P_list[0],data_file3,data_file3);
                
                pclose(gp);



                //ヒートマップ
            // data_file4="Map_I.dat";data_file5="Map_E.dat";data_file6="Map_S.dat";
            // data4=fopen(data_file4,"w");data5=fopen(data_file5,"w");data6=fopen(data_file6,"w");
            // for(l=0;l<N;l++){
            //     for(m=0;m<N;m++){
            //         if(con[l][m]==3){
            //             fprintf(data4,"%d\t%d\n",l,m);
            //         }
            //         else if(con[l][m]==2){
            //             fprintf(data5,"%d\t%d\n",l,m);
            //         }
            //         else if(con[l][m]==1){
            //             fprintf(data6,"%d\t%d\n",l,m);
            //         }
            //     }

               
            // }
                
            // fclose(data4);fclose(data5);fclose(data6);
            // gp=popen("gnuplot -persist","w");
            // fprintf(gp,"set terminal png\n");
                
            // fprintf(gp,"set output 'Map_%.1f_%3f.png'\n",P_list[i],F_list[h]);
            

                
            // fprintf(gp,"set xrange [0:%d]\n",N);
            // fprintf(gp,"set yrange [0:%d]\n",N); //5.0
            // fprintf(gp,"plot \'%s\' with points ps 0.4 pointtype 5 linecolor rgb 'red',\'%s\' with points ps 0.4 pointtype 5 linecolor rgb 'yellow',\'%s\' with points ps 0.4 pointtype 5 linecolor rgb 'blue'\n",data_file4,data_file5,data_file6);    
            
            // pclose(gp);

                

            }


        


        
        }
    }
    //data記載
    data_file5="SISdot.dat";
    data5=fopen(data_file5,"w");
    for(l=0;l<Rsize;l++){
        fprintf(data5,"%f\t%f\n", A_list[l],y2[l]);
    }
    fclose(data5);

    //閾値の線
    r=1.0;
    data_file4="SIS.dat";
    data4=fopen(data_file4,"w");
    for(l=0;l<10000;l++){
        r=(double)1.0+(double)0.05*l;
        pr=(z*(z*(1/r)*(1-2*ep)-ep*(z-1)))/(ep*z*(z-1)*(1/r)*(1/r)+(z*z*(1-ep)-z*(1-2*ep))*(1/r)-ep*(z-1)*(z-1));
        def=r*((z+1)*pow(1/r,2.0)+(-2*(1-ep)+2*z*(1-ep)+pow(z,2.0))*(1/r)+(-4*ep*(z-1)+2*z*(1-ep)-pow(z,2.0)*(1-2*ep))+pow(pow((1/r)+z-2,2.0)*(pow((z+1)*(1/r),2.0)+(4*ep-2*z+2*z*z*(3-2*ep))*(1/r)+(4*ep*ep+4*ep*z*(1-2*ep)+z*z*pow(1-2*ep,2.0))),0.5))/(2*(-(1/r)+(1-ep)*(z-1)*(z-2)));
        fprintf(data4,"%f\t%f\t%f\n", r,pr,def);
    }
    fclose(data4);printf("%f\n",r);
    gp=popen("gnuplot -persist","w");
        fprintf(gp,"set terminal png\n");
        fprintf(gp,"set logscale\n");
        // fprintf(gp,"set output 'Addefunction_f_%2f_P_%2f.png'\n",F_list[0],P_list[0]);
        fprintf(gp,"set output 'SIS.png'\n");

        
        fprintf(gp,"set xrange [1:%d]\n",500);
        fprintf(gp,"set xlabel 'alpha'\n");
        fprintf(gp,"set yrange [1:%d]\n",500);//5.0
        fprintf(gp,"set ylabel 'beta'\n");
        fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 linecolor rgb 'red',\'%s\' using 1:3 with lines linetype 1 linecolor rgb 'red',\'%s\' using 1:2 with points pointtype 1\n",data_file4,data_file4,data_file5);

    // for(i=0;i<T;i++){
    //     y1[i]=y1[i]/Q;y2[i]=y2[i]/Q;y3[i]=y3[i]/Q;y4[i]=y4[i]/Q;y5[i]=y5[i]/Q;y6[i]=y6[i]/Q;y7[i]=y7[i]/Q;
    // }

    //図の描画
        // data_file1="Linear.dat";
        // // data_file1="outadde_f.dat";
        // data1=fopen(data_file1,"w");
        // for(l=0;l<4;l++){
        //     fprintf(data1,"%f\t%f\n",y2[l]*c,y2[l]);
        // }
        // fclose(data1);
        // gp=popen("gnuplot -persist","w");
        // fprintf(gp,"set terminal png\n");
        // fprintf(gp,"set logscale\n");
        // // fprintf(gp,"set output 'Addefunction_f_%2f_P_%2f.png'\n",F_list[0],P_list[0]);
        // fprintf(gp,"set output 'Linear.png'\n");

        
        // fprintf(gp,"set xrange [0:%d]\n",100);
        // fprintf(gp,"set xlabel 'alpha'\n");
        // fprintf(gp,"set yrange [0:%d]\n",500);//5.0
        // fprintf(gp,"set ylabel 'beta'\n");
        // fprintf(gp,"f(x)=(8/7)*(1+x)\n");

        // fprintf(gp,"plot \'%s\' using 1:2 with points pointtype 1,f(x) with lines linetype 1 linecolor rgb 'red'\n",data_file1);
        // fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"f=%f \"",data_file1,F_list[0]);
        // for(l=1;l<6;l++){
        //     fprintf(gp,",\'%s\' using 1:%d with lines linetype %d title \"f=%f \"",data_file1,l+2,l+1,F_list[l]);
        // }
        // fprintf(gp,",\'%s\' using 1:%d with lines linetype %d title \"f=%f \"\n",data_file1,9,8,F_list[6]);
        
        
        
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
    
    data_file1="SIS_beta.dat";
        // data_file1="outadde_f.dat";
        data1=fopen(data_file1,"w");
        for(l=0;l<T;l++){
            fprintf(data1,"%d\t%f\n",l,y[l]);
        }
        fclose(data1);
        gp=popen("gnuplot -persist","w");
        fprintf(gp,"set terminal png\n");
        // fprintf(gp,"set logscale\n");
        // fprintf(gp,"set output 'Addefunction_f_%2f_P_%2f.png'\n",F_list[0],P_list[0]);
        fprintf(gp,"set output 'SIS_beta_%2f.png'\n",r);

        
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



