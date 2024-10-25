#include <stdio.h>
#include <stdlib.h>
#include <time.h>


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

void infection(int con[][N],double con_a[][N],double con_b[][N],double P_list[],double F_list[],double r,int P_size,int f_size,int x[],double y[],int S[],int I[],int O[],int g,double a);
void sum(int con[][N],int l,int m,int *i,int *s,int *o);

int main(){
    double con_a[N][N],P_list[2]={0.0},con_b[N][N],F_list[7]={0.5};//{0.1,0.5,1.0,2.0,3.0,10.0,100.0}
    double R_list[15]={3.0};
    int i,j,k,l,m,n,h,q,con[N][N],Rsize,count;//{2,3,4,5,7,9,13,17,20,30,50,60,80,100,140,180,200,300};
    double ramda,a,b,pr,def,c,q_ss,q_is,q_os,q_si,q_ii,q_oi,q_so,q_io,q_oo,q_sj,q_ij,q_jj,q_oj;
    int ns,ni,ne,no,nj,i_s,i_i,i_o,s_s,s_i,s_o,o_s,o_i,o_o,j_o,j_s,j_i,j_j;
    ns=0;ni=0;ne=0;i_s=0;i_i=0;i_o=0;s_s=0;s_i=0;s_o=0;o_s=0;o_i=0;o_o=0;

    a=1.0;c=0.0;
    b=5.0;Rsize=1;//24;
    FILE *gp,*data1,*data2,*data3,*data4,*data5,*data6;
    char *data_file1,*data_file2,*data_file3,*data_file4,*data_file5,*data_file6;

    int x[T],t[15000],I[15000],S[15000],E[15000],O[15000];
    double y1[T],y2[T],y3[T],y4[T],y5[T],y6[T],y7[T],y8[T],y9[T],y10[T],y11[T],a_low[T],a_hi[T];

    srand((unsigned int)time(NULL));

    printf("%d",Rsize);

    data_file2="Disbef.dat";
    
    for(q=0;q<Q;q++){
        for(h=0;h<Rsize;h++){
            printf("R=%f\n",R_list[h]); 
            for(i=0;i<20;i++){
                //b=400;
                count=0;
                y1[i]=0.01+0.05*i;
                a=0.01+0.05*i;
                b=3*a;
                // printf("P=%f\n",P_list[0]);
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
                        int n;
                        n =(int)rand()%3;
                        if(n==0){
                            con[j][k]=0;
                            con_a[j][k]=0;
                            con_b[j][k]=0;
                        }
                        else if(n==1){
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
                //         int n;
                //         n =(int)rand()%2;
                //         if(n==0){
                //             con[j][k]=0;
                //             con_a[j][k]=0;
                //             con_b[j][k]=0;
                //         }
                //         else if(n==1){
                //             con[j][k]=1;
                //             con_a[j][k]=0;
                //             con_b[j][k]=0;

                //         }
                //     }
                // }
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
            //printf("yes\n");
            infection(con,con_a,con_b,P_list,F_list,R_list[h],0,0,x,y1,S,I,O,0,a);
            //printf("yes\n");
            int ns,ni,ne,no,i_s,i_i,i_o,s_s,s_i,s_o,o_s,o_i,o_o;
            ns=0;ni=0;ne=0;no=0;i_s=0;i_i=0;i_o=0;s_s=0;s_i=0;s_o=0;o_s=0;o_i=0;o_o=0;
            for(l=0;l<N;l++){
                for(m=0;m<N;m++){
                    //i_s=0;i_i=0;i_o=0;s_s=0;s_i=0;s_o=0;o_s=0;o_i=0;o_o=0;
                    if(con[l][m]==2){
                        ni=ni+1;
                        sum(con,l,m,&i_i,&i_s,&i_o);

                    }else if(con[l][m]==1){
                        ns=ns+1;
                        
                    }else{
                        no=no+1;
                    }
                    
                        
                }
            }
            printf("i\n");
            printf("ns=%d,no=%d,ni=%d\n",ns,no,ni); 
            y3[i]=b;y4[i]=ni;
            if(ni==0){
                q_si=0.0;
            }else{q_si=(double)i_s/(4*ni);}
            printf("q_si=%f\n",q_si);
            /*変異株導入*/
            for(l=1;l<100;l++){ 
                nj=0;j_o=0;j_s=0;j_i=0;j_j=0;ramda=0;q_sj=0.0;
                for(j=0;j<N;j++){
                    for(k=0;k<N;k++){
                        if(con[j][k]==2){
                            if((double)rand()/RAND_MAX < 0.1){
                                //con[j][k]=3;
                                con_a[j][k]=l*0.01;
                                con_b[j][k]=3*con_a[j][k];
                                sum(con,j,k,&j_i,&j_s,&j_o);nj=nj+1;//printf("j_s=%d\n",j_s);
                                // if((double)rand()/RAND_MAX < 0.5){
                                //     // con_a[j][k]=a+v_muu;
                                //     con_b[j][k]=con_b[j][k]+v_muu;
                                // }
                                // else{
                                //     con_b[j][k]=con_b[j][k]-v_muu;
                                //     if(con_b[j][k]<b_min){
                                //         // con_a[j][k]=a;
                                //         con_b[j][k]=b_min;
                                //     }
                                // }

                            }
                            // else{
                            //     con_a[j][k]=a;
                            //     con_b[j][k]=b;
                            // }
                        }
                    }
                }
                if(nj==0){
                q_sj=0.0;
                }else{q_sj=(double)j_s/(4*nj);}
                ramda=(d+a)/b-(d+l*0.01)/(3*l*0.01)+(1-P_list[0])*(q_sj-q_si);//printf("nj=%d,q_sj=%f,ramda=%f\n",nj,q_sj,ramda);
                if(ramda>0){
                    if(count==0){
                        a_low[i]=l*0.01;
                        a_hi[i]=l*0.01;
                        count=count+1;
                    }
                    else{
                        a_hi[i]=l*0.01;
                    }
                }


            }

            

            





                /*メインの平衡状態を作る*/
                



                //ヒートマップ
            

                

        }
         
        

        



           

        


        
        }
    

    // for(i=0;i<T;i++){
    //     y1[i]=y1[i]/Q;y2[i]=y2[i]/Q;y3[i]=y3[i]/Q;y4[i]=y4[i]/Q;y5[i]=y5[i]/Q;y6[i]=y6[i]/Q;y7[i]=y7[i]/Q;
    // }

    //図の描画
        data_file1="PIPs.dat";
        // data_file1="outadde_f.dat";
        data1=fopen(data_file1,"w");
        for(l=0;l<20;l++){
            fprintf(data1,"%f\t%f\t%f\n", y1[l],a_low[l],a_hi[l]);
        }
        fclose(data1);
        gp=popen("gnuplot -persist","w");
        fprintf(gp,"set terminal png\n");
        //fprintf(gp,"set logscale\n");
        // fprintf(gp,"set output 'Addefunction_f_%2f_P_%2f.png'\n",F_list[0],P_list[0]);
        fprintf(gp,"set output 'PIPs.png'\n");

        
        fprintf(gp,"set xrange [0:%f]\n",1.0);
        fprintf(gp,"set xlabel 'alpha_I'\n");
        fprintf(gp,"set yrange [0:%f]\n",1.0);//5.0
        fprintf(gp,"set ylabel 'alpha_J'\n");
        // fprintf(gp,"f(x)=2*x/(x-1)\n");

        fprintf(gp,"plot \'%s\'using 1:2 with lines linetype 1 linecolor rgb 'red',\'%s\' using 1:3 with lines linetype 1 linecolor rgb 'red'\n",data_file1,data_file1);
        // fprintf(gp,"plot \'%s\' using 1:2 with lines linetype 1 title \"f=%f \"",data_file1,F_list[0]);
        // for(l=1;l<6;l++){
        //     fprintf(gp,",\'%s\' using 1:%d with lines linetype %d title \"f=%f \"",data_file1,l+2,l+1,F_list[l]);
        // }
        // fprintf(gp,",\'%s\' using 1:%d with lines linetype %d title \"f=%f \"\n",data_file1,9,8,F_list[6]);
        
        
        
        pclose(gp);



    }
}







// void infection(int con[][N],double con_a[][N],double con_b[][N],double P_list[],double F_list[],int P_size,int f_size,int x[],double y[],int S[],int E[],int I[],int O[],int a);
void infection(int con[][N],double con_a[][N],double con_b[][N],double P_list[],double F_list[],double r,int P_size,int f_size,int x[],double y[],int S[],int I[],int O[],int g,double a){
    int i,j,k,m,l,MM;
    double v_mu,pr;
    double vir;
    int ni,ne,ns,xx,yy;
    int ss,ee,ii,oo;
    
    FILE *gp,*data1,*data2,*data3,*data4,*data5,*data6;
    char *data_file1,*data_file2,*data_file3,*data_file4,*data_file5,*data_file6;
    srand((unsigned int)time(NULL));

    if(g==0){v_mu=0;MM=100;}else{v_mu=v_muu;MM=M;}//MM=100
    for(j=0;j<T;j++){ //M
        //printf("OK");
        for(k=0;k<MM*N*N;k++){
            int x,y,n;
            x=(int)rand()%N,y=(int)rand()%N;
            n=(int)rand()%4;pr=(double)rand()/RAND_MAX;
            if(con[x][y]==1){
                            
                if(pr < r*dt){
                    if((double)rand()/RAND_MAX<P_list[P_size]){
                        xx=(int)rand()%N,yy=(int)rand()%N;
                        if(con[xx][yy]==0){
                            con[xx][yy]=1;
                        }
                    }
                    else{
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
                else if(pr < r*dt+d*dt){
                    con[x][y]=0;
                }
            }
            
            else if(con[x][y]==2){
                
                xx=(int)rand()%N,yy=(int)rand()%N;
                pr=(double)rand()/RAND_MAX;
                            
                if(pr < con_b[x][y]*dt){
                    if( (double)rand()/RAND_MAX < P_list[P_size]){
                        if(con[xx][yy]==1){
                            con[xx][yy]=2;con_a[xx][yy]=a;
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
                                //con_b[xx][yy]=3*con_a[xx][yy];
                                con_b[xx][yy]=con_b[x][y];
                            }
                            else{
                                //con_a[xx][yy]=a;
                                con_b[xx][yy]=con_b[x][y];
                            }
                                        
                        }
                    }
                    else{
                        if(n==0){
                            if(con[(x-1)%N][y]==1){
                                con[(x-1)%N][y]=2;con_a[(x-1)%N][y]=a;
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
                                    con_b[(x-1)%N][y]=con_b[x][y];
                                }
                                else{
                                    con_b[(x-1)%N][y]=con_b[x][y];
                                        // con_b[(x-1)%N][y]=3*con_a[(x-1)%N][y];
                                }
                            }
                        }
                        else if(n==1){
                            if(con[(x+1)%N][y]==1){
                                con[(x+1)%N][y]=2;con_a[(x+1)%N][y]=a;
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
                                    con_b[(x+1)%N][y]=con_b[x][y];
                                    }
                                    else{
                                        con_b[(x+1)%N][y]=con_b[x][y];
                                        // con_b[(x+1)%N][y]=3*con_a[(x-1)%N][y];
                                    }
                            }
                        }
                            
                        else if(n==2){
                            if(con[x][(y-1)%N]==1){
                                con[x][(y-1)%N]=2;con_a[x][(y-1)%N]=a;
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
                                    con_b[x][(y-1)%N]=con_b[x][y];
                                }
                                else{
                                    con_b[x][(y-1)%N]=con_b[x][y];
                                    // con_b[x][(y-1)%N]=3*con_a[x][(y-1)%N];
                                }
                            }
                        }
                        else{
                            if(con[x][(y+1)%N]==1){
                                con[x][(y+1)%N]=2;con_a[x][(y+1)%N]=a;
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
                                    //con_b[x][(y-1)%N]=con_b[x][(y-1)%N];
                                    con_b[x][(y+1)%N]=con_b[x][y];
                                }
                                else{
                                    con_b[x][(y+1)%N]=con_b[x][y];
                                    // con_b[x][(y-1)%N]=3*con_a[x][(y-1)%N];
                                }
                            }
                        }
                    }
                }
                else if(pr < con_b[x][y]*dt+(d+con_a[x][y])*dt){
                    con[x][y]=0;
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
        ss=0,ii=0,ee=0,oo=0;vir=0.0;
        for(l=0;l<N;l++){
            for(m=0;m<N;m++){
                if(con[l][m]==2){
                    ii=ii+1;
                    vir=vir+con_b[l][m];
                }
                else if(con[l][m]==1){
                    ss=ss+1;
                }
                else if(con[l][m]==0){
                    oo=oo+1;
                }
                

            }

        }
        S[j]=ss;I[j]=ii;O[j]=oo;x[j]=j;
        // if(ss==0 && ii==0){
        //      return 0; 
        // }            
        // else{vir=(double)vir/ii;}//vir=(double)vir/ii;
        // y[j]=vir;
    
                    




    //     vir=0.0,ni=0,ne=0;
    //     for(l=0;l<N;l++){
    //         for(m=0;m<N;m++){
    //             if(con[l][m]==2){
    //                 ni=ni+1;
    //                 vir=vir+con_b[l][m];
                                
    //             }
    //             else if(con[l][m]==1){
    //                 ns=ns+1;
    //             }
                            

    //         }

    //     }
                    
    //     if(ni==0){
    //          return 0; 
    //     }            
    //     else {vir=(double)vir/ni;}
    //     //vir=(double)vir/ni;
    // x[j]=j;
    // y[j]=vir;

    
    // if(f_size==0){
    //     y1[j]=y1[j]+vir;
                    
    // }
    // else if(f_size==1) {y2[j]=y2[j]+vir;}
    // else if(f_size==2){y3[j]=y3[j]+vir;}else if(f_size==3){y4[j]=y4[j]+vir;}else if(f_size==4){y5[j]=y5[j]+vir;}else if(f_size==5){y6[j]=y6[j]+vir;}else if(f_size==6){y7[j]=y7[j]+vir;}//else if(i==7){y8[j]=ns;}
                

    }
    
        
        
       




    // pr=0.0;
    // for(l=T-50;l<T;l++){
    //     pr=pr+y[l];
    // }
    // pr=pr/50.0;return pr;
}


void sum(int con[][N],int l,int m,int *i,int *s,int *o){
    int n;                    
                    
    for(n=0;n<4;n++){
        if(n==0){
            if(con[(l-1)%N][m]==2){
                *i=*i+1;
            }
            else if(con[(l-1)%N][m]==1){
                *s=*s+1;
            }
            else{
                *o=*o+1;
            }
        }else if(n==1){
            if(con[(l+1)%N][m]==2){
                *i=*i+1;
            }
            else if(con[(l+1)%N][m]==1){
                *s=*s+1;
            }
            else{
                *o=*o+1;
            }
        }else if(n==2){
            if(con[l][(m-1)%N]==2){
                *i=*i+1;
            }
            else if(con[l][(m-1)%N]==1){
                *s=*s+1;
            }
            else{
                *o=*o+1;
            }
        }else{
            if(con[l][(m+1)%N]==2){
                *i=*i+1;
            }
            else if(con[l][(m+1)%N]==1){
                *s=*s+1;
            }
            else{
                *o=*o+1;
            }
        }
    }



}