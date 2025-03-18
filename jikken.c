#include <stdio.h>
#include <stdlib.h>
#include <time.h>
void sum(int *i,int *s,int *o);
int main(){
    int i=1,s=0,o=2;
    
    sum(&i,&s,&o);
    printf("i=%d",i);
}
void sum(int *i,int *s,int *o){
    *i=*i+1;
}