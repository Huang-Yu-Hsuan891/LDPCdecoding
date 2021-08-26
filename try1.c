#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define SNR 1.322513
#define SIGMA 0.869561
unsigned long long SEED = 3122891; // have value
//unsigned long long SEED = 3881111387891;
unsigned long long RANV;
int RANI = 0;
#define N 8000
#define RC 4000
#define DV 3
#define DC 6

int M[N][DV];
int L[RC][DC];
double q[RC][DC], u[RC][DC];
double  Lj[N];
//double log_table[8] = {0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, 0};
double Ranq1();
void normal(double sig, double *n1, double *n2);
/*double sgn(double L);
double minabs(double L1, double L2);
double triangle(double L1, double L2);*/
double bottom_up(int i, int j);
double top_down(int i, int j);
double app(int j);
double CHK(double L1, double L2);

int n, rc, dv, dc;

int main(){
    // for declaration
    int i, j, k, m;                  //for counting
    //int n, rc;              // n is column and rc is row
    //int dv,dc;              // dv: column have #1 and dc: row have #1
    int a;                  // no need
    int *codarray;          //codeword, Dynamic memory allocation, length = codarraylen
    int codarraylen = n;
    int *codearray;         //0->1; 1->-1,  Dynamic memory allocation, length = codearraylen
    int codearraylen = n;
    double *outp;           // codeword + noise, Dynamic memory allocation, length = outparray
    int outparray = n;
    int *output;            // result of interative algorithm decoding, Dynamic memory allocation, length = outputarray
    int outputarray = n;
    double *qj;
    int qjlen = n;
    int *checkbit;
    int checkbitlen = rc;
    int step;               // test 6 EbN0(dB) result
    int s = 0;              // receive 100 error block
    long long num = 0;            // do compute block
    int error;              // error block has the number of error
    int totalerror=0;       // after receive 100 error block, how many total error 
    int restart = 0;        // compute check node = 0
    //double sigma;
    //double ebn0;
    //double ebn0s[3];
    //double bers[3];
    //double berscompare[3];
    double x, y;
    int stp;
    int noteq = 0;
    //double SIGMA;
    //double SNR;
    // declaration end
    //for (i = 0; i < 3; i++) bers[i] = 0;
    
    /*ebn0s[0] = 0.8279;
    ebn0s[1] = 1.214;
    ebn0s[2] = 1.584;
    berscompare[0] = 0.081;
    berscompare[1] = 0.01171;
    berscompare[2] = 1.696 * pow(10,-5); */
    // open file
    //ebn0 = ebn0s[2];
    //SIGMA = sqrt(1.0 / (/*log2(rc) **/ pow(10, ebn0/10)));
    //SNR = pow(10, ebn0/10);
    FILE *fpr;
    fpr=fopen("parity8000.txt","r");
    fscanf(fpr,"%d",&n);
    fscanf(fpr,"%d",&rc);
    printf("column = %d\n", n);
    printf("row = %d\n", rc);
    fscanf(fpr,"%d",&dv);
    fscanf(fpr,"%d",&dc);
    printf("dv = %d\n", dv);
    printf("dc = %d\n", dc); 
    for (i = 0; i < n; i++) fscanf(fpr,"%d",&a);
    for (i = 0; i < rc; i++) fscanf(fpr,"%d",&a);
    for (j = 0; j < n; j++) {
        for (i = 0; i < dv; i++) {
            fscanf(fpr,"%d",&M[j][i]);
            M[j][i] = M[j][i] -1;
            //printf("M[%d][%d] = %d; ", j, i, M[j][i]);
        }
        //printf("\n");
    }
    for (i = 0; i < rc; i++) {
        for (j = 0; j < dc; j++) {
            fscanf(fpr,"%d",&L[i][j]);
            L[i][j] = L[i][j] -1;
            //printf("L[%d][%d] = %d; ", i, j, L[i][j]);
        }
        //printf("\n");
    }
    fclose(fpr);
    // close file

    // for CODE part
    codarraylen = n;
    codarray = (int *)malloc(codarraylen * sizeof(int));
    if( codarray == NULL ) {
        // 無法取得記憶體空間
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    codearraylen = n;
    codearray = (int *)malloc( codearraylen * sizeof(int) );
    if( codearray == NULL ) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    outparray = n;
    outp = (double *)malloc( outparray * sizeof(double) );
    if (outp == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    outputarray = n;
    output = (int *)malloc( outputarray * sizeof(int) );
    if (output == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    qjlen = n;          //qjlen = n = 816
    qj = (double *)malloc(qjlen * sizeof(double));
     
    checkbitlen = rc;   // checkbitlen = rc = 408
    checkbit = (int *)malloc(checkbitlen * sizeof(int));

    //printf("computlen = %d; comput1len = %d; comput2len = %d\n" ,computlen,comput1len,comput2len);
    /*for (step =1; step < 2; step++) {*/
    s = 0;
    num = 0;
    totalerror = 0;
    while (s < 1/*s < 100*/) {
        num++;                                  // compute the number of transmit block 
        //printf("cod = \n");                     // pretend encoder
        for (i = 0; i < n; i++) {
            codarray[i] = 0;
            //printf("%d ", codarray[i]);
        }
        //printf("code = \n");                    //input to AWGN channel normalized to +-1
        for (i = 0; i < n; i++) {
            if (codarray[i] == 0) codearray[i] = 1;
            else codearray[i] = -1;
            //printf("%d ", codearray[i]);
        }
            
        for(i = 0; i < n; i = i + 2) {
            normal(SIGMA, &x, &y);
            outp[i] = codearray[i] + x;
            outp[i + 1] = codearray[i + 1] + y;
        }
        for(i = 0; i < n; i++) {
            if (outp[i] >= 0) output[i] = 0;
            else output[i] = 1;
        }
        int row_sum;
        for(i=0;i<rc;i++){
            row_sum = 0;
            for(j=0;j<dc;j++){
                row_sum = row_sum + output[L[i][j]];
            }
            if((row_sum%2)==1){
                noteq = 1;
                row_sum = 0;
                break;
            }
        }
        //printf("Lj = ");
        for(i = 0; i < n; i++) {
            Lj[i] = 2 * SNR * outp[i];     //  0.5 * 1.2544 = Es/N0
                // printf("Lj[%d] = %g   ", i, Lj[i]);
        }
            // the interative decoding algotrithm
        for (j = 0; j < rc; j++) {                               // initialization
            for (i = 0; i < dc; i++) {
                    q[j][i] = Lj[L[j][i]];
            }  
        }
        for (k = 0; k < 100/*k < 100*/ && noteq!=0; k++) {         // message passing, for predetermined threshold = 100
            restart = 0;   
            noteq = 0;
            // bottom_up
            for (i = 0; i < rc; i++) {
                for (j = 0; j < dc; j++) {
                    u[i][j] = bottom_up(i,j);
                    if (i == 3999) printf("u[%d][%d] = %g; \n", i, j, u[i][j]);
                }
                //printf("\n");
            }

            // top-down
            for (j = 0; j < rc; j++) {
                for (i = 0; i < dc; i++) {
                    q[j][i] = top_down(j,i);
                    //printf("q[%d][%d] = %g; ", j, i, q[j][i]);
                }
                //printf("\n");
            }

            // decision
             //printf("output = ");
             for (j = 0; j < n; j++) {
                qj[j] = app(j);
                   
            }
            for (j = 0; j < n; j++) {
                if(qj[j] >= 0) output[j] = 0;
                else output[j] = 1;
            }

            // to check Hx=0
            //printf("checkbit = ");          
            for (i = 0; i < rc; i++) {
                row_sum = 0;
                 for (j = 0; j < 6; j++) {
                    row_sum = row_sum + output[L[i][j]];
                }
                if((row_sum%2) != 0){
                for(int k=0;k<n;k++)
                    if(output[k] != codarray[k]) noteq = noteq + 1;
                    break;
                }        
            }
        }
        printf("s = %d; k[%lld] = %d\n", s, num, k);    
        totalerror += noteq;

        if (noteq != 0) {
            s++;
        }
    }
    double ber;
    ber = (double)totalerror / (num * n);
    printf("totalerror = %d\n", totalerror);
    printf("BER = %g\n", ber);
 
    /*}*/

    // CODE  end

    // open write file
    /*FILE *outfp;
    
    outfp = fopen("result8000.txt","w");
    for (i = 0; i < 3; i++) {
         fprintf(outfp,"%g ",ebn0s[i]);
         fprintf(outfp,"%g ",berscompare[i]);
         fprintf(outfp,"%g ",bers[i]);
         fprintf(outfp,"\n");
    }
    fclose(outfp);*/
    // close write file

    // free dynametic
    free(codarray);
    free(codearray);
    free(outp);
    free(output);
    return 0;
}
void normal(double sig, double *n1, double *n2)
{   
    double x1,x2;
    double s;
    //printf("sigma = %g\n", sig);
    do{
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);
    *n1 = sig * x1 * sqrt((-2.0 * log(s))/ s);
    *n2 = sig * x2 * sqrt((-2.0 * log(s))/ s);
    
}

double Ranq1() {
    if ( RANI == 0 ){
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;

    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

/*double sgn (double L){
    if (L > 0) return 1.0;
    else if (L == 0) return 0.0;
    else return -1.0;
}

double minabs(double L1, double L2) {
    if(L1 <= 0) L1 = (-1) * L1;
    else L1 = L1;
    if(L2 <= 0) L2 = (-1) * L2;
    else L2 = L2;
    if(L1>=L2) return L2;
    else return L1;
}

double triangle(double L1, double L2) {
    double temp1, temp2;
    double ope1,ope2;
    double answer;
    ope1 = L1 + L2;
    ope2 = L1 - L2;
    if (ope1 <= 0) ope1 = ope1;
    else ope1 = (-1) * ope1;
    if (ope2 <= 0) ope2 = ope2;
    else ope2 = (-1) * ope2;
    temp1 = 1 + exp(ope1);
    temp2 = 1 + exp(ope2);
    answer = log(temp1 / temp2);
    return answer;
}*/
double bottom_up(int a, int b)
{
    double u = 0;
    for(int j=1; j<dc-1; j++){
        if(j == 1) {
            u = CHK(q[a][(b+j)%dc] , q[a][(b+j+1)%dc]);
        }
        else u = CHK(u, q[a][(b+j+1)%dc]);
    }
    return u;
}
double top_down(int a, int b)
{
    double q = 0;

    for(int i=0; i<dv; i++){
        for(int j=0; j<dc; j++){
            if(L[M[L[a][b]][i]][j] == L[a][b]){
                q = q + u[M[L[a][b]][i]][j];
            }
        }
    }
    q = q - u[a][b] + Lj[L[a][b]];
    return q;
}
double app(int col)
{
    double q = 0;

    for(int i=0; i<dv; i++){
        for(int j=0; j<dc; j++){
            if(L[M[col][i]][j] == col){
                q = q + u[M[col][i]][j];
            }
        }
    }
    q = q + Lj[col];
    return q;
}


double CHK(double L1, double L2) 
{   
    double sgn1, sgn2, min;

    if(L1>0) sgn1 = 1;
    else if(L1 == 0) sgn1 = 0;
    else sgn1 = -1;
    if(L2>0) sgn2 = 1;
    else if(L2 == 0) sgn2 = 0;
    else sgn2 = -1;
    if(fabs(L1) >= fabs(L2)) min = fabs(L2);
    else min = fabs(L1);
    double temp1, temp2;
    double ope1,ope2;
    double answer;
    ope1 = L1 + L2;
    ope2 = L1 - L2;
    if (ope1 <= 0) ope1 = ope1;
    else ope1 = (-1) * ope1;
    if (ope2 <= 0) ope2 = ope2;
    else ope2 = (-1) * ope2;
    temp1 = 1 + exp(ope1);
    temp2 = 1 + exp(ope2);
    answer = log(temp1 / temp2);

    return sgn1 * sgn2 * min + answer;
}