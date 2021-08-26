#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 8000
#define RC 4000
#define DV 3
#define DC 6

#define SNR 1.322513
#define SIGMA 0.869561

#define ERROR_BLOCK 50

int n, rc, dv, dc;
int M[N][DV];
int L[RC][DC];
double q[RC][DC], u[RC][DC];
double y[N], result[N], LLR[N];
int x[N];


double log_table[8] = {0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05, 0};

double bottom_up(int i, int j);
double top_down(int i, int j);
double app(int j);
double CHK(double L1, double L2);
double Ranq1();
void normal(double* n1, double* n2, double sigma);
unsigned long long SEED = 17;
//SEED must be an unsigned integer smaller than 4101842887655102017
unsigned long long RANV;
int RANI = 0;

double *n1, *n2;

int main(void){
    int skip;
    int round = 0;
    int noteq = 0;
    int all_not_equal = 0;
    int error_times = 0;
    int row_sum = 0;
    long long all_times = 0;
    int first_codeword = 1;
    //for(int i=0; i<N; i++) A_col[i] = -1;

    FILE *fin;
    FILE* outfile;
    outfile = fopen("out_decoding_8000.txt", "w");

    double v1, v2;
    v1 = 0;
    v2 = 0;
    n1 = &v1;
    n2 = &v2;

    //fin
    fin = fopen("parity8000.txt", "rt");
    if(fin == NULL){
        printf("Fail To Open File in1.txt!!\n");
        return 1;
    }
    fscanf(fin, "%d", &n);
    fscanf(fin, "%d", &rc);
    fscanf(fin, "%d", &dv);
    fscanf(fin, "%d", &dc);

    for(int i=0; i<n+rc; i++) fscanf(fin, "%d", &skip);

    for(int i=0; i<n; i++){
        for(int j=0; j<dv;j++){
            fscanf(fin, "%d", &M[i][j]);
            M[i][j] = M[i][j] -1;
        }
    }

    for(int i=0; i<rc; i++){
        for(int j=0; j<dc;j++){
            fscanf(fin, "%d", &L[i][j]);
            L[i][j] = L[i][j] -1;
        }
    }

    int codeword[N];

    while(error_times < ERROR_BLOCK){
        round = 0;
        noteq = 0;
        for(int j=0; j<N; j++) codeword[j] = 0;
        //進入檢測

        //noise
        for(int j=0; j<n; j=j+2){
            normal(n1, n2, SIGMA);
            if(codeword[j] == 0) y[j] = 1 + *n1;
            else y[j] = -1 + *n1;
            if(codeword[j+1] == 0) y[j+1] = 1 + *n2;
            else y[j+1] = -1 + *n2;
        }

        for(int j=0; j<n; j++){
            if(y[j] >= 0) x[j] = 0;
            else x[j] = 1;
        }

        for(int i=0;i<rc;i++){
            row_sum = 0;
            for(int j=0;j<dc;j++){
                row_sum = row_sum + x[L[i][j]];
            }
            if((row_sum%2)==1){
                noteq = 1;
                row_sum = 0;
                break;
            }
        }

        //decoding
        //initialization
        for(int j=0; j<n; j++) LLR[j] = 2 * SNR * y[j];

        for(int i=0; i<rc; i++){
            for(int j=0; j<dc;j++){
                q[i][j] = LLR[L[i][j]];
            }
        }

        //message passing
         while(noteq!=0 && round<100){
            noteq = 0;

            for(int i=0; i<rc; i++){
                for(int j=0; j<dc;j++){
                    u[i][j] = bottom_up(i,j);
                }
            }

            for(int i=0; i<rc; i++){
                for(int j=0; j<dc;j++){
                    q[i][j] = top_down(i,j);
                }
            }

            for(int j=0; j<n; j++) result[j] = app(j);

            //decision
            for(int j=0; j<n; j++){
                if(result[j]>=0) x[j] = 0;
                else x[j] = 1;
            }

            for(int i=0; i<rc; i++){
                row_sum = 0;
                for(int j=0;j<dc;j++){
                    row_sum = row_sum + x[L[i][j]];
                }
                if((row_sum%2) != 0){
                    for(int k=0;k<n;k++)
                        if(x[k] != codeword[k]) noteq = noteq + 1;
                    break;
                }
            }

            round++;
        }

        all_not_equal = all_not_equal + noteq;

        if(noteq != 0){
            error_times ++;
            printf("%d ", noteq);
        }

        all_times ++;
        first_codeword = 0;
    }

    printf("all not equal = %d ", all_not_equal);
    fprintf(outfile, "all not equal = %d\n", all_not_equal);

    printf("all times = %lld\n", all_times);
    fprintf(outfile, "all times = %lld\n", all_times);

    //善後
    fclose(outfile);
    fclose(fin);


    return 0;
}

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
    q = q - u[a][b] + LLR[L[a][b]];
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
    q = q + LLR[col];
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

    /*double x1 = fabs(L1 + L2);
    double x2 = fabs(L1 - L2);
    double log1, log2;

    if((0 <= x1) && (x1 < 0.196)) log1 = log_table[0];
    else if((0.196 <= x1) && (x1 < 0.433)) log1 = log_table[1];
    else if((0.433 <= x1) && (x1 < 0.71)) log1 = log_table[2];
    else if((0.71 <= x1) && (x1 < 1.05)) log1 = log_table[3];
    else if((1.05 <= x1) && (x1 < 1.508)) log1 = log_table[4];
    else if((1.508 <= x1) && (x1 < 2.252)) log1 = log_table[5];
    else if((2.252 <= x1) && (x1 < 4.5)) log1 = log_table[6];
    else if(4.5 <= x1) log1 = log_table[7];

    if((0 <= x2) && (x2 < 0.196)) log2 = log_table[0];
    else if((0.196 <= x2) && (x2 < 0.433)) log2 = log_table[1];
    else if((0.433 <= x2) && (x2 < 0.71)) log2 = log_table[2];
    else if((0.71 <= x2) && (x2 < 1.05)) log2 = log_table[3];
    else if((1.05 <= x2) && (x2 < 1.508)) log2 = log_table[4];
    else if((1.508 <= x2) && (x2 < 2.252)) log2 = log_table[5];
    else if((2.252 <= x2) && (x2 < 4.5)) log2 = log_table[6];
    else if(4.5 <= x2) log2 = log_table[7];*/
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

    return sgn1 * sgn2 * min + (answer);
}

void normal(double* n1,double* n2,double sigma)
{
    double x1, x2, s;
    do{
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2*x1 - 1;
        x2 = 2*x2 - 1;
        s = x1*x1 + x2*x2;
    }while(s>=1.0);
    *n1 = sigma * x1 * sqrt(-2*log(s)/s);
    *n2 = sigma * x2 * sqrt(-2*log(s)/s);
}

double Ranq1()
{
    if(RANI == 0){
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
