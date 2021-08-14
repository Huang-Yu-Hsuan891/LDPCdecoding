#include <stdio.h>
#include <math.h>
#include <stdlib.h>

unsigned long long SEED = 3122891; // have value
unsigned long long RANV;
int RANI = 0;

double Ranq1();
void normal(double sig, double *n1, double *n2);
int sgn(double L);
double minabs(double L1, double L2);
double triangle(double L1, double L2);

int main() {
    double x,y;             // for normal random variable
    int i, j, k, m;         // for counting
    int step;
    int cod[7];           // codeword
    int code[7];
    double outp[7];       // codeword+noise
    int output[7];     // out of channel
    int num = 0;                // do compute block

    double Lj[7] = {0};         // LLR
    double qij[3][7]= {0};   // from down to top
    //double uij[3][8000];
    double uij[7][3] = {0};   // from top to down
    
    int L[7][3];          // check L(i)
    int M[7][3];          // check M(j)
    int n,rc;   // n is column and rc is row
    int dv,dc;  // dv: column have #1 and dc: row have #1
    int a;      // no need
    int input;
    double tempqij[2]; 
    double tempuij;
    double temp1uij[3];
    double temp1qij;   
    double qj[7];
    int checkbit[7];
    int s = 0;      // receive 100 error block
    int restart = 0;
    int totalerror = 0; 
    int error;
    double sigma;
    double ebn0 = 0.9844;
    double ebn0s[1];
    //double bers[6];
    //double berscompare[6];
    int stp;
    int valL;
    int valL2; 
    int comput[7] = {0};
    int comput1[7] = {0};
    double app;
    double app1;  
    ebn0s[0] = 0.9844;

    FILE *fpr;
    fpr=fopen("parity77.txt","r");
    fscanf(fpr,"%d",&n);
    fscanf(fpr,"%d",&rc);
    //printf("column = %d\n", n);
    //printf("row = %d\n", rc);
    fscanf(fpr,"%d",&dv);
    fscanf(fpr,"%d",&dc);
    //printf("dv = %d\n", dv);
    //printf("dc = %d\n", dc); 
    for (i = 0; i < 7; i++) fscanf(fpr,"%d",&a);
    for (i = 0; i < 7; i++) fscanf(fpr,"%d",&a);

    for (j = 0; j < 7; j++) {
        fscanf(fpr,"%d",&M[j][0]);
        fscanf(fpr,"%d",&M[j][1]);
        fscanf(fpr,"%d",&M[j][2]);
        //printf("%d\n", M[j][2]);
    }
    
    for (i = 0; i < 408; i++) {
        fscanf(fpr,"%d",&L[i][0]);
        fscanf(fpr,"%d",&L[i][1]);
        fscanf(fpr,"%d",&L[i][2]);
    }
    fclose(fpr);

    while (s < 1) {
        num++;   // into this while loop compute add 1
        for (i = 0; i < 7; i++) cod[i] = 0;
        for (i = 0; i < 7; i++) {
            if (cod[i] == 0) code[i] = 1;
            else code[i]= -1;
        }
        ebn0 = ebn0s[0];
        sigma = sqrt(pow(pow(10, ebn0/10), -1));
        //printf("%g\n", sigma);
        for (i = 0; i < 7; i++) {
            normal(sigma, &x, &y);
            printf("x[%d] = %g; y[%d] = %g\n", i, x, i, y);
            outp[i] = code[i] + x;
            printf("outp[%d] = %g\n",i, outp[i]);
        }
        ebn0 = pow(10, ebn0/10);
        printf("ebno = %g\n",ebn0);
        printf("%g\n", log2(408));
        for(i = 0; i < 7; i++) {
            Lj[i] = log2(408) * 4 * 0.5 * ebn0 * outp[i];
            printf("Lj[%d] = %g\n", i, Lj[i]);
        }
        for (j = 0; j < 7; j++) {
            for (i = 0; i < 3; i++) {
                qij[i][j] = Lj[j];
                printf("qij[%d][%d] = %g\n", i, j, qij[i][j]);
            }
        }
        for (k = 0; k < 1/*100*/ && restart != 7; k++) {
            restart = 0;
            comput[7] = {0};
            comput1[7] = {0};
            for (i = 0; i < 2; i++) tempqij[i] = 0.0;
            for (i = 0; i < 7; i++) {
                    for (j = 0; j < 3; j++) {
                        for (m = 0; m < 2; m++) {
                            if (m < j) {
                                valL = L[i][m]-1;
                                tempqij[m] = qij[comput[valL]][valL];
                                //printf("qij[%d][%d] = %g\n",comput[valL], valL, qij[comput[valL]][valL]);
                                //printf("tempqij[%d] = %g \n", m, tempqij[m]);
                            } 
                            else if (m >= j) {
                                valL = L[i][m+1]-1;
                                tempqij[m] = qij[comput[valL]][valL];
                                //printf("qij[%d][%d] = %g\n",comput[valL], valL, qij[comput[valL]][valL]);
                                //printf("tempqij[%d] = %g \n", m, tempqij[m]);
                            }
                        }
                        tempuij = tempqij[0];
                        for(m = 1; m < 2; m++) {
                            //printf("tempuij = %g\n", tempuij);
                            //printf("tempqij[%d] = %g\n", m, tempqij[m]);
                            app = sgn(tempuij) * sgn(tempqij[m]) * minabs(tempuij,tempqij[m]);
                            app1 = triangle(tempuij,tempqij[m]);
                            //printf("app = %g\n", app);
                            //printf("app1 = %g\n", app1);
                            //tempuij = ((sgn(tempuij) * sgn(tempqij[m]) * minabs(tempuij,tempqij[m])) + triangle(tempuij,tempqij[m]));
                            tempuij = app + app1;
                            //printf("tempqij[%d] = %g\n", m, tempqij[m]);
                            //printf("sgn(tempqij[%d]) = %d\n", m, sgn(tempqij[m]));
                            //printf("minabs(tempuij,tempqij[m]) = %g\n", minabs(tempuij,tempqij[m]));
                            //printf("triangle(tempuij,tempqij[m]) = %g\n", triangle(tempuij,tempqij[m]));
                            //printf("tempuij[%d] = %g\n", m, tempuij); 
                        }
                        //valL2 = L[i][j]-1;
                        //uij[i][valL2] = tempuij;
                        uij[i][j] = tempuij;
                        printf("uij[%d][%d] = %g\n",i,j,uij[i][j]);
                    }
                    for (m = 0; m < 3; m++) {
                            comput[L[i][m] - 1] += 1;
                    }
            }
            for (i = 0; i < 7; i++) {
                for(j = 0; j < 3; j++) {
                    printf("uij[%d][%d] = %g\n",i,j,uij[i][j]);
                }
            }
            for(i = 0; i < 3; i++) temp1uij[i] = 0.0;
            for (j = 0; j < 7; j++) {
                for (i = 0; i < 3; i++) {
                    for (m = 0; m < 2; m++) {
                        if (m < i) { 
                            valL = M[j][m] - 1;
                            printf("%d\n", M[j][m] - 1);
                            temp1uij[m] = uij[valL][comput1[valL]]; 
                            printf("temp1uij[%d] = %g; uij[%d][%d] = %g; \n",m, temp1uij[m],valL, comput1[valL], uij[valL][comput1[valL]]);
                        }
                        else if (m >= i) {
                            valL = M[j][m + 1] - 1;
                            printf("%d\n", M[j][m + 1] - 1);
                            temp1uij[m] = uij[valL][comput1[valL]];
                            printf("temp1uij[%d] = %g; uij[%d][%d] = %g; \n",m, temp1uij[m],valL, comput1[valL], uij[valL][comput1[valL]]);
                        }
                    }
                    temp1uij[2] = Lj[j];
                    //valL2 = M[j][i] - 1;
                    qij[i][j] = temp1uij[0] + temp1uij[1] + temp1uij[2];
                }
                for (m = 0; m < 3; m++) {
                    comput1[M[j][m] - 1] += 1;
                }
            }
        }
        s++;
    }
}

void normal(double sig, double *n1, double *n2)
{   
    double x1,x2;
    double s;

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

int sgn (double L){
    if (L >= 0) return 1;
    else return -1;
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
}