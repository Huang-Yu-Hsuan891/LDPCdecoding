#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define SNR 1.440124
#define SIGMA 0.833297
unsigned long long SEED = 3122891; // have value
//unsigned long long SEED = 3881111387891;
unsigned long long RANV;
int RANI = 0;

double Ranq1();
void normal(double sig, double *n1, double *n2);
double CHK(double L1, double L2);
int main(){
    // for declaration
    int i, j, k, m;                  //for counting
    int n, rc;              // n is column and rc is row
    int dv,dc;              // dv: column have #1 and dc: row have #1
    int a;                  // no need
    int **L = NULL;                 // check node connect 6 variable nodes
    int Llenrow = rc;        // rc = 408
    int Llencolumn = dc;     // dc = 6
    int **M = NULL;
    int Mlenrow = n;       // n = 816
    int Mlencolumn = dv;   //dv = 3
    int *codarray;          //codeword, Dynamic memory allocation, length = codarraylen
    int codarraylen = n;
    int *codearray;         //0->1; 1->-1,  Dynamic memory allocation, length = codearraylen
    int codearraylen = n;
    double *outp;           // codeword + noise, Dynamic memory allocation, length = outparray
    int outparray = n;
    int *output;            // result of interative algorithm decoding, Dynamic memory allocation, length = outputarray
    int outputarray = n;
    double *Lj;             // LLR
    int Ljlen = n;
    double **qij = NULL;    // from down to top
    int qijrow = dv;        // qijrow = dv = 3
    int qijcolumn = n;      // qijcolumn = n = 816
    double **uij = NULL;    // from top to down
    int uijrow = rc;        // uijrow = rc = 408
    int uijcolumn = dc;     // uijcloumn = dc 3
    double tempqij[5]; 
    double tempuij;
    double temp1uij[3];
    double temp1qij; 
    int *comput;
    int computlen = n;      // computlen = 816
    int *comput1;
    int comput1len = rc;    // comput1len = 408
    int *comput2;
    int comput2len = rc;    //  comput2len = rc = 408
    double *qj;
    int qjlen = n;
    int *checkbit;
    int checkbitlen = rc;


    int step;               // test 6 EbN0(dB) result
    int s = 0;              // receive 100 error block
    int num = 0;            // do compute block
    int error;              // error block has the number of error
    int totalerror=0;       // after receive 100 error block, how many total error 
    int restart = 0;        // compute check node = 0
    double sigma;
    double ebn0;
    double ebn0s[6];
    double bers[6];
    double berscompare[6];
    double x, y;
    int valL;
    int valL2;
    double app;
    double app1;
    int stp;
    
    // declaration end
    
    /*ebn0s[0] = 0.9844;
    ebn0s[1] = 1.289;
    ebn0s[2] = 1.584;
    ebn0s[3] = 1.868;
    ebn0s[4] = 2.411;
    ebn0s[5] = 3.046;
    berscompare[0] = 0.0583;
    berscompare[1] = 0.03052;
    berscompare[2] = 0.01011;
    berscompare[3] = 0.003;
    berscompare[4] = 8.34 * pow(10,-5);
    berscompare[5] = 2.799 * pow(10,-7);*/

    for (i = 0; i < 6; i++) bers[i] = 0;
    
    ebn0s[0] = 0.8279;
    ebn0s[1] = 1.214;
    ebn0s[2] = 1.584;
    berscompare[0] = 0.081;
    berscompare[1] = 0.01171;
    berscompare[2] = 1.696 * pow(10,-5); 
    // open file
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
    Llenrow = rc;
    Llencolumn = dc;
    Mlenrow = n; 
    Mlencolumn = dv;
    M = (int **)malloc(Mlenrow * sizeof(int *));
    for (i = 0; i < Mlenrow; i++) M[i] = (int *)malloc(Mlencolumn * sizeof(int));
    L = (int **)malloc(Llenrow * sizeof(int *));
    for (i = 0; i < Llenrow; i++) L[i] = (int *)malloc(Llencolumn * sizeof(int));
    //printf("Mlenrow = %d; Mlencolumn = %d\n", Mlenrow, Mlencolumn);
    for (j = 0; j < Mlenrow; j++) {
        for (i = 0; i < Mlencolumn; i++) {
            fscanf(fpr,"%d",&M[j][i]);
            //printf("M[%d][%d] = %d; ", j, i, M[j][i]);
        }
        //printf("\n");
    }
    int temp;
    for (j = 0; j < n; j++) {
        for (i = 0; i < dv; i++) {
            for (m = i; m < dv; m++) {
                if (M[j][m] < M[j][i]) {
                    temp = M[j][m];
                    M[j][m] = M[j][i];
                    M[j][i] = temp;
                }
            }
            printf("M[%d][%d] = %d; ", j, i, M[j][i]);
        }
        printf("\n");
    }
    for (i = 0; i < Llenrow; i++) {
        for (j = 0; j < dc; j++) {
            fscanf(fpr,"%d",&L[i][j]);
            //printf("L[%d][%d] = %d; ", i, j, L[i][j]);
        }
        //printf("\n");
    }
    for (i = 0; i < rc; i++) {
        for (j = 0; j < dc; j++) {
            for (m = j; m < dc; m++) {
                if (L[i][m] < L[i][j]) {
                    temp = L[i][m];
                    L[i][m] =L[i][j];
                    L[i][j] = temp;
                }
            }
            printf("L[%d][%d] = %d; ", i, j, M[i][j]);
        }
        printf("\n");
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

    Ljlen = n;
    Lj = (double *)malloc(Ljlen * sizeof(double));
    if (Lj == NULL) {
        fprintf(stderr, "Error: unable to allocate required memory\n");
        return 1;
    }
    for (i = 0; i < Ljlen; i++) Lj[i] = 0;

    qijrow = dv;        // qijrow = dv = 3
    qijcolumn = n;      // qijcolumn = n = 816
    qij = (double **)malloc(qijrow * sizeof(double *));
    for (i = 0; i < qijrow; i++) qij[i] = (double *)malloc(qijcolumn * sizeof(double));
    
    uijrow = rc;        // uijrow = rc = 408
    uijcolumn = dc;     // uijcolumn = dc = 6 
    uij = (double **)malloc(uijrow * sizeof(double *));
    for (i = 0; i < uijrow; i++) uij[i] = (double *)malloc(uijcolumn * sizeof(double));
    computlen = n;      // computlen = n = 816 
    comput = (int *)malloc(computlen * sizeof(int));
    comput1len = rc;    // comput1len = rc = 408
    comput1 = (int *)malloc(comput1len * sizeof(int));
    comput2len = rc;    // comput2len = rc = 408
    comput2 = (int *)malloc(comput2len * sizeof(int));
    qjlen = n;          //qjlen = n = 816
    qj = (double *)malloc(qjlen * sizeof(double));
     
    checkbitlen = rc;   // checkbitlen = rc = 408
    checkbit = (int *)malloc(checkbitlen * sizeof(int));

    printf("computlen = %d; comput1len = %d; comput2len = %d\n" ,computlen,comput1len,comput2len);
    for (step = 0; step < 1; step++) {
        s = 0;
        num = 0;
        totalerror = 0;
        //ebn0 = ebn0s[step];
        //sigma = sqrt(1.0 / (/*log2(rc) **/ pow(10, ebn0/10)));
        //ebn0 = pow(10, ebn0/10);
        while (s < 2/*s < 100*/) {
            //s++;
            num++;                                  // compute the number of transmit block 
            //printf("cod = \n");                     // pretend encoder
            for (i = 0; i < n; i++) {
                codarray[i] = 0;
                //printf("%d ", codarray[i]);
            }
            //printf("\n");
            //printf("code = \n");                    //input to AWGN channel normalized to +-1
            for (i = 0; i < n; i++) {
                if (codarray[i] == 0) codearray[i] = 1;
                else codearray[i] = -1;
                //printf("%d ", codearray[i]);
            }
            for(i = 0; i < n; i=i+2) {
                normal(SIGMA, &x, &y);
                outp[i] = codearray[i] + x;
                outp[i + 1] = codearray[i + 1] + y;
            }
            
            for(i = 0; i < n; i++) {
                Lj[i] =2 * SNR * outp[i]; 
            }
            // the interative decoding algotrithm
            for (j = 0; j < n; j++) {                               // initialization
                for (i = 0; i < dv; i++) {
                        qij[i][j] = Lj[j];
                }  
            }
            for (i = 0; i < rc; i++) {
                for (j = 0; j < 6; j++) uij[i][j] = 0;
            }
            for (k = 0; k < 100/*k < 100*/ && restart != rc; k++) {         // message passing, for predetermined threshold = 100
                restart = 0;
                app = 0;
                app1 = 0;    
                for (i = 0; i < 5; i++) tempqij[i] = 0.0;                        // bottom-up  
                for (i = 0; i < n; i++) comput[i] = 0;
                for (i = 0; i < rc; i++) {
                    for (j = 0; j < 6; j++) {
                        for (m = 0; m < 5; m++) {
                            if (m < j) {
                                valL = L[i][m] - 1;
                                tempqij[m] = qij[comput[valL]][valL];
                            } 
                            else if (m >= j) {
                                valL = L[i][m + 1] - 1;
                                tempqij[m] = qij[comput[valL]][valL];
                            }
                        }
                        tempuij = tempqij[0];
                        for(m = 1; m < 5; m++) {
                            tempuij = CHK(tempuij, tempqij[m]);
                        
                        }
                        uij[i][j] = tempuij;
                        //if (i == 3999) printf("uij[%d][%d] = %g; \n", i, j, uij[i][j]);
                    }
                    //printf("\n");
                    for (m = 0; m < 6; m++) comput[L[i][m] - 1] += 1;
                }

                // top-down
                for(i = 0; i < 3; i++) temp1uij[i] = 0.0;  
                for (i = 0; i < comput1len; i++) comput1[i] = 0;
                for (j = 0; j < n; j++) {
                    for (i = 0; i < 3; i++) {
                        for (m = 0; m < 2; m++) {
                            if (m < i) { 
                                valL = M[j][m] - 1;
                                temp1uij[m] = uij[valL][comput1[valL]]; 
                            }
                            else if (m >= i) {
                                valL = M[j][m + 1] - 1;
                                temp1uij[m] = uij[valL][comput1[valL]];
                            }
                        }
                        temp1uij[2] = Lj[j];
                        qij[i][j] = temp1uij[0] + temp1uij[1] + temp1uij[2];
                    }
                    for (m = 0; m < 3; m++) comput1[M[j][m] - 1] += 1;
                }
                // decision
                //printf("output = ");
                for (i = 0; i < comput2len; i++) comput2[i] = 0;
                for (j = 0; j < n; j++) {
                    qj[j] = Lj[j];
                    for (i = 0; i < 3; i++) {
                        valL = M[j][i] - 1;
                        qj[j] += uij[valL][comput2[valL]];
                    }
                    //if (j == 7999) printf("qj[%d] = %g;\n",j,qj[j]);
                    if (qj[j] >= 0) output[j] = 0;
                    else if (qj[j] < 0) output[j] = 1;
                    //if (j == 7999 && output[j] == 1) printf("%d \n", output[j]);
                    for (i = 0; i < 3; i++) {
                        comput2[M[j][i] - 1] += 1;
                    }
                }
                //printf("\n");

                // to check Hx=0
                //printf("checkbit = ");          
                for (i = 0; i < rc; i++) {
                    checkbit[i] = 0;
                    for (j = 0; j < 6; j++) {
                        checkbit[i] += output[L[i][j] - 1];
                    }
                    checkbit[i] = checkbit[i] % 2;
                    if (checkbit[i] == 0) restart += 1;
                }
                stp = 0;
                if (k == 99/*99*/ && restart != rc) {
                    //printf("failure\n");
                    stp = 1;
                    s++;
                }
            }
            /*if(k == 100) */printf("s = %d; k[%d] = %d\n", s, num, k);
            error = 0;
            for(i = 0; i < n; i++) {
                if (output[i] != codarray[i]) {
                    error += 1;
                    //printf("ouput[%d] = %d",i, output[i]);
                }
            }
            //printf("\n");
            if (error != 0 && stp == 0) s++;
            restart = 0;
            printf("error = %d\n", error);       
            totalerror += error;
            //s++;
        }
        double ber;
        ber = (double)totalerror / (num * n);
        printf("totalerror = %d\n", totalerror);
        printf("BER = %g\n", ber);
        bers[step] = ber;
    }

    for(step = 0; step < 3; step++) {
        printf("enb0s[%d] = %g\n",step, ebn0s[step]);
        printf("bers[%d] = %g\n",step, bers[step]);
        printf("berscompare[%d] = %g\n\n", step, berscompare[step]);
    }
    // CODE  end

    // open write file
    FILE *outfp;
    
    outfp = fopen("result8000.txt","w");
    for (i = 0; i < 3; i++) {
         fprintf(outfp,"%g ",ebn0s[i]);
         fprintf(outfp,"%g ",berscompare[i]);
         fprintf(outfp,"%g ",bers[i]);
         fprintf(outfp,"\n");
    }
    fclose(outfp);
    // close write file

    // free dynametic
    free(codarray);
    free(codearray);
    free(outp);
    free(output);
    for (i = 0; i < Llenrow; i++) free(L[i]);
    free(L);
    for (i = 0; i < Mlenrow; i++) free(M[i]);
    free(M);
    for (i = 0; i < qijrow; i++) free(qij[i]);
    free(qij);
    for (i = 0; i < uijrow; i++) free(uij[i]);
    free(uij);
    free(Lj);
    free(comput);
    free(comput1);
    free(comput2);

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