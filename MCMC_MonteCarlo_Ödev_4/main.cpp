#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXLEN 300000

char B[4] = {'A','G','T','C'};

int idx(char c)
{

    int ret = -1;

    switch(c)
    {
        case 'A': {ret =0;  break;}
        case 'G': {ret =1;  break;}
        case 'T': {ret =2;  break;}
        case 'C': {ret =3;  break;}
        default : {break;}
    }
    return ret;
}

char base(int i){ return B[i]; }

double urand()
{
    double randomNumber = (double)(rand() / (double)RAND_MAX);

    return randomNumber;
}

int categorical(double *p)
{
    double r = urand();
    double c = 0;

    for(int i=0;i<4;i++)
    {
        c += p[i];
        if(r < c) 
            return i;
    }
    return 3;
}

void compute_pmf(long counts[4], long total, double pmf[4]) {
    for (int i = 0; i < 4; i++) {
        if (total == 0)
            pmf[i] = 0.0;
        else
            pmf[i] = (double)counts[i] / (double)total;
    }
}

void print_pmf(double pmf[4], const char *label) {
    printf("\n=== %s PMF ===\n", label);
    printf("P(A) = %.5f\n", pmf[0]);
    printf("P(G) = %.5f\n", pmf[1]);
    printf("P(T) = %.5f\n", pmf[2]);
    printf("P(C) = %.5f\n", pmf[3]);
}

int main()
{
    srand(time(NULL));

    FILE *fc = fopen("gene.txt","r");
    FILE *fr = fopen("gene_rev.txt","r");

    if(fc == NULL || fr == NULL)
    {
        if (fc == NULL)
        {
            printf("gene.txt not found \r\n");
        }
        else
        {
            printf("gene_rev.txt not found \r\n");
        }
        system("pause");
        return 1;
    }

    char clean[MAXLEN];
    char noisy[MAXLEN];

    int n_clean = 0, n_rev = 0, ch;

    while((ch = fgetc(fc)) != EOF)
    {
        if(ch=='A'||ch=='G'||ch=='T'||ch=='C')
            clean[n_clean++] = ch;
    }

    clean[n_clean] = '\0';
    fclose(fc);

    while((ch = fgetc(fr)) != EOF)
    {
        if(ch=='A'||ch=='G'||ch=='T'||ch=='C'||ch=='X')
            noisy[n_rev++] = ch;
    }

    noisy[n_rev] = '\0';
    fclose(fr);

    // ---------------------------
    // 1) Transition matrix
    // ---------------------------

    int count[4][4] = {0};

    for(int i=0;i<n_clean-1;i++)
    {
        int a = idx(clean[i]);
        int b = idx(clean[i+1]);
        if(a>=0 && b>=0) count[a][b]++;
    }

    double M[4][4];
    for(int i=0;i<4;i++){
        double s = 0;
        for(int j=0;j<4;j++) s += count[i][j];
        for(int j=0;j<4;j++){
            if(s==0) M[i][j] = 0.25;
            else M[i][j] = count[i][j] / s;
        }
    }

    /*print Transition Matrix*/
    printf("=== Transition Matrix (M) ===\n");
    for(int k = 0; k < 4; k++)
    {
        for(int t = 0; t < 4; t++)
        {
            printf("%.3f ", M[k][t] );
            
        }
        printf("\n\r");
    }
    // ---------------------------
    // 2) X’LERY TAHMYN ETME
    // ---------------------------

    char det[MAXLEN];
    char mc[MAXLEN];

    memcpy(det, noisy, n_rev);
    memcpy(mc,  noisy, n_rev);

    long hist_det[4]={0};
    long hist_mc[4]={0};

    int last_det = -1;
    int last_mc  = -1;

    for(int i=0;i<n_rev;i++){
        int c = idx(noisy[i]);

        if(c >= 0){
            // gerçek nükleotid
            last_det = c;
            last_mc  = c;

            det[i] = base(c);
            mc[i]  = base(c);

            hist_det[c]++;
            hist_mc[c]++;
            continue;
        }

        // X bulundu
        if(last_det < 0) last_det = 0;
        if(last_mc  < 0) last_mc  = 0;

        // --- DETERMINISTIC ---
        double bestP = -1;
        int bestS = 0;
        for(int s=0;s<4;s++){
            double p = M[last_det][s];
            if(p > bestP){
                bestP = p;
                bestS = s;
            }
        }
        det[i] = base(bestS);
        hist_det[bestS]++;
        last_det = bestS;

        // --- MONTE CARLO ---
        double p4[4], sum = 0;
        for(int s=0;s<4;s++){
            p4[s] = M[last_mc][s];
            sum += p4[s];
        }
        for(int s=0;s<4;s++) p4[s] /= sum;

        int sampled = categorical(p4);
        mc[i] = base(sampled);
        hist_mc[sampled]++;
        last_mc = sampled;
    }

    // ---------------------------
    // 3) FASTA YAZDIR
    // ---------------------------

    FILE *fd = fopen("deterministic_corrected.fasta","w");
    FILE *fm = fopen("montecarlo_corrected.fasta","w");

    fprintf(fd,">Deterministic\n");
    fprintf(fm,">MonteCarlo\n");

    int col = 0;
    for(int i=0;i<n_rev;i++){
        fputc(det[i], fd);
        col++; if(col==70){ fputc('\n',fd); col=0; }
    }
    fputc('\n',fd);

    col = 0;
    for(int i=0;i<n_rev;i++){
        fputc(mc[i], fm);
        col++; if(col==70){ fputc('\n',fm); col=0; }
    }
    fputc('\n',fm);

    fclose(fd);
    fclose(fm);

    // ---------------------------
    // 4) HISTOGRAM
    // ---------------------------
    printf("=== DETERMINISTIC HISTOGRAM ===\n");
    printf("A: %ld\nG: %ld\nT: %ld\nC: %ld\n",hist_det[0], hist_det[1], hist_det[2], hist_det[3]);

    printf("\n=== MONTE CARLO HISTOGRAM ===\n");
    printf("A: %ld\nG: %ld\nT: %ld\nC: %ld\n",hist_mc[0], hist_mc[1], hist_mc[2], hist_mc[3]);


    // ---- PMF Hesaplama ----
    long total_det = hist_det[0] + hist_det[1] + hist_det[2] + hist_det[3];
    long total_mc  = hist_mc[0]  + hist_mc[1]  + hist_mc[2]  + hist_mc[3];


    double pmf_det[4], pmf_mc[4];

    compute_pmf(hist_det, total_det, pmf_det);
    compute_pmf(hist_mc,  total_mc,  pmf_mc);

    // Ekrana yazdır
    print_pmf(pmf_det, "Deterministic");
    print_pmf(pmf_mc,  "Monte Carlo");

    system("pause");
    return 0;
}


