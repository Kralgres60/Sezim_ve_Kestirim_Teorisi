#include <stdio.h>
#include <string.h>

#define LEN 1000

/*Function Prototypes*/
char Base(double v[]);
void ExtractProb(const char *s, double M[4][4]);
/*Function Prototypes*/


int main() 
{
    // Gene Squences FASTA Code from NCBI Database
    FILE *fptr = fopen("gene_rev.txt", "r");

    if (fptr == NULL)
    {
       return 1;
    }


    double M[4][4] = 
    {
        {0.364, 0.209, 0.295, 0.133},
        {0.329, 0.209, 0.283, 0.179},
        {0.250, 0.216, 0.365, 0.169},
        {0.382, 0.038, 0.384, 0.197}
    };
    
    double v[4];
   
    //First Distribution
    v[0] = 1;   // A
    v[1] = 0;   // T
    v[2] = 0;   // G
    v[3] = 0;   // C

    char ch = 0;


    // Karakter bazinda okuma ve siniflama
    while ((ch = fgetc(fptr)) != EOF)
    {
        if (ch == 'A')
        {
            v[0]=1;v[1]=v[2]=v[3]=0;
        }
        else if (ch == 'T')
        {
           v[1]=1; v[0]=v[2]=v[3]=0;
        }
        else if (ch == 'G')
        {
            v[2]=1; v[0]=v[1]=v[3]=0;
        }
        else if (ch == 'C')
        {
            v[3]=1; v[0]=v[1]=v[2]=0;
        }
        else if (ch == 'X')
        {
            double x = (v[0]*M[0][0])+(v[1]*M[1][0])+(v[2]*M[2][0])+(v[3]*M[3][0]);
            double y = (v[0]*M[0][1])+(v[1]*M[1][1])+(v[2]*M[2][1])+(v[3]*M[3][1]);
            double z = (v[0]*M[0][2])+(v[1]*M[1][2])+(v[2]*M[2][2])+(v[3]*M[3][2]);
            double w = (v[0]*M[0][3])+(v[1]*M[1][3])+(v[2]*M[2][3])+(v[3]*M[3][3]);

            double vector[4] = {x,y,z,w};
            char corrected_word = Base(vector);

            printf("predicted word = %c \r\n",corrected_word);
        }
    }

    fclose(fptr);
    return 0;
}


/*chosing the greater probabilty number*/
char Base(double v[]) 
{
    int max_index = 0;
    double max_val = v[0];

    for (int i = 1; i < 4; i++) {
        if (v[i] > max_val) {
            max_val = v[i];
            max_index = i;
        }
    }

    if (max_index == 0) return 'A';
    if (max_index == 1) return 'T';
    if (max_index == 2) return 'G';
    return 'C';
}

// Extract transition matrix from DNA
void ExtractProb(const char *s, double M[4][4]) 
{
    int len = strlen(s);
    
    /*reset matrix*/
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            M[i][j] = 0;

    double Ta=0, Tt=0, Tg=0, Tc=0;

    for (int i = 0; i < len - 1; i++) {
        char DI1 = s[i];
        char DI2 = s[i+1];

        int r, c;

        if (DI1=='A') { r=0; Ta++; }
        else if (DI1=='T') { r=1; Tt++; }
        else if (DI1=='G') { r=2; Tg++; }
        else { r=3; Tc++; }

        if (DI2=='A') c=0;
        else if (DI2=='T') c=1;
        else if (DI2=='G') c=2;
        else c=3;

        M[r][c] += 1;
    }

    // normalizing
    for (int i = 0; i < 4; i++) {
        double total = (i==0?Ta : i==1?Tt : i==2?Tg : Tc);

        if (total == 0) continue;

        for (int j = 0; j < 4; j++)
            M[i][j] = M[i][j] / total;
    }
}

