#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

//Nucleic Acids
const char NucleicAcids[4] = {'A','G','T','C'};

const int NumberOfNA = 4;

const char *NADef[]={
    "Ade",
    "Gua",
    "Thy",
    "Cyt",
};

/* 
Nucleic Acid Codes in FASTA File
A   Ade Adenine
G   Gua Guanine
T   Thy Thymine
C   Cyt Cytosine
*/

const SliceLength = 31; 
float TrMatrix[4][4] = {{0.0f,0.0f,0.0f,0.0f}, 
                         {0.0f,0.0f,0.0f,0.0f},
                         {0.0f,0.0f,0.0f,0.0f},
                         {0.0f,0.0f,0.0f,0.0f}};
float NAHist[4];
float StateHist[4][4];
float TranSum;
int StPrDist[4]; //monte carlo sonrasý duragan olasýlýk dagýlýmý
const int IntState = 2; // T
int CrtState;

char ch, tmp;
char str[10];


unsigned int i,k,j,t,size,cnt=0;
float t_cnt, s_cnt;

                         
float get_random(void){
    float r_num;
        r_num =(float)rand()/RAND_MAX;
    return r_num;
}

/** 
* Ters donusum ornekleyicisi
*/
int it_sampler(int NACode){
    
    int m, tr_value = 0;
    float u = 0.0f;
    float lc,hc;

    u = get_random(); //0-1 arasi rastgele sayi
    
    if(u <= TrMatrix[NACode][0]){ //A
        tr_value = 0;   
    }
    
    lc = TrMatrix[NACode][0];
    hc = TrMatrix[NACode][0] + TrMatrix[NACode][1];
    if(lc<=u && u<hc){ //G
        tr_value = 1;   
    }
    
    lc = TrMatrix[NACode][0] + TrMatrix[NACode][1];
    hc =TrMatrix[NACode][0] + TrMatrix[NACode][1] + TrMatrix[NACode][2];
    if(lc<=u && u<hc){ // T
        tr_value = 2;   
    }
    
    lc = TrMatrix[NACode][0] + TrMatrix[NACode][1] + TrMatrix[NACode][2];
    hc = TrMatrix[NACode][0] + TrMatrix[NACode][1] + TrMatrix[NACode][2] + TrMatrix[NACode][3];
    if(lc<=u && u<hc){ // C
        tr_value = 3;   
    }
    
    return tr_value;
}

int main() {
    
        
    for(k = 0; k < NumberOfNA; k++){
        NAHist[k] = 0.0f;
    }
    
    for(k = 0; k < 4; k++){
            
        for(t = 0; t < 4; t++){
            
            StateHist[k][t] = 0.0f;
            
        }
    }
    
    
    printf("*************************Fasta Okumalari*******************************\n\r");
    
    // Gene Squences FASTA Code from NCBI Database
    FILE *fptr = fopen("gene.txt", "r");

        
    // Karakter bazinda okuma ve siniflama
    while ((ch = fgetc(fptr)) != EOF){ 
        
        printf("%c", ch);
        
        if(cnt == 0){ //ilk karaterde durum gecisi tespiti yapilmiyor.
            
            for(i = 0; i < NumberOfNA; i++){
            
                if(ch==NucleicAcids[i]){
                    NAHist[i] += 1.0f;
                    break;
                }

            }
            tmp = ch;
            cnt++; 
        }
        else{
            
            for(i = 0; i < NumberOfNA; i++){
            
                if(ch==NucleicAcids[i]){
                    NAHist[i] += 1.0f;
                    break;
                }

            }
            
            //tmp --> ch 
        
            for(j = 0; j < 4; j++){
                    if(tmp == NucleicAcids[j]){
                        t = j;
                    }
                }
            
            for(k = 0; k < 4; k++){
                
                if(ch == NucleicAcids[k]){
                    StateHist[t][k] += 1.0f; 
                }
                
            }
            tmp = ch;

        }       

    }
    
    // dosya kapama
    fclose(fptr);
    
    cnt = 0;

   
    //satýrlardaki olasýlýk deðerlerinin normalizasyonu
    bool no_trs = false;
    float residue = 0.0f;
    int indx;
    
    for(k = 0; k < 4; k++){
            
        for(t = 0; t < 4; t++){

            if(StateHist[k][t] == 0){
             no_trs = true;
             indx = t;
            }
            TranSum += StateHist[k][t];
        }
        
        if(TranSum != NAHist[k]){
            if(no_trs == true){
                residue = (NAHist[k] - TranSum)/3.0f;
            }
            else{
                residue = (NAHist[k] - TranSum)/4.0f;
            }
            
            for(j = 0; j < 4; j++){
            
                if((j == indx) && (no_trs == true)){
                    StateHist[k][j] = 0.0f;     
                }
                else{
                    StateHist[k][j] = StateHist[k][j] + residue; 
                }
            
            }
        }
            
        no_trs = false;
        TranSum = 0.0f;
        
    }
    
     
    //Tr matriks hesabý
    printf("\n\r");
    printf("Gecis Matrisi : \n\r");
    
    for(k = 0; k < 4; k++){
            
        for(t = 0; t < 4; t++){
            
            t_cnt = StateHist[k][t];
            s_cnt = NAHist[k];
            TrMatrix[k][t] = t_cnt/s_cnt;
            printf("%.3f ", TrMatrix[k][t] );
            
        }
        printf("\n\r");
    }
 

    printf("\n\r");
    printf("*****************************Sekans Sonu*******************************\n\r");

    printf("\n\r");
    printf("Nukleik Asit Histogram Sonuclari : \n\r");
    

    fptr = fopen("gene_histg.txt", "w+"); //proje klasorunde
    
    for(j = 0; j < NumberOfNA; j++){
        
        printf("%s %.3f\n\r", NADef[j], NAHist[j]);
        sprintf(str,"%s %.3f\n", NADef[j], NAHist[j]);  
        size = strlen(str);
        
        for(t = 0; t < 10; t++){
            if(str[t]!=0){
                putc(str[t],fptr);
            }
            else{
                break;
            }
        }
    }
  
       
    // dosya kapama
    fclose(fptr);
    
    //MH sim
    CrtState = IntState; //index T'yi temsil ediyor.
    StPrDist[IntState] += 1;
    
    int s, m, n, acc;
    float mu, crnt;
    
    for(s = 0; s < 100; s++){
        
        m = it_sampler(CrtState);
        StPrDist[m] += 1;
        CrtState = m;
/*      
        //duragan dagýlýma ulasýldý mý
        if(s > 10){
            acc = StPrDist[0]+StPrDist[1]+StPrDist[2]+StPrDist[3];
            mu  = (float)(acc) / s;
            crnt = ((float)acc * 0.98f) / s;
            if(mu > crnt){
                break;
            }
        }
*/
        
    }
    
    printf("%d iter \n\r", s);
    
    printf("\n\r");
    printf("Monte Carlo Histogram Sonuclari : \n\r");   
    
    for(k = 0; k < 4; k++){
        printf("%s  = %d\n\r", NADef[k], StPrDist[k]);
    }

        

    
    return 0;
}
