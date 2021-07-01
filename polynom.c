#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

typedef struct polynome{
    int* val;
    int degre;
    int modulo;
} polynom;

int add(int a, int b, int modulo){ // H : a,b déja modulo

    int res = a + b;

    if(res >= modulo){res -= modulo;}

    return res;

}

int sous(int a, int b, int modulo){ //H : a,b déja modulo

    int c = a-b;

    if (c<0){
        c += modulo;
    }

    return c;
}

int mult(int a, int b, int modulo){//H : a,b déja modulo

    int res = a*b;

    while (res >= modulo)
    {
        res -=modulo;
    }
    
    return res;

}

int maximum(int n,int m){

    if(n < m ){return m;}
    else{return n;}
    
}

int minimum(int n,int m){

    if(n < m ){return n;}
    else{return m;}

}

polynom* initialisation(int degre,int modulo){//créer un polynom

    polynom* res = (polynom*)malloc(sizeof(polynom));
    res->val = calloc(degre+1, sizeof(int));
    res->degre = degre;
    res->modulo = modulo;

    return res;

}

void libererPoynome(polynom* p){

    if(p){
        free(p->val);
        free(p);
    }

}

int degreP(polynom* p1){// met à jour le degré du polynom

    if(p1->degre == 0){
        return 0;
    }

    int res = p1->degre;

    while (p1->val[res] == 0 && res > 0)
    {
        res--;
    }
    
    return res;

}

int pgcd(int a, int b){

    if(a == 0){
        return b;
    }else if(b == 0){
        return a;
    }

    if(a==b){
        return a;        
    }     
    else
    {
        if(a>b)
           return pgcd(a-b, b);
        else
           return pgcd(a, b-a);
    }
}

int euclide_etendu(int a, int b, int* u, int* v) { // a,b entiers naturels
    int r1 = a;
    int r2 = b; 
    int u1 =1;
    int v1 = 0;
    int u2 = 0;
    int v2 = 1;

    while(r2 != 0){
        int q = r1/r2; //(division entiere)
        int r3 = r1;
        int u3 = u1; 
        int v3 = v1;
        r1 = r2;
        u1 = u2;
        v1 = v2;
        r2 = r3 - q*r2;
        u2 = u3 - q*u2;
        v2 = v3 - q*v2;        
    }
    
    *u = u1;
    *v = v1;

    return r1; //(r1 entier naturel et u1, v1 entiers relatifs)
}



void affichage(polynom* p1){ //affiche le polynome

    for(int i = p1->degre; i >= 0; i--){



        printf("+ %dx^%d ",p1->val[i],i);
    }
    
    printf("\n");

}

int inverse(int a, int modulo,int* u, int* v){//trouve l'inverse modulaire (code à surement revoir)

    if(pgcd(a,modulo) != 1){return -1;}
    else{

        euclide_etendu(a,modulo,u,v);

    }

    while(*u < 0){*u = *u + modulo;}
    return *u;
}

int* random_tab(int deg,int modulo){//rends un tab de valeurs aléatoire à un certain modulo (positif)

    int i = 0;

    int* res = (int*)malloc(sizeof(int)*(deg+1));

    while (i < (deg + 1))
    {
        res[i] = (rand()%(modulo));
        i++;
    }

    return res;
        
}

polynom* randomP(int degre, int modulo){//rend un random polynom (pour les tests)

    polynom* res = initialisation(degre,modulo);

    res->degre = degre;
    res->modulo = modulo;
    res->val = random_tab(degre,modulo);
    res->degre = degreP(res);

    return res;

}

polynom* addition(polynom* p1, polynom* p2){//addition modulaire de deux polynomes

    int degre = maximum(p1->degre,p2->degre);
    int min = minimum(p1->degre,p2->degre);

    polynom* res = initialisation(degre,p1->modulo);

    int i = 0;

    while (i <= min)
    {
        res->val[i] = add(p1->val[i],p2->val[i],p1->modulo);
        i++;
    }
    
    if(p1->degre == p2->degre){
        res->degre = degreP(res);
        return res;

    }
    else if(p1->degre > p2->degre){

        while (i <= degre)
        {
            res->val[i] = p1->val[i];
            i++;
        }
        res->degre = degreP(res);
        return res;

    }
    else{

        while (i <= degre)
        {
            res->val[i] = p2->val[i];
            i++;
        }
        res->degre = degreP(res);
        return res;

    }
    
}

polynom* soustraction(polynom* p1, polynom* p2){// soustraction modulaire de deux polynomes

    // t1 - t2

    if(p1->modulo != p2->modulo){
        printf("Les polynomes n'ont pas des coefficient d'un meme corps \n");
        return NULL;
    }

    int modulo = p1->modulo;

    int degre = maximum(p1->degre,p2->degre);

    int taille = degre + 1;

    polynom* res = initialisation(degre,modulo);

    int min = minimum(p1->degre,p2->degre);

    for(int i=0; i<=min; i++){

        int p2val = p2->val[i];
        int val = sous(p1->val[i],p2val,modulo);
        res->val[i] = val;
    
    }
    
    if(p1->degre == p2->degre){

        res->degre = degreP(res);
        return res;
    }
    else{    
        if(min == p1->degre){

            for(int i=0; i<=taille; i++){
                res->val[i] = sous(0,p2->val[i],modulo);
            }
            

        }
        else{

            for(int i=0; i<=taille; i++){
                res->val[i] = p1->val[i];
            }

        }
    
    }

    res->degre = degreP(res);
    return res;
}

polynom* multiplication(polynom* p1, polynom* p2){

    polynom* res = initialisation(p1->degre + p2->degre,p1->modulo);

    int deg = p1->degre + p2->degre + 1;

    for(int i=0; i <= p1->degre; i++){

        for(int j=0; j <= p2->degre; j++){

            res->val[i + j] = add(res->val[i+j],mult(p1->val[i],p2->val[j],p1->modulo),p1->modulo);

        }
    }

    res->degre = degreP(res);

    return res;

}

polynom* augmenter_puissance(polynom* p, int newdegre, int coeff1, int coeff2, int* u, int* v){

    polynom* res = initialisation(newdegre, p->modulo);

    int diffdeg = newdegre - p->degre;

    for(int i = 0; i <= p->degre; i++){

        int val = mult(mult(coeff1, inverse(coeff2, p->modulo, u, v), p->modulo),p->val[i], p->modulo);
        res->val[i+diffdeg] = val;
    
    }
    
    res->degre = degreP(res);

    return res;
}

int newdegre(int* t1,int taillemax){

    for(int i = taillemax; i >= 0; i--){
        if(t1[i] !=0){
            return i;
        }
    }
    
    return 0;
}

polynom* reste_division_polynome(polynom* p1, polynom* p2, int* u, int* v){ //renvoie le reste de la division euclidienne de p1/p2

    //hypothèse :  deg p1 >= deg tp2

    int deg1 = p1->degre;
    int deg2 = p2->degre;

    int degReste = deg1;

    polynom* reste = initialisation(deg1, p1->modulo);

    reste->val = p1->val;

    while (degReste >= deg2)
    {
        int coeff1 = reste->val[degReste];
        int coeff2 = p2->val[deg2];

        polynom* b = augmenter_puissance(p2, degReste, coeff1, coeff2, u, v);
        
        reste = soustraction(reste,b);
        degReste = reste->degre;

        libererPoynome(b);
        
    }
    
    return reste;
}

int estPolynomeNul(polynom* p){

    for(int i = 0; i<=p->degre; i++){
        if(p->val[i] != 0){
            return 0;
        }
    }

    return 1;

}

polynom* division_polynome(polynom* p1, polynom* p2, int* u, int* v){ //renvoie le quotient de la division euclidienne de p1/p2

    int deg1 = p1->degre;
    int deg2 = p2->degre;

    if(estPolynomeNul(p2)){
        printf("on ne peut pas diviser par le polynome nul\n");
        return NULL;
    }

    if(deg1<deg2){

        return division_polynome(p2, p1, u, v);

    }

    int degReste = deg1;
    int degQuotient = deg1 - deg2;

    polynom* reste = initialisation(deg1, p1->modulo);
    polynom* quotient = initialisation(degQuotient, p1->modulo);

    reste->val = p1->val;

    int i = degQuotient;

    while (degReste >= deg2){

        int coeff1 = reste->val[degReste];
        int coeff2 = p2->val[deg2];

        int coeffQuotient = mult(coeff1, inverse(coeff2, p1->modulo, u, v), p1->modulo);
        quotient->val[i] = coeffQuotient;

        polynom* b = augmenter_puissance(p2, degReste, coeff1, coeff2, u, v);
        reste = soustraction(reste,b);
        degReste = reste->degre;

        libererPoynome(b);

        i = degReste - deg2;
    }
    
    libererPoynome(reste);
    quotient->degre = degreP(quotient);

    return quotient;
}

polynom* pgcd_polynome(polynom* p1,polynom* p2,int* u,int* v){

    polynom* a = initialisation(p1->degre, p1->modulo);
    polynom* b = initialisation(p2->degre, p2->modulo);

    a->val = p1->val;
    b->val = p2->val;

    int deg1 = p1->degre;
    int deg2 = p2->degre;

    while (deg2 > 0){
        polynom* temp = b;
        b = reste_division_polynome(a,b,u,v);
        a = temp;
        deg1 = deg2;
        deg2 = degreP(b);
    }

    libererPoynome(b);
    
    return a;

}

polynom* derive(polynom* p){
    if (p->degre==0){
        return initialisation(0,p->modulo);
    }

    polynom* p2 = initialisation(p->degre-1,p->modulo);

    unsigned int i=0;

    while (i<p->degre){
        p2->val[i]=mult((i+1), p->val[i+1], p->modulo);
        i++;
    }

    p2->degre = degreP(p2);
    return p2;

}

polynom** algorithme_yun(polynom* p){
    polynom** LP = malloc((p->degre)*sizeof(polynom));

    int u, v;
    int i = 1;

    polynom* G = pgcd_polynome(p, derive(p), &u, &v);
    affichage(G);
    polynom* H = division_polynome(p, G, &u, &v);
    affichage(H);
    polynom* K = division_polynome(derive(p), G, &u, &v);
    affichage(K);

    while(degreP(H) >= 1){

        polynom* R = pgcd_polynome(H, soustraction(K, derive(H)), &u, &v);
        affichage(R);
        if(!estPolynomeNul(R)){
            K = division_polynome(soustraction(K, derive(H)), R, &u, &v);
            affichage(K);
        }    
        if(!estPolynomeNul(R)){
            H = division_polynome(H, R, &u, &v);
            affichage(H);
        }
        LP[i] = R;
        i++;

    }

    return LP;

    
}

int puissance(int a, int p, int mod){

    int res = pow(a,p);

    while(res >= mod){
        res -= mod;
    }

    return res;

}

int est_racine(polynom* p, int racine){

    int res = 0;

    int i = 0;

    while (i <= p->degre)
    {
        int puis = puissance(racine,i,p->modulo);
        int pval = p->val[i];
        int multi = mult(pval,puis,p->modulo);
        res = add(multi,res,p->modulo);

        i++;

    }
    
    if(res ==0){return 1;}
    
    return 0;
}


polynom* polynom_racine(int racine, int modulo){

    polynom* res = initialisation(1,modulo);

    res->val[1] = 1;

    res->val[0] = (-1)*racine + modulo;

    return res;

}

polynom** algo_naif(polynom* p){

    polynom** res = (polynom**)malloc(sizeof(polynom*)*p->degre);
    polynom* p2 = initialisation(p->degre, p->modulo);
    p2->val = p->val;
    int* racine = malloc(sizeof(int)*p->degre);
    int compteur = 0;
    int i = 0;
    int u,v;
    int degre = p->degre;

    while (i <= p->modulo && compteur != p->degre)
    {
        while(est_racine(p2,i) == 1)
        {
            // p = p / X - racine
            res[compteur] = polynom_racine(i,p->modulo);
            racine[compteur] = i;
            p2 = division_polynome(p2,res[compteur],&u,&v);
            affichage(p2);
            compteur++;

        }
        
        i++;

    }

    if(compteur != degre){

        libererPoynome(p2);
        printf("le polynom n'est pas scindable dans notre corps\n");
        return NULL;

    }
    
    /*i = 0;

    while (i < p->degre)
    {
        res[i] = initialisation(1,p->modulo);

        res[i]->val[1] = 1;

        if(racine[i] > 0){racine[i] = (-1)*racine[i] + p->modulo;}

        res[i]->val[0] = racine[i];

        i++;

    }*/
    
    libererPoynome(p2);
    return res;

}




int* Liste_deux_n(int n){

    int* res = (int*)malloc(sizeof(int)*n);

    int i = 2;

    while (i < n+1)
    {

        res[i-2] = i;
        i++;

    }
    
    res[i-2] = -1;

    return res;

}

int* reallocation(int* Liste,int nbelem){

    int i = 0;

    int compteur = 0;

    int* res = (int*)malloc((sizeof(int))*(nbelem + 1));

    while (Liste[i] != -1)
    {
        if(Liste[i] != 0){

            res[compteur] = Liste[i];
            compteur++;

        }

        i++;

    }
    
    res[compteur] = -1;

    free(Liste);

    return res;

}

int* liste_non_diviseur(int* Liste,int n,int max){

    if(n==1){return Liste;}

    if(n < 1){return NULL;}

    int* res = (int*)calloc(max,sizeof(int));
    int i = 0;

    int compteur = 0;
   
    while (Liste[i] != -1)
    {   
        if(Liste[i]%n != 0){

            res[compteur] = Liste[i];
            compteur++;
            //printf("on a ajoute %d en indice %d qui donne dans res : %d \n",Liste[i],compteur-1,res[compteur - 1]);
        }

        i++;

    }

    res[compteur] = -1;

    free(Liste);

    res = reallocation(res,compteur);
    
    return res;

}

void affiche_liste(int* Liste){

    int i = 0;

    while (Liste[i] != -1)
    {
        printf("%d ",Liste[i]);
        i++;
    }
    
    printf("\n");

}

int* crible(int n){

    int* res = calloc(n,sizeof(int));

    int i = 0;

    int* Liste = Liste_deux_n(n);

    while (Liste[0] != -1)
    {
        res[i] = Liste[0];

        Liste = liste_non_diviseur(Liste,Liste[0],n);

        i++;

    }

    res[i] = -1;

    res = reallocation(res,i);

    return res;
     
}

int* erastothene(int n){

    int taille = (int)n/2 + 1;

    int* Prem = crible(n);
 
    int compteur = 0;

    int i = 0;

    int reste = n;

    int* res = calloc(taille,sizeof(int));

    while (reste != 1)
    {
        while ((reste % Prem[i] == 0)&&(reste !=1))
        {

            res[compteur] = Prem[i];
            compteur++;
            reste = reste/Prem[i];

        }
        
        i++;

    }
    
    res[compteur] = -1;

    free(Prem);

    res = reallocation(res,compteur);

    return res;

}



int main(int argc, char const *argv[]){

    srand(time(NULL));
    
    polynom* p1 = initialisation(5,5);
    polynom* p2 = initialisation(3,5);

    int t1[6] = {2,2,3,4,3,4};
    int t2[4] = {2,4,3,2};

    p1->val = t1;
    p2->val = t2;

    printf("\nP1 : ");
    affichage(p1);
    printf("P2 : ");    
    affichage(p2);

    polynom* p3 = addition(p1,p2);
    printf("\nAddition P1 + P2 : ");
    affichage(p3);

    polynom* p4 = multiplication(p1,p2);
    printf("Multiplication P1*P2 : ");
    affichage(p4);

    polynom* p5 = soustraction(p1,p2);
    printf("Soustraction P1 - P2 : ");
    affichage(p5);

    polynom* p8 = derive(p1);
    printf("Derivée de P1 : ");
    affichage(p8);

    polynom* p12 = derive(p2);
    printf("Derivée de P2 : ");
    affichage(p12);
    
    int u = 0;
    int v = 0;

    polynom* p6 = reste_division_polynome(p1, p2, &u, &v);
    printf("\nReste de la divison euclidienne P1/P2 : ");
    affichage(p6);

    polynom* p10 = division_polynome(p1, p2, &u, &v);
    printf("Quotient de la division euclidienne P1/P2 : ");
    affichage(p10);

    polynom* p7 = multiplication(p10, p2);
    printf("Multiplication du quotient avec P2 : ");
    affichage(p7);

    polynom* p11 = addition(p7, p6);
    printf("Puis on ajoute le reste et on trouve P1 = ");
    affichage(p11);

    polynom* p9 = pgcd_polynome(p1, p2, &u, &v);
    printf("\nPGCD : ");
    affichage(p9);
    
    polynom* p20 = initialisation(3, 5);
    int tab4[4] = {4,1,4,1};
    p20->val = tab4;
    printf("\np20 = ");
    affichage(p20);
    
    polynom** p19 = algo_naif(p20);
    for(int i=0; i<p20->degre; i++){
        printf("facteur %d : ", i);
        affichage(p19[i]);
    }

    //free(p1);
    //free(p2);
    libererPoynome(p3);
    libererPoynome(p4);
    libererPoynome(p5);
    libererPoynome(p6);
    libererPoynome(p7);
    libererPoynome(p8);
    libererPoynome(p9);
    libererPoynome(p10);
    libererPoynome(p11);
    libererPoynome(p12);
    //libererPoynome(p13);
    //libererPoynome(p14);
    //libererPoynome(p15);
    //libererPoynome(p16);
    //(p17);
    //(p18);
    //libererPoynome(p19);
    //libererPoynome(p20);

    return 0;
}


/*

int main(int argc, char const *argv[])
{
    
    int* test1 = Liste_deux_n(100);

    int* test2 = liste_non_diviseur(test1,2,100);

    affiche_liste(test2);
    
    int* test3 = liste_non_diviseur(test2,3,100);

    affiche_liste(test3);
    printf("pd 1 \n");
    int* test4 = liste_non_diviseur(test3,5,100);
    printf("pd 2 \n");
    affiche_liste(test4);


    free(test1);
    free(test2);
    free(test3);
    free(test4);

    //int* gabriel = crible(150000);

    //affiche_liste(gabriel);

    int* yesin = erastothene(1234567);

    affiche_liste(yesin);

    free(yesin);

    //free(gabriel);

    return 0;
}
*/