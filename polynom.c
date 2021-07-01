#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

typedef struct polynome{
    int* val;
    int degre;
    int modulo;
} polynom;

// Fonctions calculs de base sur entiers modulaires

int mod_add(int a, int b, int modulo){ // H : a,b déja modulo

    int res = a + b;

    if(res >= modulo)
        res -= modulo;

    return res;
}

int mod_sous(int a, int b, int modulo){ //H : a,b déja modulo

    int c = a-b;

    if (c<0)
        c += modulo;

    return c;
}

int mod_mult(int a, int b, int modulo){//H : a,b déja modulo

    int res = a*b;
    int tmp = res/modulo;
    res -= tmp*modulo;

    return res;
}

int int_maximum(int n,int m){

    if(n < m)
        return m;
    
    return n;
}

int int_minimum(int n,int m){

    if(n < m)
        return n;
    
    return m;
}

int mod_puissance2(int a, int p, int mod){

    unsigned long long res = pow(a,p);
    int tmp = res/mod;
    res -= tmp*mod;

    return res;
}

int mod_puissance(int a, int p, int mod){

    int res = 1;

    for(int i=0; i<p; i++){
        res = mod_mult(res, a, mod);
    }

    return res;
}

// Fonctions un peu plus avancées sur les entiers

int int_pgcd(int a, int b){

    if(a == 0){
        return b;
    }else if(b == 0){
        return a;
    }

    if(a == b){
        return a;        
    }else if(a>b){
        return int_pgcd(a-b, b);
    }

    return int_pgcd(a, b-a);
}

int int_euclide_etendu(int a, int b, int* u, int* v) { // a,b entiers naturels
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

int mod_inverse(int a, int modulo,int* u, int* v){ //trouve l'mod_inverse modulaire

    if(int_pgcd(a,modulo) != 1){
        return -1;
    }
    else{
        int_euclide_etendu(a,modulo,u,v);
    }

    while(*u < 0){
        *u += modulo;
    }

    return *u;
}

int* mod_tab_random(int deg,int modulo){ //rends un tab de valeurs aléatoire à un certain modulo (positif)

    int i = 0;

    int* res = (int*)malloc(sizeof(int)*(deg+1));

    while (i < (deg + 1)){
        res[i] = (rand()%(modulo));
        i++;
    }

    return res;  
}

int* Liste_deux_A_n(int n){

    int* res = (int*)malloc(sizeof(int)*n);

    int i;

    for(i = 2; i<n+1; i++){
        res[i-2] = i;
    }
    
    res[i-2] = -1;

    return res;
}

int* liste_realocation(int* Liste,int nbelem){

    int i = 0;

    int compteur = 0;

    int* res = (int*)malloc((sizeof(int))*(nbelem + 1));

    while (Liste[i] != -1){

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

    if(n==1)
        return Liste;

    if(n < 1)
        return NULL;

    int* res = (int*)calloc(max,sizeof(int));

    int i = 0;
    int compteur = 0;
   
    while (Liste[i] != -1){  

        if(Liste[i]%n != 0){
            res[compteur] = Liste[i];
            compteur++;
        }

        i++;
    }

    res[compteur] = -1;
    free(Liste);
    res = liste_realocation(res,compteur);
    
    return res;
}

void liste_affichage(int* Liste){

    int i = 0;

    while (Liste[i] != -1){
        printf("%d ",Liste[i]);
        i++;
    }
    
    printf("\n");
}

int* int_crible(int n){

    int* res = calloc(n,sizeof(int));
    int* Liste = Liste_deux_A_n(n);

    int i = 0;

    while (Liste[0] != -1){
        res[i] = Liste[0];
        Liste = liste_non_diviseur(Liste,Liste[0],n);
        i++;
    }

    res[i] = -1;
    res = liste_realocation(res,i);

    return res;
}

int* erastothene(int n){ // la fonction erastothene ne nous a pas servi dans le projet mais elle nous semblait utile pour la décompostion en facteur premiers d'un entier

    int taille = (int)n/2 + 1;
    
    int* res = calloc(taille,sizeof(int));
    int* Prem = int_crible(n);
    
    int compteur = 0;
    int i = 0;
    int reste = n;

    while (reste != 1){

        while ((reste % Prem[i] == 0)&&(reste !=1)){
            res[compteur] = Prem[i];
            compteur++;
            reste = reste/Prem[i];
        }
        
        i++;
    }
    
    res[compteur] = -1;
    free(Prem);
    res = liste_realocation(res,compteur);

    return res;
}

// Fonction de base sur polynome a coeff modulaire

polynom* poly_initialisation(int degre,int modulo){//créer un polynom

    polynom* res = (polynom*)malloc(sizeof(polynom));
    res->val = calloc(degre+1, sizeof(int));
    res->degre = degre;
    res->modulo = modulo;

    return res;

}

void poly_liberer(polynom* p){

    if(p){
        if(p->val)
            free(p->val);
        free(p);
    }
}

int poly_degre(polynom* p1){// met à jour le degré du polynom

    if(p1->degre == 0)
        return 0;

    int res = p1->degre;

    while (p1->val[res] == 0 && res > 0){
        res--;
    }
    
    return res;
}

polynom* poly_random(int degre, int modulo){//rend un random polynom (pour les tests)

    polynom* res = poly_initialisation(degre,modulo);

    res->degre = degre;
    res->modulo = modulo;
    res->val = mod_tab_random(degre,modulo);
    res->degre = poly_degre(res);

    return res;
}

polynom* poly_add(polynom* p1, polynom* p2){//poly_add modulaire de deux polynomes

    int degre = int_maximum(p1->degre,p2->degre);
    int min = int_minimum(p1->degre,p2->degre);

    polynom* res = poly_initialisation(degre,p1->modulo);

    int i = 0;

    while (i <= min){
        res->val[i] = mod_add(p1->val[i],p2->val[i],p1->modulo);
        i++;
    }
    
    if(p1->degre == p2->degre){
        res->degre = poly_degre(res);
        return res;
    }
    else if(p1->degre > p2->degre){

        while (i <= degre){
            res->val[i] = p1->val[i];
            i++;
        }

        res->degre = poly_degre(res);
        return res;
    }
    else{

        while (i <= degre){
            res->val[i] = p2->val[i];
            i++;
        }

        res->degre = poly_degre(res);
        return res;
    }
}

polynom* poly_sous(polynom* p1, polynom* p2){// poly_sous modulaire de deux polynomes

    // t1 - t2

    if(p1->modulo != p2->modulo){
        printf("Les polynomes n'ont pas des coefficient d'un meme corps \n");
        return NULL;
    }


    int modulo = p1->modulo;
    int degre = int_maximum(p1->degre,p2->degre);
    int min = int_minimum(p1->degre,p2->degre);

    polynom* res = poly_initialisation(degre,modulo);
    
    for(int i=0; i<=min; i++){
        int p2val = p2->val[i];
        int val = mod_sous(p1->val[i],p2val,modulo);
        res->val[i] = val;
    }
    
    if(p1->degre == p2->degre){

        res->degre = poly_degre(res);
        return res;
    }
    else{    
        if(min == p1->degre){

            for(int i=0; i<=degre; i++){
                res->val[i] = mod_sous(0,p2->val[i],modulo);
            }
        }
        else{
            for(int i=0; i<=degre; i++){
                res->val[i] = p1->val[i];
            }
        }
    }

    res->degre = poly_degre(res);

    return res;
}

polynom* poly_mult(polynom* p1, polynom* p2){

    polynom* res = poly_initialisation(p1->degre + p2->degre,p1->modulo);

    int deg = p1->degre + p2->degre + 1;

    for(int i=0; i <= p1->degre; i++){

        for(int j=0; j <= p2->degre; j++){

            res->val[i + j] = mod_add(res->val[i+j],mod_mult(p1->val[i],p2->val[j],p1->modulo),p1->modulo);

        }
    }

    res->degre = poly_degre(res);

    return res;
}

polynom* poly_puissance(polynom* p, int puissance){

    polynom* res = poly_initialisation(p->degre, p->modulo);
    res->val = p->val;

    for(int i = 1; i<puissance; i++){
        res = poly_mult(res, p);
    }

    return res;
}

int poly_estNul(polynom* p){

    for(int i = 0; i<=p->degre; i++){
        if(p->val[i] != 0)
            return 0;
    }

    return 1;
}

int poly_estUn(polynom* p){

    if(p->val[0] != 1)
        return 0;

    for(int i = 1; i<p->degre+1; i++){
        if(p->val[i] != 0)
            return 0;
    }
    
    return 1;
}

void poly_affichage(polynom* p1){ //affiche le polynome (très mignonne)

    if(p1->degre == 0 && p1->val[0] == 0){
        printf("0\n");
        return;
    }

    for(int i = p1->degre; i >= 0; i--){

        if(p1->val[i] == 0){
        }
        else if(p1->val[i] == 1){
            if(i == p1->degre && p1->degre == 1){
                printf("x ");
            }
            else if(i == p1->degre && p1->degre >= 1){
                printf("x^%d ",i);
            }else if(i == 0){
                printf("+ 1");
            }else if(i == 1){
                printf("+ x ");
            }
            else{       
                printf("+ %dx^%d ",p1->val[i],i);
            }
        }
        else{
            if(i == p1->degre && p1->degre >= 1){
                printf("%dx^%d ",p1->val[i],i);
            }else if(i == 0){
                printf("+ %d",p1->val[i]);
            }else if(i == 1){
                printf("+ %dx ",p1->val[i]);
            }
            else{
                printf("+ %dx^%d ",p1->val[i],i);
            }
        }
    }

    printf("\n");
}

// Fonctions un peu plus avancées sur les polynomes

polynom* poly_ajuster_puissance(polynom* p, int newdegre, int coeff1, int coeff2, int* u, int* v){

    polynom* res = poly_initialisation(newdegre, p->modulo);

    int diffdeg = newdegre - p->degre;

    int inv = mod_inverse(coeff2, p->modulo, u, v);
    int multi = mod_mult(coeff1, inv, p->modulo);

    for(int i = 0; i <= p->degre; i++){
        int val = mod_mult(multi,p->val[i], p->modulo);
        res->val[i+diffdeg] = val;
    }
    
    res->degre = poly_degre(res);

    return res;
}

polynom* poly_reste_division(polynom* p1, polynom* p2, int* u, int* v){ //renvoie le reste de la division euclidienne de p1/p2

    //hypothèse :  deg p1 >= deg tp2

    polynom* reste = poly_initialisation(p1->degre, p1->modulo);
    reste->val = p1->val;

    int degreste = p1->degre;
    int deg2 = p2->degre;

    int coeff2 = p2->val[deg2];

    while (degreste >= deg2){

        int coeff1 = reste->val[degreste];
        
        polynom* b = poly_ajuster_puissance(p2, degreste, coeff1, coeff2, u, v);
        
        reste = poly_sous(reste,b);
        degreste = reste->degre;

        poly_liberer(b);
    }
    
    return reste;
}

polynom* poly_quotient_division(polynom* p1, polynom* p2, int* u, int* v){ //renvoie le quotient de la division euclidienne de p1/p2

    if(poly_estNul(p2)){
        printf("on ne peut pas diviser par le polynome nul\n");
        return NULL;
    }

    if(p1->degre < p2->degre){
        printf("on ne peut pas diviser par un polynome de degre supérieur\n");
        return NULL;
    }

    int degreste = p1->degre;
    int deg2 = p2->degre;
    int degQuotient = degreste - deg2;

    polynom* reste = poly_initialisation(degreste, p1->modulo);
    polynom* quotient = poly_initialisation(degQuotient, p1->modulo);

    reste->val = p1->val;

    int i = degQuotient;
    int cpt = 0;

    int coeff2 = p2->val[deg2];
    int invCoeff2 = mod_inverse(coeff2, p1->modulo, u, v);

    while (degreste >= deg2){

        int coeff1 = reste->val[degreste];

        int coeffQuotient = mod_mult(coeff1, invCoeff2, p1->modulo);
        quotient->val[i] = coeffQuotient;

        polynom* b = poly_ajuster_puissance(p2, degreste, coeff1, coeff2, u, v);
        reste = poly_sous(reste,b);
        degreste = reste->degre;

        poly_liberer(b);
        
        if(cpt == 1){
            break;
        }
        if(degreste == 0){
            cpt++;
        }

        i = degreste - deg2;
    }
    
    poly_liberer(reste);
    quotient->degre = poly_degre(quotient);

    return quotient;
}

polynom* poly_pgcd(polynom* p1,polynom* p2,int* u,int* v){

    polynom* a = poly_initialisation(p1->degre, p1->modulo);
    polynom* b = poly_initialisation(p2->degre, p2->modulo);

    a->val = p1->val;
    b->val = p2->val;

    int deg1 = p1->degre;
    int deg2 = p2->degre;

    while (deg2 > 0){
        polynom* temp = b;
        b = poly_reste_division(a,b,u,v);
        a = temp;
        deg1 = deg2;
        deg2 = poly_degre(b);
    }

    poly_liberer(b);
    
    return a;

}

polynom* poly_derive(polynom* p){

    if (p->degre==0)
        return poly_initialisation(0,p->modulo);

    polynom* p2 = poly_initialisation(p->degre-1,p->modulo);

    for(int i = 0; i<p->degre; i++){
        p2->val[i]=mod_mult((i+1), p->val[i+1], p->modulo);
    }

    p2->degre = poly_degre(p2);

    return p2;
}

int poly_estRacine(polynom* p, int racine){

    int res = 0;

    for(int i = 0; i<=p->degre; i++){
        int puis = mod_puissance(racine,i,p->modulo);
        int pval = p->val[i];
        int mod_multi = mod_mult(pval,puis,p->modulo);
        res = mod_add(mod_multi,res,p->modulo);
    }
    
    if(res == 0)
        return 1;
    
    return 0;
}

polynom* poly_facteurRacine(int racine, int modulo){

    polynom* res = poly_initialisation(1,modulo);

    res->val[1] = 1;

    if(racine == 0){
        res->val[0] = 0;
        return res; 
    }

    res->val[0] = modulo - racine;

    return res;
}

// Algorithmes de factorisation 

polynom** algo_naif(polynom* p){ // mais pas trop quand meme

    polynom** res = (polynom**)malloc(sizeof(polynom*)*p->degre);
    polynom* p2 = poly_initialisation(p->degre, p->modulo);
    p2->val = p->val;
    int* racine = malloc(sizeof(int)*p->degre);

    int compteur = 0;
    int i = 0;
    int u,v;
    int degre = p->degre;

    while (i <= p->modulo && compteur != p->degre){
        
        while(poly_estRacine(p2,i) == 1){
            res[compteur] = poly_facteurRacine(i,p->modulo);
            racine[compteur] = i;
            p2 = poly_quotient_division(p2,res[compteur],&u,&v);
            compteur++;
        }
        
        i++;
    }

    if(compteur != degre){
        poly_liberer(p2);
        free(res);
        printf("le polynome n'est pas scindable dans notre corps\n");
        return NULL;
    }
    
    poly_liberer(p2);
    return res;
}

polynom** algo_yun(polynom* p, int* nbRacine){// H tout les ai sont supposés premiers entre eux pensez à faire le teste en début de fonction

    polynom** res = (polynom**)malloc((p->degre)*sizeof(polynom*));

    int u,v;
    int compteur = 0;
    int nbRac = 0;

    polynom* pDerive = poly_derive(p);
    polynom* A = poly_pgcd(p,pDerive,&u,&v);
    polynom* B = poly_quotient_division(p,A,&u,&v);
    polynom* C = poly_quotient_division(pDerive,A,&u,&v);
    polynom* bDerive = poly_derive(B);
    polynom* D = poly_sous(C,bDerive);

    while (B->degre > 0){ // il faut faire degP tours de boucles

        A = poly_pgcd(B,D,&u,&v);

        polynom* Abis = poly_initialisation(0, A->modulo);
        Abis->val[0] = A->val[A->degre];

        // ajout de ai dans le résultat
        res[compteur] = poly_quotient_division(A,Abis,&u,&v);;   
        nbRac++;     
        compteur++;

        poly_liberer(Abis);

        B = poly_quotient_division(B,A,&u,&v);

        if(poly_estUn(B)){
            break;
        }

        C = poly_quotient_division(D,A,&u,&v);

        bDerive = poly_derive(B);

        D = poly_sous(C,bDerive);
    }
    
    *nbRacine = nbRac;
    return res;
}

// Un joli main

int main(int argc, char const *argv[]){

    srand(time(NULL));
    
    polynom* p1 = poly_initialisation(5,5);
    polynom* p2 = poly_initialisation(3,5);

    int t1[6] = {2,2,3,4,3,4};
    int t2[4] = {2,4,3,2};

    p1->val = t1;
    p2->val = t2;

    printf("\nP1 : ");
    poly_affichage(p1);
    printf("P2 : ");    
    poly_affichage(p2);

    polynom* p3 = poly_add(p1,p2);
    printf("\naddition P1 + P2 : ");
    poly_affichage(p3);

    polynom* p4 = poly_mult(p1,p2);
    printf("multiplication P1*P2 : ");
    poly_affichage(p4);

    polynom* p5 = poly_sous(p1,p2);
    printf("soustraction P1 - P2 : ");
    poly_affichage(p5);

    polynom* p6 = poly_derive(p1);
    printf("Derivée de P1 : ");
    poly_affichage(p6);

    polynom* p7 = poly_derive(p2);
    printf("Derivée de P2 : ");
    poly_affichage(p7);
    
    int u = 0;
    int v = 0;

    polynom* p8 = poly_reste_division(p1, p2, &u, &v);
    printf("\nReste de la divison euclidienne P1/P2 : ");
    poly_affichage(p8);

    polynom* p9 = poly_quotient_division(p1, p2, &u, &v);
    printf("Quotient de la division euclidienne P1/P2 : ");
    poly_affichage(p9);

    polynom* p10 = poly_mult(p9, p2);
    printf("Multuplication du quotient avec P2 : ");
    poly_affichage(p10);

    polynom* p11 = poly_add(p10, p8);
    printf("Puis on ajoute le reste et on trouve P1 = ");
    poly_affichage(p11);

    polynom* p12 = poly_pgcd(p1, p2, &u, &v);
    printf("\nPGCD entre P1 et P2: ");
    poly_affichage(p12);

    printf("________________________________________\n\nDecompositon du polynome (modulo 7) ");
    polynom* p13 = poly_initialisation(4, 7);
    int tab1[5] = {0,0,6,0,1};
    p13->val = tab1;
    printf("P = ");
    poly_affichage(p13);
    

    printf("________________________________________\n\nAlgo NAIF : \n\n");

    polynom** tp1 = algo_naif(p13);
    for(int i=0; i<p13->degre; i++){
        printf("facteur %d : ", i);
        poly_affichage(tp1[i]);
    }

    for(int i = 0; i<p13->degre; i++){
        poly_liberer(tp1[i]);
    }

    printf("________________________________________\n\nAlgo de YUN : \n\n");

    int nbRacines;

    polynom** tp2 = algo_yun(p13, &nbRacines);
    polynom** tv1 = malloc(nbRacines*sizeof(polynom*));
    for(int i=0; i<nbRacines; i++){
        printf("A%d : ", i+1);
        poly_affichage(tp2[i]);
        tv1[i] = poly_puissance(tp2[i], i+1);
    }

    polynom* v1 = tv1[0];
    for(int i = 1; i<nbRacines; i++){
        v1 = poly_mult(v1, tv1[i]);
    }

    polynom* s1 = poly_sous(p13,v1);

    printf("\nverification : P - la multiplication des Ai = ");
    poly_affichage(s1);

    if(nbRacines > 1)
        poly_liberer(v1);

    for(int i = 0; i<nbRacines; i++){
        poly_liberer(tp2[i]);
        if(i == 0){
            free(tv1[i]);
        }
    }
    
    printf("________________________________________\n\nDecompositon du polynome (modulo 5) ");

    polynom* p14 = poly_initialisation(6, 5);
    int tab2[7] = {3,4,2,3,0,4,1};
    p14->val = tab2;
    printf("P = ");
    poly_affichage(p14);

    printf("________________________________________\n\nAlgo NAIF : \n\n");

    polynom** tp3 = algo_naif(p14);
    for(int i=0; i<p14->degre; i++){
        printf("facteur %d : ", i);
        poly_affichage(tp3[i]);
    }

    for(int i = 0; i<p14->degre; i++){
        poly_liberer(tp3[i]);
    }

    printf("________________________________________\n\nAlgo de YUN : \n\n");
 
    polynom** tp4 = algo_yun(p14, &nbRacines);
    polynom** tv2 = malloc(nbRacines*sizeof(polynom*));
    for(int i=0; i<nbRacines; i++){
        printf("A%d : ", i+1);
        poly_affichage(tp4[i]);
        tv2[i] = poly_puissance(tp4[i], i+1);
    }

    polynom* v2 = tv2[0];
    for(int i = 1; i<nbRacines; i++){
        v2 = poly_mult(v2, tv2[i]);
    }
    polynom* s2 = poly_sous(p14,v2);

    printf("\nverification : P - la multiplication des Ai = ");
    poly_affichage(s2);

    if(nbRacines > 1)
        poly_liberer(v2);

    for(int i = 0; i<nbRacines; i++){
        poly_liberer(tp4[i]);
        if(i == 0){
            free(tv2[i]);
        }    
    }

    printf("________________________________________\n\nDecompositon du polynome (modulo 13) ");

    polynom* p15 = poly_initialisation(10, 13);
    int tab3[11] = {10,1,10,3,1,11,11,11,10,4,1};
    p15->val = tab3;
    printf("P = ");
    poly_affichage(p15);

    printf("________________________________________\n\nAlgo NAIF : \n\n");

    polynom** tp5 = algo_naif(p15);
    for(int i=0; i<p15->degre; i++){
        printf("facteur %d : ", i);
        poly_affichage(tp5[i]);
    }

    for(int i = 0; i<p15->degre; i++){
        poly_liberer(tp5[i]);
    }

    printf("________________________________________\n\nAlgo de YUN : \n\n");
 
    polynom** tp6 = algo_yun(p15, &nbRacines);
    polynom** tv3 = malloc(nbRacines*sizeof(polynom*));
    for(int i=0; i<nbRacines; i++){
        printf("A%d : ", i+1);
        poly_affichage(tp6[i]);
        tv3[i] = poly_puissance(tp6[i], i+1);
    }

    polynom* v3 = tv3[0];
    for(int i = 1; i<nbRacines; i++){
        v3 = poly_mult(v3, tv3[i]);
    }

    polynom* s3 = poly_sous(p15,v3);

    printf("\nverification : P - la multiplication des Ai = ");
    poly_affichage(s3);

    if(nbRacines > 1)
        poly_liberer(v3);

    for(int i = 0; i<nbRacines; i++){
        poly_liberer(tp6[i]);
        if(i == 0){
            free(tv3[i]);
        }
    }

    printf("________________________________________\n\nDecompositon du polynome (modulo 53) ");

    polynom* p16 = poly_initialisation(15, 53);
    int tab4[16] = {0,0,0,0,0,4,36,50,42,47,10,5,48,29,30,1};
    p16->val = tab4;
    printf("P = ");
    poly_affichage(p16);

    printf("________________________________________\n\nAlgo NAIF : \n\n");

    polynom** tp7 = algo_naif(p16);
    for(int i=0; i<p16->degre; i++){
        printf("facteur %d : ", i);
        poly_affichage(tp7[i]);
    }

    for(int i = 0; i<p16->degre; i++){
        poly_liberer(tp7[i]);
    }

    printf("________________________________________\n\nAlgo de YUN : \n\n");
 
    polynom** tp8 = algo_yun(p16, &nbRacines);
    polynom** tv4 = malloc(nbRacines*sizeof(polynom*));
    for(int i=0; i<nbRacines; i++){
        printf("A%d : ", i+1);
        poly_affichage(tp8[i]);
        tv4[i] = poly_puissance(tp8[i], i+1);
    }

    polynom* v4 = tv4[0];
    for(int i = 1; i<nbRacines; i++){
        v4 = poly_mult(v4, tv4[i]);
    }

    polynom* s4 = poly_sous(p16,v4);

    printf("\nverification : P - la multiplication des Ai = ");
    poly_affichage(s4);

    if(nbRacines > 1)
        poly_liberer(v4);

    for(int i = 0; i<nbRacines; i++){
        poly_liberer(tp8[i]);
        if(i == 0)
            free(tv4[i]);
    }

    printf("\n");

    free(p1);
    free(p2);
    poly_liberer(p3);
    poly_liberer(p4);
    poly_liberer(p5);
    poly_liberer(p6);
    poly_liberer(p7);
    poly_liberer(p8);
    poly_liberer(p9);
    poly_liberer(p10);
    poly_liberer(p11);
    poly_liberer(p12);
    free(p13);
    free(p14);
    free(p15);
    free(p16);
    free(tp1);
    free(tp2);
    free(tp3);
    free(tp4);
    free(tp5);
    free(tp6);
    free(tp7);
    free(tp8);
    free(tv1);
    free(tv2);
    free(tv3);
    free(tv4);
    poly_liberer(s1);
    poly_liberer(s2);
    poly_liberer(s3);
    poly_liberer(s4);

    return 0;
}

