// Implementation of the framework as described in https://arxiv.org/pdf/1706.10268.pdf.
// It is built upon Thaler's IP protocol for matrix-matrix multiplication available at 
// http://people.cs.georgetown.edu/jthaler/Tcode.htm
// This work is licensed under CC BY-NC-SA 3.0. Refer to the licesne file for more information.

#include <cstdlib>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>

#include <sstream>
#include <vector>

#define MASK 4294967295 //2^32-1
#define PRIME 2305843009213693951 //2^61-1

typedef unsigned long long uint64;

using namespace std;


//computes x^b
uint64 myPow(uint64 x, uint64 b) 
{
    uint64 i = 1;
    for (int j = 0; j <b; j++)  i *= x;
    return i;
}

//efficient modular arithmetic function for p=2^61-1. Only works for this value of p.
//This function might
//return a number slightly greater than p (possibly by an additive factor of 8);
//It'd be cleaner to check if the answer is greater than p
//and if so subtract p, but I don't want to pay that
//efficiency hit here, so the user should just be aware of this.
uint64 myMod(uint64 x)
{
    return (x >> 61) + (x & PRIME);
}  

//efficient modular multiplication function mod 2^61-1
inline uint64 myModMult(uint64 x, uint64 y)
{
    uint64 hi_x = x >> 32; 
    uint64 hi_y = y >> 32;
    uint64 low_x = x & MASK;
    uint64 low_y = y & MASK;

    //since myMod might return something slightly large than 2^61-1,  
    //we need to multiply by 8 in two pieces to avoid overflow.
    uint64 piece1 = myMod((hi_x * hi_y)<< 3);
    uint64 z = (hi_x * low_y + hi_y * low_x);
    uint64 hi_z = z >> 32;
    uint64 low_z = z & MASK;

    //Note 2^64 mod (2^61-1) is 8
    uint64 piece2 = myMod((hi_z<<3) + myMod((low_z << 32)));
    uint64 piece3 = myMod(low_x * low_y);
    uint64 result = myMod(piece1 + piece2 + piece3);

    return result;
}

//computes b^e mod p using repeated squaring. p should be 2^61-1
uint64 myModPow(uint64 b, uint64 e)
{
    uint64 result;
    if(e==1) return b;
    if(e == 0) return 1;
    if((e & 1) == 0)
    {
        result = myModPow(b, (e>>1));
        return myModMult(result, result);
    }
    else
        return myModMult(myModPow(b, e-1), b);
}

//Performs Extended Euclidean Algorithm
//Used for computing multiplicative inverses mod p
void extEuclideanAlg(uint64 u, uint64* u1, uint64* u2, uint64* u3)
{
    *u1 = 1;
    *u2 = 0;
    *u3 = u;
    uint64 v1 = 0;
    uint64 v2 = 1;
    uint64 v3 = PRIME;
    uint64 q;
    uint64 t1;
    uint64 t2;
    uint64 t3;
    do
    {
        q = *u3 / v3;
        //t1 = *u1 + p - q * v1;
        //t2 = *u2 + p - q*v2;
        //t3 = *u3 + p - q*v3;
        t1 = myMod((*u1) + PRIME - myModMult(q, v1));
        t2 = myMod((*u2) + PRIME - myModMult(q, v2));
        t3 = myMod((*u3) + PRIME - myModMult(q, v3));
        (*u1) = v1;
        (*u2) = v2;
        (*u3) = v3;
        v1 = t1;
        v2 = t2;
        v3 = t3;
    }while(v3 != 0 && v3!= PRIME);
}


//Computes the modular multiplicative inverse of a modulo m,
//using the extended Euclidean algorithm
//only works for p=2^61-1
uint64 inv(uint64 a)
{
    uint64 u1;
    uint64 u2;
    uint64 u3;
    extEuclideanAlg(a, &u1, &u2, &u3);
    if(u3==1)
        return myMod(u1);
    else
        return 0;
}

//computes chi_v(r), where chi is the Lagrange polynomial that takes
//boolean vector v to 1 and all other boolean vectors to 0. (we view v's bits as defining
//a boolean vector. n is dimension of this vector.
//all arithmetic done mod p
uint64 chi(uint64 v, uint64* r, uint64 n)
{ 
    uint64 x=v;
    uint64 c = 1;
    for(uint64 i = 0; i <n; i++)
    {
        if( x&1 )
            c=myModMult(c, r[i]);
        else
            c=myModMult(c, 1+PRIME-r[i]);
        x=x>>1;
    }
    return c;
}

//evaluates V_i polynomial at location r.
//Here V_i is described in GKR08; it is the multi-linear extension
//of the vector of gate values at level i of the circuit
uint64 evaluate_V_i(int mi, int ni, uint64* level_i, uint64* r)
{
    uint64 ans=0;
    for(uint64 k = 0; k < ni; k++)
        ans=myMod(ans + myModMult(level_i[k], chi(k, r, mi)));
    return ans;
}

//extrapolate the polynomial implied by vector vec of length n to location r
uint64 extrap(uint64* vec, uint64 n, uint64 r)
{
    uint64 result=0;
    uint64 mult=1;
    for(uint64 i = 0; i < n; i++)
    {
        mult=1;
        for(uint64 j=0; j < n; j++)
        {
            if (i>j)
                mult=myModMult(myModMult(mult, myMod(r-j+PRIME)), inv(i-j) );
            if (i<j)
                mult=myModMult(myModMult(mult, myMod(r-j+PRIME)), inv(myMod(i+PRIME-j)) );
        }
        result=myMod(result+myModMult(mult, vec[i]));
    }
    return result;
}

//fills in high-order variable xi wih ri
void updateV(uint64* V, int num_new, uint64 ri)
{
    for(int i = 0; i < num_new; i++)
        V[i] = myMod(myModMult(V[i], 1+PRIME-ri) + myModMult(V[i+num_new], ri));
}

// evaluate the MLE of I at q and random location r in O(logn)
uint64 evaluate_I(uint64* q, uint64* r, int d)
{
    uint64 ans=1;
    for(uint64 k = 0; k < d; k++)
        ans = myModMult(ans, myMod(myModMult(q[k],r[k]) + myModMult(1+PRIME-q[k], 1+PRIME-r[k])) );
    return ans; 
}

struct runtime {
    double unverifiable;
    double prover;
    double verifier;
};

runtime update_time(runtime t, runtime nt)
{
    t.unverifiable += nt.unverifiable;
    t.prover += nt.prover;
    t.verifier += nt.verifier;
    return t;
}

runtime set_time(runtime t, double ut, double pt, double vt)
{
    t.unverifiable = ut;
    t.prover = pt;
    t.verifier = vt;
    return t;
}

void check_bias_layer(uint64* q, uint64* r, int d, uint64 n, uint64* Iin, uint64* Vin, uint64* B, uint64** F, uint64* check)
{
    //initialize Iin values
    uint64 Iin_tmp;
    uint64 steps=1;
    for (int i=0; i<d; i++)
    {
        for (int k=0; k<steps; k++)
        {
            Iin_tmp = Iin[k];

            Iin[k] = myModMult(Iin_tmp, 1+PRIME-q[i]);
            Iin[k+steps] = myModMult(Iin_tmp, q[i]);
           
        }
        steps = steps << 1;
    }

    uint64* S = (uint64*) calloc(n, sizeof(uint64));
    for (int i=0; i<n; i++)
        S[i] = myMod(Vin[i] + B[i]);

    steps=n;
    uint64 temp0; uint64 temp1; uint64 cross;
    for (int i=0; i<d; i++)
    {
        steps = steps >> 1;
        for (int k=0; k<steps; k++)
        {
            temp0 = myModMult(Iin[k], S[k]);
            temp1 = myModMult(Iin[k+steps], S[k+steps]);
            cross = myModMult(myMod(PRIME - Iin[k] + 2*Iin[k+steps]), myMod(PRIME - S[k] + 2*S[k+steps]));

            F[i][0] = myMod(F[i][0]+temp0);
            F[i][1] = myMod(F[i][1]+temp1);
            F[i][2] = myMod(F[i][2]+cross);
        }
        updateV(Iin, steps, r[d-1-i]);
	updateV(S, steps, r[d-1-i]);

        check[i] = extrap(F[i], 3, r[d-1-i]);
    }
}


runtime verify_bias(int d, int i, int L)
{
    uint64 n = myPow(2,d);
  //inputs to activation layer
    uint64* Vin = (uint64*) malloc(n*sizeof(uint64));
    for (int i=0; i<n; i++)
        Vin[i] = abs(rand()) % 100;

    //activations
    uint64* B = (uint64*) malloc(n*sizeof(uint64));
    for (int i=0; i<n; i++)
        B[i] = abs(rand()) % 100;

    //Iin values
    uint64* Iin = (uint64*) malloc(n*sizeof(uint64));
    for (int i=0; i<n; i++)
        Iin[i]=1;
    
   
    uint64** F = (uint64**) calloc(d, sizeof(uint64*));
    for (int i=0; i<d; i++)
        F[i] = (uint64*) calloc(3, sizeof(uint64));
 
    // bias layer    
    uint64* S = (uint64*) calloc(n, sizeof(uint64));
    
    // evaludate the layer
    clock_t t=clock();
    for (int i=0; i<n; i++)
        S[i] = myMod(Vin[i]+B[i]);

    t = clock()-t;
    double ut = (double)t/CLOCKS_PER_SEC;
    cout << "unverifiable time for bias = " << ut << endl;

    uint64* r = (uint64*) calloc(d, sizeof(uint64));
   // r is Vin's random coin tosses for this iteration.
    // really should choose r <--F_p, but this chooses r <-- [2^32]
    for (int i=0; i<d; i++)
        r[i] = abs(rand());

    // calculations of Fi(ri) for checking    
    uint64* check = (uint64*) calloc(d, sizeof(uint64));
   
    uint64* q = (uint64*) calloc(d, sizeof(uint64));
    for (int i=0; i<d; i++)
        q[i] = abs(rand());
    
    uint64 Ieval=0;
    uint64 Vieval=0;
    uint64 Beval=0;
    uint64 a1 = 0;          // ai-1
    uint64 a2 = 0;          // ai    

    // At the output layer, verifier evaluates a random point in the MLE of the returned matrix
    // For middle layers, this assertion is returned by prover
    clock_t otime = clock();
    a1 = evaluate_V_i(d, n, S, q);
    otime = clock()-otime;
    
    t=clock();
    check_bias_layer(q, r, d, n, Iin, Vin, B, F, check);
    t = clock() - t;
    if (i!=L-1)
        t+=otime;
    double pt = (double)t/CLOCKS_PER_SEC;
    cout << "additional prover time = " << pt << endl;

    // assertion about the input of this layer returned by the prover (output of mm mult layer)
    Vieval = evaluate_V_i(d, n, Vin, r); 

    t=clock();
    if (a1 != myMod(F[0][0]+F[0][1]))
        cout << "bias layer first check failed" << endl, exit(1);

    for (int i=1; i<d; i++)
    {
        if ((myMod(F[i][0] + F[i][1]) != check[i-1]) && (myMod(F[i][0] + F[i][1]) + PRIME != check[i-1]))
            cout<< "bias layer check " << i << " failed" << endl, exit(1);
    }

    Ieval = evaluate_I(q,r,d);
    Beval = evaluate_V_i(d, n, B, r);

    //last check
    a2 = myModMult(myMod(Vieval + Beval), Ieval);
    
    if (a2 != check[d-1])
        cout << "bias layer last check failed" << endl, exit(1);
    
    t = clock() - t;
    if (i==L-1)
        t += otime;
    double vt = (double)t/CLOCKS_PER_SEC;
    cout << "verifier time = " << vt << endl;

 
    free(Vin);
    free(Iin);
    free(B);
    for (int i=0; i<d; i++)
        free(F[i]);
    free(F);
    free(S);
    free(r);
    free(q);
    free(check);
    
    runtime bias_runtime;
    return set_time(bias_runtime, ut, pt, vt);
}


//V0 is supposed to be list of vals of matrix A in row-major order, V1 vals of matrix B.
//This function has been modified to incorporate matrix-matrix mult of size (m,n)*(n,p)
void sum_check_mm(uint64* V0, uint64* V1, int d, int e, int f, int mi, int ni, uint64*r, uint64** F, uint64* z, uint64* check)
{

    for(int i = 0; i < f+e; i++)
        r[d+i] = z[i];

    for(int i = 0; i < d; i++)
        r[i] = abs(rand()) + 3;

    
    int num_terms = mi;

    for(int round = 0; round < e; round++)
    {
        updateV(V0, num_terms >> 1, r[f+d+e-1-round]);
        num_terms = num_terms >> 1;
    }

    num_terms = ni;
    for(int round = e; round < f+e; round++)
    {
        updateV(V1, num_terms >> 1, r[f+d+e-1-round]);
        num_terms = num_terms >> 1;
    }

    uint64 temp1; uint64 temp0; uint64 cross;
    for(int round = 0; round < d; round++)
    {
        for(int i = 0; i < (num_terms >> 1); i++)
        {
            temp0 = myModMult(V0[i], V1[i]);
            temp1 = myModMult(V0[i + (num_terms>>1)], V1[i + (num_terms>>1)]);
            cross = myModMult( myMod(PRIME - V0[i] +  2*V0[i + (num_terms>>1)]), myMod(PRIME - V1[i] + 2*V1[i + (num_terms>>1)]));

            F[round][0] = myMod(F[round][0] + temp0);
            F[round][1] = myMod(F[round][1] + temp1);
            F[round][2] = myMod(F[round][2] + cross);
        }		
        updateV(V0, num_terms >> 1, r[d-1-round]);
        updateV(V1, num_terms >> 1, r[d-1-round]);
        num_terms = num_terms >> 1;

        check[round] = extrap(F[round], 3, r[d-1-round]); 
    }
}


runtime verify_mm(int e, int d, int f, int i, int L)
{
    uint64 n = myPow(2, d);
    uint64 m = myPow(2, e);
    uint64 p = myPow(2, f);

    uint64* V = (uint64*) malloc((m*n+n*p)*sizeof(uint64));
    for(int i = 0; i < m*n+n*p; i++)
        V[i] = abs(rand()) % 100;

    uint64* C = (uint64*) calloc(m*p, sizeof(uint64));

    uint64* Vcopy = (uint64*) malloc((m*n+n*p)*sizeof(uint64));

    uint64* z = (uint64*) calloc(f+d+e, sizeof(uint64));
    uint64* r = (uint64*) calloc(f+d+e, sizeof(uint64));

    for(int i = 0; i < f+e; i++)
        z[i] = abs(rand())+3;

    uint64** F = (uint64**) calloc((d), sizeof(uint64*));
    for(int i = 0; i < d; i++)
        F[i] = (uint64*) calloc(4, sizeof(uint64));


    clock_t t=clock();
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < p; j++)
        {
            for(int k = 0; k < n; k++)
            {
                C[i*p+j] = myMod( C[i*p+j] + myModMult(V[i*n+k], V[j*n+k+m*n]));
            }
        }
    }
    t=clock()-t;
    double ut = (double)t/CLOCKS_PER_SEC;
    cout << "unverifiable time for matrix-matrix mult = " << ut << endl;

    uint64 a1=0;    //ai-1
    uint64 a2=0;    //ai

    // calculations of Fi(ri) for checking    
    uint64* check = (uint64*) calloc(d, sizeof(uint64));

    for(int j = 0; j < m*n+n*p; j++)
        Vcopy[j] = V[j];


    t=clock();
    // prover evaluates the output of the mm mult layer (input to bias layer)
    a1 = evaluate_V_i(f+e, m*p, C, z);

    sum_check_mm(V, V+ m*n, d, e, f, m*n, n*p, r, F, z, check);
    t = clock()-t;
    double pt = (double)t/CLOCKS_PER_SEC;
    cout << "additional P time = " << pt << endl;

    // set the high order of values to be those of corresponding to index i, and the low order values of z to be those corresponding to index k
    for(int i = 0; i < d; i++)
        z[i] = r[i]; //set the low-order values of z
    for(int i = d; i < d+e; i++)
        z[i] = r[f+i]; //set the low-order values of z	

    // assertion about the input of this layer returned by the prover (output of sqr activation layer)
    // when reaching first layer, this is evaluated by the verifer
    clock_t itime = clock();
    uint64 Aeval = evaluate_V_i(d+e, m*n, Vcopy, z);
    itime = clock()-itime;
    
    t=clock();	
    if (a1 != myMod(F[0][0]+F[0][1]))
        cout << "matrix-matrix mult layer first check failed" << endl, exit(1);

    for (int i=1; i<d; i++)
    {
        if ((myMod(F[i][0] + F[i][1]) != check[i-1]) && (myMod(F[i][0] + F[i][1]) + PRIME != check[i-1]))
            cout<< "matrix-matrix mult layer check " << i << " failed" << endl, exit(1);
    }

    // Beval corresponds to layer weight (w), which the verifier evaluates
    uint64 Beval = evaluate_V_i(d+f, n*p, &(Vcopy[m*n]), r);

    a2 = myModMult(Aeval, Beval);

    if (a2 != check[d-1])
        cout  << "matrix-matrix mult layer last check failed" << endl, exit(1);

    t = clock()-t;

    // V evaluates the MLE of input for the first layer
    if (i==0)
        t += itime;

    double vt = (double)t/CLOCKS_PER_SEC;
    cout << "verifier time = " << vt << endl;
        
    free(V);
    free(Vcopy);
    free(C);
    free(z);
    free(r);
    for(int i = 0; i < d; i++)
        free(F[i]);
    free(F);
    free(check);

    runtime mm_runtime;
    return set_time(mm_runtime, ut, pt, vt);
}


// Protocol reduces verifying a claim that v_i-1(q)=a_i-1 to verifying that v_i(q')=a_i
void sum_check_sqr_activation(uint64* q, uint64* r, int d, uint64 n, uint64* Iin, uint64* I_t, uint64* Vin, uint64* V_t, uint64* K_t, uint64** F, uint64* check)
{
    //initialize Iin values
    uint64 Iin_tmp;
    uint64 steps=1;
    for (int i=0; i<d; i++)
    {
        for (int k=0; k<steps; k++)
        {
            Iin_tmp = Iin[k];

            Iin[k] = myModMult(Iin_tmp, 1+PRIME-q[i]);
            Iin[k+steps] = myModMult(Iin_tmp, q[i]);

        }

        steps = steps << 1;
    }

    for (int i=0; i<n; i++)
    {
        V_t[i] = Vin[i];
        I_t[i] = Iin[i];
    }

    // partial sums for calculating F at each round
    uint64* parsumV=(uint64*) calloc(4, sizeof(uint64));
    uint64* parsumI=(uint64*) calloc(4, sizeof(uint64));

    steps=n;
    int j;

    for (int i=0; i<d; i++)
    {
        steps = steps >> 1;
        for (int k=0; k<steps; k++)
        {
            j = 2*k;

            parsumV[0] = V_t[j];
            parsumV[1] = V_t[j+1];
            parsumV[2] = myMod(2*V_t[j+1] + PRIME - V_t[j]);
            parsumV[3] = myMod(3*V_t[j+1] + 2*(PRIME - V_t[j]));

            parsumI[0] = I_t[j];
            parsumI[1] = I_t[j+1];
            parsumI[2] = myMod(2*I_t[j+1] + PRIME - I_t[j]);
            parsumI[3] = myMod(3*I_t[j+1] + 2*(PRIME - I_t[j]));

            V_t[j+1]= myModMult(V_t[j+1], r[i]);
            I_t[j+1]= myModMult(I_t[j+1], r[i]);

            V_t[j]= myModMult(V_t[j], 1+PRIME-r[i]);
            I_t[j]= myModMult(I_t[j], 1+PRIME-r[i]);

            V_t[k] = myMod(V_t[j] + V_t[j+1]);
            I_t[k] = myMod(I_t[j] + I_t[j+1]);


            for (int m=0; m<4; m++)
            {
                F[i][m] = myMod(F[i][m]+myModMult(myModMult(parsumV[m],parsumV[m]),parsumI[m]));
                parsumV[m]=0;
                parsumI[m]=0;
            }
        }

        //calculate Fi(ri) 
        check[i] = extrap(F[i], 4, r[i]);
    }

    free(parsumV);
    free(parsumI);
}

runtime verify_sqr_activation(int d)
{
    uint64 n = myPow(2,d);

    // indeces of the Vs in table
    uint64* K_t = (uint64*) calloc(n, sizeof(uint64));
    uint64 k;
    for (k=0; k<n; k++)
        K_t[k] = k;

    //inputs to activation layer
    uint64* Vin = (uint64*) malloc(n*sizeof(uint64));
    for (int i=0; i<n; i++)
        Vin[i] = abs(rand()) % 100;

    // table for V_tilda holding contributions of initial Vs at each round, updated every round
    uint64* V_t = (uint64*) calloc(n, sizeof(uint64));

    //Iin values
    uint64* Iin = (uint64*) malloc(n*sizeof(uint64));
    for (int i=0; i<n; i++)
        Iin[i]=1;
    uint64* I_t = (uint64*) calloc(n, sizeof(uint64));

    uint64** F = (uint64**) calloc(d, sizeof(uint64*));
    for (int i=0; i<d; i++)
        F[i] = (uint64*) calloc(4, sizeof(uint64));

    // square activation layer    
    uint64* A = (uint64*) calloc(n, sizeof(uint64));

    clock_t t=clock();
    for (int i=0; i<n; i++)
        A[i] = myModMult(Vin[i],Vin[i]);
    t = clock()-t;
    double ut = (double)t/CLOCKS_PER_SEC;
    cout << "unverifiable time for sqr activation = " << ut << endl;

    uint64* r = (uint64*) calloc(d, sizeof(uint64));
    // r is Vin's random coin tosses for this iteration.
    // really should choose r <--F_p, but this chooses r <-- [2^32]
    for (int i=0; i<d; i++)
        r[i] = abs(rand());

    // calculations of Fi(ri) for checking    
    uint64* check = (uint64*) calloc(d, sizeof(uint64));

    uint64* q = (uint64*) calloc(d,sizeof(uint64));
    for (int i=0; i<d; i++)
        q[i] = abs(rand());

    uint64 Ieval=0;
    uint64 Vieval=0;
    uint64 a1 = 0;          // ai-1
    uint64 a2 = 0;          // ai    
    
    t=clock();
    // prover evaluates the output of the sqr activation layer (input to mm mult layer)
    a1 = evaluate_V_i(d, n, A, q);

    sum_check_sqr_activation(q, r, d, n, Iin, I_t, Vin, V_t, K_t, F, check);
    t = clock() - t;
    double pt = (double)t/CLOCKS_PER_SEC;
    cout << "additional prover time = " << pt << endl;

    // assertion about the input of this layer returned by the prover (output of bias layer)
    Vieval = evaluate_V_i(d,n,Vin, r);

    clock_t v_t=clock();
    if (a1 != myMod(F[0][0]+F[0][1]))
        cout << "square activation layer first check failed" << endl, exit(1);

    for (int i=1; i<d; i++)
    {
        if ((myMod(F[i][0] + F[i][1]) != check[i-1]) && (myMod(F[i][0] + F[i][1]) + PRIME != check[i-1]))
            cout<< "square activation layer check " << i << " failed." << endl, exit(1);
    }

    Ieval = evaluate_I(q,r,d);

    //last check
    a2 = myModMult(myModMult(Vieval, Vieval), Ieval);
    if (a2 != check[d-1])
        cout << "square activation layer last check failed" << endl, exit(1);

    v_t = clock() - v_t;
    double vt = (double)v_t/CLOCKS_PER_SEC;
    cout << "verifier time = " << vt << endl;

    free(Vin);
    free(Iin);
    free(K_t);
    free(V_t);
    free(I_t);
    for (int i=0; i<d; i++)
        free(F[i]);
    free(F);
    free(A);
    free(r);
    free(check);
    free(q);
    
    runtime sqr_runtime;
    return set_time(sqr_runtime, ut, pt, vt);
}


int main(int argc, char** argv)
{
    if (argc!=2)
        cout << "Enter the architecture file as argument." << endl, exit(1);
    ifstream archfile(argv[1]);

    string line;
    vector <int*> layers;
    
    // get batch size
    double batch;
    getline(archfile, line);
    stringstream(line) >> batch;
    batch = ceil(log2(batch));

    // read input size
    double prevl, currl;
    getline(archfile, line);
    stringstream(line) >> prevl;
    prevl = ceil(log2(prevl));

    // read layer sizes
    while (getline(archfile, line))
    {
        stringstream(line) >> currl;
        currl = ceil(log2(currl));
        int *n = new int[3];
        n[0] = batch;
        n[1] = prevl;
        n[2] = currl;
        layers.push_back(n);
        prevl = currl;
    }

    int L = layers.size();
    int e, d, f;
    runtime verify_time;
    runtime total_time;

    total_time = set_time(total_time, 0, 0, 0);

    cout << "Verifying the neural network layer by layer:" << endl;

    // Verify MM and activation layers
    for (int i=L-1; i>=0; i--)
    {
        cout << "======== Layer " << i+1 << " verification =======" << endl;
        e = layers[i][0];
        d = layers[i][1];
        f = layers[i][2];
        
        // no activation in the last layer
        if (i!=L-1)
        {
            verify_time = verify_sqr_activation(e+f);
            cout <<"\tsqr activation verification done." << endl;
            total_time = update_time(total_time, verify_time); 
        }

        verify_time = verify_bias(e+f, i, L);
        cout <<"\tbias verification done." << endl;
        total_time = update_time(total_time, verify_time); 

        verify_time = verify_mm(e, d, f, i, L);
        cout <<"\tmatrix-matrix mult verification done." << endl;
        total_time = update_time(total_time, verify_time);  
        cout << endl;
    }

    cout << "total unverifiable time = " << total_time.unverifiable << endl;
    cout << "total additional prover time = " << total_time.prover << endl;
    cout << "total verifier time = " << total_time.verifier << endl;


    archfile.close();
    for (int i=0; i<layers.size(); i++)
        delete layers[i];

    return 0;
}

