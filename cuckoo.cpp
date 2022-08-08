#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <time.h>
#include <random>
#include <algorithm>

using namespace std;

float simplebounds(float s, float Lb, float Ub)
{if(s>Ub) {s=Ub;} if(s>Lb) {s=Lb;} return s;}

// Gamma function approximation
float gamma(float x)
{return sqrt(2.0*M_PI/x)*pow(x/M_E, x);}

//Objective function definition
float fobj(float x[], int size_array)
{float t= 0;
for(int i=0; i<size_array; i++){t=t+(x[i]-1)*(x[i]-1);}

return t;}

//Generate pseudo random numbers in the range [0,1)
float ra(){return(float ) rand() /RAND_MAX;}

int main(){
 //Normal distribution generator

std::default_random_engine generator;
std::normal_distribution<double> distribution(0, 1.0);
srand(time(NULL));
const int n=25; //Population size
float pa=0.25; //Rate of discovery of other birds eggs
const int nd=15; //dimensions
float Lb[nd]; //Lower bounds
float Ub[nd]; //Upper bounds
float fitness[n]; //Fitness values
float beta=3/2, sigma, fmin, best[nd], bestnest[nd];
int K, r1, r2, r, randperm1[n], randperm2[n];
for(int i=0; i<nd; i++){Lb[i]=-5; Ub[i]=5;}
int N_IterTotal=1000; //higher value for better results, longer computation
//Random initialization solutions
float nest[n][nd], new_nest[n][nd];

for(int i=0; i<n; i++){
for(int j=0; j<nd;j++){
nest[i][j]=Lb[j]+(Ub[j]-Lb[j])*ra();}}

for(int i=0;i<n; i++){
fitness[i] = fobj(nest[i],nd);
randperm1[i] = i; randperm2[i] = i;}

K=0; fmin=fitness[K];

for(int i=1; i<n; i++){
if(fitness[i]<=fmin){
fmin=fitness[i];
K=i;}}

for(int j=0; j<nd; j++){
best[j] = nest[K][j];}

//Starting iterations
float s, u, v, step, stepsize, fnew;
for(int iter=0; iter<N_IterTotal; iter++)
{//Keep current best, generate new solutions
sigma=pow(gamma(1+beta)* sin(M_PI*beta/2)/(gamma((1+beta)/2)*beta*pow
(2,(beta - 1)/2)),(1/beta));
for(int i=0; i<n; i++)
{for(int j=0; j<nd; j++){
s=nest[i][j];
//Mantegna's algo
u=distribution(generator) * sigma;
v=distribution(generator);
step=u/pow(abs(v),(1/beta));
stepsize = 0.01*step*(s-best[j]);
s=s+stepsize*distribution(generator);
//Apply limits
new_nest[i][j] = simplebounds(s, Lb[j], Ub[j]);}

fnew=fobj(new_nest[i], nd);
if(fnew<=fitness[i]){fitness[i]=fnew;
for(int j=0; j<nd; j++){
nest[i][j] = new_nest[i][j];}}}

//Get best nest
K=min_element(fitness,fitness + nd) - fitness;
fmin=fitness[K];
for(int j=0; j<nd; j++){best[j] = nest[K][j];}
//Random permutation
for(int i=0; i<n; i++){r1=rand()%n; r2=rand()%n;
if(r1!=r2){ r=randperm1[r1]; randperm1[r1] = randperm1[r2];
randperm1[r2] = r;}
r1 = rand()%n; r2 = rand()%n;
if(r1!=r2){r=randperm2[r1];randperm2[r1] = randperm2[r2];
randperm2[r2]=r;}
}
//Discovery, randomization
for(int i=0; i<n; i++)
{
for(int j=0; j<nd; j++)
{
stepsize=ra()*(nest[randperm1[i]][j]-nest[ randperm2[i] ][j] );
if(ra()>pa)
{
new_nest[i][j]=nest[i][j]+ stepsize;
}
else
{
new_nest[i][j]=nest[i][j];
}
new_nest[i][j]=simplebounds(new_nest[i][j] , Lb[j] , Ub[j]);
}
}
//Evaluate the new set of solutions
for(int i =0; i<n; i++)
{
fnew=fobj(new_nest[i] , nd);
if(fnew<=fitness[i] )
{fitness[i]=fnew;
for(int j =0; j<nd; j++)
{nest[i][j]=new_nest[i][j];}}}
//Get the bestnest
K=min_element (fitness,fitness+nd)-fitness;
fmin=fitness[K];
for(int j =0; j<nd; j++){ best[j]=nest[K][j];}
//Find the best objevtive so far
if(fnew<fmin)
{fmin=fnew;
for(int j =0; j<nd; j++)
{bestnest[j]=best[j];}}}
//Output the best solution
cout<<"Best solution=[";
for(int i =0; i<nd; i++) { cout<<bestnest[i]<<" ";}
cout<<" ] "<<endl;
cout<<"Best objective="<<fmin<<endl;
getchar();
return 0;
}