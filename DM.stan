// Model specification for Dirichlet-Multinomial, shamelessly edited from Harrison et al. 2020 MER
data {
  int<lower=1> N; // number of treatments, if categorical?
  int<lower=1> nreps; // number of genes
  int<lower=1> notus; // number of samples

  int<lower=1> start[N]; // needs to be not numeric, is this to loop through the data, is this the first gene
  int<lower=1> end[N]; //needs to be not numeric, last gene?

  int datamatrix[nreps, notus]; // this is the whole data frame, response variable?
}

### this is where we indicate the parameters that I want to estimate. Literally, just want to estiamte the betas for all of the expression~q
parameters { // basically want to have the intercept, slope and error for all of the parameters. 
  real<lower=0> theta[N];
  simplex[notus] pi[N];
  simplex[notus] p[nreps];
}


model { //Josh has two loops here and I'm not entirely sure why. Need to understand why there isn't just one loop, that loops through all the genes?
  for(i in 1:N){
    target += beta_lpdf(theta[i] | 5,5);
    target += dirichlet_lpdf(pi[i] | rep_vector(0.0000001, notus)); // don't know what this is supposed to look like, check Harrison's code
    for(j in start[i]:end[i]){
      target += dirichlet_lpdf(p[j] | 5, 5); //this is the response variable
      target += beta_lpdf(datamatrix[j,] | p[j]); //this is the predictor variable?
    }
  }
}
