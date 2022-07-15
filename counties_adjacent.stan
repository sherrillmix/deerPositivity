
data {
  int<lower=1> nCounty;
  int counts[nCounty];
  int positives[nCounty];
  matrix<lower=0,upper=1>[nCounty,nCounty] adjacency;
}

transformed data{
  vector[nCounty] zeros;
  matrix<lower = 0>[nCounty, nCounty] D;
  {
    vector[nCounty] adj_rowsums;
    for (ii in 1:nCounty) {
      adj_rowsums[ii] = sum(adjacency[ii, ]);
    }
    D = diag_matrix(adj_rowsums);
  }
  zeros = rep_vector(0, nCounty);
}
parameters{
  real overallProp;
  vector[nCounty] countyProp;
  real<lower=0> countySd;
  real<lower=0,upper=1> alpha;
}

model {
  countySd~gamma(1,1);
  countyProp~multi_normal_prec(zeros,1/countySd*(D-alpha*adjacency));
  positives ~ binomial_logit(counts,overallProp+countyProp);
  overallProp~normal(-2,10);
}
