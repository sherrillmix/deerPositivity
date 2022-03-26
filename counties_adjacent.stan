
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
    for (i in 1:nCounty) {
      adj_rowsums[i] = sum(adjacency[i, ]);
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
transformed parameters{
  //vector[nCounty] countyProp;
  //countyProp=countyPropRaw*countySd;
}


model {
  positives ~ binomial_logit(counts,overallProp+countyProp);
  countySd~gamma(1,1);
  //https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html
  countyProp~multi_normal_prec(zeros,countySd*(D-alpha*adjacency));
}
