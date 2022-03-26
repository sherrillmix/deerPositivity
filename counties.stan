data {
  int<lower=1> nCounty;
  int counts[nCounty];
  int positives[nCounty];
}

parameters{
  real overallProp;
  vector[nCounty] countyPropRaw;
  real<lower=0> countySd;
}
transformed parameters{
  vector[nCounty] countyProp;
  countyProp=countyPropRaw*countySd;
}


model {
  positives ~ binomial_logit(counts,overallProp+countyProp);
  countySd~gamma(1,1);
  countyPropRaw~normal(0,1);
}
