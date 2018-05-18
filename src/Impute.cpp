#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec calculate_weights(arma::vec z, arma::mat X){
  Environment myEnv("package:VIPER");
  Function calculate_weights_fun = myEnv["calculate_weights"];
  Rcpp::NumericVector calculate_weights_res = wrap(calculate_weights_fun(z, X));
  return calculate_weights_res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List fitting_lasso(arma::vec y, arma::mat X, bool min, double alpha){
  Environment myEnv("package:VIPER");
  Function fitting_lasso_fun = myEnv["fitting_lasso"];
  if(min){
    Rcpp::List fitting_lasso_res = wrap(fitting_lasso_fun(y, X, "min", alpha));
    return fitting_lasso_res;
  } else {
    Rcpp::List fitting_lasso_res = wrap(fitting_lasso_fun(y, X, "1se", alpha));
    return fitting_lasso_res;
  }
}

// [[Rcpp::export]]

double log_factorial(int Y){
  double res = 0;
  for(int kk = 1; kk <= Y; ++kk){
    res += log(kk);
  }
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec log_factorial_calculated(int N){

  arma::vec values = arma::zeros<arma::vec>(N+1);

  for(int kk = 1; kk <= N; ++kk){
    values(kk) = values(kk-1) + log(kk);
  }

  return values;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double alpha_optim(double val){
  Environment myEnv("package:Imputation");
  Function update_alpha = myEnv["update_alpha"];
  Rcpp::NumericVector alpha = wrap(update_alpha(val));
  return alpha(0);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec log_poisson_likelihood_mix(arma::vec Y, double psi, double mu, int n, arma::vec log_factorial_Y){

  arma::vec likelihood = arma::zeros<arma::vec>(n);
  double common_term = -lgamma(psi) + psi*log(psi);

  for(int i = 0; i < n; ++i){
    double e_term = exp(mu);
    double psi_Yi = psi + Y(i);
    likelihood(i) = lgamma(psi_Yi) - log_factorial_Y(i) - psi_Yi*log(psi + e_term) + Y(i)*mu;
  }

  likelihood = likelihood + common_term;
  return(likelihood);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_all_mix(arma::vec Y, double psi, double mu, arma::vec posterior_y, int n1){

  arma::vec gradient = arma::zeros<arma::vec>(2);
  double exp_mu = exp(mu);
  double exp_mu_psi_inv = 1/(exp_mu + psi);
  double frac = exp_mu*exp_mu_psi_inv;

  arma::vec Y_psi = Y + psi;
  arma::vec a1 = Y - Y_psi*frac;
  gradient(0) = sum(a1%posterior_y);

  double com_term = log(psi) + 1 + log(exp_mu_psi_inv) - R::digamma(psi);

  arma::vec b1 = arma::zeros<arma::vec>(n1);
  for(int i = 0; i < n1; ++i){
    b1(i) = R::digamma(Y_psi(i)) - Y_psi(i)*exp_mu_psi_inv + com_term;
  }
  gradient(1) = sum(b1%posterior_y);

  return gradient;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_all_combined_mix(arma::vec Y, double psi, double mu, arma::vec posterior1, int n1){

  arma::vec gradient = arma::zeros<arma::vec>(2);
  double exp_mu = exp(mu);
  double exp_mu_psi_inv = 1/(exp_mu + psi);
  double frac = exp_mu*exp_mu_psi_inv;

  arma::vec Y_psi = Y + psi;
  arma::vec a1 = Y - Y_psi*frac;
  gradient(0) = sum(a1%posterior1);

  double com_term = log(psi) + 1 + log(exp_mu_psi_inv) - R::digamma(psi);

  arma::vec b1 = arma::zeros<arma::vec>(n1);
  for(int i = 0; i < n1; ++i){
    b1(i) = R::digamma(Y_psi(i)) - Y_psi(i)*exp_mu_psi_inv + com_term;
  }

  gradient(1) = sum(b1%posterior1);

  return gradient;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double test_stepsize_for_psi_mix(arma::vec Y, double gra_psi, double ll, double psi,
                                 double mu, int n1, double gamma, double down, arma::vec log_factorial_Y){

  double gra_psi2 = gra_psi*gra_psi*gamma;
  double start = sqrt(abs(psi/gra_psi))/2;

  double aa = start;
  double selected = psi;
  while(aa > 0){
    double aa2 = aa*aa;
    double psi_prime = psi + aa2*gra_psi;
    if(psi_prime > 0){
      double lpsi_prime = sum(log_poisson_likelihood_mix(Y, psi, mu, n1, log_factorial_Y));
      if(lpsi_prime - ll - aa2*gra_psi2 > 0){
        selected = psi_prime;
        break;
      }
    }
    aa = aa - start*down;
  }

  return selected;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double test_stepsize_for_mu_mix(arma::vec Y, double gra_mu, double ll, double psi, double mu,
                                int n1, double gamma, double down, arma::vec log_factorial_Y){

  double gra_mu_2 = gra_mu*gra_mu*gamma;
  double start = sqrt(abs(mu/gra_mu))/2;

  double aa = start;
  double selected = mu;
  while(aa > 0){
    double aa2 = aa*aa;
    double mu1_prime = mu + aa2*gra_mu;
    double lmu1_prime = sum(log_poisson_likelihood_mix(Y, psi, mu, n1, log_factorial_Y));
    if(lmu1_prime - ll - aa2*gra_mu_2 > 0){
      selected = mu1_prime;
      break;
    }

    aa = aa - start*down;
  }

  return selected;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_descent_mix(arma::vec Y, int n1, arma::vec mixing_weights_Y, arma::mat posterior_Y,
                               double psi, double mu, arma::vec log_factorial_Y, int steps, double gamma, double down){

  arma::vec est = arma::zeros<arma::vec>(2);
  arma::vec posterior1 = posterior_Y.col(1);
  arma::vec gradient = gradient_all_combined_mix(Y, psi, mu, posterior1, n1);
  double ll = sum(log_poisson_likelihood_mix(Y, psi, mu, n1, log_factorial_Y));

  double mu_prime = 0; double psi_prime = 0;

  for(int i = 0; i < steps; ++i){

    if(abs(gradient(0)) >= 0.00001){
      mu_prime = test_stepsize_for_mu_mix(Y,  gradient(0), ll, psi, mu, n1,  gamma, down, log_factorial_Y);
    } else {
      mu_prime = mu;
    }
    if(abs(gradient(1)) >= 0.00001){
      psi_prime = test_stepsize_for_psi_mix(Y, gradient(1), ll, psi, mu, n1, gamma, down, log_factorial_Y);
    } else {
      psi_prime = psi;
    }
    mu = mu_prime; psi = psi_prime;
    gradient = gradient_all_combined_mix(Y, psi, mu, posterior1, n1);
    ll = sum(log_poisson_likelihood_mix(Y, psi, mu, n1, log_factorial_Y));

  }
  est(0) = mu;
  est(1) = psi;

  return est;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::rowvec take_exp_weight(arma::rowvec x){

  arma::rowvec lx = x;
  arma::uvec x_ind_small = arma::find(x <= -500);
  int len = x_ind_small.n_elem;
  if(len > 0){
    for(int kk = 0; kk < len; ++kk){
      int id = x_ind_small(kk);
      x(id) = -500;
    }
  }
  if(x(0) == 1){
    x(0) = 0;
    arma::rowvec exp_x = exp(x);
    lx = exp_x/sum(exp_x);
  } else {
    arma::rowvec exp_x = exp(x);
    exp_x(0) = 0;
    lx = exp_x/sum(exp_x);
  }
  return lx;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List EM_discrete_mix(arma::vec Y, int steps, int iter, double gamma, double down, double cutoff){

  int n1 = Y.n_elem;
  arma::mat posterior_Y = arma::zeros<arma::mat>(n1, 2);

  double psi = 10;

  arma::vec calculated_values = log_factorial_calculated(Y.max());
  arma::vec log_factorial_Y = arma::zeros<arma::vec>(n1);

  for(int i = 0; i < n1; ++i){
    log_factorial_Y(i) = calculated_values(Y(i));
  }

  arma::vec mixing_weights_Y = arma::zeros<arma::vec>(2);

  mixing_weights_Y.fill(0.5);

  arma::uvec Y_ind_zero = arma::find(Y <= cutoff);

  int n_0_Y = Y_ind_zero.n_elem;

  double mu = sum(Y)/(n1 - n_0_Y);
  mu = log(mu);

  for(int j = 0; j < n_0_Y; ++j){
    int id = Y_ind_zero(j);
    posterior_Y(id, 0) = 1;
  }

  posterior_Y.col(1) = log_poisson_likelihood_mix(Y, psi, mu, n1, log_factorial_Y);

  for(int i = 0; i < n1; ++i){
    arma::rowvec prob_vec = take_exp_weight(posterior_Y.row(i));
    posterior_Y.row(i) = prob_vec;
  }

  mixing_weights_Y /= sum(mixing_weights_Y);

  arma::vec est = arma::zeros<arma::vec>(2);

  for(int ll = 0; ll < iter; ++ll){

    est = gradient_descent_mix(Y, n1, mixing_weights_Y, posterior_Y, psi, mu, log_factorial_Y, steps, gamma, down);

    double mu = est(0);
    double psi = est(1);
    arma::vec loglikelihood_Y = log_poisson_likelihood_mix(Y, psi, mu, n1, log_factorial_Y);
    posterior_Y.col(1) = exp(loglikelihood_Y)*mixing_weights_Y(1);

    for(int j = 0; j < n_0_Y; ++j){
      int id = Y_ind_zero(j);
      posterior_Y(id, 0) = 1;
    }

    for(int i = 0; i < n1; ++i){
      posterior_Y.row(i) /= sum(posterior_Y.row(i));
    }

    for(int i = 0; i < 2; ++i){
      mixing_weights_Y(i) = sum(posterior_Y.col(i));
    }

    mixing_weights_Y /= sum(mixing_weights_Y);


  }

  return Rcpp::List::create(Rcpp::Named("weights") = mixing_weights_Y,
                            Rcpp::Named("est") = est);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_descent_mix_easy(arma::vec Y, int n1, arma::vec mixing_weights_Y, arma::mat posterior_Y,
                                    double psi, double mu, arma::vec log_factorial_Y, int steps, double down){

  arma::vec est = arma::zeros<arma::vec>(2);
  arma::vec posterior1 = posterior_Y.col(1);
  arma::vec gradient = gradient_all_combined_mix(Y, psi, mu, posterior1, n1);

  double mu_prime = 0; double psi_prime = 0;

  for(int i = 0; i < steps; ++i){


    if(abs(gradient(0)) >= 0.00001){
      double aa = min(0.25, abs(down*gradient(0)));
      mu_prime = mu - gradient(0)*aa/abs(gradient(0));
    } else {
      mu_prime = mu;
    }
    if(abs(gradient(1)) >= 0.00001){
      double aa = min(0.5, abs(down*gradient(1)));
      psi_prime = psi - gradient(1)*aa/abs(gradient(1));
    } else {
      psi_prime = psi;
    }
    mu = mu_prime; psi = psi_prime;
    gradient = gradient_all_combined_mix(Y, psi, mu, posterior1, n1);
  }
  est(0) = mu;
  est(1) = psi;

  return est;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List EM_discrete_mix_easy(arma::vec Y, int steps, int iter, double down, double cutoff){

  int n1 = Y.n_elem;
  arma::mat posterior_Y = arma::zeros<arma::mat>(n1, 2);

  double psi = 10;

  arma::vec calculated_values = log_factorial_calculated(Y.max());
  arma::vec log_factorial_Y = arma::zeros<arma::vec>(n1);

  for(int i = 0; i < n1; ++i){
    log_factorial_Y(i) = calculated_values(Y(i));
  }

  arma::vec mixing_weights_Y = arma::zeros<arma::vec>(2);

  mixing_weights_Y.fill(0.5);

  arma::uvec Y_ind_zero = arma::find(Y <= cutoff);

  int n_0_Y = Y_ind_zero.n_elem;

  double mu = sum(Y)/(n1 - n_0_Y);
  mu = log(mu);

  for(int j = 0; j < n_0_Y; ++j){
    int id = Y_ind_zero(j);
    posterior_Y(id, 0) = 1;
  }

  posterior_Y.col(1) = log_poisson_likelihood_mix(Y, psi, mu, n1, log_factorial_Y);

  for(int i = 0; i < n1; ++i){
    arma::rowvec prob_vec = take_exp_weight(posterior_Y.row(i));
    posterior_Y.row(i) = prob_vec;
  }

  mixing_weights_Y /= sum(mixing_weights_Y);

  arma::vec est = arma::zeros<arma::vec>(2);

  for(int ll = 0; ll < iter; ++ll){

    est = gradient_descent_mix_easy(Y, n1, mixing_weights_Y, posterior_Y, psi, mu, log_factorial_Y, steps, down);

    //  Rcpp::Rcout << est << std::endl;
    double mu = est(0);
    double psi = est(1);
    arma::vec loglikelihood_Y = log_poisson_likelihood_mix(Y, psi, mu, n1, log_factorial_Y);

    posterior_Y.col(1) = exp(loglikelihood_Y)*mixing_weights_Y(1);

    for(int j = 0; j < n_0_Y; ++j){
      int id = Y_ind_zero(j);
      posterior_Y(id, 0) = 1;
    }

    for(int i = 0; i < n1; ++i){
      posterior_Y.row(i) /= sum(posterior_Y.row(i));
    }

    for(int i = 0; i < 2; ++i){
      mixing_weights_Y(i) = sum(posterior_Y.col(i));
    }

    mixing_weights_Y /= sum(mixing_weights_Y);


  }

  return Rcpp::List::create(Rcpp::Named("weights") = mixing_weights_Y,
                            Rcpp::Named("est") = est);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec reweighting_sum_new_posterior_expectation_C(arma::mat Ymat, arma::mat Yflagmat, arma::mat Ycountmat, arma::vec Y, arma::vec Yflag, arma::vec prior_weight){

  int k = Ymat.n_rows;
  double k_cut = 0.9*k;

  arma::vec res = Y;
  arma::uvec to_impute = arma::find(Yflag == 0);
  double to_impute_num = to_impute.n_elem;

  for(int i = 0; i < to_impute_num; ++i){

    //arma::vec Ymat_weights = prior_weight;

    int id = to_impute(i);
    arma::vec Yval = Ymat.col(id);
    arma::uvec Y_ind_zero = arma::find(Yflagmat.col(id) == 0);
    double zero_num = Y_ind_zero.n_elem;

    if(zero_num < k_cut & zero_num != 0){

      Rcpp::List res = EM_discrete_mix_easy(Ycountmat.col(id), 10, 10, 0.02, 0);

      arma::vec para = res["est"];
      arma::vec weights = res["weights"];

      double count_mean = exp(para(0));
      double pN = R::dpois(0, count_mean, FALSE);
      double non_dropout = pN/(weights(0) + pN);

      double kkk = log(count_mean + 0.1);
      double exp_val = non_dropout*kkk;
      for(int j = 0; j < zero_num; ++j){
        int zero_id = Y_ind_zero(j);
        Yval(zero_id) = exp_val;
        //    Ymat_weights(zero_id) = prior_weight(zero_id)*non_dropout;
      }

      //  Ymat_weights = Ymat_weights/sum(Ymat_weights);
    }
    // res(id) = sum(Ymat_weights%Ymat.col(id));
    res(id) = sum(prior_weight%Yval);
  }

  return res;

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec reweighting_sum_new_posterior_expectation_simplified_C(arma::mat Ymat, arma::mat Yflagmat, arma::mat Ycountmat, arma::vec Y, arma::vec Yflag, arma::vec prior_weight){

  int k = Ymat.n_rows;
  double k_cut = 0.9*k;

  arma::vec res = Y;
  arma::uvec to_impute = arma::find(Yflag == 0);
  double to_impute_num = to_impute.n_elem;

  for(int i = 0; i < to_impute_num; ++i){

    //arma::vec Ymat_weights = prior_weight;

    int id = to_impute(i);
    arma::vec Yval = Ymat.col(id);
    arma::uvec Y_ind_zero = arma::find(Yflagmat.col(id) == 0);
    double zero_num = Y_ind_zero.n_elem;

    if(zero_num < k_cut & zero_num != 0){

      double weights = zero_num/k;
      double count_mean = sum(Ycountmat.col(id))/(k - zero_num);
      double pN = R::dpois(0, count_mean, FALSE);
      double non_dropout = pN/(weights + pN);

      double kkk = log(count_mean + 0.1);
      double exp_val = non_dropout*kkk;
      for(int j = 0; j < zero_num; ++j){
        int zero_id = Y_ind_zero(j);
        Yval(zero_id) = exp_val;
        //    Ymat_weights(zero_id) = prior_weight(zero_id)*non_dropout;
      }

      //  Ymat_weights = Ymat_weights/sum(Ymat_weights);
    }
    // res(id) = sum(Ymat_weights%Ymat.col(id));
    res(id) = sum(prior_weight%Yval);
  }

  return res;

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List imputation_by_samples_posterior_expectation(arma::mat data, arma::mat selected_logxx, arma::mat logxx, arma::mat zero_matrix, arma::mat xx,
                                                       int n, int p, bool minbool, double alpha){

  arma::mat imputed = logxx;
  arma::mat t_logxx = logxx.t();
  arma::mat t_zero_matrix = zero_matrix.t();
  arma::mat t_xx = xx.t();
  arma::mat sample_weights = arma::zeros<arma::mat>(n, n);

  arma::uvec ind_list = arma::zeros<arma::uvec>(n-1);
  for(int j = 0; j < n-1; ++j){
    ind_list(j) = j;
  }

  for(int j = 0; j < n; ++j){
    arma::uvec ind_new = ind_list;
    for(int l = j; l < n-1; ++l){
      ind_new(l) += 1;
    }

    Rcpp::List res = fitting_lasso(data.col(j), data.cols(ind_new), minbool, alpha);
    arma::vec coeff = res["coeff"];
    arma::uvec selected = res["selected"] ;
    selected = selected - 1;
    // Rcpp::Rcout << j << std::endl;
    //  Rcpp::Rcout << selected.n_elem << std::endl;
    if(selected.n_elem < 3){

      sample_weights.row(j).fill(-1);

    } else {

      arma::mat selected_submat = selected_logxx.cols(ind_new);
      arma::vec prior_weight = calculate_weights(selected_logxx.col(j), selected_submat.cols(selected));
      arma::uvec nonzero_ind = arma::find(prior_weight >= 0.0001);
      arma::uvec sub_selected = selected(nonzero_ind);
      arma::uvec new_lab = ind_new(sub_selected);
      arma::vec sub_prior_weight = prior_weight(nonzero_ind);
      sub_prior_weight /= sum(sub_prior_weight);

      for(int s = 0; s < new_lab.n_elem; ++s){
        sample_weights(j, new_lab(s)) = sub_prior_weight(s);
      }

      arma::mat Ymat = t_logxx.rows(new_lab);
      arma::mat Yflagmat = t_zero_matrix.rows(new_lab);
      arma::mat Ycountmat = t_xx.rows(new_lab);

      imputed.col(j) = reweighting_sum_new_posterior_expectation_C(Ymat, Yflagmat, Ycountmat, logxx.col(j), zero_matrix.col(j), sub_prior_weight);

    }

  }

  return Rcpp::List::create(Rcpp::Named("imputed") = imputed,
                            Rcpp::Named("sample_weights") = sample_weights);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List imputation_by_samples_posterior_expectation_simplified(arma::mat data, arma::mat selected_logxx, arma::mat logxx, arma::mat zero_matrix, arma::mat xx,
                                                                  int n, int p, bool minbool, double alpha){

  arma::mat imputed = logxx;
  arma::mat t_logxx = logxx.t();
  arma::mat t_zero_matrix = zero_matrix.t();
  arma::mat t_xx = xx.t();
  arma::mat sample_weights = arma::zeros<arma::mat>(n, n);

  arma::uvec ind_list = arma::zeros<arma::uvec>(n-1);
  for(int j = 0; j < n-1; ++j){
    ind_list(j) = j;
  }

  for(int j = 0; j < n; ++j){
    arma::uvec ind_new = ind_list;
    for(int l = j; l < n-1; ++l){
      ind_new(l) += 1;
    }

    Rcpp::List res = fitting_lasso(data.col(j), data.cols(ind_new), minbool, alpha);
    arma::vec coeff = res["coeff"];
    arma::uvec selected = res["selected"] ;
    selected = selected - 1;
    // Rcpp::Rcout << j << std::endl;
    //  Rcpp::Rcout << selected.n_elem << std::endl;
    if(selected.n_elem < 3){

      sample_weights.row(j).fill(-1);

    } else {

      arma::mat selected_submat = selected_logxx.cols(ind_new);
      arma::vec prior_weight = calculate_weights(selected_logxx.col(j), selected_submat.cols(selected));
      arma::uvec nonzero_ind = arma::find(prior_weight >= 0.0001);
      arma::uvec sub_selected = selected(nonzero_ind);
      arma::uvec new_lab = ind_new(sub_selected);
      arma::vec sub_prior_weight = prior_weight(nonzero_ind);
      sub_prior_weight /= sum(sub_prior_weight);

      for(int s = 0; s < new_lab.n_elem; ++s){
        sample_weights(j, new_lab(s)) = sub_prior_weight(s);
      }

      arma::mat Ymat = t_logxx.rows(new_lab);
      arma::mat Yflagmat = t_zero_matrix.rows(new_lab);
      arma::mat Ycountmat = t_xx.rows(new_lab);

      imputed.col(j) = reweighting_sum_new_posterior_expectation_simplified_C(Ymat, Yflagmat, Ycountmat, logxx.col(j), zero_matrix.col(j), sub_prior_weight);

    }

  }

  return Rcpp::List::create(Rcpp::Named("imputed") = imputed,
                            Rcpp::Named("sample_weights") = sample_weights);
}

