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

arma::vec reweighting_sum_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec Y, arma::vec Yflag, arma::vec prior_weight, bool ImputeAll){

  int p = Ymat.n_cols;
  int k = Ymat.n_rows;
  double k_cut = 0.9*k;

  arma::vec res = Y;
  arma::mat Ymat_prod = Ymat%Yflagmat;

  if(ImputeAll)  {


    for(int i = 0; i < p; ++i){

      arma::vec Ymat_weights = arma::ones<arma::vec>(k);

      arma::uvec Y_ind_zero = arma::find(Yflagmat.col(i) == 0);
      double zero_num = Y_ind_zero.n_elem;
      double nonzero_num = k - zero_num;

      if(zero_num < k_cut & zero_num != 0){

        double mean_Y = sum(Ymat_prod.col(i))/nonzero_num;
        double var_Y = 1;
        if(nonzero_num != 1){
          var_Y = sum(Ymat_prod.col(i)%Ymat_prod.col(i))/nonzero_num - mean_Y*mean_Y;
          var_Y *= nonzero_num/(nonzero_num-1);
        }
        if(var_Y <= 0.01) {
          var_Y = 1;
        }
        double pN = R::dnorm(0, mean_Y, sqrt(var_Y), FALSE);
        double non_dropout = nonzero_num*pN/(zero_num + nonzero_num*pN);
        for(int j = 0; j < zero_num; ++j){
          int zero_id = Y_ind_zero(j);
          Ymat_weights(zero_id) = non_dropout;
        }
        Ymat_weights %= prior_weight;
        Ymat_weights /= sum(Ymat_weights);

      } else {

        Ymat_weights = prior_weight;

      }
      res(i) = sum(Ymat_weights%Ymat.col(i));
    }

  } else {

    arma::uvec to_impute = arma::find(Yflag == 0);
    double to_impute_num = to_impute.n_elem;
    //
    for(int i = 0; i < to_impute_num; ++i){

      arma::vec Ymat_weights = prior_weight;

      int id = to_impute(i);

      arma::uvec Y_ind_zero = arma::find(Yflagmat.col(id) == 0);
      double zero_num = Y_ind_zero.n_elem;
      double nonzero_num = k - zero_num;

      if(zero_num < k_cut & zero_num != 0){

        double mean_Y = sum(Ymat_prod.col(id))/nonzero_num;
        double var_Y = 1;
        if(nonzero_num != 1){
          var_Y = sum(Ymat_prod.col(i)%Ymat_prod.col(i))/nonzero_num - mean_Y*mean_Y;
          var_Y *= nonzero_num/(nonzero_num-1);
        }
        if(var_Y == 0) {
          var_Y = 1;
        }
        double pN = R::dnorm(0, mean_Y, sqrt(var_Y), FALSE);
        double non_dropout = nonzero_num*pN/(zero_num + nonzero_num*pN);

        for(int j = 0; j < zero_num; ++j){
          int zero_id = Y_ind_zero(j);
          Ymat_weights(zero_id) = prior_weight(zero_id)*non_dropout;
        }

        Ymat_weights = Ymat_weights/sum(Ymat_weights);
      }
      res(id) = sum(Ymat_weights%Ymat.col(id));

    }

  }
  return res;

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec reweighting_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec Y, arma::vec Yflag){

  int p = Ymat.n_cols;
  int k = Ymat.n_rows;
  double k_cut = 0.9*k;

  arma::vec res = Y;
  arma::mat Ymat_prod = Ymat%Yflagmat;

  for(int i = 0; i < p; ++i){

    arma::vec Ymat_weights = arma::ones<arma::vec>(k+1);
    arma::uvec Y_ind_zero = arma::find(Yflagmat.col(i) == 0);
    double zero_num = Y_ind_zero.n_elem;
    double nonzero_num = k - zero_num;
    if(Y(i) == 0) nonzero_num = nonzero_num + 1;

    if(zero_num < k_cut & zero_num != 0){

      double mean_Y = (sum(Ymat_prod.col(i)) + Y(i))/nonzero_num;
      double pN = R::dpois(0, mean_Y, FALSE);
      double non_dropout = nonzero_num*pN/(zero_num + nonzero_num*pN);
      for(int j = 0; j < zero_num; ++j){
        int zero_id = Y_ind_zero(j);
        Ymat_weights(zero_id) = non_dropout;
      }
      if(Y(i) == 0) {
        Ymat_weights(k) = non_dropout;
      }
    }
    Ymat_weights /= sum(Ymat_weights);
    arma::vec Values = arma::ones<arma::vec>(k+1);
    Values(arma::span(0, k-1)) = Ymat.col(i);
    Values(k) = Y(i);
    res(i) = sum(Ymat_weights%Values);
  }

  return res*(k+1);

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List imputation_by_samples(arma::mat data, arma::mat selected_logxx, arma::mat logxx, arma::mat zero_matrix, int n, int p, bool minbool, double alpha){

  arma::mat imputed = logxx;
  arma::mat t_logxx = logxx.t();
  arma::mat t_zero_matrix = zero_matrix.t();
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

      imputed.col(j) = reweighting_sum_C(Ymat, Yflagmat, logxx.col(j), zero_matrix.col(j), sub_prior_weight, TRUE);

    }

  }

  return Rcpp::List::create(Rcpp::Named("imputed") = imputed,
                            Rcpp::Named("sample_weights") = sample_weights);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec replacing_by_expectation_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec Y, arma::vec Yflag, arma::mat Ymat_sub, arma::mat Yflagmat_sub, arma::vec prior_weight, bool ImputeAll){

  int p = Ymat.n_cols;
  int k = Ymat.n_rows;
  double k_cut = 0.9*k;

  arma::vec res = Y;
  arma::mat Ymat_prod = Ymat%Yflagmat;

  if(ImputeAll)  {

    for(int i = 0; i < p; ++i){

      arma::uvec Y_ind_zero = arma::find(Yflagmat.col(i) == 0);
      double zero_num = Y_ind_zero.n_elem;
      double nonzero_num = k - zero_num;

      if(zero_num < k_cut & zero_num != 0){

        double mean_Y = sum(Ymat_prod.col(i))/nonzero_num;
        double var_Y = 1;
        if(nonzero_num != 1){
          var_Y = sum(Ymat_prod.col(i)%Ymat_prod.col(i))/nonzero_num - mean_Y*mean_Y;
          var_Y *= nonzero_num/(nonzero_num-1);
        }
        if(var_Y <= 0.01) {
          var_Y = 1;
        }
        double pN = R::dnorm(0, mean_Y, sqrt(var_Y), FALSE);
        double expected = (1-pN)*mean_Y;

        arma::uvec Ysub_ind_zero = arma::find(Yflagmat_sub.col(i) == 0);
        double zero_num_sub = Ysub_ind_zero.n_elem;

        for(int j = 0; j < zero_num_sub; ++j){
          int zero_id = Ysub_ind_zero(j);
          Ymat_sub(zero_id, i) = expected;
        }

      }
      res(i) = sum(prior_weight%Ymat_sub.col(i));
    }

  } else {

    arma::uvec to_impute = arma::find(Yflag == 0);
    double to_impute_num = to_impute.n_elem;
    //
    for(int i = 0; i < to_impute_num; ++i){


      int id = to_impute(i);

      arma::uvec Y_ind_zero = arma::find(Yflagmat.col(id) == 0);
      double zero_num = Y_ind_zero.n_elem;
      double nonzero_num = k - zero_num;

      if(zero_num < k_cut & zero_num != 0){

        double mean_Y = sum(Ymat_prod.col(id))/nonzero_num;
        double var_Y = 1;
        if(nonzero_num != 1){
          var_Y = sum(Ymat_prod.col(i)%Ymat_prod.col(i))/nonzero_num - mean_Y*mean_Y;
          var_Y *= nonzero_num/(nonzero_num-1);
        }
        if(var_Y == 0) {
          var_Y = 1;
        }

        double pN = R::dnorm(0, mean_Y, sqrt(var_Y), FALSE);
        double expected = (1-pN)*mean_Y;

        arma::uvec Ysub_ind_zero = arma::find(Yflagmat_sub.col(i) == 0);
        double zero_num_sub = Ysub_ind_zero.n_elem;

        for(int j = 0; j < zero_num_sub; ++j){
          int zero_id = Ysub_ind_zero(j);
          Ymat_sub(zero_id, i) = expected;
        }

      }
      res(i) = sum(prior_weight%Ymat_sub.col(i));

    }

  }
  return res;

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List imputation_by_samples_expectation(arma::mat data, arma::mat selected_logxx, arma::mat logxx, arma::mat zero_matrix, int n, int p, bool minbool, double alpha){

  arma::mat imputed = logxx;
  arma::mat t_logxx = logxx.t();
  arma::mat t_zero_matrix = zero_matrix.t();
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
      arma::vec sub_prior_weight = prior_weight(nonzero_ind);
      sub_prior_weight /= sum(sub_prior_weight);
      arma::uvec sub_selected = selected(nonzero_ind);
      arma::uvec new_lab = ind_new(sub_selected);
      for(int s = 0; s < new_lab.n_elem; ++s){
        sample_weights(j, new_lab(s)) = sub_prior_weight(s);
      }
      arma::mat Ymat_sub = t_logxx.rows(new_lab);
      arma::mat Yflagmat_sub = t_zero_matrix.rows(new_lab);

      arma::uvec new2_lab = ind_new(selected);
      arma::mat Ymat = t_logxx.rows(new2_lab);
      arma::mat Yflagmat = t_zero_matrix.rows(new2_lab);

      imputed.col(j) = replacing_by_expectation_C(Ymat, Yflagmat, logxx.col(j), zero_matrix.col(j), Ymat_sub, Yflagmat_sub, sub_prior_weight, TRUE);

    }

  }

  return Rcpp::List::create(Rcpp::Named("imputed") = imputed,
                            Rcpp::Named("sample_weights") = sample_weights);
}



