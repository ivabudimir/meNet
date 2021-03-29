// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;


double absC(double x){
  if(x<0){
    return(-x);
  }else{
    return(x);
  }
}

NumericVector sliceColumnC(NumericMatrix m, NumericVector rows, int col){
  int n = rows.size();
  NumericVector slice(n);
  for(int i=0; i<n; ++i){
    slice[i]=m(rows[i],col);
  }
  return(slice);
}

double corC(NumericVector x, NumericVector y) {
  int n = x.size();
  double x_sum = 0, y_sum = 0, xy_sum = 0, x2_sum = 0, y2_sum=0;
  for(int i=0; i<n; ++i){
    x_sum += x[i];
    y_sum += y[i];
    xy_sum += x[i]*y[i];
    x2_sum += pow(x[i],2);
    y2_sum += pow(y[i],2);
  }
  double cor = (n*xy_sum-x_sum*y_sum)/(sqrt(n*x2_sum-pow(x_sum,2))*sqrt(n*y2_sum-pow(y_sum,2)));
  return (cor);
  }

// [[Rcpp::export]]
NumericMatrix cor_meanSd(NumericMatrix data, int size, bool replace, int n_repetitions, double zero_precision=0.000001){
  int nrow = data.nrow(), ncol = data.ncol();
  NumericMatrix output_matrix(ncol*(ncol-1)/2,2); // initialized as zero matrix
  NumericVector indices_all(nrow);
  NumericVector indices_temp(size);
  int current_row;
  NumericVector dataVec_i(size), dataVec_j(size);
  double cor_ij;
  for(int i=0; i<nrow; ++i){
    indices_all[i] = i;
  }
  for(int rep=0; rep<n_repetitions; ++rep){
    indices_temp = RcppArmadillo::sample(indices_all, size, replace);
    current_row = 0;
    for(int j=1; j<ncol; ++j) // order is by columns of upper triangular matrix: (1,2) (1,3) (2,3) (1,4) (2,4) (3,4) (1,5)...
      for(int i=0; i<j; ++i){
        dataVec_i = sliceColumnC(data, indices_temp, i);
        dataVec_j = sliceColumnC(data, indices_temp, j);
        cor_ij = corC(dataVec_i, dataVec_j);
        output_matrix(current_row,0) += cor_ij;
        output_matrix(current_row,1) += pow(cor_ij,2);
        ++current_row;
      }
  }
  for(int i=0; i<output_matrix.nrow(); ++i){
    output_matrix(i,1) = output_matrix(i,1)/(n_repetitions-1)-pow(output_matrix(i,0),2)/(n_repetitions*(n_repetitions-1));
    if(absC(output_matrix(i,1))<zero_precision){
      output_matrix(i,1) = 0;
    }else{
      output_matrix(i,1) = sqrt(output_matrix(i,1));
    }
    output_matrix(i,0) = output_matrix(i,0)/n_repetitions;
  }
  return(output_matrix);
}

// [[Rcpp::export]]
NumericVector cor_prTest(NumericMatrix data, int n_repetitions, String alternative, double zero_precision=0.000001){
  int nrow = data.nrow(), ncol = data.ncol();
  NumericVector pval_vector(ncol*(ncol-1)/2); // initialized as zero vector
  NumericVector indices_vector(nrow), indices_permuted(nrow);
  for(int i=0; i<nrow; ++i){
    indices_vector[i] = i;
  }
  //NumericMatrix resampling_matrix(nrow, n_repetitions);
  //for(int rep=0; rep<n_repetitions; ++rep){
  //  indices_permuted = RcppArmadillo::sample(indices_vector, nrow, false);
  //  for(int sample=0; sample<nrow; ++sample){
  //    resampling_matrix[sample,rep] = indices_permuted[sample];
  //  }
  //}
  int current_row=0;
  NumericVector dataVec_i(nrow), dataVec_j(nrow);
  double cor_true, cor_permutation;
  for(int j=1; j<ncol; ++j) // order is by columns of upper triangular matrix: (1,2) (1,3) (2,3) (1,4) (2,4) (3,4) (1,5)...
    for(int i=0; i<j; ++i){
      dataVec_i = sliceColumnC(data, indices_vector, i);
      dataVec_j = sliceColumnC(data, indices_vector, j);
      cor_true = corC(dataVec_i, dataVec_j);
      for(int rep=0; rep<n_repetitions; ++rep){
        indices_permuted = RcppArmadillo::sample(indices_vector, nrow, false);
        dataVec_j = sliceColumnC(data, indices_permuted, j);
        cor_permutation = corC(dataVec_i, dataVec_j);
        if(alternative=="less"){
          pval_vector[current_row] += int(cor_permutation<cor_true);
        }else if(alternative=="greater"){
          pval_vector[current_row] += int(cor_permutation>cor_true);
        }else if(alternative=="two_sided"){
          pval_vector[current_row] += int(absC(cor_permutation)>absC(cor_true));
        }else{ //two_sided_signed
          if(cor_true>0){
            pval_vector[current_row] += int(cor_permutation>cor_true);
          }else if(cor_true<0){
            pval_vector[current_row] += int(cor_permutation<cor_true);
          }else{
            pval_vector[current_row] += int(cor_permutation!=cor_true);
          }
        }
      }
      pval_vector[current_row] /= double(n_repetitions);
      ++current_row;
    }
  return(pval_vector);
}

