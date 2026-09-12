// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <limits>

// [[Rcpp::export]]
Rcpp::IntegerVector find_closest_point_cpp(
    const Eigen::MatrixXd& Z,
    const Eigen::MatrixXd& Y,
    int num_threads = 1)
{
    if (Z.rows() != Y.rows()) {
        Rcpp::stop("Z and Y must have the same number of rows (dimensions)");
    }
    const int D = Z.rows();
    const int N = Z.cols();
    const int K = Y.cols();
    if (D < 1 || N < 1 || K < 1) {
        Rcpp::stop("Z and Y must be non-empty");
    }

    std::vector<int> closest(static_cast<size_t>(N), 1);

#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
#endif

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; ++i) {
        double min_dist_sq = std::numeric_limits<double>::infinity();
        int min_k = 0;
        for (int k = 0; k < K; ++k) {
            double dist_sq = 0.0;
            for (int d = 0; d < D; ++d) {
                double diff = Z(d, i) - Y(d, k);
                dist_sq += diff * diff;
            }
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                min_k = k;
            }
        }
        closest[static_cast<size_t>(i)] = min_k + 1;
    }

    return Rcpp::wrap(closest);
}
