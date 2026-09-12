// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <limits>
#include <string>
#include <vector>

static inline void set_omp_threads(int num_threads)
{
#ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
#else
    (void)num_threads;
#endif
}

static inline double theta_from_disp(double disp)
{
    if (disp > 0.0 && std::isfinite(disp)) {
        return 1.0 / disp;
    }
    return 1.0;
}

static inline double clamp_eta(double eta)
{
    if (eta > 20.0) return 20.0;
    if (eta < -20.0) return -20.0;
    return eta;
}

static inline double calc_loglik_nb(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& beta,
    double theta)
{
    const int n = X.rows();
    double ll = 0.0;
    for (int i = 0; i < n; ++i) {
        double eta = clamp_eta(X.row(i).dot(beta));
        double mu = std::exp(eta);
        double yi = y(i);
        ll += yi * (std::log(mu) - std::log(theta + mu)) +
              theta * (std::log(theta) - std::log(theta + mu));
    }
    return ll;
}

static inline void fill_mu(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& beta,
    Eigen::VectorXd& mu)
{
    const int n = X.rows();
    mu.resize(n);
    Eigen::VectorXd eta = X * beta;
    for (int i = 0; i < n; ++i) {
        mu(i) = std::exp(clamp_eta(eta(i)));
    }
}

static inline bool fit_nb_irls(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    double theta,
    Eigen::VectorXd& beta,
    Eigen::VectorXd& mu_out,
    double& final_ll,
    int max_iter = 30,
    double tol = 1e-3)
{
    const int n = X.rows();
    const int p = X.cols();
    if (n < 1 || p < 1) {
        return false;
    }

    double y_mean = y.mean();
    if (y_mean < 1e-6) {
        y_mean = 1e-6;
    }

    beta = Eigen::VectorXd::Zero(p);
    beta(0) = std::log(y_mean);
    double current_ll = calc_loglik_nb(X, y, beta, theta);
    if (!std::isfinite(current_ll)) {
        return false;
    }

    Eigen::VectorXd eta(n);
    Eigen::VectorXd mu(n);
    Eigen::VectorXd w(n);
    Eigen::VectorXd z(n);

    for (int iter = 0; iter < max_iter; ++iter) {
        eta.noalias() = X * beta;
        for (int i = 0; i < n; ++i) {
            eta(i) = clamp_eta(eta(i));
            mu(i) = std::exp(eta(i));
        }

        for (int i = 0; i < n; ++i) {
            double m = mu(i);
            w(i) = (m * theta) / (theta + m);
            if (w(i) < 1e-10) w(i) = 1e-10;
            z(i) = eta(i) + (y(i) - m) / (m > 1e-10 ? m : 1e-10);
        }

        Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
        Eigen::VectorXd XtWz = X.transpose() * w.cwiseProduct(z);
        Eigen::VectorXd beta_target = XtWX.ldlt().solve(XtWz);
        if (!beta_target.allFinite()) {
            break;
        }

        Eigen::VectorXd delta = beta_target - beta;
        double step = 1.0;
        bool improved = false;
        for (int half = 0; half < 10; ++half) {
            Eigen::VectorXd beta_cand = beta + step * delta;
            double cand_ll = calc_loglik_nb(X, y, beta_cand, theta);
            if (std::isfinite(cand_ll) && cand_ll >= current_ll - 1e-4) {
                beta = beta_cand;
                current_ll = cand_ll;
                improved = true;
                break;
            }
            step *= 0.5;
        }

        if (!improved || (step * delta).cwiseAbs().maxCoeff() < tol) {
            break;
        }
    }

    if (!beta.allFinite() || !std::isfinite(current_ll)) {
        return false;
    }

    fill_mu(X, beta, mu_out);
    final_ll = current_ll;
    return true;
}

static inline void extract_y_sparse(
    const Eigen::SparseMatrix<double>& exprs_t,
    int g,
    int n_cells,
    bool has_size_factors,
    const Rcpp::NumericVector& size_factors,
    Eigen::VectorXd& y,
    double& y_sum)
{
    y = Eigen::VectorXd::Zero(n_cells);
    y_sum = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(exprs_t, g); it; ++it) {
        int cell_idx = it.row();
        double val = it.value();
        if (has_size_factors) {
            val = val / size_factors[cell_idx];
        }
        double rounded_val = std::round(val);
        if (rounded_val < 0.0) rounded_val = 0.0;
        y(cell_idx) = rounded_val;
        y_sum += rounded_val;
    }
}

static inline void extract_y_dense(
    const Eigen::MatrixXd& exprs,
    int g,
    int n_cells,
    bool has_size_factors,
    const Rcpp::NumericVector& size_factors,
    Eigen::VectorXd& y,
    double& y_sum)
{
    y.resize(n_cells);
    y_sum = 0.0;
    for (int i = 0; i < n_cells; ++i) {
        double val = exprs(g, i);
        if (has_size_factors) {
            val = val / size_factors[i];
        }
        double rounded_val = std::round(val);
        if (rounded_val < 0.0) rounded_val = 0.0;
        y(i) = rounded_val;
        y_sum += rounded_val;
    }
}

template <typename FillY>
static Rcpp::List nb_lrt_from_y(
    int n_genes,
    int n_cells,
    const Eigen::MatrixXd& X_full,
    const Eigen::MatrixXd& X_red,
    const Rcpp::NumericVector& disp_guesses,
    const std::string& family_name,
    FillY fill_y)
{
    const int df = X_full.cols() - X_red.cols();
    std::vector<std::string> status(static_cast<size_t>(n_genes), "OK");
    std::vector<std::string> family(static_cast<size_t>(n_genes), family_name);
    std::vector<double> pvals(static_cast<size_t>(n_genes), 1.0);

    #pragma omp parallel for schedule(dynamic, 64)
    for (int g = 0; g < n_genes; ++g) {
        double theta = theta_from_disp(disp_guesses[g]);
        Eigen::VectorXd y;
        double y_sum = 0.0;
        fill_y(g, y, y_sum);

        // A gene with no counts anywhere carries no information. Report it as
        // OK with p = 1 rather than FAIL so that the BH adjustment keeps the
        // same test count as the VGAM path, which fits such genes and reports
        // status OK with p = 1.
        if (y_sum <= 0.0) {
            status[static_cast<size_t>(g)] = "OK";
            pvals[static_cast<size_t>(g)] = 1.0;
            continue;
        }

        Eigen::VectorXd beta_f, beta_r;
        Eigen::VectorXd mu_f(n_cells), mu_r(n_cells);
        double ll_f = 0.0, ll_r = 0.0;
        bool ok_f = fit_nb_irls(X_full, y, theta, beta_f, mu_f, ll_f);
        bool ok_r = fit_nb_irls(X_red, y, theta, beta_r, mu_r, ll_r);
        if (!ok_f || !ok_r) {
            status[static_cast<size_t>(g)] = "FAIL";
            pvals[static_cast<size_t>(g)] = 1.0;
            continue;
        }

        double stat = 2.0 * (ll_f - ll_r);
        if (stat < 0.0 || !std::isfinite(stat)) {
            stat = 0.0;
        }
        double pval = R::pchisq(stat, df, 0, 0);
        if (!std::isfinite(pval)) {
            pvals[static_cast<size_t>(g)] = 1.0;
            status[static_cast<size_t>(g)] = "FAIL";
        } else {
            pvals[static_cast<size_t>(g)] = pval;
            status[static_cast<size_t>(g)] = "OK";
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("status") = Rcpp::wrap(status),
        Rcpp::Named("family") = Rcpp::wrap(family),
        Rcpp::Named("pval") = Rcpp::wrap(pvals)
    );
}

template <typename FillY>
static Rcpp::List nb_fit_predict_from_y(
    int n_genes,
    int n_cells,
    const Eigen::MatrixXd& X_fit,
    const Eigen::MatrixXd& X_pred,
    const Rcpp::NumericVector& disp_guesses,
    bool want_resid,
    FillY fill_y)
{
    const int n_pred = X_pred.rows();
    std::vector<std::string> status(static_cast<size_t>(n_genes), "OK");
    Eigen::MatrixXd mu_pred(n_genes, n_pred);
    mu_pred.setConstant(NA_REAL);
    Eigen::MatrixXd resid(want_resid ? n_genes : 0, want_resid ? n_cells : 0);
    if (want_resid) {
        resid.setConstant(NA_REAL);
    }

    #pragma omp parallel for schedule(dynamic, 64)
    for (int g = 0; g < n_genes; ++g) {
        double theta = theta_from_disp(disp_guesses[g]);
        Eigen::VectorXd y;
        double y_sum = 0.0;
        fill_y(g, y, y_sum);
        // mu -> 0 is the limiting MLE for a gene with no counts; 0/0 is what
        // the VGAM path reports for these genes (up to floating point noise),
        // so return zeros rather than NA to keep the two paths interchangeable.
        if (y_sum <= 0.0) {
            if (want_resid) {
                resid.row(g).setZero();
            }
            mu_pred.row(g).setZero();
            status[static_cast<size_t>(g)] = "OK";
            continue;
        }

        Eigen::VectorXd beta, mu_fit;
        double ll = 0.0;
        if (!fit_nb_irls(X_fit, y, theta, beta, mu_fit, ll)) {
            status[static_cast<size_t>(g)] = "FAIL";
            continue;
        }

        Eigen::VectorXd mu_g;
        fill_mu(X_pred, beta, mu_g);
        mu_pred.row(g) = mu_g.transpose();
        if (want_resid) {
            resid.row(g) = (y - mu_fit).transpose();
        }
        status[static_cast<size_t>(g)] = "OK";
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("status") = Rcpp::wrap(status),
        Rcpp::Named("mu_pred") = mu_pred
    );
    if (want_resid) {
        out["resid"] = resid;
    }
    return out;
}

static void check_lrt_dims(
    int n_cells,
    int n_genes,
    const Eigen::MatrixXd& X_full,
    const Eigen::MatrixXd& X_red,
    const Rcpp::NumericVector& disp_guesses)
{
    if (n_genes < 1 || n_cells < 1) {
        Rcpp::stop("expression matrix is empty");
    }
    if (X_full.rows() != n_cells || X_red.rows() != n_cells) {
        Rcpp::stop("design matrix rows must equal the number of cells");
    }
    if (X_full.cols() < X_red.cols()) {
        Rcpp::stop("full model must have at least as many columns as the reduced model");
    }
    if (disp_guesses.size() != n_genes) {
        Rcpp::stop("disp_guesses length must equal the number of genes");
    }
}

// [[Rcpp::export]]
Rcpp::List fast_diff_test_sparse_cpp(
    const Eigen::SparseMatrix<double>& exprs_t,
    const Eigen::MatrixXd& X_full,
    const Eigen::MatrixXd& X_red,
    const Rcpp::NumericVector& disp_guesses,
    const Rcpp::NumericVector& size_factors,
    bool relative_expr = true,
    int num_threads = 1,
    std::string family_name = "negbinomial.size")
{
    const int n_cells = X_full.rows();
    const int n_genes = exprs_t.cols();
    if (exprs_t.rows() != n_cells) {
        Rcpp::stop("transposed expression matrix rows must equal the number of cells");
    }
    check_lrt_dims(n_cells, n_genes, X_full, X_red, disp_guesses);
    bool has_size_factors = (size_factors.size() == n_cells) && relative_expr;
    set_omp_threads(num_threads);

    auto fill_y = [&](int g, Eigen::VectorXd& y, double& y_sum) {
        extract_y_sparse(exprs_t, g, n_cells, has_size_factors, size_factors, y, y_sum);
    };
    return nb_lrt_from_y(n_genes, n_cells, X_full, X_red, disp_guesses, family_name, fill_y);
}

// [[Rcpp::export]]
Rcpp::List fast_diff_test_dense_cpp(
    const Eigen::MatrixXd& exprs,
    const Eigen::MatrixXd& X_full,
    const Eigen::MatrixXd& X_red,
    const Rcpp::NumericVector& disp_guesses,
    const Rcpp::NumericVector& size_factors,
    bool relative_expr = true,
    int num_threads = 1,
    std::string family_name = "negbinomial.size")
{
    const int n_cells = X_full.rows();
    const int n_genes = exprs.rows();
    if (exprs.cols() != n_cells) {
        Rcpp::stop("dense expression matrix columns must equal the number of cells");
    }
    check_lrt_dims(n_cells, n_genes, X_full, X_red, disp_guesses);
    bool has_size_factors = (size_factors.size() == n_cells) && relative_expr;
    set_omp_threads(num_threads);

    auto fill_y = [&](int g, Eigen::VectorXd& y, double& y_sum) {
        extract_y_dense(exprs, g, n_cells, has_size_factors, size_factors, y, y_sum);
    };
    return nb_lrt_from_y(n_genes, n_cells, X_full, X_red, disp_guesses, family_name, fill_y);
}

// [[Rcpp::export]]
Rcpp::List fast_nb_fit_predict_sparse_cpp(
    const Eigen::SparseMatrix<double>& exprs_t,
    const Eigen::MatrixXd& X_fit,
    const Eigen::MatrixXd& X_pred,
    const Rcpp::NumericVector& disp_guesses,
    const Rcpp::NumericVector& size_factors,
    bool relative_expr = true,
    bool want_resid = false,
    int num_threads = 1)
{
    const int n_cells = X_fit.rows();
    const int n_genes = exprs_t.cols();
    if (exprs_t.rows() != n_cells) {
        Rcpp::stop("transposed expression matrix rows must equal the number of cells");
    }
    if (X_pred.cols() != X_fit.cols()) {
        Rcpp::stop("prediction design must have the same number of columns as the fit design");
    }
    if (disp_guesses.size() != n_genes) {
        Rcpp::stop("disp_guesses length must equal the number of genes");
    }
    bool has_size_factors = (size_factors.size() == n_cells) && relative_expr;
    set_omp_threads(num_threads);

    auto fill_y = [&](int g, Eigen::VectorXd& y, double& y_sum) {
        extract_y_sparse(exprs_t, g, n_cells, has_size_factors, size_factors, y, y_sum);
    };
    return nb_fit_predict_from_y(n_genes, n_cells, X_fit, X_pred, disp_guesses, want_resid, fill_y);
}

// [[Rcpp::export]]
Rcpp::List fast_nb_fit_predict_dense_cpp(
    const Eigen::MatrixXd& exprs,
    const Eigen::MatrixXd& X_fit,
    const Eigen::MatrixXd& X_pred,
    const Rcpp::NumericVector& disp_guesses,
    const Rcpp::NumericVector& size_factors,
    bool relative_expr = true,
    bool want_resid = false,
    int num_threads = 1)
{
    const int n_cells = X_fit.rows();
    const int n_genes = exprs.rows();
    if (exprs.cols() != n_cells) {
        Rcpp::stop("dense expression matrix columns must equal the number of cells");
    }
    if (X_pred.cols() != X_fit.cols()) {
        Rcpp::stop("prediction design must have the same number of columns as the fit design");
    }
    if (disp_guesses.size() != n_genes) {
        Rcpp::stop("disp_guesses length must equal the number of genes");
    }
    bool has_size_factors = (size_factors.size() == n_cells) && relative_expr;
    set_omp_threads(num_threads);

    auto fill_y = [&](int g, Eigen::VectorXd& y, double& y_sum) {
        extract_y_dense(exprs, g, n_cells, has_size_factors, size_factors, y, y_sum);
    };
    return nb_fit_predict_from_y(n_genes, n_cells, X_fit, X_pred, disp_guesses, want_resid, fill_y);
}
