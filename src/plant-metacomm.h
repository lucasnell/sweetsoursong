# ifndef __SWEETSOURSONG_LANDSCAPE_CONSTANTF_H
# define __SWEETSOURSONG_LANDSCAPE_CONSTANTF_H



/*
 Lanscape of plants where the total number of flowers (F) at each plant is
 constant and where pollinators are not affected by flower densities.
 We don't even include F and instead model Y and B as proportions
 of total flowers with non-colonized N = 1 - Y - B.
 */

#include <RcppArmadillo.h>
#include <vector>

#include "ode.h"

using namespace Rcpp;



inline bool lanscape_constF_arg_checks(const std::vector<double>& m,
                                       const std::vector<double>& d_yp,
                                       const std::vector<double>& d_b0,
                                       const std::vector<double>& d_bp,
                                       const std::vector<double>& g_yp,
                                       const std::vector<double>& g_b0,
                                       const std::vector<double>& g_bp,
                                       const std::vector<double>& L_0,
                                       const double& u,
                                       const double& X,
                                       const std::vector<double>& Y0,
                                       const std::vector<double>& B0,
                                       const double& dt,
                                       const double& max_t) {

    size_t np = m.size();
    /*
     I can't just use 'stop()' because it causes a segfault (or similar).
     The system below is my workaround.
     */
    bool err = false;

    len_check(err, d_yp, "d_yp", np);
    len_check(err, d_b0, "d_b0", np);
    len_check(err, d_bp, "d_bp", np);
    len_check(err, g_yp, "g_yp", np);
    len_check(err, g_b0, "g_b0", np);
    len_check(err, g_bp, "g_bp", np);
    len_check(err, L_0, "L_0", np);
    len_check(err, Y0, "Y0", np);
    len_check(err, B0, "B0", np);

    min_val_check(err, d_yp, "d_yp", 0);
    min_val_check(err, d_b0, "d_b0", 0);
    min_val_check(err, d_bp, "d_bp", 0);
    min_val_check(err, g_yp, "g_yp", 0);
    min_val_check(err, g_b0, "g_b0", 0);
    min_val_check(err, g_bp, "g_bp", 0);
    min_val_check(err, L_0, "L_0", 0);
    min_val_check(err, u, "u", 0);
    min_val_check(err, X, "X", 0);

    min_val_check(err, Y0, "Y0", 0);
    min_val_check(err, B0, "B0", 0);
    max_val_check(err, Y0, "Y0", 1);
    max_val_check(err, B0, "B0", 1);

    min_val_check(err, dt, "dt", 0, false);
    min_val_check(err, max_t, "max_t", std::max(dt, 0.0), false);

    return err;

}






class LandscapeConstF
{
public:
    std::vector<double> m;
    std::vector<double> d_yp;
    std::vector<double> d_b0;
    std::vector<double> d_bp;
    std::vector<double> g_yp;
    std::vector<double> g_b0;
    std::vector<double> g_bp;
    std::vector<double> L_0;
    double u;
    double X;
    size_t n_plants;
    bool open_sys;
    // These are only used for MetaObsStochSummRep::fill_output()
    double max_t;
    double season_len;


    LandscapeConstF(const std::vector<double>& m_,
                    const std::vector<double>& d_yp_,
                    const std::vector<double>& d_b0_,
                    const std::vector<double>& d_bp_,
                    const std::vector<double>& g_yp_,
                    const std::vector<double>& g_b0_,
                    const std::vector<double>& g_bp_,
                    const std::vector<double>& L_0_,
                    const double& u_,
                    const double& X_,
                    const bool& open_sys_,
                    const double& max_t_,
                    const double& season_len_)
        : m(m_),
          d_yp(d_yp_),
          d_b0(d_b0_),
          d_bp(d_bp_),
          g_yp(g_yp_),
          g_b0(g_b0_),
          g_bp(g_bp_),
          L_0(L_0_),
          u(u_),
          X(X_),
          n_plants(m_.size()),
          open_sys(open_sys_),
          max_t(max_t_),
          season_len(season_len_),
          weights(m_.size()) {};


    LandscapeConstF(const LandscapeConstF& other)
        : m(other.m),
          d_yp(other.d_yp),
          d_b0(other.d_b0),
          d_bp(other.d_bp),
          g_yp(other.g_yp),
          g_b0(other.g_b0),
          g_bp(other.g_bp),
          L_0(other.L_0),
          u(other.u),
          X(other.X),
          n_plants(other.n_plants),
          open_sys(other.open_sys),
          max_t(other.max_t),
          season_len(other.season_len),
          weights(other.weights) {};

    LandscapeConstF& operator=(const LandscapeConstF& other) {
        m = other.m;
        d_yp = other.d_yp;
        d_b0 = other.d_b0;
        d_bp = other.d_bp;
        g_yp = other.g_yp;
        g_b0 = other.g_b0;
        g_bp = other.g_bp;
        L_0 = other.L_0;
        u = other.u;
        X = other.X;
        n_plants = other.n_plants;
        open_sys = other.open_sys;
        max_t = other.max_t;
        season_len = other.season_len;
        weights = other.weights;
        return *this;
    }


    void operator()(const MatType& x,
                    MatType& dxdt) {

        make_weights(this->weights, x);

        double Ybar = arma::mean(x.col(0));
        double Bbar = arma::mean(x.col(1));

        for (size_t i = 0; i < n_plants; i++) {
            one_plant(i, Ybar, Bbar, x, dxdt);
        }

        return;
    }

    void operator()(const MatType& x,
                  MatType& dxdt,
                  const double& t) {
        this->operator()(x, dxdt);
        return;
    }


    void make_weights(std::vector<double>& wts_vec,
                      const MatType& x) const {

        if (wts_vec.size() != n_plants) wts_vec.resize(n_plants);

        double wt_sum = 0;
        double YN; // proportion palatable nectar
        for (size_t i = 0; i < n_plants; i++) {
            YN = 1 - x(i,1);
            wts_vec[i] = std::pow(YN, u);
            wt_sum += wts_vec[i];
        }

        double mult = static_cast<double>(n_plants) / (2 * (X + wt_sum));
        for (double& w : wts_vec) w *= mult;

        return;
    }



private:

    std::vector<double> weights;

   void one_plant(const size_t& i,
                  const double& Ybar,
                  const double& Bbar,
                  const MatType& x,
                  MatType& dxdt) {

        const double& Y(x(i,0));
        const double& B(x(i,1));
        double N = 1 - Y - B;

        double& dYdt(dxdt(i,0));
        double& dBdt(dxdt(i,1));

        double P = weights[i];

        double Lambda = P / (L_0[i] + P);

        double delta, gamma;

        if (open_sys) {
            delta = (g_yp[i] + d_yp[i] * Y) * Lambda;
            gamma = g_b0[i] + (d_b0[i] * B) + (g_bp[i] + d_bp[i] * B) * Lambda;
        } else {
            delta = (g_yp[i] * Ybar + d_yp[i] * Y) * Lambda;
            gamma = g_b0[i] * Bbar + (d_b0[i] * B) + (g_bp[i] * Bbar + d_bp[i] * B) * Lambda;
        }

        dYdt = delta * N - m[i] * Y;
        dBdt = gamma * N - m[i] * B;

        return;
    }

};




/*
 ==============================================================================
 ==============================================================================
 Observer derived classes
 ==============================================================================
 ==============================================================================
 */

/*
 These are based on the Observer and ObserverBurnEvery classes (in `ode.h`)
 but include code to fill output
 */

// For deterministic simulations:
struct MetaObs : public Observer<MatType> {

    // Fill output:
    void fill_output(MatType& output,
                     const LandscapeConstF& system) const {

        const size_t& np(system.n_plants);
        size_t n_steps = this->data.size();

        std::vector<double> wts(np);
        output.set_size(n_steps * np, 5U);

        size_t i = 0;
        for (size_t t = 0; t < n_steps; t++) {
            system.make_weights(wts, this->data[t]);
            for (size_t k = 0; k < np; k++) {
                output(i,0) = this->time[t];
                output(i,1) = k;
                output(i,2) = this->data[t](k,0);
                output(i,3) = this->data[t](k,1);
                output(i,4) = wts[k];
                i++;
            }

        }

        return;

    }
};


// For stochastic simulations with potential burnin, not saving every time step,
// and reps that need recorded:
struct MetaObsStoch : public ObserverBurnEvery<MatType> {

    MetaObsStoch(const double& burnin_, const double& save_every_,
                 const std::vector<double>& remainders_)
        : ObserverBurnEvery<MatType>(burnin_, save_every_, remainders_) {};

    // Fill output for one repetition:
    void fill_output(MatType& output,
                     const LandscapeConstF& system,
                     const double& dbl_rep) const {

        const size_t& np(system.n_plants);
        size_t n_steps = this->data.size();

        if (n_steps < 1) return;

        std::vector<double> wts(np);
        output.set_size(n_steps * np, 6U);
        // colnames(output) = CharacterVector::create("rep", "t", "p", "Y", "B", "P");
        size_t i = 0;
        for (size_t t = 0; t < n_steps; t++) {
            system.make_weights(wts, this->data[t]);
            for (size_t k = 0; k < np; k++) {
                output(i,0) = dbl_rep;
                output(i,1) = this->time[t];
                output(i,2) = k;
                output(i,3) = this->data[t](k,0);
                output(i,4) = this->data[t](k,1);
                output(i,5) = wts[k];
                i++;
            }
        }

        return;

    }
};

// For stochastic simulations with potential burnin, not saving every time step,
// reps that need recorded, and summarizing by both rep and time point
struct MetaObsStochSumm : public ObserverBurnEvery<MatType> {

    MetaObsStochSumm(const double& burnin_, const double& save_every_,
                     const std::vector<double>& remainders_)
        : ObserverBurnEvery<MatType>(burnin_, save_every_, remainders_) {};

    void operator()(const MatType& x, const double& t) {
        if (t > burnin) {
            time.push_back(t);
            data.push_back(MatType(1, 4));
            data.back()(0,0) = dissimilarity(x);       // BC
            data.back()(0,1) = diversity(x);           // H
            data.back()(0,2) = arma::accu(x.col(0));   // sum Y
            data.back()(0,3) = arma::accu(x.col(1));   // sum B
        }
        return;
    }

    /*
     Fill output for one repetition.
     Note: `system` only used as argument here for template compatibility
     with MetaObsStoch
     */
    void fill_output(MatType& output,
                     const LandscapeConstF& system,
                     const double& dbl_rep) {

        size_t n_steps = this->data.size();

        if (n_steps < 1) return;

        output.set_size(n_steps, 6U);
        // colnames(output) = CharacterVector::create("rep", "t", "BC", "H", "sumY", "sumB");
        for (size_t t = 0; t < n_steps; t++) {
            output(t,0) = dbl_rep;
            output(t,1) = this->time[t];
            output(t,2) = this->data[t](0,0);
            output(t,3) = this->data[t](0,1);
            output(t,4) = this->data[t](0,2);
            output(t,5) = this->data[t](0,3);
        }

        return;

    }


protected:


    // mean community dissimilarity at one time point:
    double dissimilarity(const MatType& x) const {

        size_t n = x.n_rows;

        const arma::subview_col<arma::mat::elem_type> yeast(x.col(0));
        const arma::subview_col<arma::mat::elem_type> bact(x.col(1));

        double n_combos = 0; // only count combos where both patches aren't empty
        double bc_sum = 0;
        double min_y, min_b, denom;
        for (size_t i = 0; i < (n-1U); i++) {
            for (size_t j = i+1U; j < n; j++) {
                denom = yeast(i) + yeast(j) + bact(i) + bact(j);
                if (denom > 0) {
                    min_y = (yeast(i) < yeast(j)) ? yeast(i) : yeast(j);
                    min_b = (bact(i) < bact(j)) ? bact(i) : bact(j);
                    bc_sum += (1 - (2 * (min_y + min_b)) / denom);
                    n_combos += 1.0;
                }
            }
        }
        double bc_mean = bc_sum / n_combos;
        return bc_mean;
    }

    // mean species diversity at one time point:
    double diversity(const MatType& x,
                     double zero_threshold = 2.220446e-16) const {

        size_t n = x.n_rows;

        const arma::subview_col<arma::mat::elem_type> yeast(x.col(0));
        const arma::subview_col<arma::mat::elem_type> bact(x.col(1));

        // Do calculation while accounting for zeros:
        double p_yeast, p_bact, total, H_i, H_mean = 0;
        for (size_t i = 0; i < n; i++) {
            /*
             This if statement implicitly makes H = 0 when one or more
             species is extinct.
             */
            if (bact(i) > zero_threshold && yeast(i) > zero_threshold) {
                total = yeast(i) + bact(i);
                p_yeast = yeast(i) / total;
                p_bact = bact(i) / total;
                H_i = - p_yeast * std::log(p_yeast) - p_bact * std::log(p_bact);
                H_mean += H_i;
            }
        }
        H_mean /= static_cast<double>(n);
        return H_mean;
    }


};



// For stochastic simulations with potential burnin, reps that need recorded,
// and summarizing by just rep.
// Mostly the same as MetaObsStochSumm, except for the `fill_output` method
struct MetaObsStochSummRep : public MetaObsStochSumm {

    MetaObsStochSummRep(const double& burnin_,
                        const double& save_every_,
                        const std::vector<double>& remainders_)
        : MetaObsStochSumm(burnin_, save_every_, remainders_) {};

    /*
     Fill output for one repetition.
     Note: `system` only used as argument here for template compatibility
     with MetaObsStoch
     */
    void fill_output(MatType& output,
                     const LandscapeConstF& system,
                     const double& dbl_rep) const {

        size_t n_steps = this->data.size();

        if (n_steps < 1) return;

        double dbl_n = static_cast<double>(n_steps); // for mean calcs

        // colnames = "rep", "BC", "H",
        //            "minY", "maxY", "meanY", "meanlogY",
        //            "minB", "maxB", "meanB", "meanlogB"
        output.set_size(1U, 11U);

        output(0, 0) = dbl_rep;
        output(0, 1) = 0.0;  // mean BC
        output(0, 2) = 0.0;  // mean H
        // these are summarized for only the last season (if seasons included):
        output(0, 3) = this->data[0](0,2);   // minimum(sum(Y))
        output(0, 4) = this->data[0](0,2);   // maximum(sum(Y))
        output(0, 5) = 0.0;                  // mean(sum(Y))
        output(0, 6) = 0.0;                  // mean(log(sum(Y)))
        output(0, 7) = this->data[0](0,3);   // minimum(sum(B))
        output(0, 8) = this->data[0](0,3);   // maximum(sum(B))
        output(0, 9) = 0.0;                  // mean(sum(B))
        output(0,10) = 0.0;                  // mean(log(sum(B)))

        double dbl_ls = 0; // # of last season's time points
        // time point directly before the last season:
        double last_season_t0 = system.max_t - system.season_len;
        for (size_t i = 0; i < this->data.size(); i++) {
            const MatType& m(this->data[i]);
            // metrics
            output(0,1) += m(0,0);
            output(0,2) += m(0,1);
            // yeast and bacteria
            if (this->time[i] > last_season_t0) {
                add_YB__(output, m);
                dbl_ls += 1.0;
            }
        }

        output(0, 1) /= dbl_n;
        output(0, 2) /= dbl_n;
        output(0, 5) /= dbl_ls;
        output(0, 6) /= dbl_ls;
        output(0, 9) /= dbl_ls;
        output(0,10) /= dbl_ls;

        return;

    }

    // Add summary info for yeast and bacteria
    inline void add_YB__(MatType& output, const MatType& m) const {
        // yeast
        const double& Yi(m(0,2));
        if (Yi < output(0,3)) output(0,3) = Yi;
        if (Yi > output(0,4)) output(0,4) = Yi;
        output(0, 5) += Yi;
        output(0, 6) += std::log(Yi + 1.0);
        // bacteria
        const double& Bi(m(0,3));
        if (Bi < output(0,7)) output(0,7) = Bi;
        if (Bi > output(0,8)) output(0,8) = Bi;
        output(0, 9) += Bi;
        output(0,10) += std::log(Bi + 1.0);
        return;
    }

};




#endif
