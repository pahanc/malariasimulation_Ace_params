/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "mosquito_ode.h"
#include "mosquito_ode.h"
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;


integration_function_t create_ode(MosquitoModel& model) {
    return [&model](const state_t& x, state_t& dxdt, double t) {
        auto K = carrying_capacity(
            t,
            model.model_seasonality,
            model.g0,
            model.g,
            model.h,
            model.K0,
            model.R_bar
        );
        auto beta = eggs_laid(model.beta, model.mum, model.f);
        auto n_larvae = x[get_idx(ODEState::E)] + x[get_idx(ODEState::L)];

	double K_max;
	double K_day[365];
	for (int i=0;i<365;i++){
		K_day[i] = carrying_capacity(
            i,
            model.model_seasonality,
            model.g0,
            model.g,
            model.h,
            model.K0,
            model.R_bar
        );
	}
	K_max=*std::max_element(K_day,K_day+365);//different for each species
	//if (int(100*t) % 500==0) Rcpp::Rcout <<"K_max " << K_max <<" t " << t << endl;
	//if (int(100*t) % 500==0) { 
	//	for (int i=0;i<365;i++) Rcpp::Rcout << K_day[i] << " ";
	//Rcpp::Rcout << endl;
	//}
	
	K=K+0.000065*K_max;


        if (model.use_Ace_mosq){
	dxdt[get_idx(ODEState::E)] = beta * (model.total_M) //new eggs
            - x[get_idx(ODEState::E)] / model.de //growth to late larval stage
	    - x[get_idx(ODEState::E)] * model.mue; //early larval deaths

	dxdt[get_idx(ODEState::L)] = x[get_idx(ODEState::E)] / model.de //growth from early larval
            - x[get_idx(ODEState::L)] / model.dl //growth to pupal
            - x[get_idx(ODEState::L)] * model.mul; //late larval deaths
	}
	if (!model.use_Ace_mosq){

	dxdt[get_idx(ODEState::E)] = beta * (model.total_M) //new eggs
             - x[get_idx(ODEState::E)] / model.de //growth to late larval stage
            - x[get_idx(ODEState::E)] * model.mue * (1 + model.gamma * n_larvae / K); //early larval deaths
	
	dxdt[get_idx(ODEState::L)] = x[get_idx(ODEState::E)] / model.de //growth from early larval
            - x[get_idx(ODEState::L)] / model.dl //growth to pupal
            - x[get_idx(ODEState::L)] * model.mul * (1 + model.gamma * n_larvae / K); //late larval deaths
        
        dxdt[get_idx(ODEState::P)] = x[get_idx(ODEState::L)] / model.dl //growth to pupae
            - x[get_idx(ODEState::P)] / model.dp //growth to adult
            - x[get_idx(ODEState::P)] * model.mup; // death of pupae
	}
    };
}

MosquitoModel::MosquitoModel(
    double beta,
    double de,
    double mue,
    double K0,
    double gamma,
    double dl,
    double mul,
    double dp,
    double mup,
    bool use_Ace_mosq,
    double total_M_orig,
    size_t total_M,
    bool model_seasonality,
    double g0,
    std::vector<double> g,
    std::vector<double> h,
    double R_bar
    ):
    beta(beta),
    de(de),
    mue(mue),
    K0(K0),
    gamma(gamma),
    dl(dl),
    mul(mul),
    dp(dp),
    mup(mup),
    use_Ace_mosq(use_Ace_mosq),
    total_M_orig(total_M_orig),
    total_M(total_M),
    model_seasonality(model_seasonality),
    g0(g0),
    g(g),
    h(h),
    R_bar(R_bar)
    {}



//[[Rcpp::export]]
Rcpp::XPtr<MosquitoModel> create_mosquito_model(
    double beta,
    double de,
    double mue,
    double K0,
    double gamma,
    double dl,
    double mul,
    double dp,
    double mup,
    bool use_Ace_mosq,
    double total_M_orig,
    size_t total_M,
    bool model_seasonality,
    double g0,
    std::vector<double> g,
    std::vector<double> h,
    double R_bar
    ) {
    auto model = new MosquitoModel(
        beta,
        de,
        mue,
        K0,
        gamma,
        dl,
        mul,
        dp,
        mup,
	use_Ace_mosq,
	total_M_orig,
        total_M,
        model_seasonality,
        g0,
        g,
        h,
        R_bar
    );
    return Rcpp::XPtr<MosquitoModel>(model, true);
}

//[[Rcpp::export]]
void mosquito_model_update(
    Rcpp::XPtr<MosquitoModel> model,
    size_t total_M,
    double f,
    double mum
    ) {
    model->total_M = total_M;
    model->f = f;
    model->mum = mum;
}


//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_solver(
    Rcpp::XPtr<MosquitoModel> model,
    std::vector<double> init,
    double r_tol,
    double a_tol
    ) {
    return Rcpp::XPtr<Solver>(
        new Solver(init, create_ode(*model), r_tol, a_tol),
        true
    );
}
