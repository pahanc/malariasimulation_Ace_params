/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "adult_mosquito_ode.h"
#include <sstream>
#include <fstream>
#include <iostream>
using namespace std;


AdultMosquitoModel::AdultMosquitoModel(
    MosquitoModel growth_model,
    double mu,
    double tau,
    std::vector<double> mosq_suppression,
    std::vector<double> mosq_seasonality,
    std::vector<double> emergence,
    std::vector<double> total_M_orig,
    bool use_Ace_mosq,
    bool dens_indep,
    double incubating
    ) : growth_model(growth_model), mu(mu), tau(tau)
{
    for (auto i = 0u; i < tau; ++i) {
        lagged_incubating.push(incubating);
    }
}

integration_function_t create_ode(AdultMosquitoModel& model) {
    //create original ode
    auto growth_ode = create_ode(model.growth_model);
    return [&model, growth_ode](const state_t& x, state_t& dxdt, double t) {
        //set submodel total_M
        model.growth_model.total_M =
            x[get_idx(AdultODEState::S)] +
            x[get_idx(AdultODEState::E)] +
            x[get_idx(AdultODEState::I)];

        //run the growth ode
        growth_ode(x, dxdt, t);

        //run the adult ode
        auto incubation_survival = exp(-model.mu * model.tau);

	int t_day;
	t_day = t;
	int day_of_year;
	day_of_year=t-365*(t_day/365);
	//if (t>365) Rcpp::Rcout << "day of year " << day_of_year << " t " << t << " t_day/365 "<<t_day/365<<  endl;

	std::ofstream emerg_out;

        if (model.use_Ace_mosq){

		if (model.dens_indep){//this not working right now

			/*dxdt[get_idx(AdultODEState::SN)] =.5 * model.mosq_suppression[t_day] * model.mosq_seasonality[t_day] * x[get_idx(ODEState::P)] * model.growth_model.total_M/model.total_M_orig[day_of_year] / model.growth_model.dp; //growth to adult female
			dxdt[get_idx(AdultODEState::S)] =
            .5 * model.mosq_suppression[t_day] * model.mosq_seasonality[t_day] * x[get_idx(ODEState::P)] * model.growth_model.total_M/model.total_M_orig[day_of_year] / model.growth_model.dp //growth to adult female
            - x[get_idx(AdultODEState::S)] * model.foim //infections
            - x[get_idx(AdultODEState::S)] * model.mu; //deaths   


		if (t>0 && t<2) Rcpp::Rcout << " total_M " << model.growth_model.total_M << " ratio " << model.growth_model.total_M/model.total_M_orig[day_of_year] << endl; */

		}

		if (!model.dens_indep){
			dxdt[get_idx(AdultODEState::SN)] =.5 * model.mosq_suppression[t_day] * model.mosq_seasonality[t_day] * x[get_idx(ODEState::P)] / model.growth_model.dp; //growth to adult female
			dxdt[get_idx(AdultODEState::S)] =
            .5 * model.mosq_suppression[t_day] * model.mosq_seasonality[t_day] * x[get_idx(ODEState::P)]/ model.growth_model.dp //growth to adult female
            - x[get_idx(AdultODEState::S)] * model.foim //infections
            - x[get_idx(AdultODEState::S)] * model.mu; //deaths  
		}
	    //if (t==2) Rcpp::Rcout << t << " " << dxdt[get_idx(AdultODEState::SN)] << endl;
	    
	    //if (t>1 && int(10*t)%10==0){
		    //model.emergence.push_back(dxdt[get_idx(AdultODEState::SN)]);
		    //Rcpp::Rcout << "t "<< t << " emergence " << model.emergence[t] << endl;
	    	    //Rcpp::Rcout << "t " <<t << " mosq seasonality " << model.mosq_seasonality[t_day] << endl;

	    //}


	}

	if (!model.use_Ace_mosq){
	dxdt[get_idx(AdultODEState::SN)] =
            .5 * x[get_idx(ODEState::P)] / model.growth_model.dp; //growth to adult female
	dxdt[get_idx(AdultODEState::S)] =
            .5 * x[get_idx(ODEState::P)] / model.growth_model.dp //growth to adult female
            - x[get_idx(AdultODEState::S)] * model.foim //infections
            - x[get_idx(AdultODEState::S)] * model.mu; //deaths   
	}


        dxdt[get_idx(AdultODEState::E)] =
            x[get_idx(AdultODEState::S)] * model.foim  //infections
            - model.lagged_incubating.front() * incubation_survival //survived incubation period
            - x[get_idx(AdultODEState::E)] * model.mu; // deaths

        dxdt[get_idx(AdultODEState::I)] = model.lagged_incubating.front() * incubation_survival //survived incubation period
            - x[get_idx(AdultODEState::I)] * model.mu; // deaths


	//if (t==365 && model.use_Ace_mosq) {
		//Rcpp::Rcout << " emergence size " << model.emergence.size() << endl;
		//emerg_out.open("emergence.txt",ios::out | ios::app | ios::binary);
		//seems like this resets emergence for each new solver? if you could use it like
		//lagged_incubating that might work. But its not a solved quantity...
		//for (std::vector<double>::iterator i=model.emergence.begin(); i<model.emergence.end(); i++){
		//	emerg_out << *i << endl;
		//}
	//}

    };
}

//[[Rcpp::export]]
Rcpp::XPtr<AdultMosquitoModel> create_adult_mosquito_model(
    Rcpp::XPtr<MosquitoModel> growth_model,
    double mu,
    double tau,
    std::vector<double> mosq_suppression,
    std::vector<double> mosq_seasonality,
    std::vector<double> emergence,
    std::vector<double> total_M_orig,
    bool use_Ace_mosq,
    bool dens_indep,
    double susceptible
    ) {
    auto model = new AdultMosquitoModel(*growth_model, mu, tau, mosq_suppression, mosq_seasonality, emergence, total_M_orig, use_Ace_mosq, dens_indep, susceptible);
    return Rcpp::XPtr<AdultMosquitoModel>(model, true);
}

//[[Rcpp::export]]
void adult_mosquito_model_update(
    Rcpp::XPtr<AdultMosquitoModel> model,
    double mu,
    double foim,
    std::vector<double> mosq_suppression,
    std::vector<double> mosq_seasonality,
    std::vector<double> emergence,
    std::vector<double> total_M_orig,
    bool use_Ace_mosq,
    bool dens_indep,
    double susceptible,
    double f
    ) {
    model->mu = mu;
    model->foim = foim;
    model->growth_model.f = f;
    model->growth_model.mum = mu;
    model->mosq_suppression = mosq_suppression;
    model->mosq_seasonality = mosq_seasonality;
    model->emergence = emergence;
    model->total_M_orig = total_M_orig;
    model->use_Ace_mosq = use_Ace_mosq;
    model->dens_indep = dens_indep;
    model->lagged_incubating.push(susceptible * foim);
    if (model->lagged_incubating.size() > 0) {
        model->lagged_incubating.pop();
    }
}

//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_adult_solver(
    Rcpp::XPtr<AdultMosquitoModel> model,
    std::vector<double> init,
    double r_tol,
    double a_tol
    ) {
    return Rcpp::XPtr<Solver>(
        new Solver(init, create_ode(*model), r_tol, a_tol),
        true
    );
}
