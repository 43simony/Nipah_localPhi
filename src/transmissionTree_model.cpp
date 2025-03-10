//
//  transmissionNetwork_model.cpp
//  
//
//  Created by Brandon Simony on 2/12/25.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <chrono>
#include <queue>
#include <random>
#include <sstream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "gsl/gsl_randist.h"

/*///////////////////////
//// GSL probability ////
////  Distributions  ////
///////////////////////*/

//For initiating and seeding gsl random number generator.
class Gsl_rng_wrapper
{
    gsl_rng* r;
    public:
        Gsl_rng_wrapper()
        {
            std::random_device rng_dev; //To generate a safe seed.
            long seed = time(NULL)*rng_dev();
            const gsl_rng_type* rng_type = gsl_rng_default;
            r = gsl_rng_alloc(rng_type);
            gsl_rng_set(r, seed);
        }
        ~Gsl_rng_wrapper() { gsl_rng_free(r); }
        gsl_rng* get_r() { return r; }
};


//Uniform RV using gsl.
double draw_uniform_gsl(double lo, double hi)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_flat(r, lo, hi);
}


//Binomial RV using gsl.
int draw_binom_gsl(int N, double prob)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_binomial(r, prob, N);
}


// Poisson RV using gsl. parameterized by mean
int draw_poisson_gsl(double lambda)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_poisson(r, lambda);
}


//Beta RV using gsl.
double draw_beta_gsl(double a, double b)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_beta(r, a, b);
}


//Exponential RV using gsl. parameterized by a mean
double draw_exp_gsl(double mu)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_exponential(r, mu);
}


//Gamma RV using gsl. k, theta parameterization
double draw_gamma_gsl(double shape, double scale)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_gamma(r, shape, scale);
}


//Beta RV using gsl.
double draw_beta_mix(std::vector<double> a, std::vector<double> b, std::vector<double> weights)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    
    
    int idx = 0;
    double val = 0;
    double U1 = draw_uniform_gsl(0.0, 1.0);

    while(val == 0){
        if(U1 < weights[idx]){ val = draw_beta_gsl(a[idx], b[idx]); }
        idx++;
    }

    return val;
}


int draw_negBinom_gsl(double p, double n)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_negative_binomial(r, p, n);
}


// draw from poisson mixture
int draw_pois_mix(double R0, double R0_SS, double prob_SS)
{
    double U1 = draw_uniform_gsl(0.0, 1.0);
    int N;
    if(U1 < prob_SS){
        N = draw_poisson_gsl(R0_SS); // draw from super spreader offspring distribution
    }else{
        N = draw_poisson_gsl(R0); // draw from average offspring distribution
    }
    return(N);
}


// draw from poisson-gamma mixture
int draw_poisGamma_mix(double k, double theta)
{
    double v = draw_gamma_gsl(k, theta);
    int N = draw_poisson_gsl(v);
    return(N);
}


/*///////////////////////
//// Model Parameter ////
////    Structure    ////
///////////////////////*/

struct Parameters
{

    // parameters for offspring distributions
    int model_type = 0; // toggle for offspring distribution type. 0: NB, 1: pois_mix, 2: poisGamma_mix -- 0
    double R0 = 0.0; // reproductive number for NB model or non-super spreader R0 for pois_mix model -- 1
    double k = 0.0; // shape parameter for poisGamma_mix model or dispersion in NB model -- 2
    double R0_SS = 0.0; // super spreader reproductive number for pois_mix model -- 3
    double prob_SS = 0.0; // proportion of super spreaders for pois_mix model -- 4
    
    
    // network model type
    int data_type = 0; // switch indicator 0: simulation, 1: reconstructed from data -- 5
    int reps = 1; // number of replicate contact networks to simulate [for each spillover in dataset] -- 6
    int N_threshold = 10000; // total number of individuals allowed in a contact network -- 7
    std::string data_path = ""; // file name and path from working directory to data for reconstrected network. only used if data_type == 1 -- 8
    
    std::vector<std::string> conf_v;

    //Constructor either reads parameters from standard input (if no argument is given),
    //or from file (argument 1(.
    Parameters(int argc, char* argv[])
    {
        conf_v.reserve(200);
        std::stringstream buffer;
        if(argc == 3) //Config file
        {
            std::ifstream f(argv[2]);
            if(f.is_open())
            {
                buffer << f.rdbuf();
            }
            else
            {
                std::cout << "Failed to read config file \"" << argv[2] << "\"" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            for(int i=2; i<argc; ++i)
            {
                buffer << argv[i] << std::endl;
            }
        }
        while(!buffer.eof())
        { // until end of the stream
            std::string line = "";
            std::stringstream line_ss;
            // First get a whole line from the config file
            std::getline(buffer, line);
            // Put that in a stringstream (required for getline) and use getline again
            // to only read up until the comment character (delimiter).
            line_ss << line;
            std::getline(line_ss, line, '#');
            // If there is whitespace between the value of interest and the comment character
            // this will be read as part of the value. Therefore strip it of whitespace first.
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line.size() != 0)
            {
                conf_v.push_back(line);
            }
        }

        size_t conf_length = 9; //number of model input parameters
        if (conf_v.size()!=conf_length)
        {
            std::cout << "Expected configuration file with " << conf_length << " options, loaded file with "
                      << conf_v.size() << " lines." << std::endl;
            exit(EXIT_FAILURE);
        }

        // parameters for offspring distributions
        model_type = std::stoi(conf_v[0]); // toggle for offspring distribution type. 0: NB, 1: pois_mix, 2: poisGamma_mix -- 0
        R0 = std::stod(conf_v[1]); // reproductive number for NB model or non-super spreader R0 for pois_mix model -- 1
        k = std::stod(conf_v[2]); // shape parameter for poisGamma_mix model or dispersion in NB model -- 2
        R0_SS = std::stod(conf_v[3]); // super spreader reproductive number for pois_mix model -- 3
        prob_SS = std::stod(conf_v[4]); // proportion of super spreaders for pois_mix model -- 4
        
        
        // network model type
        data_type = std::stoi(conf_v[5]); // switch indicator 0: simulation, 1: reconstructed from data -- 5
        reps = std::stoi(conf_v[6]); // number of replicate contact networks to simulate [for each spillover in dataset] -- 6
        N_threshold = std::stoi(conf_v[7]); // total number of individuals allowed in a contact network -- 7
        data_path = conf_v[8]; // path to data for reconstructed network data -- 8
    }
};



struct Edge {
    int from;
    int to;
};

// Function for simulating a transmission tree
std::vector<Edge> generate_transmission_tree( const Parameters& p) {
    
    // int initial_cases = 1, int max_generations = 5, int max_infections = 3
    std::vector<Edge> edges;
    std::queue<std::pair<int, int>> queue; // (case_num, generation)
    int case_counter = 1;
    edges.push_back({1, 1});
    
    queue.push({1, 0}); // Initial case
    
    
    while (!queue.empty()) {
        // auto [current_case, generation] = queue.front();
        std::pair<int, int> current = queue.front();
        int current_case = current.first;
        int generation = current.second;
        queue.pop();
        
        // empty remaining nodes from the queue. Censor offspring distributions, as these individuals are otherwise counted as excess draws of zero
        if (case_counter >= p.N_threshold) {
            edges.push_back({current_case, -1}); // assign out node value of -1 as marker for removal
            continue;
        }
        
        int new_infections;
        
        // draw from specified offspring distribution
        if (p.model_type == 0) {
            double prob = 1.0 / (1.0 + (p.R0/p.k));
            
            new_infections = draw_negBinom_gsl(prob, p.k); // negative binomial with mean R0 defined in standard form instead of NB(mean, dispersion)
            // std:cout << prob << "; " << new_infections << std::endl;
            
        } else if (p.model_type == 1) {
            new_infections = draw_pois_mix(p.R0, p.R0_SS, p.prob_SS);
        } else {
            new_infections = draw_poisGamma_mix(p.k, p.R0/p.k); // stochastic draw from poisson-gamma mixture. Equivalent to NB model under same starting parameters. Gamma uses shape-scale parameterization.
        }
        // end draw from offspring distribution
        
        for (int i = 0; i < new_infections; ++i) {
            case_counter += 1;
            edges.push_back({current_case, case_counter});
            queue.push({case_counter, generation + 1});
        }
        
    } // end while loop
    
    return edges;
}



std::vector<Edge> transmission_tree_reconstruct(const Parameters& p) {
    
    std::vector<Edge> edges;
    std::queue<std::pair<int, int>> queue; // (case_num, generation)
    std::vector<int> non_spreaders; // Track individuals who initially produced 0 cases
    int case_counter = 1;
    
    edges.push_back({1, 1}); // Root node
    queue.push({1, 0}); // Initial case

    int maxLoops = 1e9; // Prevent infinite loops
    int iter = 0;

    while (case_counter < p.N_threshold) {
        
        iter++;
        if (iter > maxLoops) {
            edges.push_back({-1, -1}); // Indicator of incomplete network
            return edges;
        }

        if (queue.empty()) {
            // If the queue is empty but we need more edges, use non-spreaders
            if (!non_spreaders.empty()) {
                
                int random_index = floor(draw_uniform_gsl(0, non_spreaders.size() - 1)); // draw random index for revived case
                
                int revived_case = non_spreaders[random_index];
                non_spreaders.erase(non_spreaders.begin() + random_index); // Remove selected case from the list
                queue.push({revived_case, 0}); // Reintroduce for spreading
            } else {
                break; // No more candidates to spread, terminate loop
            }
        }

        std::pair<int, int> current = queue.front();
        int current_case = current.first;
        int generation = current.second;
        queue.pop(); // Remove current case from queue

        int new_infections;

        // Determine new infections based on the model type
        if (p.model_type == 0) {
            double prob = 1.0 / (1.0 + (p.R0/p.k));
            new_infections = draw_negBinom_gsl(p.k, prob); // negative binomial with mean R0 defined in standard form instead of NB(mean, dispersion)
        } else if (p.model_type == 1) {
            new_infections = draw_pois_mix(p.R0, p.R0_SS, p.prob_SS);
        } else {
            new_infections = draw_poisGamma_mix(p.k, p.R0/p.k); // stochastic draw from poisson-gamma mixture. Equivalent to NB model under same starting parameters. Gamma uses shape-scale parameterization.
        }

        // Adjust new infections to ensure exactly N_threshold edges
        if (case_counter + new_infections > p.N_threshold) {
            new_infections = p.N_threshold - case_counter; // Trim excess infections
        }

        if (new_infections == 0) {
            non_spreaders.push_back(current_case); // Track for potential reinfection
        }

        for (int i = 0; i < new_infections; ++i) {
            case_counter += 1;
            edges.push_back({current_case, case_counter});
            queue.push({case_counter, generation + 1});

            // Stop adding if we've reached the exact number of edges
            if (case_counter >= p.N_threshold) break;
        }
    }

    return edges;
}



int main(int argc, char* argv[]){
    
    /*
     If you run this from the command line, you should provide two arguments - the
     number of replicates and a config file (argc==3). In that case the results will
     be output to a file. The alternative is to provide all the parameters in a single
     long string instead of the config file. This is done by the R-wrapper. In that
     case the results will be written to a .txt file and saved to a table in R.
     */
    
   
    if(argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <n replicates> <config file / string of parameter values>" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // initialize parameters read from command line input
    int rep = std::stoi(argv[1]);
    Parameters p(argc, argv);
    
    
    std::vector<Edge> tree;
    if(p.data_type == 1){
        tree = transmission_tree_reconstruct(p);
    }else{
        tree = generate_transmission_tree(p);
    }
    
    
    std::cout << "from;" << "to" << std::endl;
    for (const auto& edge : tree) {
        std::cout << edge.from << ";" << edge.to << std::endl;
    }
    
    return 0;
}
