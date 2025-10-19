#include "solver.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <tuple>
#include <random>
#include <unordered_map>
#include <set>
#include <climits>
#include <algorithm>
#include "algorithm.h"
#include "greedy.h"
// #include "Simulated_Annealing.h"
#include "utils.h"

using namespace std;
using namespace std::chrono;


Solution solve(const ProblemData& problem) {
    auto start_time = steady_clock::now();
    double total_time_limit_seconds = problem.time_limit_minutes * 59.0;

    precomputeVillageDistances(problem);

    Solution greedy_solution = greedy(problem, 0.40 * (double)total_time_limit_seconds);
    int greedy_cost = cost(greedy_solution, problem);
    if(!isValid(greedy_solution, problem)) {
        greedy_cost = INT_MIN;
    }
    cout << "Greedy solution cost: " << greedy_cost << endl;
    double elapsed_seconds = chrono::duration_cast<chrono::duration<double>>(steady_clock::now() - start_time).count();
    total_time_limit_seconds -= elapsed_seconds;    // --- Time Allocation (10% HC, 90% SA) ---
    double hill_climbing_time_limit = (total_time_limit_seconds) * 0.10;

    auto before_greedy = steady_clock::now();
    // --- Phase 1: Hill Climbing ---
    cout << "Starting Phase 1: Hill Climbing (Time limit: " << hill_climbing_time_limit << "s)" << endl;
    Solution hc_solution = hillClimbingWithRestarts(problem, hill_climbing_time_limit);
    int hc_cost = cost(hc_solution, problem);
    cout << "Hill Climbing finished. Initial best cost: " << hc_cost << endl;

    // --- Phase 2: Simulated Annealing with Restarts ---
    auto after_greedy = steady_clock::now();
    double elapsed_hc_time = chrono::duration_cast<chrono::duration<double>>(after_greedy - before_greedy).count();
    double remaining_time = (total_time_limit_seconds - elapsed_hc_time) * 0.95;
    Solution final_solution = hc_solution; // Start SA with the best known solution from HC

    if (remaining_time > 2.0) {
        cout << "Starting Phase 2: Simulated Annealing with Restarts (Time limit: " << remaining_time << "s)" << endl;
        
        // Run SA with restarts, but keep track of the best solution found so far (including the HC result)
        Solution sa_best_solution = simulated_annealing_with_restarts(problem, remaining_time);
        
        // Compare final HC result with final SA result
        int sa_best_cost = cost(sa_best_solution, problem);
        if(sa_best_cost > hc_cost && sa_best_cost > greedy_cost) {
            final_solution = sa_best_solution;
            cout << "SA improved the solution. Final best cost: " << sa_best_cost << endl;
        } else if(hc_cost > greedy_cost) {
            final_solution = hc_solution;
            cout << "HC provided the best solution. Final best cost: " << hc_cost << endl;
        } else {
            final_solution = greedy_solution;
            cout << "Greedy provided the best solution. Final best cost: " << greedy_cost << endl;
        }
    } else {
        cout << "Skipping SA phase, insufficient time remaining." << endl;
    }

    return final_solution;
}