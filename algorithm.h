#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "structures.h"
#include "utils.h"
#include "neighbours.h"

#include <chrono>

void precomputeVillageDistances(const ProblemData& problem);

//Random State Generators
vector<int> clusterVillages(const ProblemData& problem, int time_limit_seconds);
vector<Trip> planTripsForHelicopters(const ProblemData& problem, int helicopter_idx, const vector<int>& village_assignment, int time_limit_seconds);
Solution random_state(const ProblemData& problem, int time_limit_seconds);

//Neighbourhood functions
pair<int, Solution> neighbouringFunction(Solution& current_solution, int current_cost, const ProblemData& problem, int time_limit_seconds);
pair<int, Solution> neighbouringFunctionSA(Solution current_solution, int current_cost, const ProblemData& problem, int time_limit_seconds);

//Hill Climbing Algorithms
Solution hillClimbing(Solution& current_solution, const ProblemData& problem, int time_limit_seconds, int num_of_steps);
Solution hillClimbingWithRestarts(const ProblemData& problem, double time_limit_seconds);
Solution hillClimbingLocalSearch(Solution current_solution, const ProblemData& problem, double time_limit_seconds);

Solution simulated_annealing_with_restarts(const ProblemData& problem, double time_limit_seconds);
Solution simulated_annealing(Solution current_solution, const ProblemData& problem, double time_limit_seconds);

#endif // ALGORITHM_H