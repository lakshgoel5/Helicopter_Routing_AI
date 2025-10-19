#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "structures.h"
#include "utils.h"

#include <chrono>
#include <random>
#include <utility>
#include <algorithm>

using namespace std::chrono;
using namespace std;

//Allotment strategies
OptimalSupplies greedyAllocationByRatio(const ProblemData& problem, double remaining_weight_capacity, int village_global_idx);
OptimalSupplies greedyAllocationByValue(const ProblemData& problem, double remaining_weight_capacity, int village_global_idx);

bool swapVillagesInDifferentHelicopters(Solution& solution, const ProblemData& problem, int& best_cost);
bool swapVillagesinDifferentTripsinSameHelicopter(Solution& solution, const ProblemData& problem, int& best_cost);
bool adjustSuppliesInTrip(Solution& solution, const ProblemData& problem, int& best_cost); //Use Linear Optimization methods
bool removeVillageFromTrip(Solution& solution, const ProblemData& problem);
void addVillageToTrip(Solution& solution, const ProblemData& problem); //Do random Sampling
bool reverseVillagesListInSameTrip(Solution& solution, const ProblemData& problem, int& best_cost);
bool SwapVillagesInSameTrip(Solution& solution, const ProblemData& problem, int& best_cost);
bool mergeTripsInSameHelicopter(Solution& solution, const ProblemData& problem, int& best_cost);
bool addNewTripToHelicopter(Solution& solution, const ProblemData& problem, int& best_cost);
bool moveVillageBetweenTrips(Solution& solution, const ProblemData& problem, int& best_cost);

bool swapVillagesInDifferentHelicoptersSA(Solution& solution, const ProblemData& problem, int& best_cost);
bool swapVillagesinDifferentTripsinSameHelicopterSA(Solution& solution, const ProblemData& problem, int& best_cost);
bool adjustSuppliesInTripSA(Solution& solution, const ProblemData& problem, int& best_cost); //Use Linear Optimization methods
bool removeVillageFromTripSA(Solution& solution, const ProblemData& problem);
void addVillageToTripSA(Solution& solution, const ProblemData& problem); //Do random Sampling
bool reverseVillagesListInSameTripSA(Solution& solution, const ProblemData& problem, int& best_cost);
bool SwapVillagesInSameTripSA(Solution& solution, const ProblemData& problem, int& best_cost);
bool mergeTripsInSameHelicopterSA(Solution& solution, const ProblemData& problem, int& best_cost);
bool addNewTripToHelicopterSA(Solution& solution, const ProblemData& problem, int& best_cost);
bool moveVillageBetweenTripsSA(Solution& solution, const ProblemData& problem, int& best_cost);



#endif // NEIGHBOURS_H