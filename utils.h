#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include "structures.h"
#include <set>

using namespace std;

struct OptimalSupplies {
    int dry_food;
    int perishable_food;
    int other_supplies;
    double total_weight;
    double total_value;
    bool feasible;
};

struct VillageState {
    double food_delivered = 0.0;      // total food (dry + perishable) delivered so far
    double supplies_delivered = 0.0;  // total other supplies delivered so far
    double max_food_needed = 0.0;     // 9 * population
    double max_supplies_needed = 0.0; // 1 * population
};

struct HelicopterStats {
    double total_distance = 0.0; // sum of distances of all its trips
};

extern std::vector<VillageState> village_states;
extern std::vector<HelicopterStats> helicopter_stats;

extern std::vector<vector<double>> village_distances; // Precomputed distances between villages
extern std::vector<vector<double>> home_to_village; // Precomputed distances from home to villages

double tripDistance(const Trip& trip, const ProblemData& problem, int helicopter_id);
double cost(const Solution& solution, const ProblemData& problem);
bool isValid(Solution& solution, const ProblemData& problem);
double getTotalHelicopterDistance(const HelicopterPlan& plan, const ProblemData& problem);
double getTripWeight(const Trip& trip, const ProblemData& problem);

double calculateTripDistance(const Trip& trip, int heli_idx);
double calculateTripWeight(const Trip& trip, const ProblemData& prob);
double calculateTripCost(const Trip& trip, int heli_idx, const ProblemData& prob);




#endif // UTILS_H