#include "utils.h"

std::vector<std::vector<double>> village_distances;
std::vector<std::vector<double>> home_to_village;

std::vector<VillageState> village_states;
std::vector<HelicopterStats> helicopter_stats;


void precomputeVillageDistances(const ProblemData& problem) {
    int num_villages = problem.villages.size();
    village_distances.assign(num_villages, vector<double>(num_villages, 0.0));
    home_to_village.assign(problem.helicopters.size(), vector<double>(num_villages, 0.0));

    for (int i = 0; i < num_villages; i++) {
        for (int j = i + 1; j < num_villages; j++) {
            double d = distance(problem.villages[i].coords, problem.villages[j].coords);
            village_distances[i][j] = d;
            village_distances[j][i] = d;
        }
    }

    for (int h = 0; h < (int)problem.helicopters.size(); h++) {
        Point home_pos = problem.cities[problem.helicopters[h].home_city_id - 1];
        for (int v = 0; v < num_villages; v++) {
            home_to_village[h][v] = distance(home_pos, problem.villages[v].coords);
        }
    }
}

double cost(const Solution& solution, const ProblemData& problem) {
    double total_travel_cost = 0;
    double total_supply_value = 0;
    for(const auto& helicopter_plan: solution){
        const auto& heli = problem.helicopters[helicopter_plan.helicopter_id-1];
        double total_helicopter_distance = 0;
        int num_trips = helicopter_plan.trips.size();

        for(const auto& trip: helicopter_plan.trips){
            total_helicopter_distance += tripDistance(trip, problem, helicopter_plan.helicopter_id);
            for(const auto& drop: trip.drops){
                total_supply_value += drop.dry_food * problem.packages[0].value;
                total_supply_value += drop.perishable_food * problem.packages[1].value;
                total_supply_value += drop.other_supplies * problem.packages[2].value;
            }
        }
        
        if (total_helicopter_distance > 0) {
            total_travel_cost += (num_trips * heli.fixed_cost) + (heli.alpha * total_helicopter_distance);
        }
    }
    return total_supply_value - total_travel_cost;
}

double tripDistance(const Trip& trip, const ProblemData& problem, int helicopter_id) {
    if (trip.drops.empty()) return 0.0;
    const auto& helicopter = problem.helicopters[helicopter_id - 1];
    Point home_pos = problem.cities[helicopter.home_city_id - 1];

    double total_distance = 0.0;
    const Point* current_pos = &home_pos;

    for (const auto& drop : trip.drops) {
        const auto& village = problem.villages[drop.village_id - 1];
        total_distance += distance(*current_pos, village.coords);
        current_pos = &village.coords;
    }

    total_distance += distance(*current_pos, home_pos);

    return total_distance;
}

// bool isValid(Solution& solution, const ProblemData& problem) {
//     const double TOL = 1e-9;
//     for (const auto& plan : solution) {
//         const auto& heli = problem.helicopters[plan.helicopter_id - 1];
        
//         double total_helicopter_distance = getTotalHelicopterDistance(plan, problem);
//         if (total_helicopter_distance > problem.d_max + TOL) {
//             return false;
//         }

//         for (const auto& trip : plan.trips) {
//             if (getTripWeight(trip, problem) > heli.weight_capacity + TOL) {
//                 return false;
//             }
            
//             if (tripDistance(trip, problem, plan.helicopter_id) > heli.distance_capacity + TOL) {
//                 return false;
//             }
//         }
//     }
//     return true;
// }

bool isValid(Solution& solution, const ProblemData& problem) {
    for (int h = 0; h < (int)solution.size(); h++) {
        const auto& plan = solution[h];
        const auto& heli = problem.helicopters[plan.helicopter_id - 1];
        
        double total_helicopter_distance = 0.0;
        
        for (const auto& trip : plan.trips) {
            double trip_distance = tripDistance(trip, problem, plan.helicopter_id);
            total_helicopter_distance += trip_distance;
            
            double trip_weight = getTripWeight(trip, problem);
            if (trip_weight > heli.weight_capacity) {
                return false;
            }
            
            if (trip_distance > heli.distance_capacity) {
                return false;
            }
            
            // Check village uniqueness
            set<int> served_villages;
            for (const auto& drop : trip.drops) {
                if (served_villages.count(drop.village_id)) {
                    return false; // Village served more than once
                }
                served_villages.insert(drop.village_id);
            }
        }
        
        if (total_helicopter_distance > problem.d_max) {
            return false;
        }
    }
    
    return true;
}

// bool State::isValid() {
//     const double TOL = 1e-9;
//     vector<double> heli_total_dist(prob->helicopters.size(), 0.0);

//     for (size_t i = 0; i < solution.size(); ++i) {
//         const auto& helicopter = prob->helicopters[i];
//         for (const auto& trip : solution[i].trips) {
//             // if weight delivered exceeds wcap then invalid
//             if (calculateTripWeight(trip) > helicopter.weight_capacity + TOL) return false;
//             // if distance of trip exceeds dcap then invalid
//             double trip_dist = calculateTripDistance(trip, i);
//             if (trip_dist > helicopter.distance_capacity + TOL) return false;
//             heli_total_dist[i] += trip_dist;

//             int d_drop = 0, p_drop = 0, o_drop = 0;
//             set<int> visited_villages;
//             for (const auto& drop : trip.drops) {
//                 d_drop += drop.dry_food;
//                 p_drop += drop.perishable_food;
//                 o_drop += drop.other_supplies;
//                 // if village visited more than once in the same trip then invalid
//                 if (visited_villages.count(drop.village_id)) return false;
//                 visited_villages.insert(drop.village_id);
//             }
//             // if more supplies are dropped than picked up then invalid
//             if (d_drop > trip.dry_food_pickup || p_drop > trip.perishable_food_pickup || o_drop > trip.other_supplies_pickup) return false;
//         }
//         // if total distance of helicopter exceeds dmax then invalid
//         if (heli_total_dist[i] > prob->d_max + TOL) return false;
//     }
//     return true;
// }

double getTripWeight(const Trip& trip, const ProblemData& problem) {
    double total_weight = 0;
    for (const auto& drop : trip.drops) {
        total_weight += drop.dry_food * problem.packages[0].weight;
        total_weight += drop.perishable_food * problem.packages[1].weight;
        total_weight += drop.other_supplies * problem.packages[2].weight;
    }
    return total_weight;
}

double getTotalHelicopterDistance(const HelicopterPlan& plan, const ProblemData& problem) {
    double total_distance = 0;
    for (const auto& trip : plan.trips) {
        total_distance += tripDistance(trip, problem, plan.helicopter_id);
    }
    return total_distance;
}

double calculateTripDistance(const Trip& trip, int heli_idx) {
    if (trip.drops.empty()) return 0.0;

    double total_dist = 0.0;
    int last_village_idx = -1; // means "not started yet"

    for (int i = 0; i < (int)trip.drops.size(); i++) {
        int village_idx = trip.drops[i].village_id - 1; // 0-based

        if (i == 0) {
            // From home city → first village
            total_dist += home_to_village[heli_idx][village_idx];
        } else {
            // From previous village → current village
            total_dist += village_distances[last_village_idx][village_idx];
        }

        last_village_idx = village_idx;
    }

    // From last village → home city
    total_dist += home_to_village[heli_idx][last_village_idx];

    return total_dist;
}


// sum of pickups O(1)
double calculateTripWeight(const Trip& trip, const ProblemData& prob) {
    return (trip.dry_food_pickup * prob.packages[0].weight) +
           (trip.perishable_food_pickup * prob.packages[1].weight) +
           (trip.other_supplies_pickup * prob.packages[2].weight);
}

double calculateTripCost(const Trip& trip, int heli_idx, const ProblemData& prob) {
    if (trip.drops.empty()) return 0.0;
    const auto& helicopter = prob.helicopters[heli_idx];
    double distance = calculateTripDistance(trip, heli_idx);
    return helicopter.fixed_cost + (helicopter.alpha * distance);
}