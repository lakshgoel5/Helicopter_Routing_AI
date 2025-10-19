#include "neighbours.h"

#include <iostream>

using namespace std;
using namespace std::chrono;


//change values of demand in food_demand and other_supplies_demand, so that updated values can be used later
OptimalSupplies greedyAllocationByRatio(const ProblemData& problem, double remaining_weight_capacity, int village_global_idx) {
    OptimalSupplies supplies = {0, 0, 0, 0.0, 0.0, false};

    // Calculate remaining demand for the specific village from global state.
    // Use max(0.0, ...) to avoid negative demand if precision errors occur elsewhere.
    double food_needed = max(0.0, village_states[village_global_idx].max_food_needed - village_states[village_global_idx].food_delivered);
    double other_needed = max(0.0, village_states[village_global_idx].max_supplies_needed - village_states[village_global_idx].supplies_delivered);

    // Cast to int for allocation logic, ensuring we attempt to fulfill small fractional demands.
    int current_food_demand = static_cast<int>(ceil(food_needed));
    int current_other_demand = static_cast<int>(ceil(other_needed));

    if (current_food_demand <= 0 && current_other_demand <= 0) return supplies;

    // Get package data
    double w_dry = problem.packages[0].weight;
    double v_dry = problem.packages[0].value;
    double w_perishable = problem.packages[1].weight;
    double v_perishable = problem.packages[1].value;
    double w_other = problem.packages[2].weight;
    double v_other = problem.packages[2].value;

    // Calculate value-to-weight ratios
    double ratio_dry = (w_dry > 0) ? v_dry / w_dry : 0;
    double ratio_perishable = (w_perishable > 0) ? v_perishable / w_perishable : 0;
    double ratio_other = (w_other > 0) ? v_other / w_other : 0;

    int max_perishable_demand = current_food_demand;
    int max_other_demand = current_other_demand;

    vector<pair<double, int>> items; // {ratio, type} where type: 0=dry, 1=perishable, 2=other
    items.push_back({ratio_perishable, 1}); // Prioritize perishable first for tie-breaking if ratios are equal
    items.push_back({ratio_dry, 0});
    items.push_back({ratio_other, 2});

    sort(items.rbegin(), items.rend()); // Sort descending by ratio

    double remaining_capacity = remaining_weight_capacity;
    int dry_allocated = 0, perishable_allocated = 0, other_allocated = 0;

    for (auto& item : items) {
        int type = item.second;
        if (remaining_capacity <= 0) break;

        if (type == 0) { // Dry food
            int demand_to_fill = max_perishable_demand - perishable_allocated; // Fill remaining food demand
            int max_possible = (w_dry > 0) ? static_cast<int>(remaining_capacity / w_dry) : 0;
            dry_allocated = min(max_possible, demand_to_fill);
            if (dry_allocated < 0) dry_allocated = 0;
            remaining_capacity -= (dry_allocated * w_dry);
        } else if (type == 1) { // Perishable food
            int demand_to_fill = max_perishable_demand;
            int max_possible = (w_perishable > 0) ? static_cast<int>(remaining_capacity / w_perishable) : 0;
            perishable_allocated = min(max_possible, demand_to_fill);
            if (perishable_allocated < 0) perishable_allocated = 0;
            remaining_capacity -= (perishable_allocated * w_perishable);
        } else if (type == 2) { // Other supplies
            int demand_to_fill = max_other_demand;
            int max_possible = (w_other > 0) ? static_cast<int>(remaining_capacity / w_other) : 0;
            other_allocated = min(max_possible, demand_to_fill);
            if (other_allocated < 0) other_allocated = 0;
            remaining_capacity -= (other_allocated * w_other);
        }
    }
    
    // Final check to ensure total food doesn't exceed demand due to independent calculation
    int food_allocated = dry_allocated + perishable_allocated;
    if (food_allocated > current_food_demand) {
        int excess = food_allocated - current_food_demand;
        if (ratio_dry < ratio_perishable) { // Remove less valuable items first
            int reduce_dry = min(dry_allocated, excess);
            dry_allocated -= reduce_dry;
            excess -= reduce_dry;
            perishable_allocated = max(0, perishable_allocated - excess);
        } else {
            int reduce_perishable = min(perishable_allocated, excess);
            perishable_allocated -= reduce_perishable;
            excess -= reduce_perishable;
            dry_allocated = max(0, dry_allocated - excess);
        }
    }

    supplies.dry_food = dry_allocated;
    supplies.perishable_food = perishable_allocated;
    supplies.other_supplies = other_allocated;
    supplies.total_weight = dry_allocated * w_dry + perishable_allocated * w_perishable + other_allocated * w_other;
    supplies.total_value = dry_allocated * v_dry + perishable_allocated * v_perishable + other_allocated * v_other;
    supplies.feasible = (supplies.total_weight > 0 && supplies.total_weight <= remaining_weight_capacity + 1e-9); // Add tolerance

    if (supplies.feasible) {
        village_states[village_global_idx].food_delivered += supplies.dry_food + supplies.perishable_food;
        village_states[village_global_idx].supplies_delivered += supplies.other_supplies;
    }

    return supplies;
}

OptimalSupplies greedyAllocationByValue(const ProblemData& problem, double remaining_weight_capacity, int village_global_idx) {
    OptimalSupplies supplies = {0, 0, 0, 0.0, 0.0, false};

    double food_needed = max(0.0, village_states[village_global_idx].max_food_needed - village_states[village_global_idx].food_delivered);
    double other_needed = max(0.0, village_states[village_global_idx].max_supplies_needed - village_states[village_global_idx].supplies_delivered);
    int current_food_demand = static_cast<int>(ceil(food_needed));
    int current_other_demand = static_cast<int>(ceil(other_needed));

    if (current_food_demand <= 0 && current_other_demand <= 0) return supplies;

    double w_dry = problem.packages[0].weight;
    double v_dry = problem.packages[0].value;
    double w_perishable = problem.packages[1].weight;
    double v_perishable = problem.packages[1].value;
    double w_other = problem.packages[2].weight;
    double v_other = problem.packages[2].value;

    vector<tuple<double, int, double>> items; // {value, type, weight} type: 0=dry, 1=perishable, 2=other
    items.emplace_back(v_perishable, 1, w_perishable);
    items.emplace_back(v_dry, 0, w_dry);
    items.emplace_back(v_other, 2, w_other);

    sort(items.rbegin(), items.rend());

    double remaining_capacity = remaining_weight_capacity;
    int dry_allocated = 0, perishable_allocated = 0, other_allocated = 0;

    for (const auto& item : items) {
        int type = get<1>(item);
        double weight = get<2>(item);
        if (remaining_capacity <= 0 || weight <= 0) continue;

        if (type == 0) { // Dry food
            int demand_to_fill = current_food_demand - perishable_allocated; // Fill remaining food demand
            int max_possible = static_cast<int>(remaining_capacity / weight);
            dry_allocated = min(max_possible, demand_to_fill);
            if (dry_allocated < 0) dry_allocated = 0;
            remaining_capacity -= (dry_allocated * weight);
        } else if (type == 1) { // Perishable food
            int demand_to_fill = current_food_demand;
            int max_possible = static_cast<int>(remaining_capacity / weight);
            perishable_allocated = min(max_possible, demand_to_fill);
            if (perishable_allocated < 0) perishable_allocated = 0;
            remaining_capacity -= (perishable_allocated * weight);
        } else if (type == 2) { // Other supplies
            int demand_to_fill = current_other_demand;
            int max_possible = static_cast<int>(remaining_capacity / weight);
            other_allocated = min(max_possible, demand_to_fill);
            if (other_allocated < 0) other_allocated = 0;
            remaining_capacity -= (other_allocated * weight);
        }
    }

    int food_allocated = dry_allocated + perishable_allocated;
    if (food_allocated > current_food_demand) {
        int excess = food_allocated - current_food_demand;
        if (v_dry < v_perishable) { // Remove less valuable items first based on value sort order
            int reduce_dry = min(dry_allocated, excess);
            dry_allocated -= reduce_dry;
            excess -= reduce_dry;
            perishable_allocated = max(0, perishable_allocated - excess);
        } else {
            int reduce_perishable = min(perishable_allocated, excess);
            perishable_allocated -= reduce_perishable;
            excess -= reduce_perishable;
            dry_allocated = max(0, dry_allocated - excess);
        }
    }

    supplies.dry_food = dry_allocated;
    supplies.perishable_food = perishable_allocated;
    supplies.other_supplies = other_allocated;
    supplies.total_weight = dry_allocated * w_dry + perishable_allocated * w_perishable + other_allocated * w_other;
    supplies.total_value = dry_allocated * v_dry + perishable_allocated * v_perishable + other_allocated * v_other;
    supplies.feasible = (supplies.total_weight > 0 && supplies.total_weight <= remaining_weight_capacity + 1e-9);

    if (supplies.feasible) {
        village_states[village_global_idx].food_delivered += supplies.dry_food + supplies.perishable_food;
        village_states[village_global_idx].supplies_delivered += supplies.other_supplies;
    }

    return supplies;
}



bool swapVillagesInDifferentHelicopters(Solution& solution, const ProblemData& problem, int& best_cost) {
    static random_device rd;
    static mt19937 gen(rd());
    
    uniform_int_distribution<> h_dist(0, solution.size() - 1);
    int h1 = h_dist(gen);
    int h2;
        do {
            h2 = h_dist(gen);
        } while (h1 == h2);    
    if(solution[h1].trips.size() > 0 && solution[h2].trips.size() > 0) {
        int t1 = uniform_int_distribution<>(0, solution[h1].trips.size() - 1)(gen);
        int t2 = uniform_int_distribution<>(0, solution[h2].trips.size() - 1)(gen);
        
        if(solution[h1].trips[t1].drops.size() > 0 && solution[h2].trips[t2].drops.size() > 0) {
            int d1 = uniform_int_distribution<>(0, solution[h1].trips[t1].drops.size() - 1)(gen);
            int d2 = uniform_int_distribution<>(0, solution[h2].trips[t2].drops.size() - 1)(gen);

            // Create neighbor by swapping villages
            Solution neighbor_solution = solution;
            
            // Perform the swap in the neighbor
            swap(neighbor_solution[h1].trips[t1].drops[d1], neighbor_solution[h2].trips[t2].drops[d2]);
            
            // Check if neighbor is valid and calculate its cost
            if (isValid(neighbor_solution, problem)) {
                //Since only swap, therefore no change in village state
                //There is a change in trip pickup amounts, so we need to rebalance them
                neighbor_solution[h1].trips[t1].dry_food_pickup = neighbor_solution[h1].trips[t1].dry_food_pickup - solution[h1].trips[t1].drops[d1].dry_food + neighbor_solution[h1].trips[t1].drops[d1].dry_food;
                neighbor_solution[h1].trips[t1].perishable_food_pickup = neighbor_solution[h1].trips[t1].perishable_food_pickup - solution[h1].trips[t1].drops[d1].perishable_food + neighbor_solution[h1].trips[t1].drops[d1].perishable_food;
                neighbor_solution[h1].trips[t1].other_supplies_pickup = neighbor_solution[h1].trips[t1].other_supplies_pickup - solution[h1].trips[t1].drops[d1].other_supplies + neighbor_solution[h1].trips[t1].drops[d1].other_supplies;
                
                neighbor_solution[h2].trips[t2].dry_food_pickup = neighbor_solution[h2].trips[t2].dry_food_pickup - solution[h2].trips[t2].drops[d2].dry_food + neighbor_solution[h2].trips[t2].drops[d2].dry_food;
                neighbor_solution[h2].trips[t2].perishable_food_pickup = neighbor_solution[h2].trips[t2].perishable_food_pickup - solution[h2].trips[t2].drops[d2].perishable_food + neighbor_solution[h2].trips[t2].drops[d2].perishable_food;
                neighbor_solution[h2].trips[t2].other_supplies_pickup = neighbor_solution[h2].trips[t2].other_supplies_pickup - solution[h2].trips[t2].drops[d2].other_supplies + neighbor_solution[h2].trips[t2].drops[d2].other_supplies;
                //also change helicopter stats i.e. total distance
                int trip_distance_h1_old = tripDistance(solution[h1].trips[t1], problem, solution[h1].helicopter_id);
                int trip_distance_h1_new = tripDistance(neighbor_solution[h1].trips[t1], problem, neighbor_solution[h1].helicopter_id);
                helicopter_stats[h1].total_distance = helicopter_stats[h1].total_distance - trip_distance_h1_old + trip_distance_h1_new;
                int trip_distance_h2_old = tripDistance(solution[h2].trips[t2], problem, solution[h2].helicopter_id);
                int trip_distance_h2_new = tripDistance(neighbor_solution[h2].trips[t2], problem, neighbor_solution[h2].helicopter_id);
                helicopter_stats[h2].total_distance = helicopter_stats[h2].total_distance - trip_distance_h2_old + trip_distance_h2_new;

                int neighbor_new_cost = cost(neighbor_solution, problem);
                if (neighbor_new_cost > best_cost) {
                    best_cost = neighbor_new_cost;
                    solution = neighbor_solution;
                    return true;
                }
            }
        }
    }
    return false;
}

bool SwapVillagesInSameTrip(Solution& solution, const ProblemData& problem, int& best_cost) {
    static random_device rd;
    static mt19937 gen(rd());
    
    if (solution.empty()) return false;
    int h = uniform_int_distribution<>(0, solution.size() - 1)(gen);
    
    if(solution[h].trips.size() > 0) {
        int t = uniform_int_distribution<>(0, solution[h].trips.size() - 1)(gen);
        
        if(solution[h].trips[t].drops.size() > 1) {
            int d1 = uniform_int_distribution<>(0, solution[h].trips[t].drops.size() - 1)(gen);
            int d2;
            do {
                d2 = uniform_int_distribution<>(0, solution[h].trips[t].drops.size() - 1)(gen);
            } while (d1 == d2);

            // Create neighbor by swapping villages within the same trip
            Solution neighbor_solution = solution;
            auto& neighbor_drops = neighbor_solution[h].trips[t].drops;
            
            // Perform the swap in the neighbor
            swap(neighbor_drops[d1], neighbor_drops[d2]);
            
            // Check if neighbor is valid and calculate its cost
            if (isValid(neighbor_solution, problem)) {
                //No need of updating village states
                //Need to update helicopter stats
                int trip_distance_old = tripDistance(solution[h].trips[t], problem, solution[h].helicopter_id);
                int trip_distance_new = tripDistance(neighbor_solution[h].trips[t], problem, neighbor_solution[h].helicopter_id);
                helicopter_stats[h].total_distance = helicopter_stats[h].total_distance - trip_distance_old + trip_distance_new;

                int neighbor_new_cost = cost(neighbor_solution, problem);
                if (neighbor_new_cost > best_cost) {
                    best_cost = neighbor_new_cost;
                    solution = neighbor_solution;
                    return true;
                }
            }
        }
    }
    return false;
}

bool reverseVillagesListInSameTrip(Solution& solution, const ProblemData& problem, int& best_cost) {
    static random_device rd;
    static mt19937 gen(rd());
    
    if (solution.empty()) return false;
    int h = uniform_int_distribution<>(0, solution.size() - 1)(gen);
    
    if(solution[h].trips.size() > 0) {
        int t = uniform_int_distribution<>(0, solution[h].trips.size() - 1)(gen);
        
        if(solution[h].trips[t].drops.size() > 3) { // Need at least 4 villages for meaningful 2-opt
            int i = uniform_int_distribution<>(0, solution[h].trips[t].drops.size() - 3)(gen);
            int j = uniform_int_distribution<>(i + 2, solution[h].trips[t].drops.size() - 1)(gen);
            
            // Create neighbor by applying 2-opt (reverse the segment between i+1 and j)
            Solution neighbor_solution = solution;
            
            // Reverse the segment from i+1 to j (inclusive)
            reverse(neighbor_solution[h].trips[t].drops.begin() + i + 1, neighbor_solution[h].trips[t].drops.begin() + j + 1);
            
            // Check if neighbor is valid and calculate its cost
            if (isValid(neighbor_solution, problem)) {
                //No need of updating village states
                //Need to update helicopter stats
                int trip_distance_old = tripDistance(solution[h].trips[t], problem, solution[h].helicopter_id);
                int trip_distance_new = tripDistance(neighbor_solution[h].trips[t], problem, neighbor_solution[h].helicopter_id);
                helicopter_stats[h].total_distance = helicopter_stats[h].total_distance - trip_distance_old + trip_distance_new;
                int neighbor_new_cost = cost(neighbor_solution, problem);
                if (neighbor_new_cost > best_cost) {
                    best_cost = neighbor_new_cost;
                    solution = neighbor_solution;
                    return true;
                }
            }
        }
        
    }
    return false;
}

bool swapVillagesinDifferentTripsinSameHelicopter(Solution& solution, const ProblemData& problem, int& best_cost) {
    static random_device rd;
    static mt19937 gen(rd());

    // Also try swap villages within same helicopter but different trips (Type 3 from original)
    int h = uniform_int_distribution<>(0, solution.size() - 1)(gen);
    if (solution[h].trips.size() > 1) {
        int t1 = uniform_int_distribution<>(0, solution[h].trips.size() - 1)(gen);
        int t2;
        do {
            t2 = uniform_int_distribution<>(0, solution[h].trips.size() - 1)(gen);
        } while(t1 == t2);

        if (solution[h].trips[t1].drops.size() > 0 && solution[h].trips[t2].drops.size() > 0) {
            int d1 = uniform_int_distribution<>(0, solution[h].trips[t1].drops.size() - 1)(gen);
            int d2 = uniform_int_distribution<>(0, solution[h].trips[t2].drops.size() - 1)(gen);

            // Create neighbor by swapping villages between trips of same helicopter
            Solution neighbor_solution = solution;
            
            // Perform the swap in the neighbor
            swap(neighbor_solution[h].trips[t1].drops[d1], neighbor_solution[h].trips[t2].drops[d2]);
            
            // Check if neighbor is valid and calculate its cost
            if (isValid(neighbor_solution, problem)) {
                //Rebalance trip pickups
                neighbor_solution[h].trips[t1].dry_food_pickup = neighbor_solution[h].trips[t1].dry_food_pickup - solution[h].trips[t1].drops[d1].dry_food + neighbor_solution[h].trips[t1].drops[d1].dry_food;
                neighbor_solution[h].trips[t1].perishable_food_pickup = neighbor_solution[h].trips[t1].perishable_food_pickup - solution[h].trips[t1].drops[d1].perishable_food + neighbor_solution[h].trips[t1].drops[d1].perishable_food;
                neighbor_solution[h].trips[t1].other_supplies_pickup = neighbor_solution[h].trips[t1].other_supplies_pickup - solution[h].trips[t1].drops[d1].other_supplies + neighbor_solution[h].trips[t1].drops[d1].other_supplies;

                neighbor_solution[h].trips[t2].dry_food_pickup = neighbor_solution[h].trips[t2].dry_food_pickup - solution[h].trips[t2].drops[d2].dry_food + neighbor_solution[h].trips[t2].drops[d2].dry_food;
                neighbor_solution[h].trips[t2].perishable_food_pickup = neighbor_solution[h].trips[t2].perishable_food_pickup - solution[h].trips[t2].drops[d2].perishable_food + neighbor_solution[h].trips[t2].drops[d2].perishable_food;
                neighbor_solution[h].trips[t2].other_supplies_pickup = neighbor_solution[h].trips[t2].other_supplies_pickup - solution[h].trips[t2].drops[d2].other_supplies + neighbor_solution[h].trips[t2].drops[d2].other_supplies;

                //Need to update helicopter stats
                int trip_distance_t1_old = tripDistance(solution[h].trips[t1], problem, solution[h].helicopter_id);
                int trip_distance_t1_new = tripDistance(neighbor_solution[h].trips[t1], problem, neighbor_solution[h].helicopter_id);
                helicopter_stats[h].total_distance = helicopter_stats[h].total_distance - trip_distance_t1_old + trip_distance_t1_new;
                int trip_distance_t2_old = tripDistance(solution[h].trips[t2], problem, solution[h].helicopter_id);
                int trip_distance_t2_new = tripDistance(neighbor_solution[h].trips[t2], problem, neighbor_solution[h].helicopter_id);
                helicopter_stats[h].total_distance = helicopter_stats[h].total_distance - trip_distance_t2_old + trip_distance_t2_new;

                int neighbor_new_cost = cost(neighbor_solution, problem);
                if (neighbor_new_cost > best_cost) {
                    best_cost = neighbor_new_cost;
                    solution = neighbor_solution;
                    return true;
                }
            }
            
        }
    }
    return false;
}


void addVillageToTrip(Solution& solution, const ProblemData& problem) {
    vector<int> candidate_village_indices;
    for (size_t i = 0; i < problem.villages.size(); ++i) {
        if (village_states[i].food_delivered < village_states[i].max_food_needed ||
            village_states[i].supplies_delivered < village_states[i].max_supplies_needed) {
            candidate_village_indices.push_back(i);
        }
    }
    if (candidate_village_indices.empty()) return; // All villages fully supplied

    static mt19937 gen(random_device{}());
    uniform_int_distribution<> village_dist(0, candidate_village_indices.size() - 1);
    int village_global_idx = candidate_village_indices[village_dist(gen)];
    int village_id = problem.villages[village_global_idx].id;

    uniform_int_distribution<> heli_dist(0, solution.size() - 1);
    int h_idx = heli_dist(gen);
    const auto& helicopter = problem.helicopters[h_idx];

    // Decide whether to add to existing trip or create new one. For simplicity, let's pick existing or create if none exist.
    int t_idx;
    if (solution[h_idx].trips.empty()) {
        solution[h_idx].trips.emplace_back();
        t_idx = 0;
    } else {
        uniform_int_distribution<> trip_dist(0, solution[h_idx].trips.size() - 1);
        t_idx = trip_dist(gen);
    }
    Trip& trip = solution[h_idx].trips[t_idx];

    double current_trip_weight = calculateTripWeight(trip, problem); // Assuming calculateTripWeight exists from utils.cpp
    double remaining_weight_capacity = helicopter.weight_capacity - current_trip_weight;
    if (remaining_weight_capacity <= 0) return; // No capacity in this trip

    // Use a greedy or DP allocation strategy. Using greedy from your existing code for consistency.
    OptimalSupplies new_supplies = greedyAllocationByRatio(problem, remaining_weight_capacity, village_global_idx);
    if (!new_supplies.feasible || new_supplies.total_weight <= 0) return;

    Drop new_drop;
    new_drop.village_id = village_id;
    new_drop.dry_food = new_supplies.dry_food;
    new_drop.perishable_food = new_supplies.perishable_food;
    new_drop.other_supplies = new_supplies.other_supplies;

    int best_pos = -1;
    double least_distance_increase = numeric_limits<double>::max();
    double original_trip_distance = calculateTripDistance(trip, h_idx); // From utils.cpp

    for (int pos = 0; pos <= (int)trip.drops.size(); ++pos) {
        Trip temp_trip = trip;
        temp_trip.drops.insert(temp_trip.drops.begin() + pos, new_drop);
        double new_distance = calculateTripDistance(temp_trip, h_idx);
        double distance_increase = new_distance - original_trip_distance;

        if (distance_increase < least_distance_increase) {
            least_distance_increase = distance_increase;
            best_pos = pos;
        }
    }

    double final_new_distance = original_trip_distance + least_distance_increase;
    if (final_new_distance <= helicopter.distance_capacity &&
        helicopter_stats[h_idx].total_distance + least_distance_increase <= problem.d_max) {
        
        // Commit changes to solution
        trip.drops.insert(trip.drops.begin() + best_pos, new_drop);

        // Update trip pickups (from utils.h: rebalanceTripPickups or manual update)
        trip.dry_food_pickup += new_drop.dry_food;
        trip.perishable_food_pickup += new_drop.perishable_food;
        trip.other_supplies_pickup += new_drop.other_supplies;

        // Update global helicopter stats
        helicopter_stats[h_idx].total_distance += least_distance_increase;

        // Update global village states
        village_states[village_global_idx].food_delivered += new_drop.dry_food + new_drop.perishable_food;
        village_states[village_global_idx].supplies_delivered += new_drop.other_supplies;
    }
}

bool removeVillageFromTrip(Solution& solution, const ProblemData& problem) {
    (void)problem; // Unused parameter
    static mt19937 gen(random_device{}());
    std::vector<std::tuple<int, int, int>> all_drops; // Stores <helicopter_idx, trip_idx, drop_idx>

    for (size_t h = 0; h < solution.size(); ++h) {
        for (size_t t = 0; t < solution[h].trips.size(); ++t) {
            for (size_t d = 0; d < solution[h].trips[t].drops.size(); ++d) {
                all_drops.emplace_back(h, t, d);
            }
        }
    }

    if (all_drops.empty()) {
        return false; // No villages in any trip to remove.
    }

    std::uniform_int_distribution<> dist(0, all_drops.size() - 1);
    auto [helicopter_idx, trip_idx, drop_idx] = all_drops[dist(gen)];

    Trip& trip = solution[helicopter_idx].trips[trip_idx];
    Drop removed_drop = trip.drops[drop_idx];
    
    double original_distance = calculateTripDistance(trip, helicopter_idx);

    trip.drops.erase(trip.drops.begin() + drop_idx);

    double new_distance = calculateTripDistance(trip, helicopter_idx);
    double distance_change = new_distance - original_distance;

    // Update trip pickups
    trip.dry_food_pickup -= removed_drop.dry_food;
    trip.perishable_food_pickup -= removed_drop.perishable_food;
    trip.other_supplies_pickup -= removed_drop.other_supplies;

    // Update global helicopter stats (assuming helicopter_stats[helicopter_idx] corresponds to solution[helicopter_idx])
    helicopter_stats[helicopter_idx].total_distance += distance_change;

    // Update global village states
    int village_global_idx = removed_drop.village_id - 1; // Assuming village IDs are 1-indexed
    village_states[village_global_idx].food_delivered -= removed_drop.dry_food + removed_drop.perishable_food;
    village_states[village_global_idx].supplies_delivered -= removed_drop.other_supplies;
    
    return true;
}

bool adjustSuppliesInTrip(Solution& solution, const ProblemData& problem, int& best_cost) {
    static mt19937 gen(random_device{}());
    vector<pair<int, int>> non_empty_trips;
    for (size_t h = 0; h < solution.size(); ++h) {
        for (size_t t = 0; t < solution[h].trips.size(); ++t) {
            if (!solution[h].trips[t].drops.empty()) {
                non_empty_trips.push_back({h, t});
            }
        }
    }
    if (non_empty_trips.empty()) return false;
    uniform_int_distribution<> dist(0, non_empty_trips.size() - 1);
    pair<int, int> selection = non_empty_trips[dist(gen)];
    int h_idx = selection.first;
    int t_idx = selection.second;

    Trip original_trip = solution[h_idx].trips[t_idx]; // Keep a copy in case of failure/worse score
    Trip& trip_to_modify = solution[h_idx].trips[t_idx];
    const auto& helicopter = problem.helicopters[h_idx];

    // First, "return" all supplies from this trip to update global village state accurately.
    for (const auto& drop : trip_to_modify.drops) {
        int v_idx = drop.village_id - 1;
        village_states[v_idx].food_delivered -= (drop.dry_food + drop.perishable_food);
        village_states[v_idx].supplies_delivered -= drop.other_supplies;
    }

    // Reset drops in the trip we are modifying
    for (auto& drop : trip_to_modify.drops) {
        drop.dry_food = 0;
        drop.perishable_food = 0;
        drop.other_supplies = 0;
    }

    double remaining_capacity = helicopter.weight_capacity;
    double w_dry = problem.packages[0].weight;
    double v_dry = problem.packages[0].value;
    double w_perishable = problem.packages[1].weight;
    double v_perishable = problem.packages[1].value;
    double w_other = problem.packages[2].weight;
    double v_other = problem.packages[2].value;

    while (true) {
        double best_ratio = -1.0;
        int best_drop_idx = -1; // Index within trip_to_modify.drops
        int best_supply_type = -1; // 0:dry, 1:perishable, 2:other

        // Find best possible addition across all villages in the trip
        for (size_t i = 0; i < trip_to_modify.drops.size(); ++i) {
            int v_global_idx = trip_to_modify.drops[i].village_id - 1;
            double food_needed = village_states[v_global_idx].max_food_needed - village_states[v_global_idx].food_delivered;
            double other_needed = village_states[v_global_idx].max_supplies_needed - village_states[v_global_idx].supplies_delivered;

            // Check dry food addition
            if (food_needed > 0 && remaining_capacity >= w_dry) {
                if (v_dry / w_dry > best_ratio) { best_ratio = v_dry / w_dry; best_drop_idx = i; best_supply_type = 0; }
            }
            // Check perishable food addition
            if (food_needed > 0 && remaining_capacity >= w_perishable) {
                if (v_perishable / w_perishable > best_ratio) { best_ratio = v_perishable / w_perishable; best_drop_idx = i; best_supply_type = 1; }
            }
            // Check other supply addition
            if (other_needed > 0 && remaining_capacity >= w_other) {
                if (v_other / w_other > best_ratio) { best_ratio = v_other / w_other; best_drop_idx = i; best_supply_type = 2; }
            }
        }

        if (best_supply_type == -1) break; // No more items can be added (due to capacity or demand)

        // Add the best found item unit
        Drop& target_drop = trip_to_modify.drops[best_drop_idx];
        int v_global_idx = target_drop.village_id - 1;

        if (best_supply_type == 0) {
            target_drop.dry_food++;
            village_states[v_global_idx].food_delivered++;
            remaining_capacity -= w_dry;
        } else if (best_supply_type == 1) {
            target_drop.perishable_food++;
            village_states[v_global_idx].food_delivered++;
            remaining_capacity -= w_perishable;
        } else {
            target_drop.other_supplies++;
            village_states[v_global_idx].supplies_delivered++;
            remaining_capacity -= w_other;
        }
    }

    // Rebalance trip pickups based on new allocations
    trip_to_modify.dry_food_pickup = 0;
    trip_to_modify.perishable_food_pickup = 0;
    trip_to_modify.other_supplies_pickup = 0;
    for (const auto& drop : trip_to_modify.drops) {
        trip_to_modify.dry_food_pickup += drop.dry_food;
        trip_to_modify.perishable_food_pickup += drop.perishable_food;
        trip_to_modify.other_supplies_pickup += drop.other_supplies;
    }

    // Check if new solution is better.
    double new_cost = cost(solution, problem);
    if (new_cost > best_cost) {
        best_cost = new_cost;
        return true; // Improvement found
    } else {
        // Revert changes because cost did not improve
        solution[h_idx].trips[t_idx] = original_trip; // Restore original trip data

        // Restore global village states by re-adding original supplies
        for (const auto& drop : original_trip.drops) {
            int v_idx = drop.village_id - 1;
            village_states[v_idx].food_delivered += (drop.dry_food + drop.perishable_food);
            village_states[v_idx].supplies_delivered += drop.other_supplies;
        }
        return false; // No improvement
    }
}


bool moveVillageBetweenTrips(Solution& solution, const ProblemData& problem, int& best_cost) {
    static random_device rd;
    static mt19937 gen(rd());
    if (solution.empty()) return false;
    // --- 1. Select Source Village ---
    int h1_idx = uniform_int_distribution<>(0, solution.size() - 1)(gen);
    if (solution[h1_idx].trips.empty()) return false;
    int t1_idx = uniform_int_distribution<>(0, solution[h1_idx].trips.size() - 1)(gen);
    if (solution[h1_idx].trips[t1_idx].drops.empty()) return false;
    int d1_idx = uniform_int_distribution<>(0, solution[h1_idx].trips[t1_idx].drops.size() - 1)(gen);

    int h2_idx = uniform_int_distribution<>(0, solution.size() - 1)(gen);
    int t2_idx = uniform_int_distribution<>(0, solution[h2_idx].trips.size())(gen);

    Solution neighbor_solution = solution;
    Drop village_to_move = neighbor_solution[h1_idx].trips[t1_idx].drops[d1_idx];

    neighbor_solution[h1_idx].trips[t1_idx].drops.erase(neighbor_solution[h1_idx].trips[t1_idx].drops.begin() + d1_idx);

    neighbor_solution[h1_idx].trips[t1_idx].dry_food_pickup -= village_to_move.dry_food;
    neighbor_solution[h1_idx].trips[t1_idx].perishable_food_pickup -= village_to_move.perishable_food;
    neighbor_solution[h1_idx].trips[t1_idx].other_supplies_pickup -= village_to_move.other_supplies;

    if (t2_idx == (int)solution[h2_idx].trips.size()) {
        // Create a new trip for the destination helicopter
        Trip new_trip;
        new_trip.drops.push_back(village_to_move);
        // Set initial pickups for the new trip
        new_trip.dry_food_pickup = village_to_move.dry_food;
        new_trip.perishable_food_pickup = village_to_move.perishable_food;
        new_trip.other_supplies_pickup = village_to_move.other_supplies;
        neighbor_solution[h2_idx].trips.push_back(new_trip);
    } else {
        // Insert into an existing trip at a random position
        int d2_pos = uniform_int_distribution<>(0, neighbor_solution[h2_idx].trips[t2_idx].drops.size())(gen);
        neighbor_solution[h2_idx].trips[t2_idx].drops.insert(neighbor_solution[h2_idx].trips[t2_idx].drops.begin() + d2_pos, village_to_move);
        
        // --- Incremental Update for Destination Trip (Add supplies) ---
        neighbor_solution[h2_idx].trips[t2_idx].dry_food_pickup += village_to_move.dry_food;
        neighbor_solution[h2_idx].trips[t2_idx].perishable_food_pickup += village_to_move.perishable_food;
        neighbor_solution[h2_idx].trips[t2_idx].other_supplies_pickup += village_to_move.other_supplies;
    }

    if (isValid(neighbor_solution, problem)) {
        // Update helicopter stats incrementally
        double trip_dist_old_h1 = tripDistance(solution[h1_idx].trips[t1_idx], problem, solution[h1_idx].helicopter_id);
        double trip_dist_new_h1 = tripDistance(neighbor_solution[h1_idx].trips[t1_idx], problem, neighbor_solution[h1_idx].helicopter_id);
        helicopter_stats[h1_idx].total_distance += (trip_dist_new_h1 - trip_dist_old_h1);

        if (t2_idx == (int)solution[h2_idx].trips.size()) {
            double trip_dist_new_h2 = tripDistance(neighbor_solution[h2_idx].trips.back(), problem, neighbor_solution[h2_idx].helicopter_id);
            helicopter_stats[h2_idx].total_distance += trip_dist_new_h2;
        } else {
            double trip_dist_old_h2 = tripDistance(solution[h2_idx].trips[t2_idx], problem, solution[h2_idx].helicopter_id);
            double trip_dist_new_h2 = tripDistance(neighbor_solution[h2_idx].trips[t2_idx], problem, neighbor_solution[h2_idx].helicopter_id);
            helicopter_stats[h2_idx].total_distance += (trip_dist_new_h2 - trip_dist_old_h2);
        }

        int neighbor_new_cost = cost(neighbor_solution, problem);
        if (neighbor_new_cost > best_cost) {
            best_cost = neighbor_new_cost;
            solution = neighbor_solution;
            // Clean up empty trips created by removal
            for(auto& plan : solution) {
                plan.trips.erase(remove_if(plan.trips.begin(), plan.trips.end(), [](const Trip& t){ return t.drops.empty(); }), plan.trips.end());
            }
            return true;
        }
    }
    return false;
}

bool mergeTripsInSameHelicopter(Solution& solution, const ProblemData& problem, int& best_cost) {
    static random_device rd;
    static mt19937 gen(rd());

    int h_idx = uniform_int_distribution<>(0, solution.size() - 1)(gen);
    if (solution[h_idx].trips.size() < 2) return false;

    int t1_idx = uniform_int_distribution<>(0, solution[h_idx].trips.size() - 1)(gen);
    int t2_idx;
    do {
        t2_idx = uniform_int_distribution<>(0, solution[h_idx].trips.size() - 1)(gen);
    } while (t1_idx == t2_idx);

    Solution neighbor_solution = solution;
    const Trip& trip2_original = solution[h_idx].trips[t2_idx]; // Get data from original before modification

    neighbor_solution[h_idx].trips[t1_idx].dry_food_pickup += trip2_original.dry_food_pickup;
    neighbor_solution[h_idx].trips[t1_idx].perishable_food_pickup += trip2_original.perishable_food_pickup;
    neighbor_solution[h_idx].trips[t1_idx].other_supplies_pickup += trip2_original.other_supplies_pickup;
    
    neighbor_solution[h_idx].trips[t1_idx].drops.insert(
        neighbor_solution[h_idx].trips[t1_idx].drops.end(),
        trip2_original.drops.begin(), 
        trip2_original.drops.end()
    );
    
    neighbor_solution[h_idx].trips[t2_idx].drops.clear();

    if (isValid(neighbor_solution, problem)) {
        // Get old distances before merge
        double trip_dist_old_t1 = tripDistance(solution[h_idx].trips[t1_idx], problem, solution[h_idx].helicopter_id);
        double trip_dist_old_t2 = tripDistance(solution[h_idx].trips[t2_idx], problem, solution[h_idx].helicopter_id);
        
        // Get new distance of the merged trip (t1)
        double trip_dist_new_merged = tripDistance(neighbor_solution[h_idx].trips[t1_idx], problem, solution[h_idx].helicopter_id);

        // Update total distance: subtract old two distances, add new single distance.
        helicopter_stats[h_idx].total_distance -= (trip_dist_old_t1 + trip_dist_old_t2);
        helicopter_stats[h_idx].total_distance += trip_dist_new_merged;
        
        // Recalculate cost. Merging saves one fixed cost.
        int neighbor_new_cost = cost(neighbor_solution, problem);
        if (neighbor_new_cost > best_cost) {
            best_cost = neighbor_new_cost;
            solution = neighbor_solution;
            // Clean up empty trip (trip2)
            solution[h_idx].trips.erase(remove_if(solution[h_idx].trips.begin(), solution[h_idx].trips.end(), [](const Trip& t){ return t.drops.empty(); }), solution[h_idx].trips.end());
            return true;
        }
    }
    return false;
}


bool addNewTripToHelicopter(Solution& solution, const ProblemData& problem, int& best_cost) {
    static random_device rd;
    static mt19937 gen(rd());

    vector<int> unserved_village_indices;
    for (size_t i = 0; i < problem.villages.size(); ++i) {
        const auto& vs = village_states[i];
        if (vs.food_delivered < vs.max_food_needed || vs.supplies_delivered < vs.max_supplies_needed) {
            unserved_village_indices.push_back(i);
        }
    }
    if (unserved_village_indices.empty()) return false; // No work to do.

    int village_global_idx = unserved_village_indices[uniform_int_distribution<>(0, unserved_village_indices.size() - 1)(gen)];
    int h_idx = uniform_int_distribution<>(0, solution.size() - 1)(gen);
    const auto& helicopter = problem.helicopters[h_idx];

    OptimalSupplies supplies = greedyAllocationByRatio(problem, helicopter.weight_capacity, village_global_idx);
    if (!supplies.feasible) return false; // Allocation failed (e.g., village demand already met by a concurrent process, or capacity too low)

    Trip new_trip;
    Drop drop;
    drop.village_id = problem.villages[village_global_idx].id;
    drop.dry_food = supplies.dry_food;
    drop.perishable_food = supplies.perishable_food;
    drop.other_supplies = supplies.other_supplies;
    new_trip.drops.push_back(drop);
    new_trip.dry_food_pickup = supplies.dry_food;
    new_trip.perishable_food_pickup = supplies.perishable_food;
    new_trip.other_supplies_pickup = supplies.other_supplies;

    Solution neighbor_solution = solution;
    neighbor_solution[h_idx].trips.push_back(new_trip);

    if (isValid(neighbor_solution, problem)) {
        double new_trip_distance = tripDistance(new_trip, problem, helicopter.id); // Calculate distance for stats update

        // Temporarily update helicopter stats for cost calculation
        double old_total_dist = helicopter_stats[h_idx].total_distance;
        helicopter_stats[h_idx].total_distance += new_trip_distance;
        int neighbor_new_cost = cost(neighbor_solution, problem);

        if (neighbor_new_cost > best_cost) {
            best_cost = neighbor_new_cost;
            solution = neighbor_solution; // Commit change
            return true;
        } else {
            helicopter_stats[h_idx].total_distance = old_total_dist; // Revert stats update
        }
    }

    // If invalid or cost not improved, revert global state change made by greedyAllocationByRatio.
    village_states[village_global_idx].food_delivered -= supplies.dry_food + supplies.perishable_food;
    village_states[village_global_idx].supplies_delivered -= supplies.other_supplies;
    return false;
}


pair<int, Solution> generate_random_neighbor(Solution current_solution, const ProblemData& problem) {
    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<> prob(0.0, 1.0);

    double choice = prob(gen);

    Solution neighbor_solution = current_solution;
    
    int dummy_cost = 0; // The neighborhood functions update cost by reference, we ignore it here.
    if (choice < 0.5) {
        moveVillageBetweenTrips(neighbor_solution, problem, dummy_cost);
    } else {
        SwapVillagesInSameTrip(neighbor_solution, problem, dummy_cost);
    }

    int new_cost = cost(neighbor_solution, problem);
    return make_pair(new_cost, neighbor_solution);
}
