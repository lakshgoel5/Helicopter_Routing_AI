#include "algorithm.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <set>
#include <numeric>

using namespace std::chrono;
using namespace std;
vector<int> clusters;
// vector<int> villages_left_unassigned;


vector<int> clusterVillages(const ProblemData& problem, int time_limit_seconds) {
    // Initialize all villages to -1 (unassigned).
    vector<int> village_assignment(problem.villages.size(), -1);

    vector<pair<double, int>> helicopter_efficiency;
    for (int i = 0; i < (int)problem.helicopters.size(); i++) {
        double efficiency = problem.helicopters[i].fixed_cost + 10 * problem.helicopters[i].alpha;
        helicopter_efficiency.push_back({efficiency, i});
    }
    sort(helicopter_efficiency.begin(), helicopter_efficiency.end());
    auto start_time = high_resolution_clock::now();

    for (int v = 0; v < (int)problem.villages.size(); v++) {
        auto current_time = high_resolution_clock::now();
        auto elapsed = duration_cast<seconds>(current_time - start_time).count();
        if (elapsed >= time_limit_seconds - 1) break;

        double best_score = std::numeric_limits<double>::max();
        int best_heli_idx = -1;

        for (int i = 0; i < (int)helicopter_efficiency.size(); ++i) {
            int h_idx = helicopter_efficiency[i].second;
            const Helicopter& current_helicopter = problem.helicopters[h_idx];
            int city_index = current_helicopter.home_city_id - 1;

            Point village_pos = problem.villages[v].coords;
            Point heli_home = problem.cities[city_index];

            double dist = distance(village_pos, heli_home);            
            double score = dist * (1 + (double)i / helicopter_efficiency.size());
            
            if (score < best_score) {
                best_score = score;
                best_heli_idx = h_idx;
            }
        }
        if (best_heli_idx != -1) {
            // Store assignment directly by village index
            village_assignment[v] = best_heli_idx;
        }
    }
    return village_assignment;
}

vector<Trip> planTripsForHelicopters(const ProblemData& problem, int helicopter_idx, const vector<int>& village_assignment, int time_limit_seconds) {
    (void)time_limit_seconds;
    vector<Trip> trips;
    const auto& helicopter = problem.helicopters[helicopter_idx];
    static mt19937 gen(std::random_device{}());
    double total_helicopter_distance = 0.0;

    // Keep track of villages serviced by *this helicopter* to avoid re-processing them in subsequent trips.
    // Size is total number of villages for direct global index access.
    vector<bool> serviced_by_this_helicopter(problem.villages.size(), false);

    //Make a trip
    while(true){
        Trip current_trip;
        vector<OptimalSupplies> trip_village_supplies;
        vector<int> trip_village_global_indices; // Stores the sequence of global indices for this trip.
        double distance_in_this_trip = 0.0;
        double weight_in_this_trip = 0.0;

        // --- Refactored: Find available starting villages ---
        vector<int> available_start_villages;
        for(size_t global_idx = 0; global_idx < problem.villages.size(); ++global_idx) {
            // Condition 1: Village must be assigned to this helicopter.
            if (village_assignment[global_idx] != helicopter_idx) continue;
            
            // Condition 2: Village must not have been fully serviced in a previous trip by this helicopter.
            if (serviced_by_this_helicopter[global_idx]) continue;

            const auto& vs = village_states[global_idx];
            if (vs.food_delivered < vs.max_food_needed || vs.supplies_delivered < vs.max_supplies_needed) {
                available_start_villages.push_back(global_idx);
            } else {
                serviced_by_this_helicopter[global_idx] = true; // Mark as done if already serviced.
            }
        }
        if(available_start_villages.empty()) break; // All assigned villages are serviced.

        uniform_int_distribution<> start_dist(0, available_start_villages.size() - 1);
        int current_global_idx = available_start_villages[start_dist(gen)];

        // visited_in_trip tracks villages added to the *current trip being built*.
        // Size is total number of villages for direct global index access.
        vector<bool> visited_in_trip(problem.villages.size(), false);

        while(true){
            visited_in_trip[current_global_idx] = true;
            double remaining_capacity = helicopter.weight_capacity - weight_in_this_trip;

            OptimalSupplies supplies;
            //Probabilistic choice between ratio-based and value-based allocation
            if(uniform_real_distribution<>(0.0, 1.0)(gen) < 0.5) 
                supplies = greedyAllocationByValue(problem, remaining_capacity, current_global_idx);
            else 
                supplies = greedyAllocationByRatio(problem, remaining_capacity, current_global_idx);

            double dist_to_village;
            if(trip_village_global_indices.empty()) { // First village in trip
                dist_to_village = home_to_village[helicopter_idx][current_global_idx];
            } else {
                int last_global_idx = trip_village_global_indices.back();
                dist_to_village = village_distances[last_global_idx][current_global_idx];
            }

            double dist_to_home_from_here = home_to_village[helicopter_idx][current_global_idx];
            double trip_dist_if_added = distance_in_this_trip + dist_to_village + dist_to_home_from_here;
            double helicopter_dist_if_added = total_helicopter_distance + trip_dist_if_added;

            bool constraints_violated = trip_dist_if_added > helicopter.distance_capacity || helicopter_dist_if_added > problem.d_max;

            // Add to trip if feasible and constraints pass
            if (supplies.feasible && !constraints_violated) {
                trip_village_global_indices.push_back(current_global_idx);
                trip_village_supplies.push_back(supplies);
                weight_in_this_trip += supplies.total_weight;
                distance_in_this_trip += dist_to_village;

                // Check if this delivery fully services the village to update state for next trip.
                const auto& vs = village_states[current_global_idx];
                if (vs.food_delivered >= vs.max_food_needed && vs.supplies_delivered >= vs.max_supplies_needed) {
                    serviced_by_this_helicopter[current_global_idx] = true;
                }
            } else if (constraints_violated && trip_village_global_indices.empty()) {
                serviced_by_this_helicopter[current_global_idx] = true;
            }
            
            if (!supplies.feasible || constraints_violated) {
                // Stop adding more villages to this trip.
                break;
            }

            double min_dist = std::numeric_limits<double>::max();
            int next_global_idx = -1;

            for(size_t neighbor_global_idx = 0; neighbor_global_idx < problem.villages.size(); ++neighbor_global_idx) {
                if (village_assignment[neighbor_global_idx] != helicopter_idx) continue;
                if (visited_in_trip[neighbor_global_idx]) continue;
                const auto& vs = village_states[neighbor_global_idx];
                if (vs.food_delivered < vs.max_food_needed || vs.supplies_delivered < vs.max_supplies_needed) {
                    double dist = village_distances[current_global_idx][neighbor_global_idx];
                    if(dist < min_dist) {
                        min_dist = dist;
                        next_global_idx = neighbor_global_idx;
                    }
                }
            }

            if(next_global_idx == -1) break; // No more valid neighbors found.
            current_global_idx = next_global_idx;
        }

        if(trip_village_global_indices.empty()) {
            // No villages were successfully added to a trip this iteration.
            // If available_start_villages was not empty, this might indicate an issue or tight constraints.
            // Break to avoid infinite loop if start village selection logic keeps picking a failing start point.
            break; 
        }

        double return_distance = home_to_village[helicopter_idx][trip_village_global_indices.back()];
        total_helicopter_distance += distance_in_this_trip + return_distance;

        for(size_t i = 0; i < trip_village_global_indices.size(); i++) {
            Drop drop;
            drop.village_id = problem.villages[trip_village_global_indices[i]].id; // Get 1-based ID for output
            drop.dry_food = trip_village_supplies[i].dry_food;
            drop.perishable_food = trip_village_supplies[i].perishable_food;
            drop.other_supplies = trip_village_supplies[i].other_supplies;
            current_trip.drops.push_back(drop);
        }
        current_trip.dry_food_pickup = accumulate(trip_village_supplies.begin(), trip_village_supplies.end(), 0, [](int sum, const OptimalSupplies& os){ return sum + os.dry_food; });
        current_trip.perishable_food_pickup = accumulate(trip_village_supplies.begin(), trip_village_supplies.end(), 0, [](int sum, const OptimalSupplies& os){ return sum + os.perishable_food; });
        current_trip.other_supplies_pickup = accumulate(trip_village_supplies.begin(), trip_village_supplies.end(), 0, [](int sum, const OptimalSupplies& os){ return sum + os.other_supplies; });
        trips.push_back(current_trip);
    }

    helicopter_stats[helicopter_idx].total_distance = total_helicopter_distance;
    return trips;
}


Solution random_state(const ProblemData& problem, int time_limit_seconds) {
    //village states reset
    //assign it a length of number of villages
    village_states.clear();
    village_states.resize(problem.villages.size());
    for(int i = 0; i < (int)problem.villages.size(); i++) {
        village_states[i].max_food_needed = 9.0 * problem.villages[i].population;
        village_states[i].max_supplies_needed = 1.0 * problem.villages[i].population;
        village_states[i].food_delivered = 0.0;
        village_states[i].supplies_delivered = 0.0;
    }

    //helicopter stats reset
    helicopter_stats.clear();
    helicopter_stats.resize(problem.helicopters.size());
    for(int i = 0; i < (int)problem.helicopters.size(); i++) {
        helicopter_stats[i] = {0.0};
    }
    
    Solution state;

    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<> prob(0.0, 1.0);

    double p1 = 0.1, p2 = 0.7;
    double choice = prob(gen);

    if(choice < p1) {
        auto start_time = chrono::steady_clock::now();
        
        vector<int> village_assignment = clusterVillages(problem, time_limit_seconds);

        for(int i = 0; i < (int)problem.helicopters.size(); i++){
            auto now = chrono::steady_clock::now();
            chrono::duration<double> elapsed = now - start_time;
            if(elapsed.count() >= time_limit_seconds) break;
            
            HelicopterPlan plan;
            plan.helicopter_id = problem.helicopters[i].id;
            
            // Pass the entire assignment map. The function will filter for helicopter 'i'.
            plan.trips = planTripsForHelicopters(problem, i, village_assignment, time_limit_seconds - elapsed.count());
            state.push_back(plan);
        }
    }
    else if(choice<p1+p2){
        //empty State
        for(int i = 0; i < (int)problem.helicopters.size(); i++){
            HelicopterPlan plan;
            plan.helicopter_id = problem.helicopters[i].id;
            state.push_back(plan);
        }
    }
    else{
        //completely random, but valid state
        for(int i = 0; i < (int)problem.helicopters.size(); i++){
            HelicopterPlan plan;
            plan.helicopter_id = problem.helicopters[i].id;
            state.push_back(plan);
        }
    }

    return state;
}

pair<int, Solution> neighbouringFunction(Solution& current_solution, int current_cost, const ProblemData& problem, int time_limit_seconds) {
    (void)time_limit_seconds;
    Solution best_neighbor = current_solution;
    int best_cost = current_cost;

    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<> prob(0.0, 1.0);

    int total_trips = 0;
    for (const auto& plan : current_solution) {
        total_trips += plan.trips.size();
    }

    // Probabilities for operators: [AddNewTrip, MoveVillage, MergeTrips, SwapSameTrip, ReverseSegment]
    double p_add_new_trip = 0.1;
    double p_move_village = 0.3;
    double p_merge_trips = 0.2;
    double p_swap_same_trip = 0.3;
    // p_reverse_segment makes up the rest

    // Adjust probabilities based on sparsity
    if (total_trips < (int)problem.helicopters.size()) {
        // Solution is very sparse, prioritize adding new trips.
        p_add_new_trip = 0.8;
        p_move_village = 0.1;
        p_merge_trips = 0.0;
        p_swap_same_trip = 0.1;
    } else if (total_trips > problem.villages.size() * 0.5) {
        // Solution is dense, prioritize optimization over adding new trips.
        p_add_new_trip = 0.05;
        p_move_village = 0.35;
        p_merge_trips = 0.2;
        p_swap_same_trip = 0.3;
    }

    double choice = prob(gen);

    if (choice < p_add_new_trip) {
        addNewTripToHelicopter(best_neighbor, problem, best_cost);
    } else if (choice < p_add_new_trip + p_move_village) {
        moveVillageBetweenTrips(best_neighbor, problem, best_cost);
    } else if (choice < p_add_new_trip + p_move_village + p_merge_trips) {
        mergeTripsInSameHelicopter(best_neighbor, problem, best_cost);
    } else if (choice < p_add_new_trip + p_move_village + p_merge_trips + p_swap_same_trip) {
        SwapVillagesInSameTrip(best_neighbor, problem, best_cost);
    } else {
        reverseVillagesListInSameTrip(best_neighbor, problem, best_cost);
    }
    
    // In case no improvement was found in the selected operator, return original cost and solution.
    // The best_neighbor variable is only updated inside the operators if cost improves.
    return make_pair(best_cost, best_neighbor);
}


pair<int, Solution> neighbouringFunctionSA(Solution current_solution, int current_cost, const ProblemData& problem, int time_limit_seconds) {
    (void)time_limit_seconds;
    int best_cost = current_cost;

    static random_device rd;
    static mt19937 gen(rd());
    uniform_real_distribution<> prob(0.0, 1.0);

    int total_trips = 0;
    for (const auto& plan : current_solution) {
        total_trips += plan.trips.size();
    }

    // Probabilities for operators: [AddNewTrip, MoveVillage, MergeTrips, SwapSameTrip, ReverseSegment]
    double p_add_new_trip = 0.1;
    double p_move_village = 0.3;
    double p_merge_trips = 0.2;
    double p_swap_same_trip = 0.3;
    // p_reverse_segment makes up the rest

    // Adjust probabilities based on sparsity
    if (total_trips < (int)problem.helicopters.size()) {
        // Solution is very sparse, prioritize adding new trips.
        p_add_new_trip = 0.8;
        p_move_village = 0.1;
        p_merge_trips = 0.0;
        p_swap_same_trip = 0.1;
    } else if (total_trips > problem.villages.size() * 0.5) {
        // Solution is dense, prioritize optimization over adding new trips.
        p_add_new_trip = 0.05;
        p_move_village = 0.35;
        p_merge_trips = 0.2;
        p_swap_same_trip = 0.3;
    }

    double choice = prob(gen);

    if (choice < p_add_new_trip) {
        addNewTripToHelicopter(current_solution, problem, best_cost);
    } else if (choice < p_add_new_trip + p_move_village) {
        moveVillageBetweenTrips(current_solution, problem, best_cost);
    } else if (choice < p_add_new_trip + p_move_village + p_merge_trips) {
        mergeTripsInSameHelicopter(current_solution, problem, best_cost);
    } else if (choice < p_add_new_trip + p_move_village + p_merge_trips + p_swap_same_trip) {
        SwapVillagesInSameTrip(current_solution, problem, best_cost);
    } else {
        reverseVillagesListInSameTrip(current_solution, problem, best_cost);
    }
    
    // In case no improvement was found in the selected operator, return original cost and solution.
    // The current_solution variable is only updated inside the operators if cost improves.
    return make_pair(best_cost, current_solution);
}

Solution hillClimbing(Solution& current_solution, const ProblemData& problem, int time_limit_seconds, int num_of_steps=10000) {
    auto start_time = chrono::steady_clock::now();
    int current_cost = cost(current_solution, problem);
    bool improved = true;
    int iterations = 0;

    while(iterations!=num_of_steps && improved){
        //Break based on time limit
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed = now - start_time;
        if(elapsed.count() >= time_limit_seconds) break;

        int new_cost;
        Solution new_solution;
        tie(new_cost, new_solution) = neighbouringFunction(current_solution, current_cost, problem, time_limit_seconds - elapsed.count());

        if(new_cost > current_cost){
            current_cost = new_cost;
            current_solution = new_solution;
        }
        else{
            improved = false;
        }
        iterations++;
    }
    return current_solution;
}

Solution hillClimbingWithRestarts(const ProblemData& problem, double time_limit_seconds) {
    auto start_time = chrono::steady_clock::now();

    Solution best_solution;
    int best_value = std::numeric_limits<int>::min();

    while(true){
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed = now - start_time;
        if(elapsed.count() >= time_limit_seconds) break;

        Solution current_solution = random_state(problem, time_limit_seconds - elapsed.count());

        current_solution = hillClimbing(current_solution, problem, time_limit_seconds - elapsed.count());

        int current_value = cost(current_solution, problem);
        if (current_value > best_value) {
            best_value = current_value;
            best_solution = current_solution;
        }

    }
    return best_solution;
}

Solution hillClimbingLocalSearch(Solution current_solution, const ProblemData& problem, double time_limit_seconds) {
    auto start_time = chrono::steady_clock::now();
    int current_cost = cost(current_solution, problem);
    
    while(true){
        auto now = chrono::steady_clock::now();
        chrono::duration<double> elapsed = now - start_time;
        if(elapsed.count() >= time_limit_seconds) break;

        int new_cost;
        Solution new_solution;
        tie(new_cost, new_solution) = neighbouringFunction(current_solution, current_cost, problem, time_limit_seconds - elapsed.count());

        if(new_cost > current_cost){
            current_cost = new_cost;
            current_solution = new_solution;
        }
    }
    return current_solution;
}


Solution simulated_annealing(Solution current_solution, const ProblemData& problem, double time_limit_seconds) {
    auto start_time = steady_clock::now();
    
    const double INITIAL_TEMPERATURE = 100000.0;
    const double COOLING_RATE = 0.9995; // Slow cooling rate for thorough search
    const double FINAL_TEMPERATURE = 1e-4;

    double current_temperature = INITIAL_TEMPERATURE;
    int current_cost = cost(current_solution, problem);
    Solution best_solution = current_solution;
    int best_cost = current_cost;

    static mt19937 gen(random_device{}());
    uniform_real_distribution<> prob_dist(0.0, 1.0);

    while (current_temperature > FINAL_TEMPERATURE) {
        // Time check
        auto now = steady_clock::now();
        chrono::duration<double> elapsed = now - start_time;
        if (elapsed.count() >= time_limit_seconds) break;

        // 1. Generate a neighbor candidate
        // For SA, we generate a neighbor based on the *current* solution, not best_solution

        // std::vector<std::vector<double>> village_distances_copy = village_distances;
        // std::vector<std::vector<double>> home_to_village_copy = home_to_village;


        pair<int, Solution> neighbor_data = neighbouringFunctionSA(current_solution, current_cost, problem, time_limit_seconds - elapsed.count());
        int neighbor_cost = neighbor_data.first;
        Solution neighbor_solution = neighbor_data.second;

        // 2. Decide whether to accept the neighbor state
        if (neighbor_cost > current_cost) {
            // Improvement: always accept
            current_solution = neighbor_solution;
            current_cost = neighbor_cost;
            if (current_cost > best_cost) {
                best_solution = current_solution;
                best_cost = current_cost;
            }
        } else {
            // Worse solution: accept probabilistically
            double cost_difference = neighbor_cost - current_cost; // This value is negative
            double acceptance_probability = exp(cost_difference / current_temperature);
            
            if (prob_dist(gen) < acceptance_probability) {
                current_solution = neighbor_solution;
                current_cost = neighbor_cost;
            }
            // else{
                // village_distances = village_distances_copy;
                // home_to_village = home_to_village_copy;
            // }
        }

        // 3. Cool down
        current_temperature *= COOLING_RATE;
    }
    return best_solution;
}

Solution simulated_annealing_with_restarts(const ProblemData& problem, double time_limit_seconds) {
    auto start_time = steady_clock::now();
    Solution best_solution_overall;
    int best_cost_overall = numeric_limits<int>::min();

    while (true) {
        auto now = steady_clock::now();
        chrono::duration<double> elapsed = now - start_time;
        if (elapsed.count() >= time_limit_seconds) break;

        Solution initial_solution = random_state(problem, time_limit_seconds - elapsed.count());

        // now = steady_clock::now();
        // elapsed = now - start_time;
        // if (elapsed.count() >= time_limit_seconds) break;

        Solution sa_run_solution = simulated_annealing(initial_solution, problem, time_limit_seconds - elapsed.count());

        int sa_run_cost = cost(sa_run_solution, problem);
        if (sa_run_cost > best_cost_overall) {
            best_cost_overall = sa_run_cost;
            best_solution_overall = sa_run_solution;
        }
    }
    return best_solution_overall;
}
