#include "solver.h"
#include <iostream>
#include <chrono>
#include <limits> // for numeric_limits
#include <random>
#include <queue>
#include <algorithm>
#include <set>
#include <fstream>
using namespace std;
const double INF = numeric_limits<double>::infinity();

std::chrono::steady_clock::time_point GLOBAL_DEADLINE;

vector<double> d_clos; // dist of closest village for each helicopter , index using h.id
vector<Helicopter> helicopters;
vector<Village> villages;
vector<Point> cities;          
vector<PackageInfo> packages;    
double d_max;                     
double time_limit_minutes;     
vector<int> food_demand;
vector<int> other_supplies_demand;
int best_package; //index of package with best value/weight ratio
vector<vector<vector<int>>> village_dist_wrt_city; // stores k least extra dist villages for each village wrt each city
vector<vector<int>> city_dist; // stores k closest villages for each city
int num_villages;
int num_cities;
int num_neighbors;  // k
int num_helis;
vector<double> rem_dist; // index using h.id - 1
vector<vector<int>> best_pickups; // index using h.id - 1

int get_best_package () {
    int ans = -1;
    int best_ratio = 0;
    for(int i = 0; i<3; i++){
        int curr_ratio = packages[i].value/packages[i].weight;
        if(curr_ratio>best_ratio){
            ans = i;
            best_ratio = curr_ratio;
        }
    }
    return ans;
}

void compute_village_dist_wrt_city() {
    for(int i = 0; i<num_villages; i++) {
        vector<vector<int>> v (num_cities);
        for(int j = 0; j<num_cities; j++) {
            priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
            for(int k = 0; k<num_villages; k++) {
                if(k == i) continue;
                double extra_dist = distance(villages[i].coords, villages[k].coords) + distance(villages[k].coords, cities[j]) - distance(villages[i].coords, cities[j]);
                pq.push({extra_dist, k});
            }
            for(int k = 0; k<num_neighbors; k++) {
                v[j].push_back(pq.top().second);
                pq.pop();
            }
        }
        village_dist_wrt_city.push_back(v);
    }
}

void compute_city_dist() {
    city_dist.assign(num_cities, vector<int>(num_villages, 0.0));
    for(int i = 0; i<num_cities; i++) {
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
        for(int j = 0; j<num_villages; j++) {
            pq.push({distance(cities[i], villages[j].coords), j});
        }
        for(int j = 0; j<num_neighbors; j++) {
            city_dist[i][j] = pq.top().second;
            pq.pop();
        }
    }
}

void compute_d_clos() {
    d_clos.assign(helicopters.size() + 1, 0.0);
    for (const auto &heli : helicopters) {
        double best_dist = INF;
        Point home = cities[heli.home_city_id-1];

        for (const auto &v : villages) {
            double dist = distance(home, v.coords);
            if (dist < best_dist) {
                best_dist = dist;
            }
        }
        d_clos[heli.id] = best_dist;
    }
}

void discard_helicopters () {
    vector<Helicopter> refined;
    double wt_best_package = packages[best_package].weight;
    double val_best_package = packages[best_package].value;
    for(auto &heli : helicopters){
        if(d_max < 2*d_clos[heli.id]) continue;
        double min_cost = heli.fixed_cost + 2*heli.alpha*d_clos[heli.id];
        double max_value_delivered = (ceil(heli.weight_capacity/wt_best_package)*val_best_package);
        if(max_value_delivered<=min_cost) continue;
        refined.push_back(heli);
    }
    helicopters = refined;
    return;
}

vector<int> calc_best_pickup (Helicopter& heli) {
    vector<int> max_units (3);
    for(int i = 0; i<3; i++) {
        if(i==best_package){
            max_units[i] = heli.weight_capacity / packages[i].weight;
            continue;
        }
        int high = heli.weight_capacity / packages[i].weight;
        int low = 0;
        while(low < high) {
            int mid = (low + high+1)/2;
            int x = mid*(packages[i].weight)/packages[best_package].weight;
            int y = mid*(packages[i].value)/packages[best_package].value;
            if(x!=y) high = mid - 1;
            else low = mid;
        }
        max_units[i] = low;
    }
    double best_value = 0;
    vector<int> best_pickup (3, 0);
    int greatest = 0;
    int least = 1e9;
    int greatest_index = -1;
    int least_index = -1;
    for(int i = 0; i<3; i++){
        if(max_units[i]>greatest){
            greatest = max_units[i];
            greatest_index = i;
        }
    }
    for(int i = 0; i<3; i++){
        if (i==greatest_index) continue;
        if(max_units[i]<least){
            least = max_units[i];
            least_index = i;
        }
    }   
    int middle_index = 3 - greatest_index - least_index;
    for(int i = 0; i<=least; i++){
        for(int j = 0; j<=max_units[middle_index]; j++){
            int k = (heli.weight_capacity - i*packages[least_index].weight - j*packages[middle_index].weight)/packages[greatest_index].weight;
            if(k<=0) continue;
            double curr_value = i*packages[least_index].value + j*packages[middle_index].value + k*packages[greatest_index].value;
            if(curr_value>best_value){
                best_value = curr_value;
                best_pickup[least_index] = i;
                best_pickup[middle_index] = j;
                best_pickup[greatest_index] = k;
            }
        }
    }
    return best_pickup;
}

int sample(const vector<double> &prob_weights) {
    static random_device rd;
    static mt19937 gen(rd());
    discrete_distribution<int> dist(prob_weights.begin(), prob_weights.end());
    return dist(gen);
}
double calc_value (int d, int p, int o){
    return d*packages[0].value + p*packages[1].value + o*packages[2].value;
}

int find_next_village (Helicopter &heli, int curr_village, double& dist_travelled, set<int> &visited) {
    vector<double> prob_weights(num_neighbors, 0.0);
    for(int i = 0; i<num_neighbors; i++){
        int v = village_dist_wrt_city[curr_village][heli.home_city_id - 1][i];

        if(visited.count(v)) {
            prob_weights[i] = 0;
            continue;
        }
        double extra_dist = distance(villages[curr_village].coords, villages[v].coords) + distance(villages[v].coords, cities[heli.home_city_id - 1]) - distance(villages[curr_village].coords, cities[heli.home_city_id - 1]);
        if(dist_travelled + extra_dist > min(rem_dist[heli.id - 1], heli.distance_capacity)){
            prob_weights[i] = 0;
            continue;
        }
        if(food_demand[v]==0 && other_supplies_demand[v]==0){
            prob_weights[i] = 0;
            continue;
        }
        prob_weights[i] = (extra_dist>0) ?  1/(extra_dist*extra_dist) : 1e9;
    }
    bool all_zero = true;
    for(auto &x : prob_weights) {
        if(x > 0) {
            all_zero = false;
            break;
        }
    }
    if(all_zero) return -1;
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::discrete_distribution<int> dist(prob_weights.begin(), prob_weights.end());
    int chosen_idx = dist(gen);
    int chosen_village = village_dist_wrt_city[curr_village][heli.home_city_id - 1][chosen_idx];
    return chosen_village;               
}

int find_first_village (Helicopter &heli){
    vector<double> prob_weights(num_neighbors, 0.0);
    int home_city = heli.home_city_id - 1;
    for(int i = 0; i<num_neighbors; i++){
        int v = city_dist[home_city][i];
        double d = distance(cities[home_city], villages[v].coords);
        if(2*d > min(rem_dist[heli.id - 1], heli.distance_capacity)){
            prob_weights[i] = 0;
            continue;
        }
        if(food_demand[v]==0 && other_supplies_demand[v]==0){
            prob_weights[i] = 0;
            continue;
        }
        prob_weights[i] = 1/(d*d);
    }
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::discrete_distribution<int> dist(prob_weights.begin(), prob_weights.end());

    int chosen_idx = dist(gen);               // index into city_dist[home_city]
    int chosen_village = city_dist[home_city][chosen_idx]; // actual village id

    return chosen_village;
}

void add_random_trip (Solution &solution) {
    vector<double> prob_weights(num_helis, 0.0);
    for(int i = 0; i<num_helis; i++){
        if (rem_dist[helicopters[i].id - 1] > d_clos[helicopters[i].id] * 2){
            prob_weights[i] = rem_dist[helicopters[i].id - 1]/ (helicopters[i].fixed_cost + 2*helicopters[i].alpha*d_clos[helicopters[i].id]);
        }
        else prob_weights[i] = 0.0;
    }
    double value_delivered;
    int perishable_pickup;
    int dry_pickup;
    int other_pickup;
    vector<double> prob = {0.25,0.3,0.45}; // best_pack, only perishable, 7:2:2
    int kind = sample(prob);
    int h_idx = sample(prob_weights);
    double f = helicopters[h_idx].fixed_cost;
    double alpha = helicopters[h_idx].alpha;
    if(kind==0){
        dry_pickup = best_pickups[helicopters[h_idx].id - 1][0];
        perishable_pickup = best_pickups[helicopters[h_idx].id - 1][1];
        other_pickup = best_pickups[helicopters[h_idx].id - 1][2];
    }
    else if (kind==1){
        dry_pickup = 0;
        perishable_pickup = helicopters[h_idx].weight_capacity/packages[1].weight;
        other_pickup = 0;
    }
    else{
        double denominator = 7*packages[1].weight + 2*packages[0].weight + 2*packages[2].weight;
        int x = helicopters[h_idx].weight_capacity/denominator;
        dry_pickup = 2*x;
        perishable_pickup = 7*x;
        other_pickup = 2*x;
    }
    value_delivered = calc_value(dry_pickup, perishable_pickup, other_pickup);
    Trip trip;
    double dist_travelled = 0.0;
    set<int> visited;
    double max_dist = min(rem_dist[helicopters[h_idx].id - 1], helicopters[h_idx].distance_capacity);
    max_dist = min(max_dist, (value_delivered - helicopters[h_idx].fixed_cost)/helicopters[h_idx].alpha);
    int first_village = find_first_village(helicopters[h_idx]);
    vector<int> path;
    Drop first_drop;
    first_drop.village_id = villages[first_village].id;
    first_drop.perishable_food = min(food_demand[first_village], perishable_pickup);
    first_drop.dry_food = min(food_demand[first_village] - first_drop.perishable_food, dry_pickup);
    first_drop.other_supplies = min(other_supplies_demand[first_village], other_pickup);
    trip.drops.push_back(first_drop);
    perishable_pickup-=first_drop.perishable_food;
    dry_pickup-=first_drop.dry_food;
    other_pickup-=first_drop.other_supplies;
    double value_achieved = calc_value(first_drop.dry_food, first_drop.perishable_food, first_drop.other_supplies);
    path.push_back(first_village);
    visited.insert(first_village);
    dist_travelled+=(2*distance(cities[helicopters[h_idx].home_city_id-1], villages[first_village].coords));
    if(dist_travelled > max_dist) return;
    while (true){
        if (std::chrono::steady_clock::now() >= GLOBAL_DEADLINE) break;
        double cost = value_achieved - (f + 2*alpha*dist_travelled);
        vector<double> bernoulli(2);
        bernoulli[0] = value_achieved; // prob of stopping
        bernoulli[1] = value_delivered - value_achieved; // prob of continuing
        if(cost<=0){
            bernoulli[0] = 0;
            bernoulli[1] = 1;
        }
        int stop = sample(bernoulli);
        if(stop==0) break;
        int next_village = find_next_village(helicopters[h_idx], path.back(), dist_travelled, visited);
        if(next_village==-1) break;
        int prev_village = path.back();
        path.push_back(next_village);
        visited.insert(next_village);
        double extra_dist = distance(villages[prev_village].coords, villages[next_village].coords) + distance(villages[next_village].coords, cities[helicopters[h_idx].home_city_id- 1]) - distance(villages[prev_village].coords, cities[helicopters[h_idx].home_city_id - 1]);\
        dist_travelled+=extra_dist;
        Drop next_drop;
        next_drop.village_id = villages[next_village].id;
        next_drop.perishable_food = min(food_demand[next_village], perishable_pickup);
        next_drop.dry_food = min(food_demand[next_village] - next_drop.perishable_food, dry_pickup);
        next_drop.other_supplies = min(other_supplies_demand[next_village], other_pickup);
        trip.drops.push_back(next_drop);
        perishable_pickup-=next_drop.perishable_food;
        dry_pickup-=next_drop.dry_food;
        other_pickup-=next_drop.other_supplies;
        value_achieved+=calc_value(next_drop.dry_food, next_drop.perishable_food, next_drop.other_supplies);
    }
    if(value_achieved > f + 2*alpha*dist_travelled){
        trip.dry_food_pickup = 0;
        trip.perishable_food_pickup = 0;
        trip.other_supplies_pickup = 0;
        for(auto drop : trip.drops){
            trip.dry_food_pickup+=drop.dry_food;
            trip.perishable_food_pickup+=drop.perishable_food;
            trip.other_supplies_pickup+=drop.other_supplies;
            
            int curr_village = drop.village_id - 1;
            food_demand[curr_village]-=drop.dry_food;
            food_demand[curr_village]-=drop.perishable_food;
            other_supplies_demand[curr_village] -= drop.other_supplies;
        }
        rem_dist[helicopters[h_idx].id - 1]-=dist_travelled;
        solution[helicopters[h_idx].id - 1].helicopter_id = helicopters[h_idx].id;
        solution[helicopters[h_idx].id - 1].trips.push_back(trip);
    }
}

void add_initial_trips (Solution &solution) {
    if (std::chrono::steady_clock::now() >= GLOBAL_DEADLINE) return;

    vector<pair<double,pair<int,int>>> village_heli_cost;
    for(int i = 0; i<num_villages; i++){
        for(int j = 0; j<num_helis; j++){
            double cost = helicopters[j].fixed_cost + 2*helicopters[j].alpha*distance(cities[helicopters[j].home_city_id-1], villages[i].coords);
            village_heli_cost.push_back({cost, {i,j}});
        }
    }   
    sort(village_heli_cost.begin(), village_heli_cost.end());
    for(auto &p : village_heli_cost){
        if (std::chrono::steady_clock::now() >= GLOBAL_DEADLINE) return;
        double perishable_wt = packages[1].weight;
        double perishable_val = packages[1].value;
        double other_wt = packages[2].weight;
        double other_val = packages[2].value;
        int v = p.second.first;
        int h = p.second.second;
        double d = distance(cities[helicopters[h].home_city_id-1], villages[v].coords);
        if(2*d > helicopters[h].distance_capacity) continue;
        int f = food_demand[v];
        int o = other_supplies_demand[v];
        if(f==0 && o==0) continue;
        int x = (helicopters[h].weight_capacity)/(9*perishable_wt + other_wt);
        if(x==0) continue;
        double trip_value = (9*x*perishable_val + x*other_val);
        if(trip_value <= p.first) continue;
        int food_units_per_trip = 9*x;
        int num_trips = f/food_units_per_trip;
        num_trips = min(num_trips, int((rem_dist[helicopters[h].id - 1])/(2*d)));
        for(int i = 0; i<num_trips; i++){
            Trip trip;
            trip.dry_food_pickup = 0;
            trip.perishable_food_pickup = food_units_per_trip;
            trip.other_supplies_pickup = x;
            Drop drop;
            drop.village_id = villages[v].id;
            drop.dry_food = trip.dry_food_pickup;
            drop.perishable_food = trip.perishable_food_pickup;
            drop.other_supplies = trip.other_supplies_pickup;
            trip.drops.push_back(drop);
            food_demand[v] -= (trip.dry_food_pickup + trip.perishable_food_pickup);
            other_supplies_demand[v] -= trip.other_supplies_pickup;
            rem_dist[helicopters[h].id - 1] -= 2*d;
            HelicopterPlan &plan = solution[helicopters[h].id - 1];
            plan.helicopter_id = helicopters[h].id;
            plan.trips.push_back(trip);
        }
    }
    for(auto &p : village_heli_cost){
        if (std::chrono::steady_clock::now() >= GLOBAL_DEADLINE) return;
        int v = p.second.first;
        int h = p.second.second;
        double d = distance(cities[helicopters[h].home_city_id-1], villages[v].coords);
        if(2*d  > helicopters[h].distance_capacity) continue;
        int food_units_per_trip = 0;
        for(int i = 0; i< (int)packages.size() - 1; i++) food_units_per_trip += best_pickups[h][i];
        int other_units_per_trip = best_pickups[h][packages.size() - 1];
        int dry_food_pickup = best_pickups[h][0];
        int perishable_food_pickup = best_pickups[h][1];
        int other_pickup = best_pickups[h][2];
        double dry_food_value = packages[0].value;
        double perishable_value = packages[1].value;
        double other_value = packages[2].value;
        double trip_value = dry_food_pickup*dry_food_value + perishable_food_pickup*perishable_value + other_pickup*other_value;
        if(trip_value <= p.first) continue;
        int num_trips = (food_units_per_trip==0) ? 1e9 : food_demand[v]/food_units_per_trip;
        num_trips = (other_units_per_trip == 0) ? num_trips : min(num_trips, other_supplies_demand[v]/other_units_per_trip);
        num_trips = min(num_trips, int((rem_dist[helicopters[h].id - 1])/(2*d)));
        for(int i = 0; i<num_trips; i++){
            Trip trip;
            trip.dry_food_pickup = dry_food_pickup;
            trip.perishable_food_pickup = perishable_food_pickup;
            trip.other_supplies_pickup = other_pickup;
            Drop drop;
            drop.village_id = villages[v].id;
            drop.dry_food = trip.dry_food_pickup;
            drop.perishable_food = trip.perishable_food_pickup;
            drop.other_supplies = trip.other_supplies_pickup;
            trip.drops.push_back(drop);
            food_demand[v] -= (trip.dry_food_pickup + trip.perishable_food_pickup);
            other_supplies_demand[v] -= trip.other_supplies_pickup;
            rem_dist[helicopters[h].id - 1] -= 2*d;
            HelicopterPlan &plan = solution[helicopters[h].id - 1];
            plan.helicopter_id = helicopters[h].id;
            plan.trips.push_back(trip);
        }
    }
}


Solution greedy (const ProblemData& problem, double time_limit) {
    time_limit_minutes = time_limit;
    auto start = std::chrono::steady_clock::now();
    GLOBAL_DEADLINE = start + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(time_limit_minutes * 0.95));

    helicopters = problem.helicopters;
    villages = problem.villages;
    cities = problem.cities;
    packages = problem.packages;
    d_max = problem.d_max;
    food_demand.resize(villages.size());
    for(int i = 0; i<(int)villages.size(); i++) food_demand[i] = 9*villages[i].population;
    other_supplies_demand.resize(villages.size());
    for(int i = 0; i<(int)villages.size(); i++) other_supplies_demand[i] = villages[i].population;
    num_villages = villages.size();
    num_cities = cities.size();
    num_neighbors = max(100, num_villages/100);
    num_neighbors = min(num_neighbors, num_villages-1);
    if(packages[0].weight >= packages[1].weight) packages[0].weight = INF;
    best_pickups.resize(helicopters.size());
    for(int i = 0; i<(int)helicopters.size(); i++){
        best_pickups[i] = calc_best_pickup(helicopters[i]);
    }
    rem_dist.assign(helicopters.size(), d_max);
    
    compute_d_clos();
    discard_helicopters();
    compute_village_dist_wrt_city();
    compute_city_dist();
    best_package = get_best_package();
    num_helis = helicopters.size();

    Solution solution;
    for (const auto &heli : problem.helicopters) {
        HelicopterPlan plan;
        plan.helicopter_id = heli.id;
        plan.trips = {};
        solution.push_back(plan);
    }
    add_initial_trips(solution);
    while(true){
        if (std::chrono::steady_clock::now() >= GLOBAL_DEADLINE) break;
        add_random_trip(solution);
    }
    return solution;
}