<div align="center">

  <h1 align="center">Helicopter Routing AI</h1>
  <p align="center">
    An intelligent search-based system for optimally routing helicopters to deliver relief supplies to flood-affected villages, balancing capacity, range, and cost constraints to maximize aid value.
  </p>

</div>

## Authors

- **Laksh Goel** (2023CS10848)
- **Rachit Bhalani** (2023CS10961)

## Table of Contents
- [Overview](#overview)
- [Problem Statement](#problem-statement)
- [Requirements](#requirements)
- [How to Run](#how-to-run)
  - [Compilation](#compilation)
  - [Execution](#execution)
  - [Input/Output Format](#inputoutput-format)
- [Algorithm Explanations](#algorithm-explanations)
  - [Greedy Algorithm](#greedy-algorithm)
  - [Hill Climbing with Restarts](#hill-climbing-with-restarts)
  - [Simulated Annealing with Restarts](#simulated-annealing-with-restarts)
  - [Neighborhood Functions](#neighborhood-functions)
- [Key Learnings](#key-learnings)

---

## Overview

This project tackles a complex AI search problem where helicopters must be optimally routed to deliver relief supplies (dry food, perishable food, and other supplies) to flood-affected villages. The solution must maximize the total value of aid delivered while respecting multiple constraints:

- **Helicopter constraints**: Weight capacity, distance capacity per trip, fixed operating costs, and per-distance costs
- **Global constraint**: Maximum total distance allowed across all helicopters
- **Village needs**: Each village has specific population-based requirements for food (9 units per person) and other supplies (1 unit per person)

The objective is to maximize:
```
Total Value Delivered - Total Helicopter Operating Costs
```

This is formulated as a **search problem** and solved using local search algorithms with sophisticated neighborhood exploration strategies.

---

## Problem Statement

The complete problem statement can be found in the [A1.pdf](A1.pdf) document included in this repository. The problem involves:

1. **N villages** with known locations and population
2. **M helicopters** with varying capacities and costs, each based in a specific city
3. **Three types of supplies** with different weights and values
4. **Routing challenge**: Design trips for each helicopter to maximize net value delivery

The problem is NP-hard and requires intelligent heuristic search strategies to find high-quality solutions within the time limit.

---

## Requirements

### Build Dependencies
- **C++ Compiler**: g++ with C++17 support
- **Make**: GNU Make for build automation
- **Operating System**: Linux/Unix (tested on Ubuntu)

### Development Tools
- Git (for version control)
- Address Sanitizer support (included in compilation flags for debugging)

---

## How to Run

### Compilation

The project uses a Makefile for easy compilation. To build the project:

```bash
# Clean previous builds
make clean

# Build the main solver
make

# Build the format checker (optional)
make checker
```

The compilation produces an executable named `main` in the project directory.

### Execution

To run the helicopter routing solver:

```bash
./main <input_filename> <output_filename>
```

**Example:**
```bash
./main SampleInputOutput/input1.txt SampleInputOutput/output1.txt
```

The solver will:
1. Read the problem instance from the input file
2. Apply the hybrid algorithm (Greedy + Hill Climbing + Simulated Annealing)
3. Write the solution to the output file
4. Display progress and final cost on the console

### Input/Output Format

**Input Format:**
```
<time_limit_minutes>
<d_max>
<package_weight_dry> <package_value_dry> <package_weight_perishable> <package_value_perishable> <package_weight_other> <package_value_other>
<num_cities> <city1_x> <city1_y> ... <cityN_x> <cityN_y>
<num_villages> <village1_id> <village1_x> <village1_y> <village1_population> ... <villageN_id> <villageN_x> <villageN_y> <villageN_population>
<num_helicopters> <heli1_id> <heli1_home_city> <heli1_weight_capacity> <heli1_distance_capacity> <heli1_fixed_cost> <heli1_alpha> ...
```

**Output Format:**
```
<num_trips_helicopter_1>
<trip1_dry_food_pickup> <trip1_perishable_pickup> <trip1_other_pickup> <num_drops>
<village_id> <dry_food_delivered> <perishable_delivered> <other_delivered>
...
```

See `SampleInputOutput/` directory for examples.

---

## Algorithm Explanations

The solver employs a **three-phase hybrid approach** combining greedy construction with two metaheuristic search algorithms:

### Greedy Algorithm

**Purpose**: Generate a high-quality initial solution quickly (uses ~40% of total time budget)

**Key Strategies**:
1. **Helicopter Filtering**: Discards helicopters that cannot deliver value exceeding their minimum operating cost
2. **Best Package Identification**: Identifies packages with the best value-to-weight ratio
3. **Probabilistic Trip Construction**:
   - Selects starting villages probabilistically based on distance from helicopter home cities
   - Iteratively adds nearby villages to trips using distance-weighted probability distributions
   - Uses multiple pickup strategies (best-ratio packages, perishable-only, balanced 7:2:2 ratio)
4. **Initial Trip Assignment**: Creates direct single-village trips where cost-effective

**Advantages**: Fast, produces feasible solutions with good baseline quality

### Hill Climbing with Restarts

**Purpose**: Local optimization phase (uses ~10% of remaining time)

**Algorithm**:
```
1. Generate random initial solution
2. Repeat until no improvement:
   a. Generate neighboring solution using random operator
   b. If neighbor has better cost, accept it
   c. Otherwise, stop climbing
3. Restart from new random solution
4. Return best solution found across all restarts
```

**Characteristics**:
- **Exploitation-focused**: Quickly finds local optima
- **Multiple restarts**: Explores different regions of search space
- **Deterministic acceptance**: Only accepts improvements

### Simulated Annealing with Restarts

**Purpose**: Global optimization phase with escape capability (uses ~90% of remaining time)

**Algorithm**:
```
1. Start with solution from Hill Climbing
2. Set initial temperature T_0 = 100,000
3. While temperature > T_min:
   a. Generate neighboring solution
   b. If neighbor is better, accept it
   c. If neighbor is worse, accept with probability exp(ΔE/T)
   d. Cool temperature: T = T × 0.9995
4. Restart from new random solution periodically
5. Return best solution found overall
```

**Parameters**:
- Initial temperature: 100,000
- Cooling rate: 0.9995 (slow cooling for thorough exploration)
- Final temperature: 0.0001
- Acceptance probability: `exp((new_cost - current_cost) / temperature)`

**Characteristics**:
- **Balanced exploration/exploitation**: Accepts worse solutions probabilistically
- **Escapes local optima**: Temperature allows uphill moves early on
- **Gradual refinement**: As temperature decreases, becomes more selective

### Neighborhood Functions

The solution quality heavily depends on effective neighborhood operators that explore diverse states:

#### 1. **Add New Trip** (Probability: 10% normally, 80% if solution is sparse)
- Randomly selects an unserved village
- Creates a new trip for a randomly selected helicopter
- Allocates supplies using greedy strategy

#### 2. **Move Village Between Trips** (Probability: 30%)
- Removes a village from one trip
- Inserts it into another trip (possibly different helicopter)
- Can create new trips or merge into existing ones
- Updates pickup amounts and distances incrementally

#### 3. **Merge Trips** (Probability: 20%)
- Combines two trips of the same helicopter into one
- Saves one fixed cost
- Useful when trips are short and can be consolidated

#### 4. **Swap Villages in Same Trip** (Probability: 30%)
- Changes the order of village visits within a trip
- Can reduce travel distance through better routing
- Implements local search on trip sequences

#### 5. **Reverse Segment (2-opt)** (Probability: 10%)
- Reverses a segment of villages in a trip
- Classic TSP optimization technique
- Eliminates crossing paths

**Adaptive Probability Adjustment**:
```cpp
if (total_trips < num_helicopters) {
    // Sparse solution: prioritize adding trips
    p_add_new_trip = 0.8
} else if (total_trips > 0.5 × num_villages) {
    // Dense solution: prioritize optimization
    p_add_new_trip = 0.05
}
```

---

## Key Learnings

### 1. **Importance of Comprehensive Neighborhood Functions**
The most critical factor in solution quality is having neighborhood operators that can explore **most of the search space**. Our implementation uses five different operators with adaptive probabilities to ensure:
- **Diversification**: Add/remove villages, create/merge trips
- **Intensification**: Optimize village ordering within trips
- **Flexibility**: Move villages between helicopters and trips

Without comprehensive exploration, the search gets trapped in poor local optima.

### 2. **Limitations of Pure Greedy Approaches**
Greedy algorithms perform well for initial construction but have severe limitations:
- **Premature search space cutoff**: Greedy choices eliminate potentially optimal branches
- **Inability to backtrack**: Once a decision is made, it cannot be undone
- **Local perspective**: Optimizes immediate decisions without global view
- **Example**: A greedy algorithm might assign the closest village to a helicopter, but this might block a better overall configuration where that village is served by a different helicopter

### 3. **Power of Hybrid Approaches**
Combining multiple algorithms leverages their complementary strengths:
- **Greedy**: Fast baseline solution
- **Hill Climbing**: Quick local refinement
- **Simulated Annealing**: Escape local optima and global exploration

### 4. **Balance Between Exploration and Exploitation**
- Early in the search: Accept more random moves (high temperature in SA, diverse operators)
- Later in the search: Focus on refinement (low temperature, optimize existing trips)
- Adaptive operator probabilities based on solution density ensure appropriate focus

### 5. **State Management is Critical**
Maintaining global state (village demands satisfied, helicopter distances used) across neighborhood operations requires careful bookkeeping to ensure feasibility and avoid recomputation.


