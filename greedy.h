#ifndef GREEDY_H
#define GREEDY_H

#include <iostream>
#include <chrono>
#include <limits> // for numeric_limits
#include <random>
#include <queue>
#include <algorithm>
#include <set>
#include "structures.h"
using namespace std;

Solution greedy(const ProblemData& problem, double time_limit_minutes);

#endif // GREEDY_H