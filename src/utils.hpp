#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include <fstream>
#include <random>
#include "graph.hpp"

void readGraph(FILE* input, Graph& tau, bool labeled);
void skipGraph(FILE* input);

std::vector<Graph> readInput(FILE* input, bool labeled);

inline long double logsumexp(std::vector<long double>& nums) {
    long double max_exp = nums[0], sum = 0.0;
    for (int i = 1 ; i < nums.size() ; i++)
        if (nums[i] > max_exp)
            max_exp = nums[i];
    for (int i = 0; i < nums.size() ; i++)
        sum += exp(nums[i] - max_exp);
    return log(sum) + max_exp;
}

inline long double logsumexp(long double num1, long double num2) {
    long double max_exp = num2;
    if (num1 > max_exp)
        max_exp = num1;
    long double sum = exp(num1 - max_exp) + exp(num2 - max_exp);
    return log(sum) + max_exp;
}

inline long double logsumexp2(long double num1, long double num2) {
    long double max_exp = num2;
    if (num1 > max_exp)
        max_exp = num1;
    long double sum = exp2(num1 - max_exp) + exp2(num2 - max_exp);
    return log2(sum) + max_exp;
}


#endif
