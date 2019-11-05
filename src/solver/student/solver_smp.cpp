#include "..\solver_smp.h"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <random>
#include <tbb/tbb.h>


static double best_fintess = std::numeric_limits<double>::max();
static size_t iteration = 0;

static double step = 0.95;

std::vector<std::vector<double>> matrix;


struct point {
	std::vector<double> position;
	double fitness;
};

double randomPointInBounds(double lower_bound, double upper_bound) {
	//std::cout << lower_bound << ", " << upper_bound << "\n";
	double result = ((double)(rand()) / (double)(RAND_MAX)) - 0.5;
	static std::mt19937_64 generator(123);
	std::uniform_real_distribution<double> distribution(lower_bound, upper_bound);
	result = distribution(generator);
	return result;
}

/*
double fitness(double x, double y)
{
	return (4 - 2.1 * pow(x, 2) + pow(x, 4) / 3) * pow(x, 2) + x * y + (-4 + 4 * pow(y, 2)) * pow(y, 2);
}

*/

void initializePoints(std::vector<point>& points, solver::TSolver_Setup setup) {
	for (size_t i = 0; i < points.size(); i++)
	{
		points[i].position = std::vector<double>(setup.problem_size);
		for (size_t j = 0; j < setup.problem_size; j++)
		{
			points[i].position[j] = randomPointInBounds(setup.lower_bound[j], setup.upper_bound[j]);
		}

	}
}

void printPoints(std::vector<point> points) {
	for (size_t i = 0; i < points.size(); i++)
	{
		std::cout << "[";
		for (size_t j = 0; j < points[i].position.size(); j++)
		{
			std::cout << points[i].position[j] << ", ";
		}
		std::cout << "]\n";
	}
}

point* getBestPoint(std::vector<point>& points,const void* data,solver::TObjective_Function fitness) {
	point* result = &points[0];
	for (size_t i = 0; i < points.size(); i++) {
		double old_fitness = fitness(data, result->position.data());
		double new_fitness = fitness(data, points[i].position.data());
		points[i].fitness = new_fitness;
		if (old_fitness > new_fitness) {
			best_fintess = new_fitness;
			result = &points[i];
		}
	}
	return result;
}

double matrixMultiplication(std::vector<std::vector<double>>& matrix, std::vector<double>& vector,const std::vector<double>& centroid, size_t dimension) {
	double result = 0;
	std::vector<double> row = matrix[dimension];
	for (size_t i = 0; i < vector.size(); i++)
	{
		result += row[i] * (vector[i] - centroid[i]);
	}

	return result;
}

void doSPO(std::vector<point>& points, solver::TSolver_Setup& setup, solver::TSolver_Progress& progress) {
	std::vector<point> new_points = std::vector<point>(points);
	point* center = getBestPoint(points ,setup.data ,setup.objective);

	//for (int i = 0; i < points.size(); i++) {
	tbb::parallel_for((size_t)0, points.size(), [&](int i) {
		for (size_t j = 0; j < points[i].position.size(); j++) {

			new_points[i].position[j] = center->position[j] + step * matrixMultiplication(matrix, points[i].position, center->position, j);// *(points[i].position[j] - center->position[j]);

			if (new_points[i].position[j] > setup.upper_bound[j])
			{
				//new_points[i].position[j] = setup.upper_bound[j];
				new_points[i].position[j] = randomPointInBounds(setup.lower_bound[j], setup.upper_bound[j]);
			}
			else if (new_points[i].position[j] < setup.lower_bound[j])
			{
				//new_points[i].position[j] = setup.lower_bound[j];
				new_points[i].position[j] = randomPointInBounds(setup.lower_bound[j], setup.upper_bound[j]);
			}
			if (j < 2) {
				//std::cout << points[i].position[j] << ",";
			}
		}
		//std::cout << "\n";
		}
	);
	//std::cout << "\n";
	points = new_points;
}



void initMatrix(solver::TSolver_Setup setup) {
	
	matrix = std::vector<std::vector<double>>(setup.problem_size);
	for (size_t i = 0; i < matrix.size(); i++)
	{
		matrix[i] = std::vector<double>(setup.problem_size);
		for (size_t j = 0; j < matrix[i].size(); j++)
		{
			if ((i - 1) == j)
			{
				matrix[i][j] = 1;
			}
			else if (i == 0 && j == matrix.size() - 1) {
				matrix[i][j] = -1;
			}
			else {
				matrix[i][j] = 0;
			}
		}
	}
}

void printMatrix()
{
	for (size_t i = 0; i < matrix.size(); i++)
	{
		std::cout << "[";
		for (size_t j = 0; j < matrix[i].size(); j++)
		{
			std::cout << matrix[i][j] << ",";
		}
		std::cout << "]\n";
	}
}

bool shouldStop(std::vector<point>& points, solver::TSolver_Setup setup)
{
	//double test = 0;
	iteration++;
	//for (size_t i = 0; i < points.size(); i++)
	//{
	//	test += fitness(points[i].position[0], points[i].position[1]);
	//}
	//test /= points.size();
	if (/*0 < (best_fintess - test) && (best_fintess - test) < setup.tolerance ||*/ iteration > setup.max_generations)
	{
		return true;
	}

	return false;
}



HRESULT solve_smp(solver::TSolver_Setup &setup, solver::TSolver_Progress &progress) {
	iteration = 0;
	best_fintess = std::numeric_limits<double>::max();
	std::vector<point> points(setup.population_size);
	initializePoints(points, setup);
	initMatrix(setup);
	//printMatrix();
	while (!shouldStop(points, setup)) {
		doSPO(points, setup, progress);
	}

	point* result = getBestPoint(points,setup.data,setup.objective);
	std::cout << result->position[0] << " -- " << result->position[1] << " |||| " << result->fitness << "\n";
	//printPoints(points);
	for (size_t i = 0; i < setup.problem_size; i++)
	{
		setup.solution[i] = result->position[i];
	}
	
	return S_OK;
}