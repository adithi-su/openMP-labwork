/*
Author: Adithi Upadhya
Class: ECE6122
Last Date Modified: 10/10/2023
Description: Using OpenMP to Calculate the Electric Field Produced by Array of Point Charges
*/

//g++ *.cpp -O3 -lpthread -fopenmp -std=c++17 -o Lab2

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <chrono>
#include <omp.h>


const double k = 9.0e9; // Couloub's constant

void computeFieldAt(double x, double y, double z, double charge, double x_point,  double y_point, double z_point, double& Ex, double& Ey, double& Ez)
{
	//electric field calculation at the point (x,y,z) due to charge at some position
	double dx = x_point - x; //(x-x0)
	double dy = y_point - y; //(y-y0)
	double dz = z_point - z; //(z-z0)
	double r = sqrt(dx * dx + dy * dy + dz * dz); 

	if (r != 0.0) //denominator!=0
	{
		double magnitude = (k * charge) / (r * r * r);

		//calculate E_x, E_y, and E_z
		Ex = magnitude * dx;
		Ey = magnitude * dy;
		Ez = magnitude * dz;

	}
	else
	{
		// if r==0, points are at same location
		Ex = Ey = Ez = 0.0;
	}

}


void calculateE(const std::vector<std::vector<double>>& charges, double point_x, double point_y, double point_z, double separationX, double separationY, double& Ex, double& Ey, double& Ez, int NUM_THREADS)
{
	Ex = 0.0;
	Ey = 0.0;
	Ez = 0.0;

	omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for reduction(+:Ex,Ey,Ez) schedule(dynamic) 
	for (size_t i = 0; i < charges.size(); ++i)
	{   
		for (size_t j = 0; j < charges[i].size(); ++j)
		{
			double q = charges[i][j];
			double Ex_i, Ey_i, Ez_i;
			double tempx = (i - (charges.size() - 1) / 2.0) * separationX;
			double tempy = (j - (charges[i].size() - 1) / 2.0) * separationY;
			computeFieldAt(tempx, tempy, 0.0, q, point_x, point_y, point_z, Ex_i, Ey_i, Ez_i);
			Ex += Ex_i; //sum of fields from point charges along each of the directions
			Ey += Ey_i;
			Ez += Ez_i;
		}
	}
}



int main()
{

	while (true)
	{
		int N, M, NUM_THREADS;
		std::string userInput;
		double x, y, z, q, separationX, separationY;
		double Ex_final = 0.0;
		double Ey_final = 0.0;
		double Ez_final = 0.0;

		std::cout << "Please enter the number of concurrent threads to use: ";
		std::cin >> NUM_THREADS;

		//std::vector<std::thread> threads(NUM_THREADS);

	A:	std::cout << "\nPlease enter the number of rows and columns in the NxM array: ";
		std::cin >> N >> M;
		if (N <= 0 || M <= 0)
		{
			std::cout << "Invalid input for N or M\n";
			goto A;
		}

	B:	std::cout << "Please enter the x and y separation distances in meters: ";
		std::cin >> separationX >> separationY;
		if (separationX <= 0.0 || separationY <= 0.0)
		{
			std::cout << "Invalid input for separation distances\n";
			goto B;
		}

	C:	std::cout << "Please enter the common charge on the points in micro C: ";
		if (!(std::cin >> q))
		{
			std::cout << "Invalid input for charge q\n";
			goto C;
		}

		q /= 1000000.0; // to micro C

		//std::vector<std::vector<ECE_ElectricField>> field(N, std::vector<ECE_ElectricField>(M, ECE_ElectricField(0.0, 0.0, 0.0, q / 1000000.0)));
		std::vector<std::vector<double>> charges(N, std::vector<double>(M, q));
		//container for all point charges in NxM grid, divide q for micro C


	D:	std::cout << "Please enter the location in space to determine the electric field (x y z) in meters: ";
		if (!(std::cin >> x >> y >> z))
		{
			std::cout << "Invalid input for (x y z)\n";
			goto D;
		}

		auto start = std::chrono::high_resolution_clock::now(); // Start timing here
		calculateE(charges, x, y, z, separationX, separationY, Ex_final, Ey_final, Ez_final, NUM_THREADS);
		auto end = std::chrono::high_resolution_clock::now(); //end timing here
		auto time_taken = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		std::cout << "The calculation took " << time_taken.count() << " microsec!\n";

		std::cout << "The electric field at (" << x << ", " << y << ", " << z << ") in V/m is\n";
		std::cout << "Ex = " << Ex_final << "\n";
		std::cout << "Ey = " << Ey_final << "\n";
		std::cout << "Ez = " << Ez_final << "\n";
		double E_total = std::sqrt(Ex_final * Ex_final + Ey_final * Ey_final + Ez_final * Ez_final);
		std::cout << "|E| = " << E_total << "\n";

		// Ask the user if they want to enter a new location
		std::cout << "\nDo you want to enter a new location (Y/N)? ";
		std::cin >> userInput;
		if (userInput == "Y" || userInput == "y")
			continue;
		else
			break;
		/*
		// Join threads
		for (std::thread& t : threads) {
			t.join();
		}
		*/
	}

	std::cout << "Bye!\n";

	return 0;
}