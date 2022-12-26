#pragma once
#include "Vector2d.h";
#include <vector>;
#include <math.h>;

using namespace std;

class Particle
{
public:

	double m = 1;
	Vector2d summaryForce;
	Vector2d r, v;

	Particle() {
		r.resetToZero();
		v.resetToZero();
		summaryForce.resetToZero();
	}

	Particle(Vector2d r1, Vector2d v1) {
		r = r1;
		v = v1;
	}
	
};

double Force(Vector2d r1, Vector2d r2);
void ForceCalculate(vector<Particle>& particles, int i_min, int i_max);
void SpeedCalculate(vector<Particle>& particles, double dt, int i_min, int i_max);
void CoordinateCalculate(vector<Particle>& particles, double dt, int i_min, int i_max);
