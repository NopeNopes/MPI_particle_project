#include "Particle.h"

double Force(Vector2d r1, Vector2d r2) {

	double c = 0.1, d = 2;
	double F = -c * ((r2 - r1).length() - d);

	return F;
}

void ForceCalculate(vector<Particle>& particles, int i_min, int i_max) {
	for (int i = i_min; i < i_max; i++) {
		particles[i].summaryForce.resetToZero();
		for (int j = 0; j < particles.size(); j++) {
			if (i != j) {
				Vector2d dif = Vector2d(particles[i].r - particles[j].r);
				double difNorm = dif.length();
				if (dif.length() < sqrt(4)) {
					double F = Force(Vector2d(particles[i].r), Vector2d(particles[j].r));
					particles[i].summaryForce = particles[i].summaryForce + F * (1 / difNorm) * dif;
				}
			}
		}
	}
}

void SpeedCalculate(vector<Particle>& particles, double dt, int i_min, int i_max) {
	for (int i = i_min; i < i_max; i++) {
		particles[i].v = particles[i].v + dt * particles[i].summaryForce / particles[i].m;
	}
}

void CoordinateCalculate(vector<Particle>& particles, double dt, int i_min, int i_max) {
	for (int i = i_min; i < i_max; i++) {
		particles[i].r = particles[i].r + dt * particles[i].v;
	}
}