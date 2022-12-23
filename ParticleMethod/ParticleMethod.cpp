#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include "Particle.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;


int main(int argc, char* argv[])
{   
    Particle part1, part2;
    vector<Particle> particles_mesh;
    int world_rank, world_size;

    //mpi
    int root = 0;

    // time
    int time = 80;
    double dt = 0.01; // time step
    double startwtime = 0.0;
    double endwtime;

    //number of particles
    int num_particle_x = 20;
    int num_particle_y = 20;
    int num_impactor = 4;
    int size_mesh = 0;

    // mesh
    double d = 2;
    double len_x = num_particle_x * d;
    double len_y = num_particle_y * d;

    MPI_Status status;
    // initialize MPI
    MPI_Init(&argc, &argv);
    // current world_rank
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // number of processors
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    if (world_rank == root) {
        // main processor, world_rank=0
        // generating mesh
        startwtime = MPI_Wtime();
        for (int i = 0; i < num_particle_x; i++) {
            double coord_x = i * len_x / num_particle_x;
            for (int j = 0; j < num_particle_y; j++) {
                if (i % 2 == 0) {
                    double coord_y = j * len_y / num_particle_y;
                    part1.r = Vector2d(coord_x, coord_y);
                    particles_mesh.push_back(part1);
                }
                else {
                    double coordY = j * len_y / num_particle_y + len_y / num_particle_y / 2;
                    part1.r = Vector2d(coord_x, coordY);
                    particles_mesh.push_back(part1);
                }
            }
        }

        // generating impactor
        for (int i = 1; i < num_impactor + 1; i++) {
            double coordX = i * len_x / num_particle_x;
            double coordY = len_y / num_particle_y;

            part2.r = Vector2d(len_x * 1.5 + coordX, coordY + len_y / 2);
            part2.v = Vector2d(-1, 0);
            particles_mesh.push_back(part2);
            
            part2.r = Vector2d(len_x * 1.5 + coordX, len_y / 2 - coordY);
            part2.v = Vector2d(-1, 0);
            particles_mesh.push_back(part2);
            
            part2.r = Vector2d(len_x * 1.5 + coordX, len_y / 2);
            part2.v = Vector2d(-1, 0);
            particles_mesh.push_back(part2);
        }

        // total number of particles, including impactor and main mesh of particles
        size_mesh = particles_mesh.size();
    }

    // broadcasting mesh to other processors
    MPI_Bcast(&size_mesh, 1, MPI_INT, root, MPI_COMM_WORLD);
    particles_mesh.resize(size_mesh);
    MPI_Bcast(&(particles_mesh[0]), size_mesh * 7, MPI_DOUBLE, root, MPI_COMM_WORLD);


    // CALCULATING
    
    // lower and upper bounds of number of particles assigned to a current processor (depending on world_rank)
    int i_min = world_rank * particles_mesh.size() / world_size;
    int i_max = (world_rank + 1) * particles_mesh.size() / world_size;
    int count = 0;

    string line;
    ofstream out_x, out_y;
    out_x.open("D:\\Uni\\11 semester\\c++\\ParticleMethod\\MPI_particle_project\\output\\particles_list_x.txt"); 
    out_y.open("D:\\Uni\\11 semester\\c++\\ParticleMethod\\MPI_particle_project\\output\\particles_list_y.txt"); 
    out_x << endl;
    out_y << endl;

    for (double t = 0; t <= time; t += dt) {
        count++;

        // define sending data and result (recieving) data
        vector<Particle> particles_result;
        vector<Particle> particles_send;

        for (int i = 0; i < particles_mesh.size(); i++) {
            particles_result.push_back(particles_mesh[i]);
        }

        // calculating
        ForceCalculate(particles_mesh, i_min, i_max);
        SpeedCalculate(particles_mesh, dt, i_min, i_max);
        CoordinateCalculate(particles_mesh, dt, i_min, i_max);

        // preparing sending data
        for (int i = i_min; i < i_max; i++) {
            particles_send.push_back(particles_mesh[i]);
        }

        int size_send = particles_send.size();
        int send_count_all = particles_result.size();
        int send_count_per_process = size_send * 7;
        particles_send.resize(size_send);
        particles_result.resize(send_count_all);

        // collecting every calculated part from all processors
        MPI_Gather(&particles_send[0], send_count_per_process, MPI_DOUBLE, &particles_result[0], send_count_per_process, MPI_DOUBLE, root, MPI_COMM_WORLD);

        // clearing and override/redefine mesh for the next time step
        particles_mesh.clear();
        if (world_rank == 0) {
            for (int i = 0; i < particles_result.size(); i++) {
                particles_mesh.push_back(particles_result[i]);
            }
            int size_mesh = particles_mesh.size();
        }

        // broadcasting mesh to other processors
        MPI_Bcast(&size_mesh, 1, MPI_INT, root, MPI_COMM_WORLD);
        particles_mesh.resize(size_mesh);
        MPI_Bcast(&(particles_mesh[0]), size_mesh * 7, MPI_DOUBLE, root, MPI_COMM_WORLD);

        // writing to a file 
        if (world_rank == root) {

            if (count % 100 == 0) {

                out_x << dt * count << ",";
                out_y << dt * count << ",";

                for (int i = 0; i < particles_mesh.size(); i++) {
                    out_x << particles_mesh[i].r.x << ",";
                    out_y << particles_mesh[i].r.y << ",";
                }

                out_x << endl;
                out_y << endl;
            }
        }

    }// end of calculating and time steps

    out_x.close();
    out_y.close(); 

    // resulting time
    if (world_rank == root) {
        endwtime = MPI_Wtime();
        cout << "time= " << (endwtime - startwtime) << endl;
    }
    
    // end working with MPI
    MPI_Finalize();
    return 0;
}

