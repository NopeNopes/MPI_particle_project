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

    //tags for mpi
    int root = 0;
    int send_initial_size_tag = 0;
    int send_particle_initial_tag = 1;
    int send_particle_size_calc_tag = 2;
    int send_particle_calc_tag = 3;

    // time
    int time = 80;
    double dt = 0.01; // time step
    double startwtime = 0.0;
    double endwtime;

    // mesh
    double len_x = 40;
    double len_y = 40;

    //number of particles
    int num_particle_x = 20;
    int num_particle_y = 20;
    int num_impactor = 5;
    int size_mesh = 0;


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
        part2.r = Vector2d(len_x * 1.5, len_y / 2);
        part2.v = Vector2d(-1, 0);
        particles_mesh.push_back(part2);

        for (int i = 1; i < num_impactor; i++) {
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
    
    // lower and upper bounds of number of particles assigned to current processor (depending on world_rank)
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

        // calculating
        ForceCalculate(particles_mesh, i_min, i_max);
        SpeedCalculate(particles_mesh, dt, i_min, i_max);
        CoordinateCalculate(particles_mesh, dt, i_min, i_max);

        // define sending data and result (recieving) data
        vector<Particle> particles_send;
        vector<Particle> particles_result;

        // preparing sending data
        for (int i = i_min; i < i_max; i++) {
            particles_send.push_back(particles_mesh[i]);
        }

        // send calculated data to other processors
        for (int world_rank_i = 0; world_rank_i < world_size; world_rank_i++) {
            if (world_rank_i != world_rank) {
                int size_send = particles_send.size();
                MPI_Send(&size_send, 1, MPI_INT, world_rank_i, send_particle_size_calc_tag, MPI_COMM_WORLD);
                MPI_Send(&particles_send[0], size_send * 7, MPI_DOUBLE, world_rank_i, send_particle_calc_tag, MPI_COMM_WORLD);
            }
        }

        // recieve result from other processors
        for (int world_rank_i = 0; world_rank_i < world_size; world_rank_i++) {
            vector<Particle> particles_recv;
            int size_recv;

            if (world_rank_i != world_rank) {
                // if recieving from other processors:
                MPI_Recv(&size_recv, 1, MPI_INT, world_rank_i, send_particle_size_calc_tag, MPI_COMM_WORLD, &status);
                particles_recv.resize(size_recv);
                MPI_Recv(&particles_recv[0], size_recv * 7, MPI_DOUBLE, world_rank_i, send_particle_calc_tag, MPI_COMM_WORLD, &status);
                
                //transfering data to result vector
                for (int j = 0; j < particles_recv.size(); j++) {
                    particles_result.push_back(particles_recv[j]);
                }
            }
            else
            {   
                // if NOT recieving from other processors:
                for (int i = i_min; i < i_max; i++) {
                    particles_result.push_back(particles_mesh[i]);
                }
            }
        }


        // clearing and override/redefine mesh for the next time step
        particles_mesh.clear();
        for (int i = 0; i < particles_result.size(); i++) {
            particles_mesh.push_back(particles_result[i]);
        }

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

