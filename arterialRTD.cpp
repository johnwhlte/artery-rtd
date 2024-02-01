/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Louis Kronberg, Christian Vorwerk, Bastian Schäffauer
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* rtdVal.cpp
Validation case for discrete particle RTD in OpenLB validated against
a previously run RTD in Ansys Fluent.
*/
#include "simfuncs/sim_funcs.h"

int main( int argc, char* argv[] )
{
    // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
     N,            // Number of lattice cells across the char length defined below
    (T) tau,          // relaxation time tau
    (T) charMinL,     // charPhysLength: reference length of simulation geometry
    (T) 2.*avgVel,    // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T) nuP,          // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T) rhoP          // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write(vtkFileName);

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  char stlDir[] = "geometry/";
  strcat(stlDir,stlFileName);
  STLreader<T> stlReader( stlDir, converter.getConversionFactorLength(), 0.001, 0, true );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getConversionFactorLength() );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = util::min(16*N, 8*singleton::mpi().getSize());
#else
  const int noOfCuboids = 2;
#endif
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(), noOfCuboids, "volume" );
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  Randomizer<T> randomizer;

  const unsigned latticeMaterial = 2; //Material number of wall
  const unsigned contactMaterial = 0; //Material identifier (only relevant for contact model)
  SolidBoundary<T,3> wall( std::make_unique<IndicInverse<T,3>>(stlReader),
                           latticeMaterial, contactMaterial );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  util::Timer<T> timer1( converter.getLatticeTime( fluidMaxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer1.start();

  prepareLattice( sLattice, converter, stlReader, superGeometry );

  SuperParticleSystem<T,PARTICLETYPE> superParticleSystem(superGeometry);

  std::vector<int> materials {4};
  SuperIndicatorMaterial<T,3> materialIndicator (superGeometry, materials);

  ParticleManager<T,DESCRIPTOR,PARTICLETYPE> particleManager(
    superParticleSystem, superGeometry, sLattice, converter);

  prepareParticles( converter, superParticleSystem, wall, materialIndicator,
    particleDynamicsSetup, randomizer);

  initializeParticleVelocity( sLattice, superGeometry, converter, superParticleSystem );

  timer1.stop();
  timer1.printSummary();

  // // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( fluidMaxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT <= converter.getLatticeTime( fluidMaxPhysT ); iT++ ) {
  //   // === 5th Step: Definition of Initial and Boundary Conditions ===
     setBoundaryValues( sLattice, converter, iT, superGeometry );

  //   // === 6th Step: Collide and Stream Execution ===
     sLattice.collideAndStream();

  //   // === 7th Step: Computation and Output of the Results ===
     getResults( sLattice, converter, iT, superGeometry, timer, stlReader, superParticleSystem);
   }

  for ( std::size_t iT = converter.getLatticeTime(fluidMaxPhysT); iT <= converter.getLatticeTime( fluidMaxPhysT + particleMaxPhysT); ++iT ) {
    // particles simulation starts after run up time is over
    particleManager.execute<
      couple_lattice_to_particles<T,DESCRIPTOR,PARTICLETYPE>,
      process_dynamics<T,PARTICLETYPE>,
      update_particle_core_distribution<T, PARTICLETYPE>
        >();

    if ( !getResults( sLattice, converter, iT,
                      superGeometry, timer,
                      stlReader, superParticleSystem )){
      break;
    }
  }

   timer.stop();
   timer.printSummary();


    return 0;
}
