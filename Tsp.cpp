#include <iostream>  // cout
#include <fstream>   // ifstream
#include <string.h>  // strncpy
#include <stdlib.h>  // rand
#include <math.h>    // sqrt, pow
#include <omp.h>     // OpenMP
#include "Timer.h"
#include "Trip.h"
#include <algorithm>

using namespace std;

// Already implemented. see the actual implementations below
void initialize( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] );
void select( Trip trip[CHROMOSOMES], Trip parents[TOP_X] );
void populate( Trip trip[CHROMOSOMES], Trip offsprings[TOP_X] );

// need to implement for your program 1
extern void evaluate( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] );
extern void crossover( Trip parents[TOP_X], Trip offsprings[TOP_X], int coordinates[CITIES][2] );
extern void mutate( Trip offsprings[TOP_X] );
//extern void complement1(CITIES);
char complement[CITIES + 1] = "9876543210ZYXWVUTSRQPONMLKJIHGFEDCBA";

/*
 * MAIN: usage: Tsp #threads
 */
int main( int argc, char* argv[] ) {
  Trip* trip = new Trip[CHROMOSOMES];       // all 50000 different trips (or chromosomes)
  Trip shortest;                // the shortest path so far
  int coordinates[CITIES][2];   // (x, y) coordinates of all 36 cities:
  int nThreads = 1;

//double distance = 0;
//unsigned int seedptr=1;  
  // verify the arguments
  if ( argc == 2 )
    nThreads = atoi( argv[1] );
  else {
    cout << "usage: Tsp #threads" << endl;
    if ( argc != 1 )
      return -1; // wrong arguments
  }
  cout << "# threads = " << nThreads << endl;

  // shortest path not yet initialized
  shortest.itinerary[CITIES] = 0;  // null path
  shortest.fitness = -1.0;         // invalid distance

  //complement1(CITIES);
  // initialize 5000 trips and 36 cities' coordinates
  initialize( trip, coordinates );


  // start a timer 
  Timer timer;
  timer.start( );

  // change # of threads
  omp_set_num_threads( nThreads );

  // find the shortest path in each generation

  for ( int generation = 0; generation < MAX_GENERATION; generation++ ) {

    // evaluate the distance of all 50000 trips
   evaluate( trip, coordinates );

   //call complement function to create complement array for 36 cities...

   // complement(coordinates);
    // just print out the progress
    if ( generation % 20 == 0 )
      cout << "generation: " << generation << endl;

    // whenever a shorter path was found, update the shortest path
    if ( shortest.fitness < 0 || shortest.fitness > trip[0].fitness ) {

      strncpy( shortest.itinerary, trip[0].itinerary, CITIES );
      shortest.fitness = trip[0].fitness;

      cout << "generation: " << generation 
	   << " shortest distance = " << shortest.fitness
	   << "\t itinerary = " << shortest.itinerary << endl;
    }

    // define TOP_X parents and offsprings.
    Trip* parents;
    Trip* offsprings;
    parents = new Trip[TOP_X];
    offsprings = new Trip[TOP_X];

    // choose TOP_X parents from trip
    select( trip, parents );

    // generates TOP_X offsprings from TOP_X parents
    crossover( parents, offsprings, coordinates );

    // mutate offsprings

    mutate(offsprings);
    // populate the next generation.
    populate( trip, offsprings );


  }

  // stop a timer
  cout << "elapsed time = " << timer.lap( ) << endl;
  return 0;
}



/*
 * Initializes trip[CHROMOSOMES] with chromosome.txt and coordiantes[CITIES][2] with cities.txt
 *
 * @param trip[CHROMOSOMES]:      50000 different trips
 * @param coordinates[CITIES][2]: (x, y) coordinates of 36 different cities: ABCDEFGHIJKLMNOPQRSTUVWXYZ
 */
void initialize( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] ) {
  // open two files to read chromosomes (i.e., trips)  and cities
  ifstream chromosome_file( "chromosome.txt" );
  ifstream cities_file( "cities.txt" );
  
  // read data from the files
  // chromosome.txt:                                                                                           
  //   T8JHFKM7BO5XWYSQ29IP04DL6NU3ERVA1CZG                                                                    
  //   FWLXU2DRSAQEVYOBCPNI608194ZHJM73GK5T                                                                    
  //   HU93YL0MWAQFIZGNJCRV12TO75BPE84S6KXD
  for ( int i = 0; i < CHROMOSOMES; i++ ) {
    chromosome_file >> trip[i].itinerary;
    trip[i].fitness = 0.0;
  }

  // cities.txt:                                                                                               
  // name    x       y                                                                                         
  // A       83      99                                                                                        
  // B       77      35                                                                                        
  // C       14      64                                                                                        
  for ( int i = 0; i < CITIES; i++ ) {
    char city;
    cities_file >> city;
    int index = ( city >= 'A' ) ? city - 'A' : city - '0' + 26;
    cities_file >> coordinates[index][0] >> coordinates[index][1];
  }

  // close the files.
  chromosome_file.close( );
  cities_file.close( );

  // just for debugging
  if ( DEBUG ) {
    for ( int i = 0; i < CHROMOSOMES; i++ )
      cout << trip[i].itinerary << endl;
    for ( int i = 0; i < CITIES; i++ )
      cout << coordinates[i][0] << "\t" << coordinates[i][1] << endl;
  }
}

//define new complement array to use in crossover

/*
 * Select the first TOP_X parents from trip[CHROMOSOMES]
 *
 * @param trip[CHROMOSOMES]: all trips
 * @param parents[TOP_X]:    the firt TOP_X parents
 */
void select( Trip trip[CHROMOSOMES], Trip parents[TOP_X] ) {
  // just copy TOP_X trips to parents
//#pragma omp parallel for
  for ( int i = 0; i < TOP_X; i++ )
    strncpy( parents[i].itinerary, trip[i].itinerary, CITIES + 1 );
}

/*
 * Replace the bottom TOP_X trips with the TOP_X offsprings
 */
void populate( Trip trip[CHROMOSOMES], Trip offsprings[TOP_X] ) {
  // just copy TOP_X offsprings to the bottom TOP_X trips.
#pragma omp parallel for 

  for ( int i = 0; i < TOP_X; i++ )
    strncpy( trip[ CHROMOSOMES - TOP_X + i ].itinerary, offsprings[i].itinerary, CITIES + 1 );

  if ( DEBUG ) {
    for ( int chrom = 0; chrom < CHROMOSOMES; chrom++ ) 
      cout << "chrom[" << chrom << "] = " << trip[chrom].itinerary 
	   << ", trip distance = " << trip[chrom].fitness << endl;
  }

}

int compare(const void * a ,const void * b)
{
	return (((Trip*)a)->fitness - ((Trip*)b)->fitness);
}
void evaluate(Trip trip[CHROMOSOMES],int coordinates[CITIES][2])
{
//	int noOfThread = omp_get_thread_num);
//	#pragma omp parallel  for
//	double totalDistanceTrip=0; 
	#pragma omp parallel for 
	for (int i = 0; i < CHROMOSOMES; i++)
	{
		double distance;
		double totalDistanceTrip = 0;
		#pragma omp parallel for private(distance) reduction(+:totalDistanceTrip) 
		for (int j = 0; j < CITIES - 1; j++)
		{
			char city1 = trip[i].itinerary[j];
			char city2 = trip[i].itinerary[j+1];

			int city1Index = (city1 >= 'A') ? city1 - 'A' : city1 - '0' + 26;
			int city2Index = (city2 >= 'A') ? city2 - 'A' : city2 - '0' + 26;
			int coordinate_x1 = coordinates[city1Index][0];
			int coordinate_y1 = coordinates[city1Index][1];
			int coordinate_x2 = coordinates[city2Index][0];
			int coordinate_y2 = coordinates[city2Index][1];

			int distanceX = (coordinate_x2 - coordinate_x1)*(coordinate_x2 - coordinate_x1);
			int distanceY = (coordinate_y2 - coordinate_y1)*(coordinate_y2 - coordinate_y1);
			distance = sqrt(distanceX + distanceY);
			totalDistanceTrip += distance;

		}
			//cout << "\n";
			//cout << i << "distance calculate"<< totalDistanceTrip;
			trip[i].fitness = totalDistanceTrip;
	}


	qsort(trip, CHROMOSOMES,sizeof(Trip),compare);
	//for (int n=0; n<CHROMOSOMES; n++)
	//cout<< "\n"<< "sorting" << trip[n].fitness;


}

char getCityChar(int index)
{
	return index < 26? 'A' + index : index -26 + '0';
}

int getCityIndex(char city)
{
	return (city >= 'A') ? city - 'A' : city - '0' + 26;
}

bool cityExists(char city, char trip[CITIES])
{
	bool exists = false;
	//	#pragma omp parallel for shared(exists) 
	for(int i=0;i<CITIES;i++)
	{
		if(trip[i] == city)
		{
			exists = true;
			break;
		}
	}

	return exists;
}

void crossover(Trip parents[TOP_X], Trip offsprings[TOP_X], int coordinates[CITIES][2])
{
//	int offspringIndex=0;
//	unsigned int seed = time(NULL);
	unsigned int seed = 0;
	int randomIndex;
	#pragma omp parallel for private(randomIndex)
	for (int i = 0; i < TOP_X; i+=2)
	{
		/*Trip parentA = parents[i];
		Trip parentB = parents[i + 1];
		Trip offspringA = offsprings[i];
		Trip offspringB = offsprings[i + 1]; */

		int offspringCounter = 0;
//		offspringIndex=offspringIndex+2;
		for(int j=0;j<CITIES;j++)
		{
			char nextTripNeighbourCity = 0;
			char tripCityA = parents[i].itinerary[j];
			char TripNeighbourCity = parents[i].itinerary[j+1];

			if(cityExists(tripCityA,offsprings[i].itinerary))
			{
				continue;
			}

			offspringCounter++;
			offsprings[i].itinerary[offspringCounter-1]=tripCityA;

			//look for a source city in cityParent1 at trip(i+1)
			for (int k=0;k<CITIES;k++)
			{
				if(parents[i+1].itinerary[k] == tripCityA && k < CITIES - 1)
				{
					nextTripNeighbourCity = parents[i+1].itinerary[k+1];
					break;
				}
			}

			/*calculate distance between tripCityA and TripNeighbourCity in trip 1 (i) and calculate
			distance between tripCityA and nextTripNeighbourCity */
			int city1Index = (tripCityA >= 'A') ? tripCityA - 'A' : tripCityA - '0' + 26;
			int city2Index = (TripNeighbourCity >= 'A') ? TripNeighbourCity - 'A' : TripNeighbourCity - '0' + 26;
			int coordinate_x1 = coordinates[city1Index][0];
			int coordinate_y1 = coordinates[city1Index][1];
			int coordinate_x2 = coordinates[city2Index][0];
			int coordinate_y2 = coordinates[city2Index][1];
			int distanceX = (coordinate_x2 - coordinate_x1)*(coordinate_x2 - coordinate_x1);
			int distanceY = (coordinate_y2 - coordinate_y1)*(coordinate_y2 - coordinate_y1);
			double tripADistance = sqrt(distanceX + distanceY);

			int tripBCityIndex = (nextTripNeighbourCity >= 'A') ? nextTripNeighbourCity - 'A' : nextTripNeighbourCity - '0' + 26;
			int tripBCoordinate_x2 = coordinates[tripBCityIndex][0];
			int tripBCoordinate_y2 = coordinates[tripBCityIndex][1];
			int tripBDistanceX = (tripBCoordinate_x2 - coordinate_x1)*(tripBCoordinate_x2 - coordinate_x1);
			int tripBDistanceY = (tripBCoordinate_y2 - coordinate_y1)*(tripBCoordinate_y2 - coordinate_y1);
			double tripBDistance = sqrt(tripBDistanceX+ tripBDistanceY);

			if(tripADistance < tripBDistance &&
				!cityExists(TripNeighbourCity,offsprings[i].itinerary))
			{
				offspringCounter++;
				offsprings[i].itinerary[offspringCounter-1]=TripNeighbourCity ;
			}
			else if(nextTripNeighbourCity != 0 &&
					!cityExists(nextTripNeighbourCity, offsprings[i].itinerary))
			{
				offspringCounter++;
				offsprings[i].itinerary[offspringCounter-1]=nextTripNeighbourCity;
			}
			else
			{
				// unsigned int seed = time(NULL);

			//	#pragma omp critical
			//	{
				while(true)
				{
			//	#pragma omp critical
				//	{
					//unsigned int seed = time(NULL);
				//	int randomIndex;
					#pragma omp critical  
					{	
					 randomIndex = rand_r(&seed)% CITIES;
					//	int  randomIndex = rand()% CITIES;
					}
					char randomCity = getCityChar(randomIndex);
					if(!cityExists(randomCity, offsprings[i].itinerary))
					{
						offspringCounter++;
						offsprings[i].itinerary[offspringCounter-1] = randomCity;
						break;
					}
				}
			//	}
			}
		}

		//complement of offspring[i]'s element. Using  child[i] to create chile[i+1]
		//creating child[i+1]
		for(int z =0;z<CITIES;z++)
		{
			char city = offsprings[i].itinerary[z];
			int cityIndex = getCityIndex(city);

			char complementCity= complement[cityIndex];
			offsprings[i+1].itinerary[z] = complementCity;
		}
	}
}



void mutate(Trip offsprings[TOP_X])
{
//	unsigned int seedptr= 1;

//	#pragma omp parallel for 
	for(int i = 0;i<TOP_X;i++)
	{
//	srand (time(NULL));
//	unsigned int seedptr=time(NULL);
	//	if( rand_r(&seedptr)%100 < MUTATE_RATE)
	//	{
	//		int c1 = rand_r(&seedptr)%36;
	//		int c2 = rand_r(&seedptr)%36;


		if( rand()%100 < MUTATE_RATE)
                {
                        int c1 = rand()%36;
                        int c2 = rand()%36;

			//swap two elements at index c1 and c2
			char temp1 =offsprings[i].itinerary[c1];
			offsprings[i].itinerary[c1] = offsprings[i].itinerary[c2];;
			offsprings[i].itinerary[c2] = temp1;
		}
	}

}




