# GeneticAlgorithmParallelComputing
Implementing multithreading using openMP library on Genetic algorithm
Aparajita Sahay
Programming 1: Genetic Algorithm using Open MP
The shortest trip in my program 442.563 – 449.54.
Performance improvement using 4 thread = 2.57 times.
Mutation rate = 25.
Parallel programming when applied only on:-
With increase in mutation rate shortest distance value increased , therefore used mutation rate as 25%.
Evaluate function: 
1)	Approach: Used omp for on outer loop of the function to split up loop iterations among the threads.
Evaluate function	Thread 1	Thread 4	Shortest distance 
Outer for loop : used omp parallel for 	16630556	13988424	442.563
			
 
2)	I also tried omp parallel for on inner loop to compute the city but the computation time increased and the shortest distance remain the same. Used distance variable as private, so that each thread can compute its own distance and declared totalDistanceTrip variable as reduction so that every thread uses totalDistanceTrip variable to add it’s computed value at the end. 
Evaluate function 	Thread 1	Thread 4 	Shortest distance
Inner loop : used omp parallel for	18988637	14492842	442.563

3)	Used omp for to split up loop iterations among the threads. I also used static scheduling clause to divide for loop iteration among all the thread. Reason: As each iteration will take equal time to compute the distance of a trip, because every trip includes equal cities. Therefore, by using static scheduling we can divide 50000 iteration into sections and assign it to threads. Because of static scheduling distribution of work among the thread will be easy. Iteration will be assigned to the thread before the loop is executed. However, in dynamic iteration will be assigned to the thread at run time which will also have some extra overhead. But as we know each iteration takes equal time as explained above we may choose static schedule clause.
Evaluate function 	Thread 1	Thread 4 	Shortest distance
Outter for loop : used omp parallel for with static schedule clause 	14277701	16946633	442.563
Outer for loop : used omp parallel for with dynamic schedule clause	16836734	14131155	442.563

So, by using both the scheduling, computation time was not drastically different but static scheduling takes minutely less time than dynamic scheduling. Therefore, by using parallel computing only on outer for loop with default scope for the variable got the maximum performance and accuracy for the shortest distance also remained as per expected requirement.

Crossover function
1)	Approach: Used omp for on outer loop of the function to split up loop iterations among the threads.
Also used critical section at location where random number is generated. The scope for the variable used to store random number i.e.randomIndex was declared private at the outer for loop because each thread can have their private random number to generate an index. To make sure that rand function is called by one thread at a time and each thread generate some new output, I used critical section. I also make sure minimum to use minimum lines of code inside the critical section to reduce serialization of the program.
The performance and accuracy of the shortest distance changes after running the program. Performance remain close to 2.57% but the shortest distance varies from 448.57 to 465.
Since performance is better than applying parallel computing only evaluate function, I implemented parallelization of crossover function. But the tradeoff for a better performance was the accuracy of the shortest distance in my program.
Parallel computing Crossover and evaluation and on population function 	Performance improvement 	Shortest distance 
	2.57%	448.57-465
	
 


	 



 
