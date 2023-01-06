#include <fstream>
#include <string>
#include <algorithm>
#include <chrono>
#include <ctime>

using namespace std;

int **Dijkstra(int **graphMatrix, int size, int destinationCount, bool *destinations)
{
    int **destinationDistances = (int**)malloc(sizeof(int*)*destinationCount);
    int *distances = (int*)malloc(sizeof(int)*size);
    bool *visited = (bool*)malloc(sizeof(bool)*size);

    for(int i = 0; i < destinationCount; i++)
        destinationDistances[i] = (int*)malloc(sizeof(int)*destinationCount);

    int currentDestination = -1;
    for(int i = 0; i < size; i++)
        //for each destination, run Dijkstra's algorithm to find shortest path to other destinations.
        if (destinations[i])
        {
            currentDestination++;

            for(int n = 0; n < size; n++) //Set up starting values.
            {
                distances[n] = graphMatrix[i][n];
                visited[n] = false;
            }

            distances[i] = 0;
            visited[i] = true;

            int currentNode = i;
            for(;;){
                    //Find next unvisited node with shortest distance.
                for(int n = 0; n < size; n++) if((!visited[n] && (visited[currentNode])) || (distances[n] != 0 && distances[n] > distances[currentNode]))
                        currentNode = n;

                if(visited[currentNode]) break;  //End the loop if no unvisited node was found.
                visited[currentNode] = true;

                    //Assign new distance values to its neighbours.
                for(int n = 0; n < size; n++)
                    if(graphMatrix[currentNode][n] > 0
                    && ((distances[n] == 0 && n != i) || (distances[n] > distances[currentNode]+graphMatrix[currentNode][n])))
                            distances[n] = distances[currentNode]+graphMatrix[currentNode][n];
            }

            //Get lengths between destinations and put them into the new graph.
            int n = currentDestination;
            for(int y = i+1; y < size; y++) if(destinations[y])
                {
                    n++;
                    destinationDistances[currentDestination][n] = distances[y];
                    destinationDistances[n][currentDestination] = destinationDistances[currentDestination][n];
                }
        }

    free(distances);
    free(visited);

    return destinationDistances;
}

void shufflePermutation(int *permutation, int size)
{
    for(int i = size - 1; i > 0; i--)
    {
        int toSwap = rand()%i;
        swap(permutation[toSwap], permutation[i]);
    }
}

int greedyAlgorithm(int **graphMatrix, int size)
{
    int *permutation = (int*)malloc(sizeof(int)*size);
    bool *visited = (bool*)malloc(sizeof(bool)*size);
    for(int i = 0; i < size; i++)
    {
        visited[i] = false;
        permutation[i] = i;
    }
    visited[0] = true;

    shufflePermutation(permutation, size);  //Shuffle the order of comparing vertices.

    int result = 0;

    int currentDestination = 0;
    int nextDestination;
    int closestDistance;

    for(int i = 1; i < size; i++)
    {
        nextDestination = 0;
        closestDistance = 0;
        for(int y = 0; y < size; y++)
            if(!visited[permutation[y]] && (graphMatrix[currentDestination][permutation[y]] < closestDistance || closestDistance == 0)) //Find closest unvisited neighbour to current node.
                nextDestination = permutation[y];

        result += graphMatrix[currentDestination][nextDestination];    //Otherwise, add nextDestination to the path and set it as next search point.
        visited[nextDestination] = true;
        currentDestination = nextDestination;
    }
    result += graphMatrix[currentDestination][0];  //Connect final point with first point to end the cycle.

    free(visited);
    free(permutation);
    return result;
}

int cycleLength(int **graphMatrix, int size, int *permutation)
{
    int result = 0;
    for(int i = size-1; i > 0; i--)
    {
        result += graphMatrix[permutation[i]][permutation[i-1]];
    }
    result += graphMatrix[permutation[0]][permutation[size-1]];
    return result;
}

int bruteForce(int **graphMatrix, int size)
{
    int *permutation = (int*)malloc(sizeof(int)*size);
    for(int i = 0 ; i < size; i++)
        permutation[i] = i;


    int result = cycleLength(graphMatrix, size, permutation);
    while(next_permutation(&permutation[1],&permutation[size]))
    {
        result = min(result, cycleLength(graphMatrix, size, permutation));
    }
    free(permutation);
    return result;
}

int swapDifference(int **graphMatrix, int size, int *permutation, int firstSegment, int secondSegment)
{
    int result = 0;
    if(secondSegment-firstSegment != 1)
    {
            //Erase connections in places where we are disconnecting the vertices.
            //Modulo operation is to prevent indices from exceeding table range.
        result -= graphMatrix[permutation[firstSegment]][permutation[(firstSegment+1)%size]];
        result -= graphMatrix[permutation[firstSegment]][permutation[firstSegment-1]];
        result -= graphMatrix[permutation[secondSegment]][permutation[(secondSegment+1)%size]];
        result -= graphMatrix[permutation[secondSegment]][permutation[secondSegment-1]];
            //Now attach vertices in their new locations.
        result += graphMatrix[permutation[secondSegment]][permutation[(firstSegment+1)%size]];
        result += graphMatrix[permutation[secondSegment]][permutation[firstSegment-1]];
        result += graphMatrix[permutation[firstSegment]][permutation[(secondSegment+1)%size]];
        result += graphMatrix[permutation[firstSegment]][permutation[secondSegment-1]];
    }
    else   //Edge case if two elements are next to each other.
    {
        result -= graphMatrix[permutation[firstSegment]][permutation[firstSegment-1]];
        result -= graphMatrix[permutation[secondSegment]][permutation[(secondSegment+1)%size]];
        result += graphMatrix[permutation[secondSegment]][permutation[firstSegment-1]];
        result += graphMatrix[permutation[firstSegment]][permutation[(secondSegment+1)%size]];
    }

    return result;
}

int localSearch(int **graphMatrix, int size)
{
    int *cycle = (int*)malloc(sizeof(int)*size);
    int *permutation = (int*)malloc(sizeof(int)*size);
    bool *visited = (bool*)malloc(sizeof(bool)*size);
    for(int i = 0; i < size; i++)
    {
        visited[i] = false;
        permutation[i] = i;
    }
    visited[0] = true;

    shufflePermutation(permutation, size);  //Shuffle the order of comparing vertices.
    cycle[0] = permutation[0];

    int currentDestination = 0;
    int nextDestination;
    int closestDistance;

        //Find starting solution using greedy heuristic.
    for(int i = 1; i < size; i++)
    {
        nextDestination = 0;
        closestDistance = 0;
        for(int y = 0; y < size; y++)
            if(!visited[permutation[y]] && (graphMatrix[currentDestination][permutation[y]] < closestDistance || closestDistance == 0)) //Find closest unvisited neighbour to current node.
                nextDestination = permutation[y];

        cycle[i] = nextDestination;
        visited[nextDestination] = true;
        currentDestination = nextDestination;
    }

    free(visited);
    free(permutation);

        //Search nearby solutions to see if there's an improvement.
    bool improved = true;
    while(improved == true)
    {
        improved = false;
        for(int i = 1; i < size; i++) for(int y = i+1; y < size; y++)
            if(swapDifference(graphMatrix, size, cycle, i, y) < 0)
            {
                swap(cycle[i], cycle[y]);
                improved = true;
            }
    }

    int result = cycleLength(graphMatrix, size, cycle);
    free(cycle);
    return result;
}

int main()
{
    int i, y, z, point, connection;

    int size;
    int destinationCount;

    int **graphMatrix;
    bool *destinations;

    string line;
    ifstream inputf;

    inputf.open("macierz.txt");
    if (!inputf)
        return 0;

        //Open the file and get size data.
    getline(inputf,line);
    size = line.length();
    destinationCount = 0;

        //Malloc tables.
    graphMatrix = (int**)malloc(sizeof(int*)*size);
    destinations = (bool*)malloc(sizeof(bool)*size);

    for(i = 0; i < size; i++)
        for(y = 0; y < size; y++)
        {
            graphMatrix[y] = (int*)malloc(sizeof(int)*size);
            for(z = 0; z < size; z++)
                graphMatrix[y][z] = 0;
        }

    //Get destination count
    for(i = 0; i < line.length(); i++)
    {
        if(line[i] == '1')
        {
            destinations[i] = true;
            destinationCount++;
        }else
            destinations[i] = false;
    }

    //Load graph
    for(point = 0; point < size; point++)
    {
        getline(inputf,line);
        connection = 0;
        for(i = 0; i < line.length(); i++)
        {
            if(line[i] != ' ')
            {
                graphMatrix[point][connection] *= 10;
                graphMatrix[point][connection] += (line[i] - '0');
            }else
            {
                graphMatrix[connection][point] = graphMatrix[point][connection];
                connection++;
            }
        }
    }
    inputf.close();

    int **destinationDistances = Dijkstra(graphMatrix, size, destinationCount, destinations);

    //From now on we'll be only working on the graph of destinations, so the original graph can be erased from memory.
    for(i = 0; i < size; i++)
        free(graphMatrix[i]);
    free(graphMatrix);
    free(destinations);

    //TSP
    ofstream outputf;
    outputf.open("wynik.txt");
    if(!outputf)
        return 0;

    int result;
    std::chrono::steady_clock::time_point startTime;
    std::chrono::steady_clock::time_point endTime;

    srand(time(NULL));  //Setting up the seed for RNG.

    startTime = std::chrono::steady_clock::now();
    result = bruteForce(destinationDistances, destinationCount);
    endTime = std::chrono::steady_clock::now();
    outputf << "Brute force: " << result << " time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << endl;

    startTime = std::chrono::steady_clock::now();
    result = greedyAlgorithm(destinationDistances, destinationCount);
    endTime = std::chrono::steady_clock::now();
    outputf << "Greedy: " << result << " time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << endl;

    startTime = std::chrono::steady_clock::now();
    result = localSearch(destinationDistances, destinationCount);
    endTime = std::chrono::steady_clock::now();
    outputf << "Local search: " << result << " time taken: " << std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() << endl;

    outputf.close();

    //Release graph
    for(i = 0; i < destinationCount; i++)
        free(destinationDistances[i]);
    free(destinationDistances);

    return 0;
}
