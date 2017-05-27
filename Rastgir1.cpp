/*
Munem Rastgir
Greedy TSP Algorithm

Main uses functions:
minDist- Creates n(n-1)/2 edges out of n countries from our file and puts them all into a vector of edges
Heap Sort: uses a nlogn sort to sort the vector of edges from least to greatest distance
minPath- Adds edges to the Adjancency List only if both vertices are size < 2 and there is no cycle formed if we add that edge
->Min Path uses functions: listSize, addEdge, isCycle, and pathReader
	listSize traverses a list of a given position until it hits NULL and returns the size of the list
	addEdge adds the 2 numbers of the 2 cities of the edge we are looking at into the first spot of the adjacency list
	isCycle is performed on a given edge to see if adding this edge causes a cycle. We start off with the 2nd city of this edge and traverse 
	the adjancency list until we hit NULL or reach the first city of the vertex.
	pathReader traverses the list starting from arr[1] and when we read the last vertex: we attach a new edge from that vertex to this first one
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <fstream>
#include<algorithm>
#include<string>
#include<chrono>
using namespace std;

//A country holds the latitude and longitude and position data from the input file
struct Country {
	int position;
	double latitude;
	double longitude;
};

//An edge holds data:
struct edge {
	int start;			//The starting vertex
	int finish;			//The ending vertex
	double distance;	//distance between 2 vertices
};

//adjListNode has it's position(number of what node it is) and pointer to next node
struct adjListNode {
	int position;
	double distance;
	adjListNode* next;
};

//adjList is just the head pointer to a node/list
struct adjList
{
	adjListNode* head;
};

double eucD(double x, double y, double a, double b)	//eucD takes 4 arguments and returns the euclidian distance of them
{
	double f = (x - a)*(x - a);	//x-a squared
	double s = (y - b)*(y - b);	//y-b squared
	return sqrt(f + s);			//sqrt of sum
}

void minDist(vector<edge>& eV, Country arr[], int size)	//minDist takes an array of countries and makes n(n-1)/2 edges out of them
{
	for (int i = 1; i < size - 1; i++)	//loops from first country thru 2nd to last country
	{
		for (int j = i + 1; j < size; j++)	//loops from 2nd country to last country
		{
			edge tmp;	//creates a new edge
			tmp.distance = eucD(arr[i].latitude, arr[i].longitude, arr[j].latitude, arr[j].longitude);	//calculates edge distance
			tmp.start = arr[i].position;	//start vertex number
			tmp.finish = arr[j].position;	//end vertex number
			eV.push_back(tmp);	//add the edge to the vector
		}
	}
}

//makeNode function creates a new node with argument pos which is it's position
adjListNode* makeNode(int pos)
{
	adjListNode* x = new adjListNode;
	x->position = pos;
	x->next = NULL;		//points to null
	return x;	//returns pointer to finished node
}


//Function calculates the size of a list (This is to make sure it's not greater than 2)
int listSize(adjList* x, int pos)	//takes an adj list and integer inputs
{
	adjListNode* tmp = x[pos].head;	//temp node points to the start of the desired list
	int i = 0;
	while (tmp != NULL)	//while you havent reached the end of the list
	{
		tmp = tmp->next;	//increment i by 1 everytime you go to the next item of the list
		i++;
	}

	delete tmp;	//delete the temp pointer
	return i;	//return size
}

//function adds edge to adjacency list. Since this is a undirected graph we have to add each vertex to both lists
void addEdge(adjList* x, int s, int f, double dist) //takes input adjList and start and finish vertices
{
	adjListNode* newnode = makeNode(s);	//adds s to f list
	newnode->distance = dist;
	newnode->next = x[f].head;
	
	x[f].head = newnode;

	newnode = makeNode(f);		//adds f to s list
	newnode->next = x[s].head;
	newnode->distance = dist;
	x[s].head = newnode;

}

//function to check if theres a cycle if we add the edge
bool isCycle(adjList* x, edge y)	//takes input adjList and edge we want to add
{
	int prev = y.finish;	//prev marks the previous vertex travelled so we dont infinitly travel through 2 nodes
	int tmp = prev;			
	int cycStart = y.start;	//cycstart marks the first node of the edge
	adjListNode* path = x[y.finish].head;	//Initialize Node* path that points to the 2nd node of the edge
	if (path == NULL)	//if the pointer is equal to null there is no cycle
	{
		return false;	//return false for no cycle
	}
	else {						//ELSE: We traverse the adjancency list from the 2nd Node of the edge until we hit NULL or cycStart
		while (path != NULL)	//loop until we reach null
		{
			if (path->position == cycStart)	//if the pointer points to cycStart: there is a cycle
			{
				return true;		//return true for cycle
			}

			if (path->position == prev) //if the pointer points to the most recent node previously visited
			{
				prev = tmp;			//update the prev value to current value
				path = path->next;	//go to the next node in the list. This is to avoid infinitely going through 2 nodes

			}
			if (path != NULL)	//If the pointer doesnt point to null
			{
				if (path->position == cycStart)	//check if it points to cycstart
				{
					return true;
				}
				else {
					prev = tmp;				//update the previous value
					tmp = path->position;	//update tmp to current value
					path = x[tmp].head;		//Travel to the next Vertex in the linked list
				}	
			}
		}

		//else theres no cycle: return false
			return false;
	}	
}

//pathReader outputs the path of our finished adjacency list starting at 1
void pathReader(adjList* x, Country Carr[])
{
	adjListNode* path = x[1].head;
	int prev = 1;	//prev marks the previous vertex travelled so we dont infinitly travel through 2 nodes
	int tmp = prev;
	double f = 0;
	cout << "1->";
	while (path != NULL)	//loop until we reach null
	{
		if (path->position == prev) //if the pointer points to the most recent node previously visited
		{
			prev = tmp;			//update the prev value to current value
			path = path->next;	//go to the next node in the list

		}

		if (path != NULL)	//If the pointer doesnt point to null
		{
			f = f + path->distance;
			cout << path->position << "->";
			prev = tmp;				//update the previous value
			tmp = path->position;	//update tmp to current value
			path = x[tmp].head;		//Travel to the next Vertex in the linked list
			
		}
	}
	f = f+ eucD(Carr[1].latitude, Carr[1].longitude, Carr[tmp].latitude, Carr[tmp].longitude);	//Adds an edge from the last vertex to the first;
	cout<< "1" << endl << "Total weight is " << f << endl;	//output first vertrx and total weight of TSP
}

//minPath is the greedy algorithm that picks the smallest edges to insert into our adjancency List:
void minPath(vector<edge> eV, int size, Country Carr[], adjList* arr)
{
	//This line of code creates our adjancency list pointer and makes every head point to NULL
	arr = new adjList[size];
	for (int d = 0; d < size; d++)
	{
		arr[d].head = NULL;
	}

	//Loop through our vector of edges: adding the minimum edge with no cycle and less than 2 connected
	for (int i = 0; i < eV.size(); i++)
	{
		if (listSize(arr, eV[i].start) < 2 && listSize(arr, eV[i].finish) < 2)	//If the size of our List for both vertices is less than 2:
		{
			if (!isCycle(arr, eV[i]))		//If there is no cycle
			{
				addEdge(arr, eV[i].start, eV[i].finish, eV[i].distance);	//Add the edge to the adjacency List
				
			}
		}
	}

		pathReader(arr, Carr);	//output the completed tour
}

//Function headers for Heap Sort:
inline int leftChild(int x);
void percDown(vector<edge> & heapList, int half, int vecSize);
void heapSort(vector<edge>& numList);

int main()
{
	using namespace chrono;
	//t1 and t2 are used to store the start and end times of an algorithm
	steady_clock::time_point t1;
	steady_clock::time_point t2;

	Country ta[10000];
	double t, v;
	int r;
	int o = 1;

	string a = "wi29.tsp";
	string b = "dj38.tsp";
	string c = "qa194.tsp";
	ifstream read;	//read reads from file
	string readStr;
	read.open(a);	//opens file (Change name for different country) This one is for Western Sahara

	//these FOR loops skip through the header of the tsp files and finds the size of the country
	int tsize = 0;
	for (int e = 0; e < 4; e++)
	{
		(getline(read, readStr));
	}
	read >> readStr >> readStr >> tsize;
	tsize++;
	for (int e = 0; e < 3; e++)
	{
		(getline(read, readStr));
	}

	//While loop insers r, t, and v as the country number and it's latitude and longitude into the array
	while (read >> r >> t >> v)
	{
		ta[o].position = r;
		ta[o].latitude = t;
		ta[o].longitude = v;
		o++;
	}
	read.close();

	//Initialize adjlist and edge vector
	adjList* countryAdjList = NULL;
	vector<edge> theEdges;

	t1 = steady_clock::now();
	minDist(theEdges, ta, tsize);	//create a list of all edges from the array into a vector
	heapSort(theEdges);				//Heap sort the edges from least to greatest distance
	minPath(theEdges, tsize, ta, countryAdjList);	//pick edges as long as they dont create a cycle and vertices less than 2
	t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "Time for our tour on input size " << tsize - 1 << " is " << time_span.count() << endl;
}


//HEAP SORT FUNCTION DEFINITIONS:

inline int leftChild(int x)	//returns position of the left child of a heap node
{
	return 2 * x + 1;
}

void percDown(vector<edge> & heapList, int half, int vecSize) //percDown takes arguments: vector of ints, middle position, and size of vector
{
	int child;	//child is for the child of the root
	edge temp;	//temp  is the midpoint of the list

	for (temp = move(heapList[half]); leftChild(half) < vecSize; half = child)	//perc down continuously, halving each time
	{
		child = leftChild(half);
		if (child != vecSize - 1 && heapList[child].distance < heapList[child + 1].distance)
			++child;
		if (temp.distance < heapList[child].distance)
			heapList[half] = move(heapList[child]);
		else break;

	}
	heapList[half] = move(temp);		//moves the element to it's proper place
}

void heapSort(vector<edge>& numList)
{
	for (int i = numList.size() / 2 - 1; i >= 0; i--)
	{
		percDown(numList, i, numList.size());			//Builds a heap with the vector
	}
	for (int j = numList.size() - 1; j > 0; j--)
	{
		swap(numList[0], numList[j]);	//swap element with position 0
		percDown(numList, 0, j);		//perform another percDown
	}

}
