/*
Munem Rastgir
Nearest Insertion Algorithm

Main uses functions: 
createMatrix- To create an adjancency matrix of the cities in the file
tourStart- Finds the smallest distance between 2 cities and creates a tour of them 
cityList- Makes a vector of cities not on the tour and their distances
countryTour- Creates a finished tour of the whole country, adding each new cities until our city list is empty
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <fstream>
#include<algorithm>
#include<string>
#include <chrono>
using namespace std;

const int matrixSize = 200;	//size of adjacency matrix: can increase or decrease to easily accomodate bigger countries

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

//Node has it's position(number of what node it is), the distance it takes to travel to the next node, and pointer to next node
struct Node {
	int position;
	double distance;
	Node* next;
};

double eucD(double x, double y, double a, double b)	//eucD takes 4 arguments and returns the euclidian distance of them
{
	double f = (x - a)*(x - a);	//x-a squared
	double s = (y - b)*(y - b);	//y-b squared
	return sqrt(f + s);			//sqrt of sum
}

void createMatrix(double adjMatrix[][matrixSize], Country arr[], int size)	//createMatrix forms an adjacency matrix from the array of countries
{
	for (int i = 1; i < size; i++)	//Nested for loops loop through every country:
	{
		for (int j = 1; j < size; j++)	
		{
			if (i == j)	//if i and j are equal
			{
				adjMatrix[i][j] = INFINITY;	//set it's edge value to infinity because we cant have an edge with 2 of the same vertex
			}
			else
			{
				adjMatrix[i][j] = eucD(arr[i].latitude, arr[i].longitude, arr[j].latitude, arr[j].longitude);	//calculates edge distance
			}
		}	
	}
}

//makeNode function creates a new node with argument pos which is it's position
Node* makeNode(int pos)
{
	Node* x = new Node;
	x->position = pos;	//position points to int argument
	x->next = NULL;		//points to null
	return x;	//returns pointer to finished node
}

//tourStart function creates a tour between the 2 vertices with the smallest distance:
void tourStart(Node* &x, double adjMatrix[][matrixSize], int size)	
{
	x = NULL;	//clears the pointer

	double min = adjMatrix[1][1];	//sets min equal to the first adjMatrix value
	int first = 2;		//first and second are the 2 vertices for our starting tour
	int second = 3;
		for (int i = 1; i < size; i++)	//Nested Loop through entire matrix to find the smallest edge
		{
			for (int j = 1; j < size; j++)
			{
				if (adjMatrix[i][j] < min)		//if the edge is less than current min
				{
					min = adjMatrix[i][j];	//update the min
					first = i;				//update the first and second vertices
					second = j;
				}
			}
		}

		//When we get the smallest edge:
		Node* newnode = makeNode(first);	//Create the first Node	
		newnode->distance = min;
		x = newnode;					//adds the first vertex to the list x

		newnode = makeNode(second);		//Create second Node
		newnode->distance = min;
		x->next = newnode;			//adds second vertex to the list x

		Node* tempptr = x->next;
		tempptr->next = x;		//we return back to the first vertex, adding it to the list
		tempptr = NULL;
		delete tempptr;
}


//puts all cities not in the tour in vector eV with the distance of its closest tour neighbor
void cityList(vector<edge>& eV, Node* &head,double adjMatrix[][matrixSize], int size)	
{
	//fr and snd are the 2 nodes in the tour
	int fr = head->position;
	Node* tmpptr = head->next;
	int snd = tmpptr->position;
	
	for (int i = 1; i < size; i++)
	{
		if (i != fr && i != snd)	//if i is not a vertex on the tour
		{
			edge temp;	//make a new edge
			temp.start = i;
			int min = adjMatrix[fr][i];		//set min equal to the distance of i and the first vertex
			if (adjMatrix[snd][i] < min)	//if it is closer to the second vertex: change min
			{
				min = adjMatrix[snd][i];
				temp.finish = snd;
			}
			else {
				temp.finish = fr;
			}
			temp.distance = min;
			eV.push_back(temp);	//pushback the completed edge into the vertex
		}
	}

	tmpptr = NULL;
	delete tmpptr;
}

int pickMinEdge(vector<edge> eV)	//Pick the smallest edge in the whole vector to be inserted next
{
	double min = eV[0].distance;	//min equals distance of first edge in vector
	int r = 0;
	for (int i = 0; i < eV.size(); i++)	//loops through vector to find min edge
	{
		if (eV[i].distance < min)
		{
			min = eV[i].distance;
			r = i;
		}
	}
	return r;		//returns r: the position in the vector
}

//pick edge returns the edge in the tour that will have a new city inserted in between it's vertices
edge pickEdge(Node* head, double adjMatrix[][matrixSize], int newNode)
{
	Node* tmp = head;

	int stopPath = tmp->position;
	edge min;
	min.start = tmp->position;

	int fp = tmp->position;
	int fd = tmp->distance;
	tmp = tmp->next;
	min.finish = tmp->position;

	//addedDist is the value that adding the new node gives us
	double addedDist = (adjMatrix[fp][newNode] + adjMatrix[tmp->position][newNode] - fd);	
	
	while (tmp->position != stopPath)	//go through the linked list to find the smallest added distance
	{
		fp = tmp->position;
		fd = tmp->distance;
		tmp = tmp->next;
		double newDist = (adjMatrix[fp][newNode] + adjMatrix[tmp->position][newNode] - fd);
		if (newDist < addedDist)	//swap distances to the smallest and update the start and finish
		{
			addedDist = newDist;
			min.start = fp;
			min.finish = tmp->position;
		}
	}

	tmp = NULL;
	delete tmp;

	return min;	//returns the min edge that we will use to insert our vertex in between
}

//addEdges adds an edge to a linked list tour. Arguments of A Node pointer, 3 vertices, and 2 distances
void addEdge(Node* head, int fVert, int sVert, int tVert, double dist, double sdist)
{
	Node* holder = head;	//holder is a temporary pointer to the 2nd vertex because the link to it is going to break
	while (holder->position != tVert)
	{
		holder = holder->next;
	}

	Node*tmp = head;		//tmp points to the first node
	while (tmp->position != fVert)
	{
		tmp = tmp->next;
	}
	tmp->next = makeNode(sVert);	//change tmp->next to a new node
	tmp->distance = dist;			//update it's distance

	tmp = tmp->next;
	tmp->distance = sdist;		//Insert the distance of the added node

	tmp->next = holder;			//the node points to the 3rd vertex that was originally attached to the first vertex

	//Delete the temp pointers to avoid memory leaks
	tmp = NULL;
	delete tmp;
	holder = NULL;
	delete holder;
}

//Update checks the most recently added Node to see if the distance to it 
//and the other cities are smaller than it's current recorded distance
void update(vector<edge>& eV,double adjMatrix[][matrixSize], int newest)
{
	for (int i = 0; i < eV.size(); i++)
	{
		if (adjMatrix[newest][eV[i].start] < eV[i].distance)
		{
			eV[i].distance = adjMatrix[newest][eV[i].start];
		}
	}
}

//CountryTour is our completed function that completes the tour on our 2 first tour vertices
void countryTour(Node* head, vector<edge> eV, double adjMatrix[][matrixSize])
{
	while (eV.size() > 0)	//while the list of countries not on the tour is greater than 0
	{
		int f = pickMinEdge(eV);	//find the smallest edge in the list of non tour countries
		edge g = pickEdge(head, adjMatrix, eV[f].start);	//create an edge g that is the edge that we will perform the insertion between
		int nextCity = eV[f].start;
		addEdge(head, g.start, nextCity, g.finish, adjMatrix[g.start][nextCity], adjMatrix[g.finish][nextCity]);	//Add the new edge to the tour
		eV.erase(eV.begin() + f);	//erase the city we inserted from the vector
		update(eV, adjMatrix, nextCity);	//update the distances of the remaining non tour cities
	}
	//Function to output the finished linked list:
		Node* tmp = head;
		double weight = head->distance;
		cout << tmp->position << "->";
		tmp = tmp->next;
		while (tmp->position != head->position)
		{
			cout << tmp->position << "->";
			weight = weight + tmp->distance;
			tmp = tmp->next;
		}
		cout << tmp->position << endl << "Distance is: " << weight << endl;
}

int main()
{
	using namespace chrono;
	//t1 and t2 are used to store the start and end times of an algorithm
	steady_clock::time_point t1;
	steady_clock::time_point t2;

	Country ta[1000];	//ta is an array of Countries that we will read from file
	double testadj[matrixSize][matrixSize];	//testadj is out adjacency matrix with a size gathered from the global variable matrixSize
	double t, v;
	int r;
	int o = 1;
	
	string a = "wi29.tsp";
	string b = "dj38.tsp";
	string c = "qa194.tsp";
	ifstream read;	//read reads from file
	string readStr;
	read.open(a);	//opens file(Change name for different country) This one is for Western Sahara

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

	//While loop insers r t v as the country number and it's latitude and longitude into the array
	while (read >> r >> t >> v)
	{
		ta[o].position = r;
		ta[o].latitude = t;
		ta[o].longitude = v;
		o++;
	}
	read.close();

	//Initialize empty edge vector and country linked list:
	vector<edge> notTour;
	Node* countryList = NULL;

	t1 = steady_clock::now();
	createMatrix(testadj, ta, tsize);				//Creates adjacency matrix,
	tourStart(countryList, testadj, tsize);			//the first tour of 2 smallest cities,
	cityList(notTour, countryList, testadj, tsize);	//a list of cities not on the tour,
	countryTour(countryList, notTour, testadj);		//and the finished tour of the whole country
	t2 = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "Time for our tour on input size " << tsize - 1 << " is " << time_span.count() << endl;
}
