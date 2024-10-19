#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm> 


using namespace std;

// Constants for Earth radius in kilometers
const double EARTH_RADIUS = 6371.0;
const int INF = numeric_limits<int>::max();  // To represent infinity

// Struct to hold the geographical coordinates of each block (latitude and longitude)
struct Block {
    string name;
    double latitude, longitude;  // Geographical coordinates

    Block() : name(""), latitude(0), longitude(0) {}
    Block(string n, double lat, double lon) : name(n), latitude(lat), longitude(lon) {}
};

// Function to convert degrees to radians
double toRadians(double degree) {
    return degree * M_PI / 180.0;
}

// Haversine formula to calculate distance between two geographical points
double haversineDistance(const Block& from, const Block& to) {
    double lat1 = toRadians(from.latitude);
    double lon1 = toRadians(from.longitude);
    double lat2 = toRadians(to.latitude);
    double lon2 = toRadians(to.longitude);

    double dLat = lat2 - lat1;
    double dLon = lon2 - lon1;

    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1) * cos(lat2) * sin(dLon / 2) * sin(dLon / 2);

    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double distance = EARTH_RADIUS * c;
    return distance;  // Distance in kilometers
}

// Function to calculate the bearing between two geographical points
double calculateBearing(const Block& from, const Block& to) {
    double lat1 = toRadians(from.latitude);
    double lon1 = toRadians(from.longitude);
    double lat2 = toRadians(to.latitude);
    double lon2 = toRadians(to.longitude);

    double dLon = lon2 - lon1;
    double x = sin(dLon) * cos(lat2);
    double y = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon);
    double initialBearing = atan2(x, y);

    // Convert bearing from radians to degrees
    double bearing = fmod((initialBearing * 180.0 / M_PI) + 360.0, 360.0);
    return bearing;  // Bearing in degrees
}

// Function to get direction based on bearing
string getDirection(double fromBearing, double toBearing) {
    double difference = toBearing - fromBearing;
    if (difference < 0) {
        difference += 360;
    }
    if (difference < 30 || difference > 330) {
        return "Go straight";
    } else if (difference < 180) {
        return "Turn right";
    } else {
        return "Turn left";
    }
}

class Graph {
    vector<vector<double>> adjMatrix;  // Adjacency matrix for distances
    vector<Block> blocks;  // Store block names and positions

public:
    // Constructor to initialize the graph with a given number of blocks
    Graph(int n) {
        adjMatrix.resize(n, vector<double>(n, INF));  // Initialize the adjacency matrix with INF
        blocks.resize(n);
    }

    // Add blocks with geographical coordinates
    void addBlock(int idx, string name, double latitude, double longitude) {
        blocks[idx] = Block(name, latitude, longitude);
    }

    // Add edges between blocks with distances (distances are automatically calculated)
    void addEdge(int u, int v) {
        double distance = haversineDistance(blocks[u], blocks[v]);
        adjMatrix[u][v] = distance;
        adjMatrix[v][u] = distance;  // Assuming undirected graph
    }

    // Dijkstra's algorithm to find the shortest path from source to destination
    void dijkstraWithPath(int start, int end) {
        int n = blocks.size();
        vector<double> distances(n, INF);  // Distance array
        vector<int> previous(n, -1);  // To track the path
        vector<bool> visited(n, false);  // To mark visited nodes

        distances[start] = 0;

        for (int i = 0; i < n; ++i) {
            int node = -1;

            // Pick the smallest distance node
            for (int j = 0; j < n; ++j) {
                if (!visited[j] && (node == -1 || distances[j] < distances[node])) {
                    node = j;
                }
            }

            if (distances[node] == INF) break;
            visited[node] = true;

            // Explore neighbors
            for (int neighbor = 0; neighbor < n; ++neighbor) {
                if (!visited[neighbor] && adjMatrix[node][neighbor] != INF) {
                    double newDist = distances[node] + adjMatrix[node][neighbor];
                    if (newDist < distances[neighbor]) {
                        distances[neighbor] = newDist;
                        previous[neighbor] = node;
                    }
                }
            }
        }

        if (distances[end] == INF) {
            cout << "No path found between " << blocks[start].name << " and " << blocks[end].name << endl;
        } else {
            cout << "Shortest distance between " << blocks[start].name << " and " << blocks[end].name << " is " << distances[end] << " kilometers." << endl;
            printPathWithDirections(previous, start, end);
        }
    }

    // Function to print the path and directions
    void printPathWithDirections(vector<int>& previous, int start, int end) {
        vector<int> path;

        for (int at = end; at != -1; at = previous[at]) {
            path.push_back(at);
        }

        if (path.back() != start) {
            cout << "No path found!" << endl;
            return;
        }

        reverse(path.begin(), path.end());

        cout << "Directions: " << endl;
        double prevBearing = 0.0;

        for (size_t i = 1; i < path.size(); ++i) {
            Block from = blocks[path[i - 1]];
            Block to = blocks[path[i]];
            double bearing = calculateBearing(from, to);

            if (i == 1) {
                cout << "Start at " << from.name << endl;
            } else {
                string direction = getDirection(prevBearing, bearing);
                cout << direction << " to " << to.name << " (" << haversineDistance(from, to) << " kilometers)" << endl;
            }

            prevBearing = bearing;
        }
    }
};

int main() {
    Graph campusGraph(5);

    // Adding blocks with latitude and longitude
    campusGraph.addBlock(0, "BlockA", 12.9716, 77.5946);
    campusGraph.addBlock(1, "BlockB", 12.9721, 77.5937);
    campusGraph.addBlock(2, "BlockC", 12.9736, 77.5952);
    campusGraph.addBlock(3, "BlockD", 12.9711, 77.5980);
    campusGraph.addBlock(4, "BlockE", 12.9745, 77.5971);

    // Adding edges (distances are automatically calculated)
    campusGraph.addEdge(0, 1);
    campusGraph.addEdge(1, 2);
    campusGraph.addEdge(0, 3);
    campusGraph.addEdge(2, 4);

    int startPoint = 0;
    int endPoint = 4;

    // Find the shortest path and directions between two blocks
    campusGraph.dijkstraWithPath(startPoint, endPoint);

    return 0;
}
