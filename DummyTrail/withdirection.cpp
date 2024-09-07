#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <climits>
#include <cmath>
#include <cfloat>  // For DBL_MAX
#include <stack>   // For std::stack

using namespace std;

// Constants for Earth radius in kilometers
const double EARTH_RADIUS = 6371.0;

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
    unordered_map<string, vector<pair<string, double>>> adjList;  // Adjacency list for graph
    unordered_map<string, Block> blockPositions;  // Store block names and positions

public:
    // Add blocks with geographical coordinates
    void addBlock(string name, double latitude, double longitude) {
        blockPositions[name] = Block(name, latitude, longitude);
    }

    // Add edges between blocks with distances (distances are automatically calculated)
    void addEdge(string u, string v) {
        double distance = haversineDistance(blockPositions[u], blockPositions[v]);
        adjList[u].push_back({v, distance});
        adjList[v].push_back({u, distance});  // Assuming undirected graph
    }

    // Dijkstra's algorithm to find the shortest path from source to destination
    void dijkstraWithPath(string start, string end) {
        unordered_map<string, double> distances;
        unordered_map<string, string> previous;  // To track the path
        for (auto& node : adjList) {
            distances[node.first] = DBL_MAX;  // Use DBL_MAX for initial distances
        }
        distances[start] = 0;

        // Priority queue to pick the smallest distance node
        priority_queue<pair<double, string>, vector<pair<double, string>>, greater<pair<double, string>>> pq;
        pq.push({0, start});

        while (!pq.empty()) {
            double dist = pq.top().first;
            string node = pq.top().second;
            pq.pop();

            if (node == end) break;

            for (auto& neighbor : adjList[node]) {
                string nextNode = neighbor.first;
                double weight = neighbor.second;

                if (dist + weight < distances[nextNode]) {
                    distances[nextNode] = dist + weight;
                    previous[nextNode] = node;
                    pq.push({distances[nextNode], nextNode});
                }
            }
        }

        if (distances[end] == DBL_MAX) {
            cout << "No path found between " << start << " and " << end << endl;
        } else {
            cout << "Shortest distance between " << start << " and " << end << " is " << distances[end] << " kilometers." << endl;
            printPathWithDirections(previous, start, end);
        }
    }

    // Function to print the path and directions
    void printPathWithDirections(unordered_map<string, string>& previous, string start, string end) {
        stack<string> path;  // Use std::stack to track the path
        string current = end;

        while (current != start) {
            path.push(current);
            current = previous[current];
        }
        path.push(start);

        cout << "Directions: " << endl;
        string prevBlock = path.top();
        path.pop();

        double prevBearing = 0.0;

        while (!path.empty()) {
            string nextBlock = path.top();
            path.pop();

            Block from = blockPositions[prevBlock];
            Block to = blockPositions[nextBlock];

            double bearing = calculateBearing(from, to);
            if (prevBlock != start) {
                string direction = getDirection(prevBearing, bearing);
                cout << direction << " to " << nextBlock << " (" << haversineDistance(from, to) << " kilometers)" << endl;
            } else {
                cout << "Start at " << nextBlock << endl;
            }

            prevBearing = bearing;
            prevBlock = nextBlock;
        }
    }
};

int main() {
    Graph campusGraph;

    // Adding blocks with latitude and longitude
    campusGraph.addBlock("BlockA", 12.9716, 77.5946);  // Example coordinates (lat, lon)
    campusGraph.addBlock("BlockB", 12.9721, 77.5937);
    campusGraph.addBlock("BlockC", 12.9736, 77.5952);
    campusGraph.addBlock("BlockD", 12.9711, 77.5980);
    campusGraph.addBlock("BlockE", 12.9745, 77.5971);

    // Adding edges (distances are automatically calculated)
    campusGraph.addEdge("BlockA", "BlockB");
    campusGraph.addEdge("BlockB", "BlockC");
    campusGraph.addEdge("BlockA", "BlockD");
    campusGraph.addEdge("BlockC", "BlockE");

    string startPoint = "BlockA";
    string endPoint = "BlockE";

    // Find the shortest path and directions between two blocks
    campusGraph.dijkstraWithPath(startPoint, endPoint);

    return 0;
}
