// application.cpp <Starter Code>
// Name: Adrian Velazquez
//
// University of Illinois at Chicago
// CS 251: Summer 2022
// Project #5 - Openstreet Maps
// Here, it the menu of our application that
// prompts the user into inputting 2
// building that 2 people are currently at
// and find center building along their paths
// to meet up.
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"
using namespace std;
using namespace tinyxml2;

class prioritize {  // you could also use a struct
 public:
  bool operator()(const pair<long long, double>& p1,
                  const pair<long long, double>& p2) const {
    return p1.second > p2.second;
  }
};

BuildingInfo searchBuilding(vector<BuildingInfo>& Buildings, string query);
BuildingInfo nearestBuilding(vector<BuildingInfo>& Buildings,
                             Coordinates midpoint, set<string>& invalidB);
long long nearestNode(vector<FootwayInfo>& Footways,
                      map<long long, Coordinates>& Nodes, BuildingInfo b);
void DijkstraShortestPath(long long startV, graph<long long, double>& G,
                          map<long long, double>& distances,
                          map<long long, long long>& previousNode);
void printPath(long long startV, long long endV,
               map<long long, long long>& previousNode);
void printInfo(BuildingInfo Location);
void printNode(long long building, map<long long, Coordinates>& Nodes);
//
// Implement your creative component application here
// add arguments
//
void creative() {}
//
// Implement your standard application here
// add a parameter for the graph you make.
//
void application(map<long long, Coordinates>& Nodes,
                 vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
                 graph<long long, double>& G) {
  string person1Building, person2Building;
  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);
  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    set<string> invalidB;  // keeps track of invalid center Buildings
    BuildingInfo Building1, Building2, centerB, tempCenter;
    Building1 = searchBuilding(Buildings, person1Building);
    Building2 = searchBuilding(Buildings, person2Building);
    if (Building1.Fullname ==
        "") {  // if Person 1's building doesn't exist in Buildings
      cout << "Person 1's building not found" << endl << endl;
      cout
          << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }
    if (Building2.Fullname ==
        "") {  // if Person 2's building doesn't exist in Buildings
      cout << "Person 2's building not found" << endl;
      cout
          << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }
    cout << endl;
    cout << "Person 1's point:" << endl;
    printInfo(Building1);
    cout << "Person 2's point:" << endl;
    printInfo(Building2);
    // Locate Center Building
    const double INF = numeric_limits<double>::max();
    Coordinates p1 = Building1.Coords;
    Coordinates p2 = Building2.Coords;
    Coordinates midpoint = centerBetween2Points(p1.Lat, p1.Lon, p2.Lat, p2.Lon);
    centerB = nearestBuilding(Buildings, midpoint, invalidB);
    if (centerB.Fullname == "") {
      cout << "Center Building not found" << endl;
    }
    cout << "Destination Building:" << endl;
    printInfo(centerB);
    cout << endl;
    // Find Nearest Footway Nodes from buildings 1, 2 & Center
    long long b1, b2, b3;
    b1 = nearestNode(Footways, Nodes, Building1);
    cout << "Nearest P1 node:" << endl;
    printNode(b1, Nodes);
    b2 = nearestNode(Footways, Nodes, Building2);
    cout << "Nearest P2 node:" << endl;
    printNode(b2, Nodes);
    b3 = nearestNode(Footways, Nodes, centerB);
    cout << "Nearest destination node:" << endl;
    printNode(b3, Nodes);
    cout << endl;
    // Run Dijkstra’s Algorithm
    // map containing the distance from starting V to shortest path
    // of each vertex.
    map<long long, double> distances1, distances2;
    map<long long, long long> previousNode1, previousNode2;
    DijkstraShortestPath(b1, G, distances1, previousNode1);
    DijkstraShortestPath(b2, G, distances2, previousNode2);
    if (distances1.at(b2) >= INF || distances2.at(b1) >= INF) {
      cout << "Sorry, destination unreachable." << endl;
      cout << endl;
      cout
          << "Enter person 1's building (partial name or abbreviation), or #> ";
      getline(cin, person1Building);
      continue;
    }
    while (distances1.at(b3) >= INF || distances2.at(b3) >= INF) {
      cout << "At least one person was unable to reach the destination "
              "building. Finding next closest building..."
           << endl;
      invalidB.insert(
          centerB
              .Fullname);  // insert current invalid center building into set.
      tempCenter = nearestBuilding(Buildings, midpoint, invalidB);
      invalidB.insert(
          tempCenter.Fullname);  // insert incase this center is also invalid
      cout << "New destination building:" << endl;
      printInfo(tempCenter);
      b3 = nearestNode(Footways, Nodes, tempCenter);
      cout << "Nearest destination node:" << endl;
      printNode(b3, Nodes);
      cout << endl;
      distances1.clear();
      distances2.clear();
      previousNode1.clear();
      previousNode2.clear();
      DijkstraShortestPath(b1, G, distances1, previousNode1);
      DijkstraShortestPath(b2, G, distances2, previousNode2);
    }
    // if both can reach center buiding
    cout << "Person 1's distance to dest: " << distances1.at(b3) << " miles"
         << endl;
    cout << "Path: " << b1;
    printPath(b1, b3,
              previousNode1);  // path between person 1 and center building
    cout << endl << endl;
    cout << "Person 2's distance to dest: " << distances2.at(b3) << " miles"
         << endl;
    cout << "Path: " << b2;
    printPath(b2, b3,
              previousNode2);  // path between person 2 and center building
    // another navigation?
    cout << endl << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}
//
// searches thru Buildings to see if query is a
// valid building entry and returns info if found
BuildingInfo searchBuilding(vector<BuildingInfo>& Buildings, string query) {
  BuildingInfo
      found;  // if building is found from query, store that building's info
  set<string> firstInstance;
  for (auto e : Buildings) {
    if (query == e.Abbrev) {  // case where query is equal to Abbreviation
      found = e;
      return found;
    } else if (e.Fullname.find(query) !=
               string::npos) {  // case where query is a substring in Buildings
                                // full name
      found = e;
      return found;
    }
  }
  return found;
}
// find the nearest building center and returns its info.
BuildingInfo nearestBuilding(vector<BuildingInfo>& Buildings,
                             Coordinates midpoint, set<string>& invalidB) {
  BuildingInfo FoundBuilding;
  Coordinates temp;
  double distance;
  const double INF = numeric_limits<double>::max();
  double min = INF;
  for (auto e : Buildings) {
    temp = e.Coords;
    distance =
        distBetween2Points(midpoint.Lat, midpoint.Lon, temp.Lat, temp.Lon);
    if (distance < min) {
      if (invalidB.count(e.Fullname)) {
        continue;
      }
      FoundBuilding = e;
      min = distance;
    }
  }
  return FoundBuilding;
}
// finds nearest footway node to the building
// and returns that ID.
long long nearestNode(vector<FootwayInfo>& Footways,
                      map<long long, Coordinates>& Nodes, BuildingInfo b) {
  long long ID;
  double dist;
  const double INF = numeric_limits<double>::max();
  double min = INF;
  for (auto e : Footways) {
    for (int i = 0; i < e.Nodes.size(); i++) {
      long long N1 = e.Nodes[i];
      Coordinates temp = Nodes[N1];
      Coordinates Build = b.Coords;
      dist = distBetween2Points(temp.Lat, temp.Lon, Build.Lat, Build.Lon);
      if (dist < min) {
        ID = e.Nodes[i];
        min = dist;
      }
    }
  }
  return ID;
}
// prints out Building's info
void printInfo(BuildingInfo Loc) {
  cout << Loc.Fullname << endl;
  cout << "(" << Loc.Coords.Lat << ", " << Loc.Coords.Lon << ")" << endl;
}
// prints out the nearest Footway Node's info
void printNode(long long Building, map<long long, Coordinates>& Nodes) {
  cout << Nodes[Building].ID << endl;
  cout << "(" << Nodes[Building].Lat << ", " << Nodes[Building].Lon << ")"
       << endl;
}
// loops and keep keeps track of the shortest
// distance from startV to each node and
// keep track of the predecessors
void DijkstraShortestPath(long long startV, graph<long long, double>& G,
                          map<long long, double>& distances,
                          map<long long, long long>& previousNode) {
  priority_queue<pair<long long, double>, vector<pair<long long, double>>,
                 prioritize>
      unvisitedQueue;
  const double INF = numeric_limits<double>::max();
  set<long long> visited;           // keeps track of nodes visited.
  pair<long long, double> curV;     // used to store nodes from unvisitedQueue
  for (auto g : G.getVertices()) {  // loop thru all vertices in graph
    distances[g] =
        INF;  // set each vertex distance to itself equal to infinity.
    previousNode[g] = 0;  // set each vertex previousNode to 0.
    unvisitedQueue.push(make_pair(g, INF));  // push curV into unvisitedQueue.
  }
  distances[startV] = 0;     // set startingV distance to itself equal to 0
  previousNode[startV] = 0;  // set startvingV's previosNode = 0
  unvisitedQueue.push(make_pair(startV, 0));
  while (!unvisitedQueue.empty()) {
    curV = unvisitedQueue.top();
    unvisitedQueue.pop();
    if (distances[curV.first] == INF) {  // if curV's distance is equal to INF..
      break;
    } else if (visited.count(curV.first)) {  // if we visited curV already...
      continue;
    } else {
      visited.insert(curV.first);  // store curV into visited set.
    }
    double adjDist = 0;
    for (auto e : G.neighbors(curV.first)) {  // loop thru curV's neighbors.
      double edgeWeight = 0;
      G.getWeight(curV.first, e, edgeWeight);
      adjDist = edgeWeight + distances.at(curV.first);
      if (adjDist < distances.at(e)) {
        distances.at(e) = adjDist;
        previousNode[e] = curV.first;
        unvisitedQueue.push(make_pair(e, adjDist));
      }
    }
  }
}
// Recursive function that prints the path
// from startV to endV if there is one
void printPath(long long startV, long long endV,
               map<long long, long long>& previousNode) {
  // check if our endV is a valid location
  if (previousNode.at(endV) == 0) {
    return;
  }
  // if our startV matches our endV which indicates we have a path.
  if (startV == endV) {
    return;
  }
  printPath(startV, previousNode.at(endV), previousNode);
  cout << "->" << endV;
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates> Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo> Footways;
  // info about each building, in no particular order
  vector<BuildingInfo> Buildings;
  XMLDocument xmldoc;
  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);
  string def_filename = "map.osm";
  string filename;
  cout << "Enter map filename> ";
  getline(cin, filename);
  if (filename == "") {
    filename = def_filename;
  }
  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }
  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);
  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);
  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);
  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;
  //
  // build the graph, output stats:
  //
  graph<long long, double> G;
  // add vertices
  for (auto vertex : Nodes) {
    G.addVertex(vertex.first);
  }
  // add edges
  for (auto e : Footways) {
    for (int i = 0; i < e.Nodes.size() - 1; i++) {
      long long N1 = e.Nodes[i];
      long long N2 = e.Nodes[i + 1];
      Coordinates C1 = Nodes[N1];
      Coordinates C2 = Nodes[N2];
      double dist = distBetween2Points(C1.Lat, C1.Lon, C2.Lat, C2.Lon);
      G.addEdge(N1, N2, dist);
      G.addEdge(N2, N1, dist);
    }
  }
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;
  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
       << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    application(Nodes, Footways, Buildings, G);
  } else if (userInput == "c") {
    creative();
  }
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
