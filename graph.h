// graph.h <Starter Code>
// Name: Adrian Velazquez
//
// Basic graph class using adjacency List representation.
// Modified using maps to keep track of vertices
// and edges.
// University of Illinois at Chicago
// Class: CS 251; Summer 2022
// Project #5 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

using namespace std;

template <typename VertexT, typename WeightT>
class graph {
 private:
  size_t numLocations;  // number of vertices in the graph
  map<VertexT, map<VertexT, WeightT> >
      table;  // map containing vertices and their edges with their weights.
 public:
  //
  //
  // Default Constructor:
  // here it sets numLocations to zero as the graph is empty
  graph() { this->numLocations = 0; }
  //
  // clear function:
  // Removes all vertices and edges
  void clear() {
    table.clear();
    this->numLocations = 0;
  }
  //
  // copy operator=
  //
  // Called when you assign one graph into another, i.e. this = other;
  //
  graph& operator=(const graph& other) {
    if (this == &other) {
      return *this;
    }
    this->table = other.table;
    this->numLocations = other.numLocations;
    return *this;
  }
  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const { return this->table.size(); }
  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    int numEdges = 0;       // used to keep track of number of edges.
    for (auto e : table) {  // traverse our table
      numEdges += e.second.size();
      // increase numEdges by the size of each value in our map.
    }
    return numEdges;  // returns number of edges.
  }
  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v) {
    if (table.count(v)) {  // if vertex already exists...
      return false;
    }
    table[v] = map<VertexT, WeightT>();
    numLocations++;  // increment numLocations as we added a vertex...
    return true;
  }
  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    // if either vertex is not found...
    if (!table.count(from) || !table.count(to)) {
      return false;
    }
    // table[to][from] = weight;
    table[from][to] = weight;
    return true;
  }
  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    // test where if from or to vertex is not found in graph
    // if to vertex is not found as an edge
    if (!table.count(from) || !table.at(from).count(to) || !table.count(to)) {
      return false;
    } else {  // if found return the weight associated with it...
      weight = table.at(from).at(to);
    }
    return true;
  }
  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> S;
    if (!table.count(v)) {  // if v cannot be found in our table...
      return S;             // return an empty set
    } else {
      for (auto V : table.at(v)) {
        S.insert(V.first);  // inserts edge that is associated with v...
      }
    }
    return S;
  }
  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    vector<VertexT> vertices;       // used to store vertices
    for (auto V : table) {          // loop thru our table
      vertices.push_back(V.first);  // push each vertice in our vector
    }
    return vertices;  // return vector vertices
  }
  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;
    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;
    output << endl;
    output << "**Vertices:" << endl;
    for (auto vertices : table) {
      output << vertices.first << endl;
    }
    output << endl;
    output << "**Edges:" << endl;
    for (auto vertice : table) {
      output << "row " << vertice.first << " ";
      for (auto edge : vertice.second) {
        output << "(" << edge.first << ", " << edge.second << ") ";
      }
      output << endl;
    }
    output << endl;
    output << "**************************************************" << endl;
  }
};
