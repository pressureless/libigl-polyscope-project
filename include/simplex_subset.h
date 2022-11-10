#pragma once

#include <set>

class SimplexSubset {

  public:
    std::set<int> vertices;
    std::set<int> edges;
    std::set<int> faces;
    std::set<int> tets;

    /* Initialize an empty SimplexSubset. */
    SimplexSubset() {}
    
    /* Initialize a SimplexSubset with the given vertices, edges, and faces. */
    SimplexSubset(const std::set<int>& V, const std::set<int>& E, const std::set<int>& F) {

        vertices = V;
        edges = E;
        faces = F;
    }
    
    /* Initialize a SimplexSubset with the given vertices, edges, and faces. */
    SimplexSubset(const std::set<int>& V, const std::set<int>& E, const std::set<int>& F, const std::set<int>& T) {

        vertices = V;
        edges = E;
        faces = F;
        tets = T;
    }

    /* Make a deep copy of the input SimplexSubset and return it as a new
     * SimplexSubset.
     */
    SimplexSubset deepCopy() const {
        std::set<int> newVertices = vertices;
        std::set<int> newEdges = edges;
        std::set<int> newFaces = faces;
        std::set<int> newTets = tets;
        return SimplexSubset(newVertices, newEdges, newFaces, newTets);
    }

    /* Add a vertex to this subset. */
    void addVertex(int index) {
        vertices.insert(index);
    }

    /* Add a set of vertices to this subset. */
    void addVertices(const std::set<int>& V) {
        for (std::set<int>::iterator it = V.begin(); it != V.end(); ++it) {
            vertices.insert(*it);
        }
    }

    /* Delete a vertex from this subset. */
    void deleteVertex(int index) {
        vertices.erase(index);
    }

    /* Delete all vertices. */
    void deleteVertices() {
        this->vertices.clear();
    }

    /* Delete a set of vertices from this subset. */
    void deleteVertices(const std::set<int>& V) {
        for (std::set<int>::iterator it = V.begin(); it != V.end(); ++it) {
            vertices.erase(*it);
        }
    }

    /* Add an edge to this subset. */
    void addEdge(int index) {
        edges.insert(index);
    }

    /* Add a set of edges to this subset. */
    void addEdges(const std::set<int>& E) {
        for (std::set<int>::iterator it = E.begin(); it != E.end(); ++it) {
            edges.insert(*it);
        }
    }

    /* Delete an edge from this subset. */
    void deleteEdge(int index) {
        edges.erase(index);
    }

    /* Delete all edges. */
    void deleteEdges() {
        this->edges.clear();
    }

    /* Delete a set of edges from this subset. */
    void deleteEdges(const std::set<int>& E) {
        for (std::set<int>::iterator it = E.begin(); it != E.end(); ++it) {
            edges.erase(*it);
        }
    }
    
    /* Add a face to this subset. */
    void addFace(int index) {
        faces.insert(index);
    }

    /* Add a set of faces to this subset. */
    void addFaces(const std::set<int>& F) {
        for (std::set<int>::iterator it = F.begin(); it != F.end(); ++it) {
            faces.insert(*it);
        }
    }

    /* Delete a face from this subset. */
    void deleteFace(int index) {
        faces.erase(index);
    }

    /* Delete all faces. */
    void deleteFaces() {
        this->faces.clear();
    }

    /* Delete a set of faces from this subset. */
    void deleteFaces(const std::set<int>& F) {
        for (std::set<int>::iterator it = F.begin(); it != F.end(); ++it) {
            faces.erase(*it);
        }
    }
    
    /* Add a tet to this subset. */
    void addTet(int index) {
        tets.insert(index);
    }

    /* Add a set of tets to this subset. */
    void addTets(const std::set<int>& T) {
        for (std::set<int>::iterator it = T.begin(); it != T.end(); ++it) {
            tets.insert(*it);
        }
    }

    /* Delete a tet from this subset. */
    void deleteTet(int index) {
        tets.erase(index);
    }

    /* Delete all tets. */
    void deleteTets() {
        this->tets.clear();
    }

    /* Delete a set of tets from this subset. */
    void deleteTets(const std::set<int>& T) {
        for (std::set<int>::iterator it = T.begin(); it != T.end(); ++it) {
            tets.erase(*it);
        }
    }

    /* Returns true if subsets are equivalent. */
    bool equals(const SimplexSubset& other) const {
        // == compares elements at each position; but std::set always orders elements upon insertion/initialization
        return (vertices == other.vertices) && (edges == other.edges) && (faces == other.faces) && (tets == other.tets);
    }

    /* Adds a subset's vertices, edges, and faces to this subset. */
    void addSubset(const SimplexSubset& other) {
        this->addVertices(other.vertices);
        this->addEdges(other.edges);
        this->addFaces(other.faces);
        this->addTets(other.tets);
    }

    /* Removes a subset's vertices, edges, and faces from this subset. */
    void deleteSubset(const SimplexSubset& other) {
        this->deleteVertices(other.vertices);
        this->deleteEdges(other.edges);
        this->deleteFaces(other.faces);
        this->deleteTets(other.tets);
    }

    /* Print vertices. */
    void printVertices() {
        std::cout << "Vertices: ";
        for (std::set<int>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
            std::cout << *it << ", ";
        }
        std::cerr << std::endl;
    }

    /* Print edges. */
    void printEdges() {
        std::cout << "Edges: ";
        for (std::set<int>::iterator it = edges.begin(); it != edges.end(); ++it) {
            std::cout << *it << ", ";
        }
        std::cerr << std::endl;
    }
    
    /* Print Faces. */
    void printFaces() {
        std::cout << "Faces: ";
        for (std::set<int>::iterator it = faces.begin(); it != faces.end(); ++it) {
            std::cout << *it << ", ";
        }
        std::cerr << std::endl;
    }
    
    /* Print Tets. */
    void printTets() {
        std::cout << "Tets: ";
        for (std::set<int>::iterator it = tets.begin(); it != tets.end(); ++it) {
            std::cout << *it << ", ";
        }
        std::cerr << std::endl;
    }

    void deleteAll(){
        deleteVertices();
        deleteEdges();
        deleteFaces();
        deleteTets();
    }
    
    /* Print All. */
    void printAll() {
        printTets();
        printFaces();
        printEdges();
        printVertices();
    }
};
