#ifndef SUBDIVIDE_H
#define SUBDIVIDE_H

#include "slVector3.H"

class Mesh {
  class Vertex {
  public:
	bool flag;
	int index;
	SlVector3 nx;
	SlVector3 x;
	HalfEdge *e;
  };

  class Face {
  public:
	HalfEdge *e;
  };
  
  class HalfEdge {
  public:
	Vertex *v;
	Face *f;
	HalfEdge *next, *prev, *pair;
  };

  Mesh(const std::vector<SlVector3> &vertices, const std::vector<SlTri> &triangles);
  void divideEdge(HalfEdge *e);
  void triangulateFace(Face *f);
  std::vector<Vertex> *vertices;
  std::vector<Face> *faces;
  std::vector<HalfEdge> *edges;
};


#endif
