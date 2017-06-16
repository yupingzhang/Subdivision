#include "slVector.H"
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <getopt.h>

class SlTri {
  unsigned int indices[3];
public:
  inline unsigned int &operator[](const unsigned int &i) { return indices[i];};
  inline unsigned int operator[](const unsigned int &i) const { return indices[i];};
};

class SlQuad {
  unsigned int indices[4];
public:
  inline unsigned int &operator[](const unsigned int &i) { return indices[i];};
  inline unsigned int operator[](const unsigned int &i) const { return indices[i];};
};

class Mesh {
  class HalfEdge;

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
	int index;
	HalfEdge *e;
  };
  
  class HalfEdge {
  public:
	int index;
	Vertex *v;
	Face *f;
	HalfEdge *next, *prev, *pair;
  };

  std::vector<Vertex *> vertices;
  std::vector<Face *> faces;
  std::vector<HalfEdge *> edges;
  void setIndices();

  bool LOOP;

public: 
  ~Mesh();
  Mesh(const std::vector<SlVector3> &vertices, const std::vector<SlTri> &triangles);
  Mesh(const std::vector<SlVector3> &vertices, const std::vector<SlQuad> &quads);
  void divideEdge(HalfEdge *e);
  void triangulateFace(Face *f);
  void splitFace(Face *f);
  void subdivide();
  void writeObj(char *fname);
  void printMesh();
};

void Mesh::printMesh() {
  for (unsigned int i=0; i<vertices.size(); i++) vertices[i]->index = i;
  for (unsigned int i=0; i<faces.size(); i++) faces[i]->index = i;
  for (unsigned int i=0; i<edges.size(); i++) edges[i]->index = i;
  
  for (unsigned int i=0; i<faces.size(); i++) {
	Face *f = faces[i];
	HalfEdge *e = f->e;
	while (true) {
	  std::cout<<e->v->index+1<<" ";
	  e = e->next;
	  if (e == f->e) break;
	}
	std::cout<<std::endl;
  }
}

struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const {
	return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};

Mesh::Mesh(const std::vector<SlVector3> &vertices, const std::vector<SlTri> &triangles) {
  std::unordered_map<std::pair<int, int>, HalfEdge *, pairhash > edgeMap;
  std::unordered_map<int, Vertex * > vertexMap;

  for (unsigned int t=0; t<triangles.size(); t++) {
	const SlTri &tri = triangles[t];
	Face *f = new Face();
	HalfEdge *fedges[3];
	for (unsigned int i=0; i<3; i++) {
	  unsigned int j = (i+1)%3;
	  HalfEdge *e = new HalfEdge();
	  this->edges.push_back(e);
	  e->f = f;
	  fedges[i] = e;
	  if (edgeMap.count(std::make_pair(tri[j],tri[i]))) {
		HalfEdge *pair = edgeMap.find(std::make_pair(tri[j],tri[i]))->second;
		e->pair = pair;
		pair->pair = e;
	  } else {
		edgeMap.insert(std::make_pair(std::make_pair(tri[i], tri[j]), e));
	  }
	  if (vertexMap.count(tri[j]) > 0) {
		e->v = vertexMap.find(tri[j])->second;
	  } else {
		Vertex *v = new Vertex();
		this->vertices.push_back(v);
		e->v = v;
		v->e = e;
		v->x = vertices[tri[j]];
		vertexMap.insert(std::make_pair(tri[j],e->v));
	  }
	}
	for (int i=0; i<3; i++) {
	  int j = (i+1)%3;
	  fedges[i]->next = fedges[j];
	  fedges[j]->prev = fedges[i];
	}
	f->e = fedges[0];  
	faces.push_back(f);
  }
  LOOP = true;
}

Mesh::Mesh(const std::vector<SlVector3> &vertices, const std::vector<SlQuad> &quads) {
  std::unordered_map<std::pair<int, int>, HalfEdge *, pairhash > edgeMap;
  std::unordered_map<int, Vertex * > vertexMap;

  for (unsigned int q=0; q<quads.size(); q++) {
	const SlQuad &quad = quads[q];
	Face *f = new Face();
	HalfEdge *fedges[4];
	for (unsigned int i=0; i<4; i++) {
	  unsigned int j = (i+1)%4;
	  HalfEdge *e = new HalfEdge();
	  this->edges.push_back(e);
	  e->f = f;
	  fedges[i] = e;
	  if (edgeMap.count(std::pair<int,int>(quad[j],quad[i]))) {
		HalfEdge *pair = edgeMap.find(std::pair<int,int>(quad[j],quad[i]))->second;
		e->pair = pair;
		pair->pair = e;
	  } else {
		edgeMap.insert(std::make_pair(std::make_pair(quad[i], quad[j]), e));
	  }
	  if (vertexMap.count(quad[j]) > 0) {
		e->v = vertexMap.find(quad[j])->second;
	  } else {
		Vertex *v = new Vertex();
		this->vertices.push_back(v);
		e->v = v;
		v->e = e;
		v->x = vertices[quad[j]];
		vertexMap.insert(std::make_pair(quad[j],e->v));
	  }
	}
	for (int i=0; i<4; i++) {
	  int j = (i+1)%4;
	  fedges[i]->next = fedges[j];
	  fedges[j]->prev = fedges[i];
	}
	f->e = fedges[0];  
	faces.push_back(f);
  }
  LOOP = false;
}

void Mesh::divideEdge(HalfEdge *e) {
  if (!e->v->flag) return;
  HalfEdge *e2, *pair, *e2pair;
  Vertex *v;
  // handle boundary case
  if (e->pair == NULL)
  {
  	  v = new Vertex();
	  vertices.push_back(v);
	  e2 = new HalfEdge();
	  edges.push_back(e2);
	  e2->v = v;
	  v->e = e2;
	  e2->f = e->f;

	  e2->prev = e->prev;
	  e2->next = e;
	  e->prev->next = e2;
	  e->prev = e2;
  	  // std::cout << "compute positions for new vertices at boundary..." << std::endl;
  	  HalfEdge *ep = e2->prev;    // e->prev  next->next
  	  SlVector3 x = 0.5 * (e->v->x + ep->v->x);
	  v->x = x;
	  v->nx = x;
	  v->v = false;
	  return;
  }

  if (!e->pair->v->flag) return;
  pair = e->pair;

  // internal mesh 
  e2 = new HalfEdge();
  e2pair = new HalfEdge();
  edges.push_back(e2);
  edges.push_back(e2pair);
  v = new Vertex();
  vertices.push_back(v);

  // v
  e2->v = v;
  e2pair->v = pair->v;
  pair->v = v;
  v->e = e2;
  e2pair->v->e = e2pair;

  // f
  e2->f = e->f;
  e2pair->f = pair->f;

  // pair
  e2->pair = e2pair;
  e2pair->pair = e2;

  // prev/next
  e2->prev = e->prev;
  e2->next = e;
  e->prev->next = e2;
  e->prev = e2;

  e2pair->next = pair->next;
  e2pair->prev = pair;
  pair->next->prev = e2pair;
  pair->next = e2pair;

  // std::cout << "compute positions for new vertices..." << std::endl;
  SlVector3 x = (3.0/8.0) * (e2pair->v->x + e->v->x);
  if (LOOP) 
  {
	if (e2pair->next->v->flag) {
	  x += (1.0/8.0) * e2pair->next->v->x;
	} else {
	  x += (1.0/8.0) * e2pair->next->next->v->x;
	}
	if (e->next->v->flag) {
	  x += (1.0/8.0) * e->next->v->x;
	} else {
	  assert (e->next->next->v->flag);
	  x += (1.0/8.0) * e->next->next->v->x;
	}
  }
  
  v->x = x;
  v->nx = x;
  v->flag = false;
}

void Mesh::triangulateFace(Face *face) {
  HalfEdge *e = face->e;

  for (int i=0; i<3; i++) {
	if (!e->v->flag) {
	  e = e->next;
	}
	HalfEdge *ne = new HalfEdge();
	HalfEdge *pair = new HalfEdge();
	Face *f = new Face();

	ne->v = e->prev->v;
	pair->v = e->next->v;

	e->f = f;
	ne->f = f;
	e->next->f = f;
	f->e = e;

	pair->f = face;
	face->e = pair;

	ne->next = e;
	ne->prev = e->next;
	ne->pair = pair;


	pair->pair = ne;
	pair->next = e->next->next;
	pair->prev = e->prev;

	e->prev->next = pair;
	e->next->next->prev = pair;
	e->next->next = ne;
	e->prev = ne;

	edges.push_back(ne);
	edges.push_back(pair);
	faces.push_back(f);

	e = pair->next;
  }
}

void Mesh::splitFace(Face *face) {
  HalfEdge *e = face->e;
  if (!e->v->flag) {
	e = e->next;
  }
  Vertex *v = new Vertex();
  vertices.push_back(v);
  v->nx = v->x = 0.25*(e->v->x + e->next->next->v->x + e->next->next->next->next->v->x +
	  e->next->next->next->next->next->next->v->x);

  HalfEdge *pe = new HalfEdge();
  edges.push_back(pe);
  HalfEdge *pp = new HalfEdge();
  edges.push_back(pp);

  face->e = pp;

  pe->v = e->prev->v;
  pp->v = v;
  v->e = pp;

  pe->next = e;
  pe->prev = pp;
  pp->next = pe;
  pp->prev = e->prev;
  pe->pair = pp;
  pp->pair = pe;
  e->prev->next = pp;
  e->prev = pe;
  
  for (int i=0; i<3; i++) {
	Face *f = new Face();
	faces.push_back(f);
	HalfEdge *ne = new HalfEdge();
	edges.push_back(ne);
	HalfEdge *pair = new HalfEdge();
	edges.push_back(pair);

	f->e = e;
	ne->v = v;
	pair->v = e->next->v;
	
	ne->pair = pair;
	pair->pair = ne;
	
	ne->prev = e->next;
	ne->next = pe;
	pair->next = e->next->next;
	pair->prev = pp;

	
	ne->next->prev = ne;
	ne->prev->next = ne;
	pair->next->prev = pair;
	pair->prev->next = pair;
	
	pe = pair;
	e = pair->next;
  }
}

void Mesh::subdivide() { 
  std::cout << "subdividing HalfEdge..." << std::endl;
  int size = vertices.size();
  std::cout << "vertex size: " << size << std::endl;
  for (unsigned int i=0; i<vertices.size(); i++) {
	Vertex *v = vertices[i];
	v->flag = true;
	if (LOOP) {
	  int count = 0;
	  SlVector3 x(0.0);
	  HalfEdge *e = v->e;
		// // e = e->next->pair;
		// if (e == v->e) break;
	 //  }

	 //  std::cout << "compute beta..." << std::endl;
	 //  double beta;
	 //  if (count > 3) {
		// beta = 3.0 / (8.0*count);
	 //  } else {
		// beta = 3.0/16.0;
	 //  }
	  // v->nx = beta * x + (1.0-count*beta) * v->x;
	  v->nx = v->x;
	}
  }
  
  std::cout << "dividing edges..." << std::endl;
  int n = edges.size();
  for (unsigned int i=0; i<n; i++) divideEdge(edges[i]);

  std::cout << "triangulate faces..." << std::endl;
  n = faces.size();
  if (LOOP) {
	for (unsigned int i=0; i<n; i++) triangulateFace(faces[i]);
  } else {
	for (unsigned int i=0; i<n; i++) splitFace(faces[i]);
  }	
  if (LOOP) {
	for (unsigned int i=0; i<vertices.size(); i++) vertices[i]->x = vertices[i]->nx;
  } else {
	for (unsigned int i=0; i<vertices.size(); i++) {
	  Vertex *v = vertices[i];
	  if (!v->flag) continue;
	  int k=0;
	  SlVector3 x(0.0), y(0.0);
	  HalfEdge *e = v->e;
	  while (true) {
		x += e->pair->v->x;
		y += e->pair->next->v->x;
		k++;
		e = e->next->pair;
		if (e == v->e) break;
	  }
	  v->x = v->nx = (6.0 * x + 1.0 * y + (4.0*k*k - 7*k) * v->x) / (4*k*k);
	}
  }
}

void Mesh::writeObj(char *fname) {
  std::cout << "Write to " << fname << std::endl;

  std::ofstream out(fname, std::ios::out);
  for (unsigned int i=0; i<vertices.size(); i++) {
	Vertex *v = vertices[i];
	v->index = i;
	out<<"v "<<v->x[0]<<" "<<v->x[1]<<" "<<v->x[2]<<std::endl;
  }
  for (unsigned int i=0; i<faces.size(); i++) {
	Face *f = faces[i];
	HalfEdge *e = f->e;
	out<<"f ";
	int count = 0;
	// while (true) {
	while (count < 3) {
	  out<<(e->v->index)+1<<" ";
	  e = e->next;
	  // if (e == f->e) break;
	  count++;
	}
	out<<std::endl;
  }
}

void readObject(char *fname, std::vector<SlVector3> &vertices, std::vector<SlTri> &triangles) {
  std::ifstream in(fname, std::ios::in);
  char c;
  SlVector3 pt;
  SlTri t;
  
  if (!in.good()) {
	std::cerr<<"Unable to open file \""<<fname<<"\""<<std::endl;
	abort();
  }
  
  while (in.good()) {
	in >> c;
	if (!in.good()) break;
	if (c == 'v') {
	  in >> pt[0] >> pt[1] >> pt[2];
	  vertices.push_back(pt);
	} else if (c == 'f') {
	  in >> t[0] >> t[1] >> t[2];
	  t[0]-=1; t[1]-=1; t[2]-=1;
	  triangles.push_back(t);
	}
  }

  std::cout << "Read file: " << fname << std::endl;
  std::cout << "Vertices: " << vertices.size() << std::endl;
  std::cout << "Triangles: " << triangles.size() << std::endl;
}

void readObject(char *fname, std::vector<SlVector3> &vertices, std::vector<SlQuad> &quads) {
  std::ifstream in(fname, std::ios::in);
  char c;
  SlVector3 pt;
  SlQuad q;
  
  if (!in.good()) {
	std::cerr<<"Unable to open file \""<<fname<<"\""<<std::endl;
	abort();
  }
  
  while (in.good()) {
	in >> c;
	if (!in.good()) break;
	if (c == 'v') {
	  in >> pt[0] >> pt[1] >> pt[2];
	  vertices.push_back(pt);
	} else if (c == 'f') {
	  in >> q[0] >> q[1] >> q[2] >> q[3];
	  q[0]-=1; q[1]-=1; q[2]-=1; q[3] -=1;
	  quads.push_back(q);
	}
  }
}

int main(int argc, char *argv[]) {
  int c;
  bool LOOP = true;
  std::vector<SlVector3> vertices;
  std::vector<SlTri> triangles;
  std::vector<SlQuad> quads;
  Mesh *m;
  while ((c = getopt(argc, argv, "c")) != -1) {
	switch(c) {
	case 'c':
	  LOOP = false;
	  break;
	default:
	  abort();
	}
  }
  if (LOOP) {
	readObject(argv[optind++], vertices, triangles);
	m = new Mesh(vertices, triangles);
	std::cout << "Created a new mesh." << std::endl;
  } else {
	readObject(argv[optind++], vertices, quads);
	m = new Mesh(vertices, quads);
  }

  std::cout << "Subdivide " << argv[optind+1] << " times." << std::endl;
  for (int i=0; i<atoi(argv[optind+1]); i++) {
  	std::cout << "subdividing " << i << std::endl;
	m->subdivide();
  }
  m->writeObj(argv[optind]);
  delete m;
}

Mesh::~Mesh() {
  for (unsigned int i=0; i<vertices.size(); i++) delete vertices[i];
  for (unsigned int i=0; i<faces.size(); i++) delete faces[i];
  for (unsigned int i=0; i<edges.size(); i++) delete edges[i];
}
