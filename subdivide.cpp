#include "slVector.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <getopt.h>
#include <cstring>
#include <iomanip> 
#include <cstdlib>
#include <experimental/filesystem>

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
	SlVector3 vel;
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
  Mesh(const std::vector<SlVector3> &vertices, const std::vector<SlTri> &triangles, const std::vector<SlVector3> &velocity);
  void divideEdge(HalfEdge *e);
  void triangulateFace(Face *f);
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

Mesh::Mesh(const std::vector<SlVector3> &vertices, const std::vector<SlTri> &triangles, const std::vector<SlVector3> &velocity) {
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
		v->vel = velocity[tri[j]];
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
	  // velocity
	  v->vel = 0.5 * (e->v->vel + ep->v->vel);
	  v->flag = false;
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
  SlVector3 x = (1.0/2.0) * (e2pair->v->x + e->v->x);
 //  SlVector3 x = (3.0/8.0) * (e2pair->v->x + e->v->x);
 //  if (LOOP) 
 //  {
	// if (e2pair->next->v->flag) {
	//   x += (1.0/8.0) * e2pair->next->v->x;
	// } else {
	//   x += (1.0/8.0) * e2pair->next->next->v->x;
	// }
	// if (e->next->v->flag) {
	//   x += (1.0/8.0) * e->next->v->x;
	// } else {
	//   assert (e->next->next->v->flag);
	//   x += (1.0/8.0) * e->next->next->v->x;
	// }
 //  }
  
  v->x = x;
  v->nx = x;
  v->vel = (1.0/2.0) * (e2pair->v->vel + e->v->vel);
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

void Mesh::subdivide() { 
  std::cout << "subdividing HalfEdge..." << std::endl;
  int size = vertices.size();
  std::cout << "vertex size: " << size << std::endl;
  
  for (unsigned int i=0; i<vertices.size(); i++) {
	Vertex *v = vertices[i];
	v->flag = true;
	v->nx = v->x;
  }
  
  std::cout << "dividing edges..." << std::endl;
  int n = edges.size();
  for (unsigned int i=0; i<n; i++) divideEdge(edges[i]);

  std::cout << "triangulate faces..." << std::endl;
  n = faces.size();
 
  for (unsigned int i=0; i<n; i++) triangulateFace(faces[i]);
  
  // for (unsigned int i=0; i<vertices.size(); i++) 
  // {
  // 	vertices[i]->x = vertices[i]->nx;
  // }
  
}

void Mesh::writeObj(char *fname) {
  std::cout << "Write to " << fname << std::endl;

  std::ofstream out(fname, std::ios::out);
  for (unsigned int i=0; i<vertices.size(); i++) {
	Vertex *v = vertices[i];
	v->index = i;
	out<<"v "<<v->x[0]<<" "<<v->x[1]<<" "<<v->x[2]<<std::endl;
	// add velocity
	out << "nv " << v->vel[0] << " " << v->vel[1] << " " << v->vel[2] << std::endl;
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

void readObject(char *fname, std::vector<SlVector3> &vertices, std::vector<SlTri> &triangles, std::vector<SlVector3> &velocity) 
{
  std::ifstream in(fname, std::ios::in);
  SlVector3 pt;
  SlVector3 pvel;
  SlTri t;
  
  if (!in.good()) {
	std::cerr<<"Unable to open file \""<<fname<<"\""<<std::endl;
	abort();
  }

  std::string line;
  while (getline(in, line))
    {        
        if (line.substr(0,2) == "v ")
        {
        	std::istringstream s(line.substr(2));
            s >> pt[0]; s >> pt[1]; s >> pt[2];
            vertices.push_back(pt);
        }
        else if (line.substr(0,2) == "f ")
        {
        	std::vector<unsigned int> tokens;
	  		char* pch = strtok (&line[2]," /");
  	  		while (pch != NULL)
 	  		{
    		  tokens.push_back(atoi(pch));
        	  pch = strtok (NULL, " /");
      		}

      		t[0] = tokens[0] - 1;
     		t[1] = tokens[2] - 1;
     		t[2] = tokens[4] - 1;
	  		triangles.push_back(t);
        }
        else if (line.substr(0,2) == "nv")
        {
        	/* parse velocity */
			std::istringstream s(line.substr(2));
            s >> pvel[0]; s >> pvel[1]; s >> pvel[2];
        	velocity.push_back(pvel);
        }

  }

  std::cout << "Read file: " << fname << std::endl;
  std::cout << "Vertices: " << vertices.size() << std::endl;
  std::cout << "Triangles: " << triangles.size() << std::endl << std::endl;
}

int main(int argc, char *argv[]) {

  if (argc < 4)
  {
  	printf("Usage: ./subdivide input_folder/ output_folder/ subdiv_num");
  	return 0;
  }

  int c;
  bool LOOP = true;

  while ((c = getopt(argc, argv, "c")) != -1) {
	switch(c) {
	case 'c':
	  LOOP = false;
	  break;
	default:
	  abort();
	}
  }

  std::string path = argv[1];
  std::string out_path = argv[2];
  for (auto &iter : std::filesystem::directory_iterator(path)) 
  {
    if( !stdfs::is_regular_file(*iter) )
    {
      // mkdir
      cout << iter << endl;
      std::string new_path = out_path + string(iter);
      system(("mkdir -p "+new_path).c_str());

      // 100 frames in each folder
      int count = 100;
      for (int i = 0; i < count; ++i)
      {
        Mesh *m;
        std::vector<SlVector3> vertices;
        std::vector<SlVector3> velocity;
        std::vector<SlTri> triangles;
        std::ostringstream ss;
        ss << std::setw( 5 ) << std::setfill( '0' ) << i + 1;
        std::string filename = path + ss.str() + "_00.obj";
        readObject(&filename[0], vertices, triangles, velocity); 
        m = new Mesh(vertices, triangles, velocity);
        // std::cout << "Created a new mesh." << std::endl;
        // std::cout << "Subdivide " << argv[3] << " times." << std::endl;
        for (int i=0; i<atoi(argv[3]); i++) {
          std::cout << "subdividing " << i+1 << std::endl;
          m->subdivide();
        }
      
        std::string output = new_path + ss.str() + "_00.obj";
        m->writeObj(&output[0]);
        delete m;
      }
    }
  }

}

Mesh::~Mesh() {
  for (unsigned int i=0; i<vertices.size(); i++) delete vertices[i];
  for (unsigned int i=0; i<faces.size(); i++) delete faces[i];
  for (unsigned int i=0; i<edges.size(); i++) delete edges[i];
}
