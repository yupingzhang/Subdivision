/*********************************************************************
 * carton.c
 * Authored by Adam Bargteil 2003
 *
 * Copyright 2003, Regents of the University of California 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *********************************************************************/

#include <slVector.H>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <sstream>
#include <vector>

using namespace std;

class SimpleTri {
public:
  unsigned int indices[3];
  inline unsigned int &operator[](const unsigned int &i) { return indices[i];};
  inline unsigned int operator[](const unsigned int &i) const { return indices[i];};
};

class SimpleEdge {
public:
  unsigned int indices[2];
  inline unsigned int &operator[](const unsigned int &i) { return indices[i];};
  inline unsigned int operator[](const unsigned int &i) const { return indices[i];};
};

void subdivide(std::vector<SlVector3> &vertices, 
	       std::vector<SimpleTri> &tris, 
	       int verbose=0);

void usage() {
  cerr<<"Usage:"<<endl;
  cerr<<"iv2rib [options] inputfile outputfile"<<endl;
}

void readObjFile(char *fname, std::vector<SlVector3> &vertices,
		 std::vector<SimpleTri> &triangles) {
  int numVertices=0, numFaces=0;
  bool normals = false, texture = false;
  int tint;
  char ch;
  int p, q, r;
  double x, y, z;
  std::vector<SlVector3>::iterator v;
  std::vector<SimpleTri>::iterator t;
  std::ifstream in1(fname, std::ios::in);
  in1.flags(in1.flags() & ~std::ios::skipws);

  while (in1>>ch) {
    if (ch == 'v') {
      in1>>ch;
      if (ch == ' ') numVertices++;
      else if (ch == 'n') normals = true;
      else if (ch == 't') texture = true;
      else std::cerr<<"error \'"<<ch<<"\'"<<std::endl;
    } else if (ch == 'f') numFaces++;
  }
  in1.close();

  vertices.resize(numVertices);
  triangles.resize(numFaces);
  v = vertices.begin();
  t = triangles.begin();

  std::ifstream in(fname, std::ios::in);
  while (in>>ch) {
    if (ch == 'v') {
      ch = in.peek();
      if (ch != 't' && ch != 'n') {
	in>>x>>y>>z;
	(*v)[0] = x;
	(*v)[1] = y;
	(*v)[2] = z;
	v++;
      }
    } else if (ch == 'f') {
      if (normals && texture) {
	in>>p>>ch>>tint>>ch>>tint>>q>>ch>>tint>>ch>>tint>>r>>ch>>tint>>ch>>tint;
      } else if (normals || texture) {
	in>>p>>ch>>tint>>q>>ch>>tint>>r>>ch>>tint;
      } else {
	in>>p>>q>>r;
      }
      (*t)[0] = p-1;
      (*t)[1] = q-1;
      (*t)[2] = r-1;
      t++;
    }
  }
  in.close();
}

bool NO_AVERAGE = false;

int main (int argc, char *argv[]) {
  extern int optind;

  unsigned int iter = 1;

  while (1) {
    switch(getopt(argc, argv, "i:nh")) {
    case 'i':
      iter = atoi(optarg);
      break;
    case 'n':
      NO_AVERAGE = true;
      break;
    case 'h':
      usage();
      exit(0);
    case '?':
      printf ("Unknown argument %s\n", argv[optind-2]);
      usage();
      exit(0);
    case -1:
      argv += optind-1;
      argc -= optind-1;
      goto end_while;
    }
  }
 end_while:

  cout<<"reading \""<<argv[1]<<"\" ..."<<endl;
  std::vector<SlVector3> vertices;
  std::vector<SimpleTri> tris;

  readObjFile(argv[1], vertices, tris);
  for (unsigned int i=0; i<iter; i++) subdivide(vertices,tris);

  std::ofstream out(argv[2], std::ios::out);

  for (unsigned int i=0; i<vertices.size(); i++) {
    out<<"v "<<vertices[i][0]<<" "<<vertices[i][1]<<" "<<vertices[i][2]<<endl;
  }
  for (unsigned int i=0; i<tris.size(); i++) {
    out<<"f "<<tris[i][0]+1<<" "<<tris[i][1]+1<<" "<<tris[i][2]+1<<endl;
  }

  return 0;
}



#ifdef WIN32
#if _MSC_VER >= 1300
#include <hash_map>
#ifndef HASH_MAP
#define HASH_MAP std::hash_map
#endif
#else
#include <map>
#ifndef HASH_MAP
#define HASH_MAP std::map
#endif
#endif
#else
#include <ext/hash_map>
#ifndef HASH_MAP
#define HASH_MAP __gnu_cxx::hash_map
#endif
#endif

#ifdef WIN32
class simpleEdgeHashSomething : public std::hash_compare<SimpleEdge> {
public:
  inline std::size_t operator()(const SimpleEdge& x) const {
    return (x[0] ^ x[1]);
  };
  inline bool operator()(const SimpleEdge& e1, const SimpleEdge& e2) const {
    int a=0,b=0,c=0,d=0;
    if (e1[0] < e1[1]) {a = e1[0]; b=e1[1];}
    else {a=e1[1]; b=e1[0];};
    if (e2[0] < e2[1]) {c=e2[0]; d=e2[1];}
    else {c=e2[1]; d=e2[0];};
    if (a<c) return true;
    if (c<a) return false;
    return (b<d);
  };
};

typedef HASH_MAP<SimpleEdge, int, simpleEdgeHashSomething> EdgeMap

#else

struct eqSimpleEdge {
  bool operator()(const SimpleEdge &e1, const SimpleEdge &e2) const {
    return ((e1[0] == e2[0] && e1[1] == e2[1]) || (e1[0] == e2[1] && e1[1] == e2[0]));
  }
};

struct hashSimpleEdge {
  size_t operator()(const SimpleEdge& x) const {
    return (x[0] ^ x[1]);
  }
};

typedef HASH_MAP<SimpleEdge, int, hashSimpleEdge, eqSimpleEdge> EdgeMap;

#endif

inline double sweight(int i) {
  if (i == 3) return 3.0/16.0;
  else return 3.0/(8.0*i);
}

void subdivide(std::vector<SlVector3> &vertices, 
	       std::vector<SimpleTri> &tris, 
	       int verbose) {
  std::vector<SlVector3> newVertices;
  std::vector<SimpleTri> newTris;
  std::vector<int> triCount;
  std::vector<SlVector3>::iterator v;
  std::vector<SimpleTri>::iterator t1, t2;
  SimpleEdge e1, e2, e3;
  int v1, v2, v3;
  
#ifdef WIN32
  HASH_MAP<SimpleEdge, int, simpleEdgeHashSomething> edgeMap;
#else
  HASH_MAP<SimpleEdge, int, hashSimpleEdge, eqSimpleEdge> edgeMap;
#endif

  HASH_MAP<int, int> vertMap;

  newVertices.reserve(vertices.size() + 3*tris.size()/2);
  triCount.reserve(vertices.size() + 3*tris.size()/2);
  newTris.resize(4*tris.size());
  t2 = newTris.begin();

  for (t1=tris.begin(); t1!=tris.end(); t1++) {
    v1 = (*t1)[0];
    v2 = (*t1)[1];
    v3 = (*t1)[2];
    e1[0] = v1;
    e1[1] = v2;
    e2[0] = v2;
    e2[1] = v3;
    e3[0] = v3;
    e3[1] = v1;
    
    if (edgeMap.count(e1) == 0) {
      edgeMap[e1] = newVertices.size();
      if (NO_AVERAGE) {
	newVertices.push_back(0.5*(vertices[v1]+vertices[v2]));
      } else {
	newVertices.push_back(0.125*(vertices[v3]) + 
			      0.1875*(vertices[v1]+vertices[v2]));
      }
      triCount.push_back(0);
    }  else {
      triCount[edgeMap[e1]]--;
      if (!NO_AVERAGE) {
	newVertices[edgeMap[e1]] += (0.125*(vertices[v3]) + 
				     0.1875*(vertices[v1]+vertices[v2]));
      }
    }

    if (edgeMap.count(e2) == 0) {
      edgeMap[e2] = newVertices.size();
      if (NO_AVERAGE) {
	newVertices.push_back(0.5*(vertices[v2]+vertices[v3]));
      } else {
	newVertices.push_back(0.125*(vertices[v1]) + 
			      0.1875*(vertices[v2]+vertices[v3]));
      }
      triCount.push_back(0);
    }  else { 
      triCount[edgeMap[e2]]--;
      if (!NO_AVERAGE) {
	newVertices[edgeMap[e2]] += (0.125*(vertices[v1]) + 
				     0.1875*(vertices[v2]+vertices[v3]));
      }
    }

    if (edgeMap.count(e3) == 0) {
      edgeMap[e3] = newVertices.size();
      if (NO_AVERAGE) {
	newVertices.push_back(0.5*(vertices[v3]+vertices[v1]));
      } else {
	newVertices.push_back(0.125*(vertices[v2]) + 
			      0.1875*(vertices[v3]+vertices[v1]));
      }
      triCount.push_back(0);
    }  else {
      triCount[edgeMap[e3]]--;
      if (!NO_AVERAGE) {
	newVertices[edgeMap[e3]] += (0.125*(vertices[v2]) + 
				     0.1875*(vertices[v3]+vertices[v1]));
      }
    }
    
    if (vertMap.count(v1) == 0) {
      vertMap[v1] = newVertices.size();
      newVertices.push_back(vertices[v1]);
      triCount.push_back(1);
    } else {
      triCount[vertMap[v1]]++;
    }
    if (vertMap.count(v2) == 0) {
      vertMap[v2] = newVertices.size();
      newVertices.push_back(vertices[v2]);
      triCount.push_back(1);
    } else {
      triCount[vertMap[v2]]++;
    }
    if (vertMap.count(v3) == 0) {
      vertMap[v3] = newVertices.size();
      newVertices.push_back(vertices[v3]);
      triCount.push_back(1);
    } else {
      triCount[vertMap[v3]]++;
    }
    
    (*t2)[0] = edgeMap[e1];
    (*t2)[1] = edgeMap[e2];
    (*t2)[2] = edgeMap[e3];
    t2++;
    
    (*t2)[0] = edgeMap[e3];
    (*t2)[1] = vertMap[v1];
    (*t2)[2] = edgeMap[e1];
    t2++;
    
    (*t2)[0] = edgeMap[e1];
    (*t2)[1] = vertMap[v2];
    (*t2)[2] = edgeMap[e2];
    t2++;
    
    (*t2)[0] = edgeMap[e2];
    (*t2)[1] = vertMap[v3];
    (*t2)[2] = edgeMap[e3];
    t2++;
  }

  if (!NO_AVERAGE) {
    for (unsigned int i=0; i<newVertices.size(); i++) {
      int &c = triCount[i];
      if (c > 0) {
	newVertices[i] *= 1.0-c*sweight(c);
      } else if (c == 0) {
	cout<<"problem a"<<endl;
      } else if (c < -1) {
	if (c == -3) {
	  cout<<"WARNING!!! four faces adjacent to edge "<<i<<endl;
	  newVertices[i] /= 2.0;
	} else if (c == -5) {
	  cout<<"WARNING!!! six faces adjacent to edge "<<i<<endl;
	  newVertices[i] /= 3.0;
	} else {
	  cout<<"problem b "<<c<<" "<<i<<endl;
	}
      }
    }
    
    for (t1=tris.begin(); t1!=tris.end(); t1++) {
      v1 = (*t1)[0];
      v2 = (*t1)[1];
      v3 = (*t1)[2];
      
      newVertices[vertMap[v1]] += 
	0.5*sweight(triCount[vertMap[v1]])*(vertices[v2]+vertices[v3]);
      
      newVertices[vertMap[v2]] += 
	0.5*sweight(triCount[vertMap[v2]])*(vertices[v3]+vertices[v1]);
      
      newVertices[vertMap[v3]] += 
	0.5*sweight(triCount[vertMap[v3]])*(vertices[v1]+vertices[v2]);
    }
  }


  vertices = newVertices;
  tris = newTris;
}
