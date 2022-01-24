#ifndef MESHSEGMENTATION_MESHREADER_H
#define MESHSEGMENTATION_MESHREADER_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "model/Mesh.h"

using namespace std;

class MeshReader
{
public:

    /** Reading <mesh>.off file and stores result in Mesh class */
    static Mesh readOff(const std::string& filename);
    static Mesh readObj(const std::string& filename);
    static std::vector<std::string> split(const std::string& s, const std::string& delimiter);
};

Mesh MeshReader::readObj(const std::string& filename) {
    ifstream file(filename);
    string line;
    vector<string> tmp;
    int v1, v2, v3;


    Mesh result;

    while(getline(file, line)) {

        if(line.empty() || line[0] == '#') {
            continue;
        }
        else if(line[0] == 'v') {
            tmp = split(line, " ");
            result.vertices.emplace_back(stod(tmp[1]), stod(tmp[2]), stod(tmp[3]));

        }
        else if(line[0] == 'f') {
            tmp = split(line, " ");
            v1 = stoi(tmp[1]); v1--;
            v2 = stoi(tmp[2]); v2--;
            v3 = stoi(tmp[3]); v3--;

            Indices indices(v1, v2, v3); 
            result.faces.emplace_back(indices, result.vertices[indices[0]].p, result.vertices[indices[1]].p, result.vertices[indices[2]].p);

        }
    }
    
    cout << "\tRead all of vertices and faces: " << result.vertices.size() << " " << result.faces.size() << endl;

    return result; 
}

Mesh MeshReader::readOff(const std::string& filename)
{
    ifstream file(filename);

    string line;
    getline(file, line);
    while(line == "OFF" || line.empty() || line[0] == '#')
        getline(file, line);

    cout << "\tRead the title" << endl;

    auto tmp = split(line, " ");
    int vSize = stoi(tmp[0]);
    int fSize = stoi(tmp[1]);

    cout << "\tRead the first line and got amount of vertices and faces: " << vSize << " " << fSize << endl;
    Mesh result;
    result.vertices.resize(vSize);
    result.faces.resize(fSize);

    for (int i = 0; i < vSize; i++)
    {
        getline(file, line);
        tmp = split(line, " ");
        result.vertices[i] = Vertex(stod(tmp[0]), stod(tmp[1]), stod(tmp[2]));
    }

    cout << "\tRead all vertices" << endl;

    for (int i = 0; i < fSize; i++)
    {
        getline(file, line);
        tmp = split(line, " ");
        Indices indices(stoi(tmp[1]), stoi(tmp[2]), stoi(tmp[3]));
        result.faces[i] = Face(indices, result.vertices[indices[0]].p, result.vertices[indices[1]].p, result.vertices[indices[2]].p);
    }

    cout << "\tRead all faces" << endl;

    return result;
}

vector<string> MeshReader::split(const string& s, const string& delimiter)
{
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != string::npos)
    {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

#endif //MESHSEGMENTATION_MESHREADER_H
