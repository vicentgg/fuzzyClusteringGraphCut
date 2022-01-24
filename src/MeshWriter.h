

#ifndef MESHSEGMENTATION_MESHWRITER_H
#define MESHSEGMENTATION_MESHWRITER_H


#include <string>
#include <iostream>
#include <fstream>
#include "model/Mesh.h"

using namespace std;

class MeshWriter
{
public:
    /** Save mesh in filename.off */
    static void writeOff(const std::string& filename, const Mesh& mesh);
};

void MeshWriter::writeOff(const string& filename, const Mesh& mesh)
{
    ofstream file(filename);

    cout << "Writing in " << filename << endl;

    file << "OFF" << endl;
    file << mesh.vertices.size() << " " << mesh.faces.size() << " 0" << endl;

    for (auto &v : mesh.vertices)
        file << v.p.x[0] << " " << v.p.x[1] << " " << v.p.x[2] << endl;

    for (auto &f : mesh.faces)
    {
        file << "3 " << f.indices.x << " " << f.indices.y << " " << f.indices.z << endl;
    }
}

#endif //MESHSEGMENTATION_MESHWRITER_H
