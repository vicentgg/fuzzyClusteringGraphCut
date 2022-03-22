
#include <iostream>
#include "MeshReader.h"
#include "MeshWriter.h"    
#include "model/Mesh.h"

using namespace std;

void obj2off(const string& s1, const string& s2) {
    cout << "Reading file " << s1 << endl;
    Mesh mesh = MeshReader::readObj(s1);
    MeshWriter::writeOff(s2, mesh);
    cout << "Writed file" << s2 << endl;
}

int main(int argc, char **argv)
{
    if (argc < 3)
        return -1;
    //可将obj转换为off
    if(argc == 4) obj2off(argv[1], argv[2]); 
    else {
        cout << "Reading file " << argv[1] << endl;
        Mesh mesh = MeshReader::readOff(argv[1]);
        //1.预处理阶段
        // mesh.getDual();
        // mesh.compDist();
        //2.区域分割
        // mesh.solve();

        mesh.ShapeIndex();

        mesh.writeOff(argv[2]);
    }
    
    //Mesh* res = graphMesh;
    //MeshWriter::writeOff(argv[2], *res);
}
