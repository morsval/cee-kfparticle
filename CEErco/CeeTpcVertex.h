#ifndef CEETPCVERTEX_H
#define CEETPCVERTEX_H

#include <iostream>
#include <vector>
#include <TVector3.h>
#include "FairTask.h"
#include "FairMultiLinkedData.h"

#include <TObject.h>
#include <limits>
#include <assert.h>
#include <cmath>
#include <functional>
#include <fstream>      // std::ifstream
#include <string>

//--------------------------

using namespace std; //for sample purposes
class CeeTpcVertex : public FairMultiLinkedData
{

 public:
    
  CeeTpcVertex();
  
  ~CeeTpcVertex();
  ClassDef(CeeTpcVertex,2)
  
  // Add methods to add values to vectors
  void AddX(double x) { fX.push_back(x); }
  void AddY(double y) { fY.push_back(y); }
  void AddZ(double z) { fZ.push_back(z); }
  void AddXErr(double x_err) { fXErr.push_back(x_err); }
  void AddYErr(double y_err) { fYErr.push_back(y_err); }
  void AddZErr(double z_err) { fZErr.push_back(z_err); }
  void AddRank(int rank) { fRank.push_back(rank); }
  
  // Getter methods
  const std::vector<double>& GetX() const { return fX; }
  const std::vector<double>& GetY() const { return fY; }
  const std::vector<double>& GetZ() const { return fZ; }
  const std::vector<double>& GetXErr() const { return fXErr; }
  const std::vector<double>& GetYErr() const { return fYErr; }
  const std::vector<double>& GetZErr() const { return fZErr; }
  const std::vector<int>& GetRank() const { return fRank; }
  
  // Get size of vectors
  size_t GetSize() const { return fX.size(); }
  
  // Track DCA methods - flat structure for all tracks
  void AddTrackDCA(int vertexIdx, int trackId, double dcaXY, double dcaZ, double dca3D);
  size_t GetNTrackDCA() const { return fTrackIds.size(); }
  const std::vector<int>& GetTrackIds() const { return fTrackIds; }
  const std::vector<int>& GetTrackVertexIdx() const { return fTrackVertexIdx; }
  const std::vector<double>& GetTrackDcaXY() const { return fTrackDcaXY; }
  const std::vector<double>& GetTrackDcaZ() const { return fTrackDcaZ; }
  const std::vector<double>& GetTrackDca3D() const { return fTrackDca3D; }
  
  // Clear all vectors
  void Clear() { 
    fX.clear(); 
    fY.clear(); 
    fZ.clear(); 
    fXErr.clear(); 
    fYErr.clear(); 
    fZErr.clear(); 
    fRank.clear();
    fTrackIds.clear();
    fTrackVertexIdx.clear();
    fTrackDcaXY.clear();
    fTrackDcaZ.clear();
    fTrackDca3D.clear();
  }

  private:

  // Private Data Members - all as vectors
  std::vector<double> fX;      // x coordinates
  std::vector<double> fY;      // y coordinates  
  std::vector<double> fZ;      // z coordinates
  std::vector<double> fXErr;   // x errors
  std::vector<double> fYErr;   // y errors
  std::vector<double> fZErr;   // z errors
  std::vector<int> fRank;      // rank values
  
  // Track DCA info - flat structure (one entry per track)
  std::vector<int> fTrackIds;        // track IDs
  std::vector<int> fTrackVertexIdx;  // vertex index for each track
  std::vector<double> fTrackDcaXY;   // DCA in XY plane
  std::vector<double> fTrackDcaZ;    // DCA in Z direction
  std::vector<double> fTrackDca3D;   // 3D DCA

};

#endif // CEETPCVERTEX_H
